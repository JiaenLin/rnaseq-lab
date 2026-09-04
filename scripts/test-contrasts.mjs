// The contrast algebra, checked as arithmetic rather than by eyeballing labels.
//
// Every contrast this app builds is a weight vector over GROUP MEANS. If the
// weights are right the statistics are right, and they can be verified here
// with no R, no webR and no counts: build the vector, apply it to known group
// means, and check the number that comes out is the quantity the label claims.
import {
  detectFactors, withinFactorContrasts, interactionContrast,
  blockedContrasts, suggestBlockFactor,
} from '../src/lib/design.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}
const close = (name, got, want) => {
  const ok = Math.abs(got - want) < 1e-9
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got ${got}  want ${want}`}`)
}

// Mirrors the weight vector the R engines build: UNIT weights, +1 per plus
// group and -1 per minus group, applied over group means.
const applyContrast = (plus, minus, means) =>
  plus.reduce((a, g) => a + means[g], 0) - minus.reduce((a, g) => a + means[g], 0)

const SAMPLES = [
  'Ctrl_Thermo_1', 'Ctrl_Thermo_2', 'Ctrl_Thermo_3',
  'Ctrl_Cold_1', 'Ctrl_Cold_2', 'Ctrl_Cold_3',
  'KO_Thermo_1', 'KO_Thermo_2', 'KO_Thermo_3',
  'KO_Cold_1', 'KO_Cold_2', 'KO_Cold_3',
]
const d = detectFactors(SAMPLES)
const REFS = ['Ctrl', 'Thermo']

// A gene with a known truth: cold raises it by 2 in controls, and the knockout
// blunts that response to 0.5. Genotype itself does nothing at thermoneutrality.
const MEANS = {
  Ctrl_Thermo: 10,
  Ctrl_Cold: 12,     // +2 cold response in control
  KO_Thermo: 10,     // no genotype effect at baseline
  KO_Cold: 10.5,     // +0.5 cold response in KO  -> interaction = 0.5 - 2 = -1.5
}

console.log('\nPAIRWISE CONTRASTS RECOVER THE RIGHT DIFFERENCES')
{
  const cs = withinFactorContrasts(d, REFS)
  const by = Object.fromEntries(cs.map(c => [c.id, c]))
  check('four contrasts', cs.length, 4)

  const g = id => applyContrast([by[id].numerator], [by[id].denominator], MEANS)
  close('cold response in controls  = +2',   g('Ctrl_Cold_vs_Ctrl_Thermo'), 2)
  close('cold response in KO        = +0.5', g('KO_Cold_vs_KO_Thermo'), 0.5)
  close('genotype at thermoneutral  =  0',   g('KO_Thermo_vs_Ctrl_Thermo'), 0)
  close('genotype under cold        = -1.5', g('KO_Cold_vs_Ctrl_Cold'), -1.5)
}

console.log('\nTHE INTERACTION IS A DIFFERENCE OF DIFFERENCES')
{
  const ix = interactionContrast(d, REFS)
  check('offered', !!ix, true)

  // Exactly the plus/minus split App.tsx builds for an interaction.
  const [fa, fb] = d.factors
  const aAlt = fa.levels.find(l => l !== REFS[0])
  const bAlt = fb.levels.find(l => l !== REFS[1])
  const cell = (a, b) => d.groups[fa.values.findIndex((v, k) => v === a && fb.values[k] === b)]
  const plus = [cell(aAlt, bAlt), cell(REFS[0], REFS[1])]
  const minus = [cell(REFS[0], bAlt), cell(aAlt, REFS[1])]

  check('plus  side', plus.sort(), ['Ctrl_Thermo', 'KO_Cold'].sort())
  check('minus side', minus.sort(), ['Ctrl_Cold', 'KO_Thermo'].sort())

  const got = applyContrast(plus, minus, MEANS)
  // (KO_Cold - KO_Thermo) - (Ctrl_Cold - Ctrl_Thermo) = 0.5 - 2 = -1.5
  close('equals (KO cold response) - (Ctrl cold response)', got, -1.5)

  // The identity that makes it an interaction rather than a fold change.
  const koResp = MEANS.KO_Cold - MEANS.KO_Thermo
  const ctrlResp = MEANS.Ctrl_Cold - MEANS.Ctrl_Thermo
  close('same as computing it by hand', got, koResp - ctrlResp)

  // And it is symmetric: the other way round gives the same magnitude.
  const other = (MEANS.KO_Cold - MEANS.Ctrl_Cold) - (MEANS.KO_Thermo - MEANS.Ctrl_Thermo)
  close('symmetric in the two factors', got, other)
}

console.log('\nNO CONTRAST MOVES TWO FACTORS AT ONCE')
{
  // This is the flaw the RNA-seq Service request form has: a KO+Cold vs
  // Ctrl+Thermo comparison is unattributable to either factor.
  const cs = withinFactorContrasts(d, REFS)
  const bad = cs.filter(c => {
    const n = c.numerator.split('_'), m = c.denominator.split('_')
    return n.filter((v, i) => v !== m[i]).length !== 1
  })
  check('none confounds the two factors', bad.map(c => c.id), [])
}

console.log('\nAN EMPTY CELL BLOCKS THE INTERACTION')
{
  // Same design with every KO_Thermo sample gone.
  const partial = SAMPLES.filter(s => !s.startsWith('KO_Thermo'))
  const dp = detectFactors(partial)
  check('grid no longer complete', dp.factorial, false)
  check('interaction withheld', interactionContrast(dp, REFS), null)
}


/* ------------------------------------------------------------------ *
 * BLOCKING — one fit per level, and never a contrast across two.
 *
 * The failure this guards against is silent: a contrast that spans two
 * blocks is still a perfectly well-formed weight vector, and R would
 * happily refuse it halfway through an eleven-fit run, or worse, find
 * both coefficients present and answer a question no single model asked.
 * ------------------------------------------------------------------ */
console.log('\nBLOCKING')

const TISSUES = ['Liver', 'Kidney', 'Heart', 'BAT', 'iWAT', 'eWAT',
  'Gastrocnemius', 'Brain', 'Hypothalamus', 'Lung', 'BoneMarrow']
const AGES = ['008w', '026w', '060w', '078w', '104w']
const ATLAS = TISSUES.flatMap(t => AGES.flatMap(a => [1, 2, 3, 4, 5].map(r => `${t}_${a}_${r}`)))
const atlas = detectFactors(ATLAS)

check('the atlas reads as 11 x 5', atlas.factors.map(f => f.levels.length), [11, 5])
check('and 55 groups over 275 samples',
  [atlas.groupLevels.length, atlas.groups.length], [55, 275])
check('the many-levelled factor is the one suggested', suggestBlockFactor(atlas), 'factor1')

const all = blockedContrasts(atlas, 'factor1', { scheme: 'all-pairs', reference: '008w' })
check('every pair within a tissue is C(5,2) x 11', all.contrasts.length, 110)
check('one fit per tissue', all.blocks.length, 11)
check('and the same 10 in each', all.perBlock, 10)

// The whole point: 110, not C(55,2).
check('the ungrouped appetite would have been 1485', (55 * 54) / 2, 1485)

const tissueOf = g => g.slice(0, g.lastIndexOf('_'))
const ageOf = g => g.slice(g.lastIndexOf('_') + 1)
check('no contrast crosses a tissue',
  all.contrasts.filter(c => tissueOf(c.numerator) !== tissueOf(c.denominator)).length, 0)
check('every contrast is older vs younger',
  all.contrasts.filter(c => AGES.indexOf(ageOf(c.numerator)) <= AGES.indexOf(ageOf(c.denominator))).length, 0)
check('every contrast names the block it is answered in',
  all.contrasts.filter(c => c.block !== tissueOf(c.numerator)).length, 0)
check('ids are unique', new Set(all.contrasts.map(c => c.id)).size, 110)

// Direction follows the LEVEL order, so reversing the ages reverses every
// contrast. That is the control a reader has over the sign of 110 tables.
const reversed = detectFactors(
  TISSUES.flatMap(t => [...AGES].reverse().flatMap(a => [1, 2].map(r => `${t}_${a}_${r}`))))
const revPlan = blockedContrasts(reversed, 'factor1', { scheme: 'all-pairs' })
check('reversing the level order reverses the direction',
  revPlan.contrasts[0].numerator, 'Liver_078w')

const ref = blockedContrasts(atlas, 'factor1', { scheme: 'vs-reference', reference: '008w' })
check('vs-reference is 4 per tissue', ref.contrasts.length, 44)
check('and every denominator is the reference age',
  ref.contrasts.filter(c => ageOf(c.denominator) !== '008w').length, 0)

const consec = blockedContrasts(atlas, 'factor1', { scheme: 'consecutive' })
check('consecutive is also 4 per tissue', consec.contrasts.length, 44)
check('and each spans one step',
  consec.contrasts.filter(c =>
    AGES.indexOf(ageOf(c.numerator)) - AGES.indexOf(ageOf(c.denominator)) !== 1).length, 0)

// The budget guards against the explosion, not against the atlas.
const squeezed = blockedContrasts(atlas, 'factor1', { scheme: 'all-pairs', reference: '008w', budget: 50 })
check('past the budget the scheme downgrades', squeezed.scheme, 'vs-reference')
check('and says what was asked for', squeezed.requested, 'all-pairs')
check('the atlas fits inside the default budget',
  blockedContrasts(atlas, 'factor1', {}).scheme, 'all-pairs')

// Blocking on the OTHER factor is legal and is a different question: five
// fits, one per age, comparing tissues inside each.
const byAge = blockedContrasts(atlas, 'factor2', { scheme: 'all-pairs' })
check('blocking by age gives 5 fits', byAge.blocks.length, 5)
// 5 x C(11,2) = 275, which is over the default budget — so this is the guard
// firing on a real alternative framing rather than on a contrived one, and it
// downgrades instead of emitting 275 tables.
check('every tissue pair in every age is over budget', byAge.scheme, 'vs-reference')
check('so it lands on 10 per age', byAge.perBlock, 10)
check('and C(11,2) is reachable when the budget allows it',
  blockedContrasts(atlas, 'factor2', { scheme: 'all-pairs', budget: 400 }).perBlock, 55)
check('none of which crosses an age',
  byAge.contrasts.filter(c => ageOf(c.numerator) !== ageOf(c.denominator)).length, 0)

// RENAMING A FACTOR MUST NOT MOVE OR LOSE THE BLOCKING.
//
// The App holds the blocking factor by INDEX for this reason. It held a name
// once, and typing "tissue" over "factor1" silently dropped blocking: the
// stored name matched nothing, and the page went from 110 within-tissue
// comparisons back to 95 that include Kidney-vs-Liver, with nothing on screen
// saying why. This checks the half that lives here — that the planner keys off
// the name it is handed, so a renamed design plans identically.
const renamed = { ...atlas, factors: atlas.factors.map((f, i) => ({ ...f, name: ['tissue', 'age'][i] })) }
const byName = blockedContrasts(renamed, 'tissue', { scheme: 'all-pairs', reference: '008w' })
check('a renamed factor still blocks the same way',
  [byName.contrasts.length, byName.blocks.length], [110, 11])
check('and produces exactly the same contrasts',
  byName.contrasts.map(c => c.id), all.contrasts.map(c => c.id))
check('the stale name matches nothing, rather than silently blocking elsewhere',
  blockedContrasts(renamed, 'factor1', {}).contrasts.length, 0)

// Degenerate shapes must return nothing rather than something wrong.
const oneFactor = detectFactors(['KO_1', 'KO_2', 'WT_1', 'WT_2'])
check('a one-factor design cannot be blocked',
  blockedContrasts(oneFactor, 'factor1', {}).contrasts.length, 0)
check('and suggests no blocking', suggestBlockFactor(oneFactor), null)
check('a 2x2 suggests no blocking either', suggestBlockFactor(d), null)
check('an unknown factor name blocks nothing',
  blockedContrasts(atlas, 'nope', {}).contrasts.length, 0)

// A ragged block — one tissue missing an age — still gets its own complete
// set rather than being dropped or padded with a group that has no samples.
const ragged = detectFactors([
  ...['008w', '026w', '060w'].flatMap(a => [1, 2].map(r => `Liver_${a}_${r}`)),
  ...['008w', '026w'].flatMap(a => [1, 2].map(r => `Heart_${a}_${r}`)),
])
const rp = blockedContrasts(ragged, 'factor1', { scheme: 'all-pairs' })
check('a ragged design blocks on what each level actually has',
  [rp.blocks.length, rp.contrasts.filter(c => c.block === 'Liver').length,
    rp.contrasts.filter(c => c.block === 'Heart').length], [2, 3, 1])
check('and names no group that has no samples',
  rp.contrasts.every(c => ragged.groupLevels.includes(c.numerator)
    && ragged.groupLevels.includes(c.denominator)), true)

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll contrast tests passed\n')
process.exit(failed ? 1 : 0)
