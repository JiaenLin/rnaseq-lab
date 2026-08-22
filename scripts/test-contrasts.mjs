// The contrast algebra, checked as arithmetic rather than by eyeballing labels.
//
// Every contrast this app builds is a weight vector over GROUP MEANS. If the
// weights are right the statistics are right, and they can be verified here
// with no R, no webR and no counts: build the vector, apply it to known group
// means, and check the number that comes out is the quantity the label claims.
import {
  detectFactors, withinFactorContrasts, interactionContrast,
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

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll contrast tests passed\n')
process.exit(failed ? 1 : 0)
