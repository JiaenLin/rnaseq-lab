// Should a factor be fitted separately? Answered from the DATA.
//
// The first version answered it from sample names — "the factor with the most
// levels" — and got two of four designs wrong, silently, by default. These
// cases are the ones it got wrong plus the one it got right.
import { detectFactors } from '../src/lib/design.ts'
import { separation, separationByFactor, suggestBlockFactor, SEPARATION_THRESHOLD } from '../src/lib/blocking.ts'
import { parseMatrix, probeFromCsv } from '../src/lib/matrix.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}
const note = (n, v) => console.log(`  ..   ${n}: ${v}`)

let seed = 42
const rnd = () => { seed = (seed * 1103515245 + 12345) & 0x7fffffff; return seed / 0x7fffffff }
const gauss = () => { let u = 0, v = 0; while (!u) u = rnd(); while (!v) v = rnd(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v) }

/** Counts where `blockEffect` is how strongly factor A shifts the transcriptome. */
function synth(levelsA, levelsB, reps, blockEffect, withinEffect, nGenes = 800) {
  const samples = []
  for (const a of levelsA) for (const b of levelsB) for (let r = 1; r <= reps; r++)
    samples.push({ name: `${a}_${b}_${r}`, a, b })
  const rows = [['gene_id', ...samples.map(s => s.name)]]
  for (let g = 0; g < nGenes; g++) {
    const base = 3 + rnd() * 5
    const aEff = Object.fromEntries(levelsA.map(a => [a, gauss() * blockEffect]))
    const bEff = Object.fromEntries(levelsB.map(b => [b, gauss() * withinEffect]))
    rows.push([`G${g}`, ...samples.map(s =>
      String(Math.max(0, Math.round(2 ** (base + aEff[s.a] + bEff[s.b] + gauss() * 0.3)))))])
  }
  return { rows, samples }
}

console.log('\nBLOCKING SUGGESTION — from the data, not the names')

// Eleven "tissues" x five "ages": the levels of A are wildly apart.
const tissues = ['Liver', 'Kidney', 'Heart', 'BAT', 'iWAT', 'eWAT', 'Gastroc', 'Brain', 'Hypo', 'Lung', 'Marrow']
const ages = ['008w', '026w', '060w', '078w', '104w']
const atlas = synth(tissues, ages, 5, 3.0, 0.15)
const pmA = parseMatrix(atlas.rows)
const dA = detectFactors(pmA.samples)
const sA = suggestBlockFactor(dA, pmA.probe)
note('atlas separations', separationByFactor(dA, pmA.probe).map(x => `${x.name} ${x.separation.toFixed(1)}`).join(', '))
check('a tissue-like factor is suggested', sA?.index, 0)
check('and it is the many-levelled, far-apart one', dA.factors[sA.index].levels.length, 11)
check('its separation clears the threshold', sA.separation > SEPARATION_THRESHOLD, true)

// Two genotypes x six timepoints: the OLD heuristic blocked by timepoint.
const tc = synth(['WT', 'KO'], ['0h', '2h', '6h', '12h', '24h', '48h'], 3, 0.4, 0.4)
const pmT = parseMatrix(tc.rows)
const dT = detectFactors(pmT.samples)
note('time-course separations', separationByFactor(dT, pmT.probe).map(x => `${x.name} ${x.separation.toFixed(1)}`).join(', '))
check('a time course is NOT blocked (the old rule blocked by timepoint)',
  suggestBlockFactor(dT, pmT.probe), null)

// Three doses x two cell lines: the OLD heuristic blocked by dose.
const dose = synth(['D0', 'D1', 'D10'], ['A549', 'HeLa'], 3, 0.3, 0.3)
const pmD = parseMatrix(dose.rows)
const dD = detectFactors(pmD.samples)
check('a dose series is NOT blocked (the old rule blocked by dose)',
  suggestBlockFactor(dD, pmD.probe), null)

// But two genuinely different cell lines SHOULD be, whichever way round.
const lines = synth(['A549', 'Jurkat'], ['DMSO', 'Drug'], 4, 3.0, 0.2)
const pmL = parseMatrix(lines.rows)
const dL = detectFactors(pmL.samples)
check('two far-apart cell lines are blocked even with only 2 levels',
  suggestBlockFactor(dL, pmL.probe)?.index, 0)

check('no probe means no suggestion, i.e. one fit', suggestBlockFactor(dA, undefined), null)
check('a one-factor design is never blocked',
  suggestBlockFactor(detectFactors(['KO_1', 'KO_2', 'WT_1', 'WT_2']), pmA.probe), null)

console.log('\nTHE STATISTIC')
check('separation of a constant factor is 0',
  separation(pmA.probe, pmA.samples.map(() => 'same')), 0)
check('a factor whose levels each hold one sample is 0',
  separation(pmA.probe, pmA.samples.map((_, i) => `s${i}`)), 0)
check('mismatched lengths give 0, not a wrong number',
  separation(pmA.probe, ['a', 'b']), 0)
const sepT = separation(pmA.probe, dA.factors[0].values)
check('separation is deterministic', separation(pmA.probe, dA.factors[0].values), sepT)

console.log('\nPROBE')
check('a probe from the canonical CSV agrees with the parsed one',
  probeFromCsv(pmA.countsCsv).nSamples, pmA.probe.nSamples)
check('and gives the same verdict',
  suggestBlockFactor(dA, probeFromCsv(pmA.countsCsv))?.index, 0)
check('the probe is thinned, not the whole matrix', pmA.probe.nGenes <= 2000, true)

console.log(failed ? `\n${failed} blocking test(s) failed\n` : '\nAll blocking tests passed\n')
process.exit(failed ? 1 : 0)
