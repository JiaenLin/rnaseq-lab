// The bundle handed to RNA-seq Studio.
//
// meta.json is the documented contract between the two apps and had no test at
// all, which is how `control` came to be read off whichever contrast happened to
// sort first instead of off the reference the reader chose.

import { detectFactors, withinFactorContrasts, interactionContrast, pairwiseContrasts }
  from '../src/lib/design.ts'
import { buildBundleFiles, referenceGroup } from '../src/lib/bundle.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}`
    + (ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`))
}

const DEG_HEADER = 'gene_id,gene_name,baseMean,log2FoldChange,lfcSE,pvalue,padj\n'

/** A design, its contrasts, and the bundle they produce. */
function build(samples, factorNames, refs) {
  const d = detectFactors(samples)
  const named = { ...d, factors: d.factors.map((f, i) => ({ ...f, name: factorNames[i] ?? f.name })) }
  const within = named.factors.length > 1 ? withinFactorContrasts(named, refs) : []
  const base = within.length ? within : pairwiseContrasts(named.groupLevels, refs[0])
  const ix = named.factors.length > 1 ? interactionContrast(named, refs) : null
  const specs = ix ? [...base, ix] : base
  const contrasts = specs.map(c => ({ ...c, plus: [c.numerator], minus: [c.denominator] }))
  const countsCsv = `gene_id,${samples.join(',')}\nA,${samples.map(() => 1).join(',')}\n`
  const input = {
    countsCsv,
    samples: samples.map((s, i) => {
      const rec = { sample: s, group: named.groups[i] }
      named.factors.forEach(f => { rec[f.name] = f.values[i] })
      return rec
    }),
    groupLevels: named.groupLevels, contrasts, method: 'limma',
  }
  const result = {
    contrasts: contrasts.map(c => ({ ...c, degCsv: DEG_HEADER, nDeg: 0 })),
    normCsv: countsCsv,
  }
  const files = buildBundleFiles(input, result, {
    project: 'p', species: 'mouse', method: 'limma',
    covariates: named.factors.map(f => f.name),
    control: referenceGroup(named, refs),
  })
  return { named, specs, meta: JSON.parse(new TextDecoder().decode(files['meta.json'])), files }
}

const SAMPLES_2x2 = [
  'Ctrl_Cold_1', 'Ctrl_Cold_2', 'Ctrl_Thermo_1', 'Ctrl_Thermo_2',
  'KO_Cold_1', 'KO_Cold_2', 'KO_Thermo_1', 'KO_Thermo_2',
]

console.log('\nTHE REFERENCE THE READER CHOSE IS THE ONE THE BUNDLE RECORDS')
{
  // The report: on a 2x2, picking Thermo as the temperature reference still
  // produced a bundle whose control was Ctrl_Cold. `control` was read off
  // whichever pairwise contrast sorted first, and the contrast list is built by
  // walking each factor's levels in the order they appear in the sample names —
  // so it was "the first level of the second factor", which is not a choice
  // anybody made.
  const warm = build(SAMPLES_2x2, ['genotype', 'temperature'], ['Ctrl', 'Thermo'])
  check('Thermo chosen -> Ctrl_Thermo', warm.meta.control, 'Ctrl_Thermo')

  const cold = build(SAMPLES_2x2, ['genotype', 'temperature'], ['Ctrl', 'Cold'])
  check('Cold chosen -> Ctrl_Cold', cold.meta.control, 'Ctrl_Cold')

  const ko = build(SAMPLES_2x2, ['genotype', 'temperature'], ['KO', 'Thermo'])
  check('the other genotype too', ko.meta.control, 'KO_Thermo')

  // The bug was invisible in the contrast list, which was right all along —
  // which is why it survived: every table in the bundle was correct and only
  // the field naming the baseline was wrong.
  //
  // A 2x2 asks FOUR pairwise questions, and only two of them are about
  // temperature. The genotype ones hold temperature fixed at each of its levels
  // in turn — that is the point of them — so Ctrl_Cold appearing as a
  // denominator is correct however the temperature reference is set, and is
  // exactly what made "the denominator of the first contrast" look plausible.
  const denoms = b => b.specs.filter(c => c.kind === 'pairwise').map(c => c.denominator).sort()
  check('Thermo as the reference: the temperature contrasts run against Thermo',
    denoms(warm), ['Ctrl_Cold', 'Ctrl_Thermo', 'Ctrl_Thermo', 'KO_Thermo'])
  check('Cold as the reference: they run against Cold',
    denoms(cold), ['Ctrl_Cold', 'Ctrl_Cold', 'Ctrl_Thermo', 'KO_Cold'])
  // Said as the property rather than the list: every contrast that VARIES a
  // factor uses that factor's chosen reference.
  const tempContrasts = b => b.specs.filter(c => /Cold vs Thermo|Thermo vs Cold/.test(c.label))
  check('every temperature contrast ends at Thermo',
    tempContrasts(warm).every(c => c.denominator.endsWith('Thermo')), true)
  check('and at Cold when Cold is chosen',
    tempContrasts(cold).every(c => c.denominator.endsWith('Cold')), true)

  // Whatever the control is, it has to be a group that exists — the studio
  // seeds its comparison bar from it, and a name no sample carries selects
  // nothing.
  for (const b of [warm, cold, ko]) {
    check(`${b.meta.control} is a real group`, b.meta.conditions.includes(b.meta.control), true)
  }
}

console.log('\nONE FACTOR, AND THE ODD SHAPES')
{
  const one = build(['WT_1', 'WT_2', 'KO_1', 'KO_2'], [], ['KO'])
  check('a single factor uses the level itself', one.meta.control, 'KO')

  // A design whose labels do not split uniformly recovers no factors at all;
  // the reference is then just the chosen group.
  const ragged = build(['A_1', 'A_2', 'B_x_1', 'B_x_2'], [], ['A'])
  check('a ragged design falls back to the chosen level', ragged.meta.control, 'A')

  // referenceGroup must never invent a name. If the combination of references
  // names no real group, the first group is a safer answer than a label
  // nothing carries.
  const d = detectFactors(SAMPLES_2x2)
  const named = { ...d, factors: d.factors.map((f, i) => ({ ...f, name: `f${i}` })) }
  check('an impossible combination falls back to a real group',
    named.groupLevels.includes(referenceGroup(named, ['Ctrl', 'NoSuchLevel'])), true)
  check('and no refs at all still gives one',
    named.groupLevels.includes(referenceGroup(named, [])), true)
}

console.log('\nTHE REST OF THE CONTRACT')
{
  const b = build(SAMPLES_2x2, ['genotype', 'temperature'], ['Ctrl', 'Thermo'])
  check('every declared contrast has its table',
    b.meta.contrasts.every(c => b.files[c.deg_file] !== undefined), true)
  check('samples.csv carries the factors as covariates',
    new TextDecoder().decode(b.files['samples.csv']).split('\n')[0],
    'sample,condition,genotype,temperature')
  check('a pairwise contrast is listed first', b.meta.contrasts[0].kind, 'pairwise')
  check('and the interaction is there', b.meta.contrasts.some(c => c.kind === 'interaction'), true)
  check('conditions are the groups', b.meta.conditions.sort(),
    ['Ctrl_Cold', 'Ctrl_Thermo', 'KO_Cold', 'KO_Thermo'])
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll bundle tests passed\n')
process.exit(failed ? 1 : 0)
