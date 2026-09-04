// Regression tests for counts-matrix parsing and design detection.
// Runs the real src/lib/*.ts via Node's built-in TypeScript type-stripping.
import { parseMatrix } from '../src/lib/matrix.ts'
import {
  detectFactors, pairwiseContrasts, withinFactorContrasts, interactionContrast, stripReplicate,
} from '../src/lib/design.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}

const rows = (...lines) => lines.map(l => l.split(','))

console.log('\nTELLING ANNOTATION COLUMNS FROM SAMPLES')
{
  // The exact shape nf-core/rnaseq writes. `gene_name` must NOT become a sample.
  const m = parseMatrix(rows(
    'gene_id,gene_name,Ctrl_Cold_1,Ctrl_Cold_2,KO_Cold_1',
    'ENSMUSG00000000001,Gnai3,1234,1102,980',
    'ENSMUSG00000000003,Pbsn,0,0,2',
  ))
  check('nf-core header -> 3 samples, not 4', m.samples, ['Ctrl_Cold_1', 'Ctrl_Cold_2', 'KO_Cold_1'])
  check('gene_id kept as key', m.geneIds, ['ENSMUSG00000000001', 'ENSMUSG00000000003'])
  check('gene_name kept separately', m.geneNames, ['Gnai3', 'Pbsn'])
  check('annotation columns reported', m.annotationColumns, ['gene_id', 'gene_name'])
  check('csv is gene_id + samples only',
    m.countsCsv.split('\n')[0], 'gene_id,Ctrl_Cold_1,Ctrl_Cold_2,KO_Cold_1')
  check('gene count', m.nGenes, 2)
}
{
  // The simple shape: one gene column, no symbols.
  const m = parseMatrix(rows('gene,WT_1,WT_2,KO_1,KO_2', 'Actb,10,12,20,22', 'Gapdh,5,6,7,8'))
  check('single id column', m.samples, ['WT_1', 'WT_2', 'KO_1', 'KO_2'])
  check('no symbol column -> null', m.geneNames, null)
}
{
  // A symbol column that is all digits must still read as annotation.
  const m = parseMatrix(rows('entrez_id,S1,S2,S3', '11461,10,12,20', '14433,5,6,7'))
  check('numeric-looking id column is annotation', m.samples, ['S1', 'S2', 'S3'])
}
{
  // Empty / NA cells are counts, not text.
  const m = parseMatrix(rows('gene,A1,A2,A3', 'X,1,,3', 'Y,NA,2,3'))
  check('NA and blank stay samples', m.samples, ['A1', 'A2', 'A3'])
  check('missing becomes 0', m.countsCsv.trim().split('\n')[1], 'X,1,0,3')
}
{
  let msg = ''
  try { parseMatrix(rows('gene_id,gene_name', 'A,B')) } catch (e) { msg = e.message }
  check('too few samples is a readable error', /at least 2 samples/.test(msg), true)
}

console.log('\nFACTORIAL DESIGN FROM SAMPLE NAMES')
const BAT = [
  'Ctrl_Cold_1', 'Ctrl_Cold_2', 'Ctrl_Cold_3', 'Ctrl_Cold_4', 'Ctrl_Cold_5', 'Ctrl_Cold_6',
  'Ctrl_Thermo_1', 'Ctrl_Thermo_2', 'Ctrl_Thermo_3',
  'KO_Cold_1', 'KO_Cold_2', 'KO_Cold_3', 'KO_Cold_4', 'KO_Cold_5', 'KO_Cold_6',
  'KO_Thermo_1', 'KO_Thermo_2', 'KO_Thermo_3', 'KO_Thermo_4',
]
{
  check('replicate suffix stripped', stripReplicate('KO_Cold_12'), 'KO_Cold')
  const d = detectFactors(BAT)
  check('two factors found', d.factors.length, 2)
  check('factor 1 levels', d.factors[0].levels, ['Ctrl', 'KO'])
  check('factor 2 levels', d.factors[1].levels, ['Cold', 'Thermo'])
  check('four groups', d.groupLevels.length, 4)
  check('grid fully populated', d.factorial, true)
  // 6/3/6/4 — a real delivery, and unbalanced because one sample was never shipped.
  check('correctly reported as unbalanced', d.balanced, false)
}
{
  const d = detectFactors(BAT)
  const cs = withinFactorContrasts(d, ['Ctrl', 'Thermo'])
  check('four within-factor contrasts', cs.length, 4)
  check('every contrast moves ONE factor', cs.every(c => {
    const n = c.numerator.split('_'), m = c.denominator.split('_')
    return n.filter((v, i) => v !== m[i]).length === 1
  }), true)
  check('contrast ids', cs.map(c => c.id).sort(), [
    'KO_Cold_vs_Ctrl_Cold', 'KO_Thermo_vs_Ctrl_Thermo',
    'Ctrl_Cold_vs_Ctrl_Thermo', 'KO_Cold_vs_KO_Thermo',
  ].sort())
}
{
  const d = detectFactors(BAT)
  const ix = interactionContrast(d, ['Ctrl', 'Thermo'])
  check('interaction offered on a populated grid', !!ix, true)
  check('interaction is marked as such', ix.kind, 'interaction')
}
{
  // One factor only: no interaction, plain pairwise.
  const d = detectFactors(['WT_1', 'WT_2', 'KO_1', 'KO_2'])
  check('single factor -> no factorial', d.factorial, false)
  check('single factor -> no interaction', interactionContrast(d, ['WT']), null)
  check('pairwise against reference',
    pairwiseContrasts(['WT', 'KO', 'DKO'], 'WT').map(c => c.id), ['KO_vs_WT', 'DKO_vs_WT'])
}
{
  // A ragged grid must NOT offer an interaction: the coefficient exists but is
  // not estimated from the cells the label implies.
  const d = detectFactors(['A_X_1', 'A_X_2', 'A_Y_1', 'A_Y_2', 'B_X_1', 'B_X_2'])
  check('missing cell -> not factorial', d.factorial, false)
  check('missing cell -> no interaction', interactionContrast(d, ['A', 'X']), null)
}
{
  // A constant token is not a factor.
  const d = detectFactors(['P1_WT_1', 'P1_WT_2', 'P1_KO_1', 'P1_KO_2'])
  check('constant token dropped', d.factors.length, 1)
  check('the varying token is the factor', d.factors[0].levels, ['WT', 'KO'])
}
{
  // A repeated split is one factor, not two.
  const d = detectFactors(['KO_KO_1', 'KO_KO_2', 'WT_WT_1', 'WT_WT_2'])
  check('duplicate partition collapsed', d.factors.length, 1)
}


/* ------------------------------------------------------------------ *
 * INPUTS THAT WOULD CORRUPT THE RUN ARE REFUSED AT THE DOOR.
 * ------------------------------------------------------------------ */
console.log('\nREFUSED INPUTS')

const throws = (name, rows, rx) => {
  let msg = null
  try { parseMatrix(rows) } catch (e) { msg = String(e.message) }
  const ok = msg !== null && rx.test(msg)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got ${JSON.stringify(msg)}`}`)
}

// Proven in R: counts[, cd$sample] resolves both "s1" columns to the first,
// so the second library's numbers are replaced with no error anywhere.
throws('duplicate sample names are refused',
  [['gene_id', 's1', 's1', 's2'], ['G1', '1', '2', '3'], ['G2', '4', '5', '6']],
  /[Dd]uplicate sample/)
// R refuses these itself, but only after the upload and without naming the file.
throws('duplicate gene ids are refused',
  [['gene_id', 'a', 'b'], ['G1', '1', '2'], ['G1', '3', '4']],
  /[Dd]uplicate gene/)
// Correctly quoted, and every naive splitter downstream still mis-reads it.
throws('a comma in a sample name is refused',
  [['gene_id', 'ctrl,rep1', 'b'], ['G1', '1', '2'], ['G2', '3', '4']],
  /comma, quote or line break/)
throws('so is a quote',
  [['gene_id', 'ct"rl', 'b'], ['G1', '1', '2'], ['G2', '3', '4']],
  /comma, quote or line break/)

check('a clean matrix still parses',
  parseMatrix([['gene_id', 'a', 'b'], ['G1', '1', '2'], ['G2', '3', '4']]).samples, ['a', 'b'])

console.log('\nPROBE')
const big = [['gene_id', 'a', 'b', 'c']]
for (let i = 0; i < 5000; i++) big.push([`G${i}`, String(i), String(i + 1), String(i + 2)])
const pb = parseMatrix(big)
check('the probe is thinned to a bound', pb.probe.nGenes <= 2000, true)
check('and spans every sample', pb.probe.nSamples, 3)
check('the full matrix is untouched', pb.nGenes, 5000)
check('probe values are the counts, not text',
  pb.probe.values[0] >= 0 && Number.isFinite(pb.probe.values[0]), true)

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll matrix/design tests passed\n')
process.exit(failed ? 1 : 0)

