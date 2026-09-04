// The pool's bookkeeping. A race here does not crash — it mislabels.
//
// Every failure this guards against is silent: a table filed under the wrong
// contrast, a sample's column read from another block, a library size computed
// from the wrong denominator. None of them throw, and all of them ship.
import { poolSize, blocksOf, splitByBlock, cpmCsv, mapLimit, MAX_WORKERS } from '../src/lib/parallel.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}

console.log('\nPOOL SIZE')
check('one block never gets a pool', poolSize(1, 12), 1)
check('never more workers than blocks', poolSize(3, 12), 3)
check('leaves the machine a core', poolSize(11, 4), 3)
check('capped whatever the machine has', poolSize(11, 64), MAX_WORKERS)
check('the atlas on a 12-core machine', poolSize(11, 12), 6)
check('a single-core machine still runs', poolSize(11, 1), 1)
check('a browser that reports nothing still runs', poolSize(11, 0), 1)

console.log('\nSPLITTING THE MATRIX')
const SAMPLES = [
  { sample: 'Lv_a', group: 'Lv_a', block: 'Liver' },
  { sample: 'Lv_b', group: 'Lv_b', block: 'Liver' },
  { sample: 'He_a', group: 'He_a', block: 'Heart' },
  { sample: 'He_b', group: 'He_b', block: 'Heart' },
]
// Column order in the matrix is deliberately NOT the order of `samples`.
const CSV = 'gene_id,He_a,Lv_a,He_b,Lv_b\nG1,11,12,13,14\nG2,21,22,23,24\nG3,0,0,0,0\n'
const CONTRASTS = [
  { id: 'lv', label: 'lv', numerator: 'Lv_b', denominator: 'Lv_a', kind: 'pairwise', plus: ['Lv_b'], minus: ['Lv_a'], block: 'Liver' },
  { id: 'he', label: 'he', numerator: 'He_b', denominator: 'He_a', kind: 'pairwise', plus: ['He_b'], minus: ['He_a'], block: 'Heart' },
]
check('blocks in first-appearance order', blocksOf(SAMPLES), ['Liver', 'Heart'])
const jobs = splitByBlock(CSV, SAMPLES, CONTRASTS)
check('one job per block', jobs.map(j => j.block), ['Liver', 'Heart'])
check('each job carries only its own samples', jobs.map(j => j.samples), [['Lv_a', 'Lv_b'], ['He_a', 'He_b']])
// The whole point of the split, and the easiest thing to get wrong: columns are
// selected by NAME, so a matrix whose column order differs from the sample list
// must still hand each block its own numbers.
check('Liver gets the Liver columns, by name not position',
  jobs[0].countsCsv, 'gene_id,Lv_a,Lv_b\nG1,12,14\nG2,22,24\nG3,0,0\n')
check('Heart gets the Heart columns',
  jobs[1].countsCsv, 'gene_id,He_a,He_b\nG1,11,13\nG2,21,23\nG3,0,0\n')
check('every job keeps every gene row', jobs.map(j => j.countsCsv.trim().split('\n').length - 1), [3, 3])
check('contrasts are filed to their block', jobs.map(j => j.contrasts.map(c => c.request.id)), [['lv'], ['he']])
check('and remember where they sat in the original list',
  jobs.map(j => j.contrasts.map(c => c.index)), [[0], [1]])

// A block nobody asked a question about is not worth a worker.
const noWork = splitByBlock(CSV, SAMPLES, [CONTRASTS[1]])
check('a block with no contrasts gets no job', noWork.map(j => j.block), ['Heart'])

console.log('\nCPM ASSEMBLED OUTSIDE R')
// Column totals over ALL genes: 11+21+0=32, 12+22=34, 13+23=36, 14+24=38.
const cpm = cpmCsv(CSV, new Set(['G1', 'G2']))
const rows = cpm.trim().split('\n')
check('header keeps every sample, in matrix order', rows[0], 'gene_id,gene_name,He_a,Lv_a,He_b,Lv_b')
check('only the kept genes appear', rows.slice(1).map(r => r.split(',')[0]), ['G1', 'G2'])
check('G1 is each column over its own library total',
  rows[1], `G1,G1,${(11 / 32 * 1e6).toFixed(0)},${Math.round(12 / 34 * 1e6 * 1000) / 1000},${Math.round(13 / 36 * 1e6 * 1000) / 1000},${Math.round(14 / 38 * 1e6 * 1000) / 1000}`)
// The denominator is the library, not the surviving genes: otherwise the same
// sample reads differently depending on how the run happened to be blocked.
const cpmAll = cpmCsv(CSV, new Set(['G1', 'G2', 'G3']))
check('a filtered gene does not change anyone else\'s value',
  cpmAll.trim().split('\n')[1], rows[1])
check('every column sums to about 1e6 over all genes',
  [2, 3, 4, 5].map(c => Math.round(cpmAll.trim().split('\n').slice(1)
    .reduce((a, r) => a + Number(r.split(',')[c]), 0) / 1000)), [1000, 1000, 1000, 1000])

console.log('\nORDER SURVIVES THE RACE')
// Finish out of order on purpose; results must still come back in input order.
const items = [0, 1, 2, 3, 4, 5, 6, 7]
const delays = [40, 5, 30, 1, 25, 2, 35, 3]
const seen = []
const out = await mapLimit(items, 3, async (n, slot) => {
  await new Promise(r => setTimeout(r, delays[n]))
  seen.push(n)
  return { n, slot }
})
check('results are in input order, not completion order', out.map(o => o.n), items)
check('they really did finish out of order', seen.join(',') === items.join(','), false)
check('no slot index exceeds the limit', out.every(o => o.slot < 3), true)
check('every item ran exactly once', seen.slice().sort((a, b) => a - b), items)

const one = await mapLimit([1, 2, 3], 1, async n => n * 2)
check('a limit of one still works', one, [2, 4, 6])
check('an empty list is fine', await mapLimit([], 4, async n => n), [])

// A worker that throws must not be swallowed into a half-built bundle.
let threw = false
try {
  await mapLimit([1, 2, 3], 2, async n => { if (n === 2) throw new Error('block failed'); return n })
} catch (e) { threw = /block failed/.test(String(e.message)) }
check('a failing job rejects rather than returning a hole', threw, true)

console.log(failed ? `\n${failed} parallel test(s) failed\n` : '\nAll parallel tests passed\n')
process.exit(failed ? 1 : 0)
