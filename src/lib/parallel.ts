// Running the blocks on more than one core.
//
// A block is a fit over its own samples that shares nothing with any other —
// no dispersion, no size factors, no correction. That was a statistical
// decision (see design.ts) and this is its dividend: eleven per-tissue fits are
// eleven independent jobs, and the browser has been doing them one at a time on
// one core because `getWebR` held a single instance.
//
// Each `new WebR()` is its own Web Worker with its own R session and its own
// wasm heap, and they genuinely run in parallel. Measured in Node on a 12-core
// machine, an identical CPU-bound R loop:
//
//     workers   2      4      6      8     12
//     speedup   1.89x  3.50x  5.06x  5.50x  6.69x
//     RSS       500MB  867MB  1.26GB 1.68GB 2.45GB
//
// Near-linear to six, then flat — twelve logical cores are not twelve physical
// ones — and about 200 MB per instance before any package is loaded. So the
// pool is capped rather than sized to the machine: past six the speedup stops
// paying for the memory.

import type { ContrastRequest } from './webr'

/** Base R costs ~200 MB per instance; DESeq2 and a fit push it higher. */
export const MAX_WORKERS = 6

/**
 * How many workers to run.
 *
 * Never more than there are blocks — a seventh worker on six blocks is 200 MB
 * for nothing. One less than the machine reports, so the tab stays responsive
 * and the browser keeps a core for itself. `hardwareConcurrency` is absent or
 * lies on some browsers, so it is treated as a hint with a floor of 1.
 */
export function poolSize(nBlocks: number, cores = navigatorCores(), cap = MAX_WORKERS): number {
  if (nBlocks <= 1) return 1
  return Math.max(1, Math.min(nBlocks, cap, Math.max(1, cores - 1)))
}

const navigatorCores = (): number => {
  const n = typeof navigator !== 'undefined' ? navigator.hardwareConcurrency : undefined
  return typeof n === 'number' && n > 0 ? n : 4
}

/** The blocks in a run, in first-appearance order. */
export function blocksOf(samples: readonly { block?: string }[]): string[] {
  const out: string[] = []
  const seen = new Set<string>()
  for (const s of samples) {
    const b = s.block ?? ''
    if (!seen.has(b)) { seen.add(b); out.push(b) }
  }
  return out
}

export interface BlockJob {
  block: string
  /** Sample names in this block, in matrix-column order. */
  samples: string[]
  countsCsv: string
  /** Contrasts for this block, carrying their index in the ORIGINAL request list. */
  contrasts: { request: ContrastRequest; index: number }[]
}

/**
 * Cut the counts matrix into one CSV per block, in a single pass.
 *
 * Per-block rather than per-worker because a block is the unit a fit is over,
 * and it is the whole point of the split: a worker fitting Liver needs Liver's
 * 25 columns, not all 275. On the 11-tissue atlas that turns one 29 MB matrix
 * into eleven of about 2.6 MB, so the `read.csv` every worker pays is a
 * fraction of the one the single instance paid — the IO shrinks rather than
 * being divided.
 *
 * ONE PASS over the file, not one per block. Splitting 34,514 lines into 275
 * fields costs about 9.5 million field reads; doing that eleven times over
 * would cost more than the parallelism saves.
 */
export function splitByBlock(
  countsCsv: string,
  samples: readonly { sample: string; block?: string }[],
  contrasts: readonly ContrastRequest[],
): BlockJob[] {
  const lines = countsCsv.split('\n')
  while (lines.length && lines[lines.length - 1] === '') lines.pop()
  if (!lines.length) return []

  const header = lines[0].split(',')
  // Column index of each sample name, from the header the matrix actually has.
  const colOf = new Map<string, number>()
  for (let c = 1; c < header.length; c++) {
    colOf.set(header[c].replace(/^"|"$/g, ''), c)
  }

  const blocks = blocksOf(samples)
  const plan = blocks.map(block => {
    const inBlock = samples.filter(s => (s.block ?? '') === block)
    return {
      block,
      samples: inBlock.map(s => s.sample),
      cols: inBlock.map(s => colOf.get(s.sample)).filter((c): c is number => c != null),
      out: [] as string[],
    }
  })

  for (const p of plan) {
    p.out.push(['gene_id', ...p.cols.map(c => header[c])].join(','))
  }
  for (let i = 1; i < lines.length; i++) {
    const cells = lines[i].split(',')
    for (const p of plan) {
      const row = new Array<string>(p.cols.length + 1)
      row[0] = cells[0]
      for (let k = 0; k < p.cols.length; k++) row[k + 1] = cells[p.cols[k]]
      p.out.push(row.join(','))
    }
  }

  return plan.map(p => ({
    block: p.block,
    samples: p.samples,
    countsCsv: p.out.join('\n') + '\n',
    // The ORIGINAL index travels with each contrast, so results can be put back
    // in the order the caller asked for them however the workers finish.
    contrasts: contrasts
      .map((request, index) => ({ request, index }))
      .filter(c => (c.request.block ?? '') === p.block),
  })).filter(j => j.contrasts.length > 0)
}

/**
 * CPM over every sample, computed here rather than in R.
 *
 * No worker holds the whole matrix any more, so no worker can write this. It
 * does not need one: CPM is a per-column quantity — a gene's share of its own
 * library — so it never mixes samples and is the same number whoever computes
 * it. That independence is also why CPM is what a blocked bundle carries:
 * median-of-ratios is defined against a reference built from the samples in one
 * fit, so eleven fits would give eleven incomparable scales.
 *
 * `keep` is the union of what the blocks actually fitted, so the matrix covers
 * genes detectable in at least one block and nothing else.
 */
export function cpmCsv(countsCsv: string, keep: ReadonlySet<string>): string {
  const lines = countsCsv.split('\n')
  while (lines.length && lines[lines.length - 1] === '') lines.pop()
  if (!lines.length) return ''
  const header = lines[0].split(',')
  const nCol = header.length - 1

  const totals = new Float64Array(nCol)
  const rows: { id: string; cells: string[] }[] = []
  for (let i = 1; i < lines.length; i++) {
    const cells = lines[i].split(',')
    const id = cells[0].replace(/^"|"$/g, '')
    // Totals come from EVERY gene, kept or not: a library's size is what was
    // sequenced, not what survived a filter. Using only kept genes would make
    // each column's denominator depend on the blocking, and the same sample
    // would read differently in two bundles of the same data.
    for (let c = 0; c < nCol; c++) totals[c] += Number(cells[c + 1]) || 0
    if (keep.has(id)) rows.push({ id, cells })
  }
  for (let c = 0; c < nCol; c++) if (totals[c] <= 0) totals[c] = 1

  const out: string[] = [['gene_id', 'gene_name', ...header.slice(1)].join(',')]
  for (const r of rows) {
    const vals = new Array<string>(nCol)
    for (let c = 0; c < nCol; c++) {
      const v = ((Number(r.cells[c + 1]) || 0) / totals[c]) * 1e6
      vals[c] = String(Math.round(v * 1000) / 1000)
    }
    out.push([r.cells[0], r.cells[0], ...vals].join(','))
  }
  return out.join('\n') + '\n'
}

/**
 * Run `task` over `items` with at most `limit` in flight.
 *
 * Results come back in the order of `items`, never the order they finished —
 * the bundle's contrast list, its meta.json and its DEG files all have to agree,
 * and a race that reordered them would be invisible until someone read the
 * wrong table under the right label.
 */
export async function mapLimit<T, R>(
  items: readonly T[], limit: number, task: (item: T, slot: number) => Promise<R>,
): Promise<R[]> {
  const out = new Array<R>(items.length)
  let next = 0
  const workers = Array.from({ length: Math.max(1, Math.min(limit, items.length)) }, (_, slot) =>
    (async () => {
      for (;;) {
        const i = next++
        if (i >= items.length) return
        out[i] = await task(items[i], slot)
      }
    })())
  await Promise.all(workers)
  return out
}
