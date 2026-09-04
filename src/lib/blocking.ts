// Should a factor's levels be fitted separately?
//
// This is a question about the DATA, and the first version of it was answered
// from sample names: "the factor with the most levels". That is wrong, and
// measurably so. Eleven tissues by five ages got tissue, correctly. But two
// genotypes across six timepoints got TIMEPOINT — six separate fits of six
// samples each, throwing away the dispersion shared across the time course for
// no reason — and three doses in two cell lines got DOSE. Level count
// correlates with nothing relevant. Blocking is for levels that are not
// COMPARABLE, not for levels that are merely numerous.
//
// What actually distinguishes them is how far apart the levels sit. Michael
// Love's own condition for splitting is that "the within-group variance is very
// different across groups" (support.bioconductor.org/p/67202), and the readable
// form of that is the ratio of between-level to within-level variance — how
// much of the transcriptome a factor explains. Zhou et al. make the same
// observation about this design in words: the samples "are grouped largely
// according to their organ identity instead of their aging stage".
//
// Measured on the real 275-sample atlas and on a simulated time course:
//
//     atlas       tissue    (11 levels)   520.2      <- fit separately
//     atlas       age       ( 5 levels)     0.2
//     timecourse  timepoint ( 6 levels)     3.0
//     timecourse  genotype  ( 2 levels)     0.4
//
// Three orders of magnitude between the one case that wants blocking and every
// case that does not, so the threshold below sits in a very wide gap rather
// than on a boundary someone tuned.

import type { Design } from './design.ts'
import type { Probe } from './matrix.ts'

/**
 * Above this, a factor's levels are different enough that one dispersion per
 * gene across them is not defensible.
 *
 * 26x below the atlas's tissue and 7x above the loudest non-blocking factor
 * seen. Deliberately not tighter: the cost of missing a suggestion is that the
 * reader ticks a box, and the cost of a false one is a silently worse analysis.
 */
export const SEPARATION_THRESHOLD = 20

/** Genes the statistic is computed on, most variable first. */
const NTOP = 1000

/**
 * Between-level variance over within-level variance, median across genes.
 *
 * On log2(CPM + 1), so a library-size difference is an offset rather than a
 * scale, and on the most variable genes, because the ratio over genes that do
 * not vary is noise divided by noise.
 *
 * The median, not the mean: a handful of markers utterly specific to one level
 * — and every tissue has them — produce ratios in the thousands that would set
 * the answer on their own.
 */
export function separation(probe: Probe, levels: readonly string[]): number {
  const { values, nGenes, nSamples } = probe
  if (nSamples !== levels.length || nGenes < 1) return 0

  const lv = [...new Set(levels)]
  if (lv.length < 2) return 0
  const members = new Map<string, number[]>(lv.map(l => [l, []]))
  levels.forEach((l, j) => members.get(l)!.push(j))
  // A level with one sample carries no within-level variance to contribute.
  const usable = lv.filter(l => members.get(l)!.length >= 2)
  if (usable.length < 2) return 0

  const total = new Float64Array(nSamples)
  for (let g = 0; g < nGenes; g++) {
    const at = g * nSamples
    for (let j = 0; j < nSamples; j++) total[j] += values[at + j]
  }
  for (let j = 0; j < nSamples; j++) if (total[j] <= 0) total[j] = 1

  const logged = new Float64Array(nGenes * nSamples)
  const vars = new Float64Array(nGenes)
  for (let g = 0; g < nGenes; g++) {
    const at = g * nSamples
    let sum = 0
    for (let j = 0; j < nSamples; j++) {
      const v = Math.log2((values[at + j] / total[j]) * 1e6 + 1)
      logged[at + j] = v
      sum += v
    }
    const mean = sum / nSamples
    let ss = 0
    for (let j = 0; j < nSamples; j++) { const d = logged[at + j] - mean; ss += d * d }
    vars[g] = ss
  }

  const order = Array.from({ length: nGenes }, (_, g) => g)
    .filter(g => vars[g] > 0)
    .sort((a, b) => vars[b] - vars[a])
    .slice(0, NTOP)
  if (!order.length) return 0

  const ratios: number[] = []
  for (const g of order) {
    const at = g * nSamples
    let grand = 0, n = 0
    for (const l of usable) for (const j of members.get(l)!) { grand += logged[at + j]; n++ }
    grand /= n
    let between = 0, within = 0, dfW = 0
    for (const l of usable) {
      const m = members.get(l)!
      let lm = 0
      for (const j of m) lm += logged[at + j]
      lm /= m.length
      between += m.length * (lm - grand) ** 2
      for (const j of m) within += (logged[at + j] - lm) ** 2
      dfW += m.length - 1
    }
    if (dfW > 0 && within > 0) ratios.push((between / (usable.length - 1)) / (within / dfW))
  }
  if (!ratios.length) return 0
  ratios.sort((a, b) => a - b)
  const h = ratios.length >> 1
  return ratios.length % 2 ? ratios[h] : (ratios[h - 1] + ratios[h]) / 2
}

export interface BlockSuggestion {
  /** Factor name, as the design currently calls it. */
  name: string
  index: number
  separation: number
}

/** How far apart each factor's levels sit, strongest first. */
export function separationByFactor(design: Design, probe: Probe): BlockSuggestion[] {
  return design.factors
    .map((f, index) => ({ name: f.name, index, separation: separation(probe, f.values) }))
    .sort((a, b) => b.separation - a.separation)
}

/**
 * The factor worth fitting separately, or null.
 *
 * Null WITHOUT a probe, deliberately. There is no honest structural answer, and
 * the version that guessed one silently degraded two of the four designs it was
 * tried on. No suggestion means one fit over everything, which is both the old
 * behaviour and the right default.
 */
export function suggestBlockFactor(
  design: Design, probe?: Probe, threshold = SEPARATION_THRESHOLD,
): BlockSuggestion | null {
  if (!probe || design.factors.length < 2) return null
  const ranked = separationByFactor(design, probe)
  const top = ranked[0]
  if (!top || top.separation < threshold) return null
  // Something must remain to compare INSIDE a block, or blocking leaves nothing
  // to test: a factor is only a block if another factor varies within it.
  const within = design.factors.some((f, i) => i !== top.index && f.levels.length >= 2)
  return within ? top : null
}
