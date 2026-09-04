// Reading an experimental design out of sample names.
//
// A counts matrix carries no design. `groups.ts` already recovers ONE factor —
// which group each sample belongs to. That is enough for "KO vs WT" and wrong
// for everything harder, because the studies that need help are factorial:
//
//     KO_Cold_1  KO_Cold_2 …  KO_Thermo_1 …  Ctrl_Cold_1 …  Ctrl_Thermo_1 …
//
// That is 2x2 — genotype by temperature — and it asks four pairwise questions
// plus one that no pairwise contrast can answer: does the effect of genotype
// DEPEND on temperature (the interaction). Collapsing it to two arbitrary
// groups throws away three quarters of the design.
//
// So: split the group labels on their separators, and keep a position as a
// FACTOR when it varies and is not a restatement of another position.

export interface Factor {
  name: string                       // "factor1" until the user renames it
  levels: string[]
  values: string[]                   // level per sample, in sample order
}

export interface Design {
  factors: Factor[]
  groups: string[]                   // full combination per sample, e.g. "KO_Cold"
  groupLevels: string[]
  factorial: boolean                 // >1 factor AND every cell of the grid present
  balanced: boolean
}

/** Replicate suffix, so "KO_Cold_1" contributes the group "KO_Cold". */
const REP_RX = /[_.\- ](?:r|rep|replicate)?\d+$/i

export const stripReplicate = (s: string) => {
  const t = s.replace(REP_RX, '').trim()
  return t || s
}

/**
 * How a column partitions the samples, independent of the level NAMES:
 * ["KO","KO","WT"] and ["a","a","b"] both give "0,0,1". Two positions carrying
 * the same partition are one factor written twice, not two factors.
 */
const partitionKey = (col: string[]) => {
  const seen = new Map<string, number>()
  return col.map(v => {
    if (!seen.has(v)) seen.set(v, seen.size)
    return seen.get(v)
  }).join(',')
}

/**
 * Recover factors from sample names.
 *
 * A position is a factor when it takes more than one value. Positions with a
 * single value are constants — a batch tag, a project code — and are dropped
 * rather than offered as a one-level factor, because a factor that cannot vary
 * cannot be tested and showing it invites someone to try.
 */
export function detectFactors(samples: string[]): Design {
  const groups = samples.map(stripReplicate)
  const groupLevels = [...new Set(groups)]

  // A design is only recoverable if every label splits the same way.
  const parts = groups.map(g => g.split(/[_.\-+ ]+/).filter(Boolean))
  const width = parts[0]?.length ?? 0
  const uniform = width > 1 && parts.every(p => p.length === width)

  const factors: Factor[] = []
  const takenKeys = new Set<string>()
  if (uniform) {
    for (let i = 0; i < width; i++) {
      const values = parts.map(p => p[i])
      const levels = [...new Set(values)]
      if (levels.length < 2) continue                  // constant: not a factor
      const key = partitionKey(values)
      if (takenKeys.has(key)) continue                 // same split as an earlier position
      takenKeys.add(key)
      factors.push({ name: `factor${factors.length + 1}`, levels, values })
    }
  }

  let factorial = factors.length > 1
  let balanced = false
  if (factorial) {
    const cells = new Map<string, number>()
    for (let s = 0; s < samples.length; s++) {
      const key = factors.map(f => f.values[s]).join(' ')
      cells.set(key, (cells.get(key) ?? 0) + 1)
    }
    const expected = factors.reduce((a, f) => a * f.levels.length, 1)
    factorial = cells.size === expected              // every combination present
    const counts = [...cells.values()]
    balanced = factorial && counts.every(c => c === counts[0])
  }

  return { factors, groups, groupLevels, factorial, balanced }
}

export interface ContrastSpec {
  id: string
  numerator: string
  denominator: string
  label: string
  kind: 'pairwise' | 'interaction'
  coef?: string                      // interaction terms: the DESeq2 coefficient
  /**
   * The block this contrast is answered inside, when the run is blocked.
   *
   * Absent means the whole dataset is one fit, which is the default and the
   * only thing that existed before. See `blockedContrasts`.
   */
  block?: string
}

export const contrastId = (num: string, den: string) =>
  `${num}_vs_${den}`.replace(/[^A-Za-z0-9._+-]+/g, '_')

/** Every level against one reference. What a one-factor design asks. */
export function pairwiseContrasts(levels: string[], reference: string): ContrastSpec[] {
  return levels.filter(l => l !== reference).map(l => ({
    id: contrastId(l, reference),
    numerator: l,
    denominator: reference,
    label: `${l} vs ${reference}`,
    kind: 'pairwise' as const,
  }))
}

/**
 * The contrasts a factorial design actually asks: each factor's effect held
 * WITHIN each level of the others, rather than one arbitrary pair.
 *
 * For 2x2 genotype x temperature that is
 *   KO vs Ctrl at Thermo,  KO vs Ctrl at Cold      (genotype effect, per temperature)
 *   Cold vs Thermo in Ctrl, Cold vs Thermo in KO   (temperature effect, per genotype)
 * — exactly the four a biologist writes down. None of them is the
 * "both factors moved at once" comparison a naive pairing produces, which is
 * unattributable to either factor.
 */
export function withinFactorContrasts(d: Design, refs: string[]): ContrastSpec[] {
  if (d.factors.length < 2) return []
  const out: ContrastSpec[] = []
  const join = (combo: string[]) => combo.join('_')

  for (let fi = 0; fi < d.factors.length; fi++) {
    const f = d.factors[fi]
    const otherIdx = d.factors.map((_, i) => i).filter(i => i !== fi)
    const combos = otherIdx.reduce<string[][]>(
      (acc, oi) => acc.flatMap(prefix => d.factors[oi].levels.map(l => [...prefix, l])), [[]])

    for (const combo of combos) {
      for (const lvl of f.levels) {
        if (lvl === refs[fi]) continue
        const num: string[] = new Array(d.factors.length)
        const den: string[] = new Array(d.factors.length)
        num[fi] = lvl
        den[fi] = refs[fi]
        otherIdx.forEach((oi, k) => { num[oi] = combo[k]; den[oi] = combo[k] })
        out.push({
          id: contrastId(join(num), join(den)),
          numerator: join(num),
          denominator: join(den),
          label: `${lvl} vs ${refs[fi]} (at ${combo.join(', ')})`,
          kind: 'pairwise',
        })
      }
    }
  }
  const seen = new Set<string>()
  return out.filter(c => (seen.has(c.id) ? false : (seen.add(c.id), true)))
}

/**
 * The interaction term — "does factor A's effect depend on factor B".
 *
 * Offered only on a populated grid. On a ragged design the coefficient still
 * exists, but it is estimated from cells that may not both be observed, and
 * reporting it as though it were the same quantity would be wrong.
 */
export function interactionContrast(d: Design, refs: string[]): ContrastSpec | null {
  if (!d.factorial || d.factors.length !== 2) return null
  const [a, b] = d.factors
  const aAlt = a.levels.find(l => l !== refs[0])
  const bAlt = b.levels.find(l => l !== refs[1])
  if (!aAlt || !bAlt) return null
  return {
    id: contrastId(`${aAlt}x${bAlt}`, 'interaction'),
    numerator: `${aAlt}:${bAlt}`,
    denominator: 'interaction',
    label: `Interaction: does the ${aAlt}-vs-${refs[0]} effect differ between ${bAlt} and ${refs[1]}?`,
    kind: 'interaction',
    coef: `${a.name}${aAlt}.${b.name}${bAlt}`,
  }
}

/* ------------------------------------------------------------------ *
 * BLOCKING — one fit per level of a factor, contrasts only within it.
 * ------------------------------------------------------------------ */

/**
 * A factor whose levels must not share a fit.
 *
 * The 2x2 above is the design this app was written for, and one fit over every
 * group is right for it. An organism-wide study is a different shape: eleven
 * tissues x five ages is 55 groups over 275 samples, and putting them in one
 * `~ 0 + grp` asks DESeq2 for a single dispersion per gene spanning brain and
 * brown fat. DESeq2 has no group-specific dispersion — Michael Love: "there is
 * not a way to have group-specific dispersion values with DESeq2" — so the
 * noisy tissues raise the estimate for the quiet ones. On a two-tissue dataset
 * with unequal within-group variance that cost the quiet tissue ~3,000 DE genes
 * down to 4, while the noisy one went ~4,000 up to ~8,000
 * (support.bioconductor.org/p/67202). Love names the condition himself:
 * "subsetting to only pairs of groups for running DESeq() can be useful when
 * the within-group variance is very different across groups."
 *
 * THIS IS NOT the thing that was reverted in `c6cd494`. Narrowing the fit to
 * whatever groups a contrast happened to name made the same comparison change
 * its p-value depending on what else was ticked beside it. A BLOCK is a
 * declared, stable partition of the samples: every contrast inside one comes
 * from that block's single fit, so the invariant survives — one fit, then
 * per-contrast extraction — it is just stated per block rather than once
 * globally. Ticking a different contrast still cannot move an existing one.
 *
 * Blocking is opt-in and off by default. A bundle with four groups must keep
 * behaving exactly as it does today.
 */

/** How the levels inside one block are paired up. */
export type ContrastScheme = 'all-pairs' | 'vs-reference' | 'consecutive'

export interface BlockedPlan {
  contrasts: ContrastSpec[]
  /** Block levels that carried at least two comparable groups, in level order. */
  blocks: string[]
  /** The scheme actually used — `requested`, unless the budget forced a downgrade. */
  scheme: ContrastScheme
  requested: ContrastScheme
  /** Contrasts per block, when every block has the same shape. */
  perBlock: number
}

/** The group label carrying a given level on each named factor, if any sample does. */
function groupWith(d: Design, assign: [number, string][]): string | null {
  for (let s = 0; s < d.groups.length; s++) {
    let ok = true
    for (const [fi, lv] of assign) {
      if (d.factors[fi]?.values[s] !== lv) { ok = false; break }
    }
    if (ok) return d.groups[s]
  }
  return null
}

/**
 * Contrasts within each level of a blocking factor.
 *
 * ORDER IS THE DIRECTION. The groups inside a block are enumerated by the
 * remaining factors' LEVEL order, and a pair is always written later-vs-earlier
 * — `104w vs 008w`, never the reverse. On a time course that makes every
 * log2 fold change read "up with age", so the sign means one thing across all
 * 110 tables. Reorder the factor's levels and the direction follows; that is
 * the control, and it is the same one the reference dropdown already uses.
 *
 * `all-pairs` is the default because a fixed young reference is not the only
 * question a time course asks. 8 weeks is a barely-mature mouse, so 008w→026w
 * carries maturation as well as ageing, and the transition a lifespan study is
 * usually about — 060w→104w — is not any of the four reference-anchored
 * contrasts. Blocking is what makes asking for all of them affordable: within a
 * block it is C(5,2) = 10, and across blocks 110, where the same appetite over
 * the ungrouped 55 groups would be C(55,2) = 1,485.
 *
 * The budget is a guard against that second number, not against the first. Past
 * it the scheme downgrades to `vs-reference` and says so, rather than quietly
 * emitting a bundle nobody can download.
 */
export function blockedContrasts(
  d: Design,
  blockFactor: string,
  opts: { scheme?: ContrastScheme; reference?: string; budget?: number } = {},
): BlockedPlan {
  const requested = opts.scheme ?? 'all-pairs'
  const budget = opts.budget ?? 150
  const empty: BlockedPlan = {
    contrasts: [], blocks: [], scheme: requested, requested, perBlock: 0,
  }

  const bi = d.factors.findIndex(f => f.name === blockFactor)
  if (bi < 0 || d.factors.length < 2) return empty

  const otherIdx = d.factors.map((_, i) => i).filter(i => i !== bi)
  // Every combination of the remaining factors, in their own level order — so
  // the enumeration is the design's order, not the order sample names happen
  // to appear in.
  const combos = otherIdx.reduce<string[][]>(
    (acc, oi) => acc.flatMap(prefix => d.factors[oi].levels.map(l => [...prefix, l])), [[]])

  /** The comparable groups inside one block, in level order. */
  const membersOf = (block: string) => {
    const seen = new Set<string>()
    const out: { group: string; within: string }[] = []
    for (const combo of combos) {
      const g = groupWith(d, [[bi, block], ...otherIdx.map((oi, k) => [oi, combo[k]] as [number, string])])
      if (!g || seen.has(g)) continue
      seen.add(g)
      out.push({ group: g, within: combo.join('_') })
    }
    return out
  }

  const blocks = d.factors[bi].levels.filter(b => membersOf(b).length >= 2)
  if (!blocks.length) return empty

  // Count before building, so a downgrade never materialises the list it is
  // meant to avoid.
  const sizes = blocks.map(b => membersOf(b).length)
  const countFor = (s: ContrastScheme) => sizes.reduce(
    (a, k) => a + (s === 'all-pairs' ? (k * (k - 1)) / 2 : k - 1), 0)
  const scheme: ContrastScheme =
    requested === 'all-pairs' && countFor('all-pairs') > budget ? 'vs-reference' : requested

  const contrasts: ContrastSpec[] = []
  for (const block of blocks) {
    const members = membersOf(block)
    const pairs: [number, number][] = []      // [laterIdx, earlierIdx]
    if (scheme === 'consecutive') {
      for (let i = 0; i + 1 < members.length; i++) pairs.push([i + 1, i])
    } else if (scheme === 'vs-reference') {
      // The reader's reference level when this block has it, else the first —
      // a block missing the reference still gets a complete set of contrasts
      // rather than none.
      const r = Math.max(0, members.findIndex(m => m.within === opts.reference))
      for (let i = 0; i < members.length; i++) if (i !== r) pairs.push([i, r])
    } else {
      for (let i = 0; i < members.length; i++) {
        for (let j = i + 1; j < members.length; j++) pairs.push([j, i])
      }
    }
    for (const [hi, lo] of pairs) {
      const num = members[hi], den = members[lo]
      contrasts.push({
        id: contrastId(num.group, den.group),
        numerator: num.group,
        denominator: den.group,
        label: `${num.within} vs ${den.within} (in ${block})`,
        kind: 'pairwise',
        block,
      })
    }
  }

  const perBlock = sizes.every(k => k === sizes[0])
    ? (scheme === 'all-pairs' ? (sizes[0] * (sizes[0] - 1)) / 2 : sizes[0] - 1)
    : 0
  return { contrasts, blocks, scheme, requested, perBlock }
}

/**
 * Which factor, if any, looks like it wants its own fits.
 *
 * Suggested rather than imposed: this returns a candidate for the UI to offer,
 * and the reader decides. The test is structural — a factor with at least three
 * levels, in a design that has another factor to compare within — because the
 * cost of a wrong guess here is a bundle fitted eleven ways when one would have
 * done, and the reader can see the design better than a heuristic can.
 */
export function suggestBlockFactor(d: Design): string | null {
  if (d.factors.length < 2) return null
  const cand = d.factors
    .filter(f => f.levels.length >= 3)
    .sort((a, b) => b.levels.length - a.levels.length)[0]
  return cand && d.factors.some(f => f !== cand && f.levels.length >= 2) ? cand.name : null
}
