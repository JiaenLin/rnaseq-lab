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
