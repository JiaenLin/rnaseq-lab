// Guessing experimental groups from sample-column names.
//
// A counts matrix carries no design — just column names like "517E2+RSL3_r1".
// Forcing the user to hand-assign 69 columns is how mistakes happen, so we infer
// the grouping and let them correct it.
//
// The trick is knowing how much of the name is a replicate suffix. Stripping too
// eagerly merges real groups: "shArf1-1" and "shArf1-2" are two different
// constructs, not replicates 1 and 2 of "shArf1". So rather than one regex, we
// try progressively more aggressive strips and keep the first that yields a
// plausible design — every group with at least 2 samples, and more than one
// group. If none does, each sample stands alone and the user assigns manually.

export interface DetectedGroup {
  name: string
  samples: string[]
}

/**
 * Strip patterns, least aggressive first.
 *
 * The separator classes are `+`, not single characters, and that is load
 * bearing. A real submission had `HFpEF__1` … `HFpEF__12` beside `HFpEF_13` …
 * `HFpEF_23` — someone renamed half the columns and the double underscore
 * survived. With a single-character class those strip to `HFpEF_` and `HFpEF`,
 * which is one condition silently split into two groups of 4 and 16, each
 * perfectly plausible on its own.
 */
const STRIPS: RegExp[] = [
  /[_.\- ]+(?:r|rep|replicate)[_.\- ]*\d+$/i,   // _r1, .rep2, -replicate 3, __rep4
  /[_.\- ]+\d+$/,                                // _1, .2, -3, __12
  /\d+$/,                                        // trailing digits with no separator
]

/**
 * The name with the replicate marker removed AND any separator it left behind.
 *
 * Trailing separators matter for the same reason: `HFpEF__1` under the last
 * strip gives `HFpEF__`, which is a different string from `HFpEF_`, so the two
 * halves of one condition still never meet.
 */
const baseName = (s: string, rx: RegExp) => {
  const t = s.replace(rx, '').replace(/[_.\- ]+$/, '').trim()
  return t || s
}

function groupBy(samples: string[], rx: RegExp | null): DetectedGroup[] {
  const m = new Map<string, string[]>()
  for (const s of samples) {
    const key = rx ? baseName(s, rx) : s
    const list = m.get(key)
    if (list) list.push(s)
    else m.set(key, [s])
  }
  return [...m.entries()].map(([name, list]) => ({ name, samples: list }))
}

/** How many samples a grouping actually places in a group of two or more. */
const replicated = (g: DetectedGroup[]) =>
  g.filter(x => x.samples.length >= 2).reduce((n, x) => n + x.samples.length, 0)

/**
 * Groups inferred from sample names, in first-appearance order.
 * Falls back to one group per sample when no strip produces a usable design.
 *
 * SCORED, NOT FIRST-PAST-THE-POST. The rule used to be "the first strip where
 * EVERY group has at least two samples", and the quantifier is the problem: one
 * unreplicated condition — or one stray numeric column that slipped in as a
 * sample — made every strip fail, and the fallback is one group per sample. A
 * 37-sample file with two obvious conditions arrived as 37 groups of one.
 *
 * So each strip is scored by how many samples it actually places in a replicated
 * group, and the best score wins. Ties go to the LEAST aggressive strip, which
 * is what kept `shArf1-1` and `shArf1-2` apart when a replicate suffix existed
 * elsewhere in the design, and is the whole reason this module is a progression
 * rather than one regex.
 *
 * A singleton is now allowed to be a singleton. It stays its own group, visible
 * and reassignable, instead of taking the other 36 samples down with it.
 */
export function detectGroups(samples: string[]): DetectedGroup[] {
  let best: DetectedGroup[] | null = null
  let bestScore = 0
  for (const rx of STRIPS) {
    const g = groupBy(samples, rx)
    if (g.length < 2 || g.length >= samples.length) continue
    // At least two real groups, or it is not a design.
    if (g.filter(x => x.samples.length >= 2).length < 2) continue
    const score = replicated(g)
    if (score > bestScore) { best = g; bestScore = score }
  }
  return best ?? groupBy(samples, null)
}

/**
 * The group for each sample, in sample order.
 *
 * `detectFactors` used to do this itself with a single regex whose replicate
 * marker was optional — `[_.\- ](?:r|rep|replicate)?\d+$` — which is the greedy
 * rule this module exists to avoid, and it was the LIVE one: nothing outside
 * this file's own test ever called `detectGroups`. Two implementations of one
 * decision, and the careful, documented, tested one was the dead one.
 *
 * So the progressive strip is now what the app runs, and the tests below cover
 * the code that ships.
 */
export function groupsFor(samples: string[]): string[] {
  const detected = detectGroups(samples)
  const byName = new Map<string, string>()
  for (const g of detected) for (const s of g.samples) byName.set(s, g.name)
  return samples.map(s => byName.get(s) ?? s)
}

/** True when detection actually found a design rather than giving up. */
export const isUsableDetection = (groups: DetectedGroup[], samples: string[]) =>
  groups.length > 1 && groups.length < samples.length

// ── Reading sample names a different way ────────────────────────────────────
//
// Inference is a guess, and on a name it has never seen it will be wrong. These
// are the ways a reader can say what the name MEANS, chosen so the common
// corrections are one click and the uncommon one is still possible.

export type NamingRule = 'auto' | 'before-first' | 'before-last' | 'strip-digits' | 'whole' | 'regex'

export const NAMING_RULES: { id: NamingRule; label: string; hint: string }[] = [
  { id: 'auto', label: 'Detect automatically', hint: 'strip a replicate suffix if one is there' },
  { id: 'before-first', label: 'Up to the first _ . - or space', hint: 'HFpEF_old_1 → HFpEF' },
  { id: 'before-last', label: 'Up to the last _ . - or space', hint: 'HFpEF_old_1 → HFpEF_old' },
  { id: 'strip-digits', label: 'Drop trailing numbers', hint: 'HFpEF__12 → HFpEF' },
  { id: 'whole', label: 'The whole name is the group', hint: 'every sample stands alone' },
  { id: 'regex', label: 'Custom pattern to remove…', hint: 'a regular expression, e.g. [-_]S\\d+$' },
]

/**
 * The group for each sample under an explicit rule.
 *
 * `regex` is the escape hatch and is deliberately forgiving: an incomplete
 * pattern, which is what a half-typed one always is, must not throw while
 * someone is still typing it. An unparseable pattern groups nothing and the UI
 * says so, rather than the page going blank.
 */
export function groupsByRule(samples: string[], rule: NamingRule, pattern = ''): string[] {
  const clean = (t: string, fallback: string) =>
    t.replace(/[_.\- ]+$/, '').trim() || fallback
  switch (rule) {
    case 'before-first':
      return samples.map(s => clean(s.split(/[_.\- ]/)[0] ?? s, s))
    case 'before-last':
      return samples.map(s => {
        const i = s.search(/[_.\- ][^_.\- ]*$/)
        return i < 0 ? s : clean(s.slice(0, i), s)
      })
    case 'strip-digits':
      return samples.map(s => clean(s.replace(/[_.\- ]*\d+$/, ''), s))
    case 'whole':
      return [...samples]
    case 'regex': {
      let rx: RegExp | null = null
      try { rx = new RegExp(pattern) } catch { rx = null }
      if (!rx || !pattern) return [...samples]
      return samples.map(s => clean(s.replace(rx!, ''), s))
    }
    default:
      return groupsFor(samples)
  }
}

/** True when a pattern is a regular expression this browser can compile. */
export const isValidPattern = (p: string) => {
  if (!p) return true
  try { new RegExp(p); return true } catch { return false }
}
