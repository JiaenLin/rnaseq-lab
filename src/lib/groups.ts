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

/** Strip patterns, least aggressive first. */
const STRIPS: RegExp[] = [
  /[_.\- ](?:r|rep|replicate)[_.\- ]?\d+$/i,   // _r1, .rep2, -replicate 3
  /[_.\- ]\d+$/,                                // _1, .2, -3
  /\d+$/,                                       // trailing digits with no separator
]

const baseName = (s: string, rx: RegExp) => {
  const t = s.replace(rx, '').trim()
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

/**
 * Groups inferred from sample names, in first-appearance order.
 * Falls back to one group per sample when no strip produces a usable design.
 */
export function detectGroups(samples: string[]): DetectedGroup[] {
  for (const rx of STRIPS) {
    const g = groupBy(samples, rx)
    if (g.length > 1 && g.length < samples.length && g.every(x => x.samples.length >= 2)) return g
  }
  return groupBy(samples, null)
}

/** True when detection actually found a design rather than giving up. */
export const isUsableDetection = (groups: DetectedGroup[], samples: string[]) =>
  groups.length > 1 && groups.length < samples.length
