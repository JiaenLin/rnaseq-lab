// Regression tests for group detection (npm test, and in CI before deploy).
// Runs the real src/lib/groups.ts via Node's built-in TypeScript type-stripping.
import {
  detectGroups, groupsFor, isUsableDetection, groupsByRule, isValidPattern,
} from '../src/lib/groups.ts'

let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}
const names = ss => detectGroups(ss).map(g => g.name)

console.log('\nREPLICATE SUFFIXES')
check('_r1 style', names(['WT_r1', 'WT_r2', 'KO_r1', 'KO_r2']), ['WT', 'KO'])
check('_1 style', names(['WT_1', 'WT_2', 'KO_1', 'KO_2']), ['WT', 'KO'])
check('.rep style', names(['WT.rep1', 'WT.rep2', 'KO.rep1', 'KO.rep2']), ['WT', 'KO'])
check('no separator', names(['WT1', 'WT2', 'KO1', 'KO2']), ['WT', 'KO'])
check('order follows first appearance',
  names(['KO_r1', 'WT_r1', 'KO_r2', 'WT_r2']), ['KO', 'WT'])

console.log('\nCOMBINATORIAL ARM NAMES')
const arms = ['517E2', '517E2+RSL3', '517E2+RSL3+Fer1', 'shArf1-1', 'shArf1-2', 'shAUTS43-2+CoQ10']
const cols = arms.flatMap(a => [1, 2, 3].map(r => `${a}_r${r}`))
check('every arm recovered', names(cols), arms)
check('"+" survives', names(cols).includes('517E2+RSL3+Fer1'), true)
// The important one: "-1"/"-2" are different constructs, not replicates.
check('shArf1-1 and shArf1-2 stay separate',
  names(cols).filter(n => n.startsWith('shArf1')), ['shArf1-1', 'shArf1-2'])

console.log('\nOVER-STRIPPING GUARD')
check('constructs kept when replicates are explicit',
  names(['shArf1-1_r1', 'shArf1-1_r2', 'shArf1-2_r1', 'shArf1-2_r2']),
  ['shArf1-1', 'shArf1-2'])

console.log('\nGIVING UP CLEANLY')
check('all-distinct names group individually', names(['a', 'b', 'c']), ['a', 'b', 'c'])
check('and are reported as unusable',
  isUsableDetection(detectGroups(['a', 'b', 'c']), ['a', 'b', 'c']), false)
check('a single group is unusable',
  isUsableDetection(detectGroups(['WT_r1', 'WT_r2']), ['WT_r1', 'WT_r2']), false)
check('a real design is usable',
  isUsableDetection(detectGroups(['WT_r1', 'WT_r2', 'KO_r1', 'KO_r2']),
    ['WT_r1', 'WT_r2', 'KO_r1', 'KO_r2']), true)
check('never returns an empty name',
  detectGroups(['1', '2', '3', '4']).every(g => g.name.length > 0), true)

console.log('\nEVERY SAMPLE IS ACCOUNTED FOR')
{
  const all = detectGroups(cols).flatMap(g => g.samples)
  check('no sample lost', all.length, cols.length)
  check('no sample duplicated', new Set(all).size, cols.length)
}


console.log('\nTHE 2026-10-02 SUBMISSION: inconsistent separators and a stray column')
{
  // Real file, `副本gene_count.csv`. Someone renamed half the HFpEF columns and
  // the double underscore survived; the export also carried gene_start and
  // gene_end, which are integers and so were read as two extra "samples".
  const real = [
    'CL_1','CL_10','CL_11','CL_12','CL_13','CL_14','CL_15','CL_16','CL_17','CL_18',
    'CL_19','CL_2','CL_20','CL_3','CL_7','CL_8','CL_9',
    'HFpEF__1','HFpEF__10','HFpEF__11','HFpEF__12','HFpEF_13','HFpEF_14','HFpEF_15',
    'HFpEF_16','HFpEF_17','HFpEF_18','HFpEF_19','HFpEF_22','HFpEF_23','HFpEF_3',
    'HFpEF_4','HFpEF_5','HFpEF_6','HFpEF_7','HFpEF_8','HFpEF_9',
  ]
  const g = detectGroups(real)
  check('two conditions, not 37 singletons', g.length, 2)
  check('and HFpEF is ONE group despite the double underscore',
    g.map(x => `${x.name}:${x.samples.length}`).sort(), ['CL:17', 'HFpEF:20'])

  // The all-or-nothing rule: one unreplicated sample used to take the rest down.
  const withSingleton = [...real, 'Sham_1']
  const g2 = detectGroups(withSingleton)
  check('a singleton stays a singleton instead of collapsing the design',
    g2.length, 3)
  check('and the replicated groups survive intact',
    g2.filter(x => x.samples.length >= 2).map(x => x.name).sort(), ['CL', 'HFpEF'])
}

console.log('\nTIES GO TO THE LEAST AGGRESSIVE STRIP')
{
  // The shArf1 case this module exists for: a real replicate suffix elsewhere
  // must not license merging two constructs that differ only by a number.
  const g = detectGroups(['shArf1_rep1', 'shArf1_rep2', 'shCtrl_rep1', 'shCtrl_rep2'])
  check('replicate suffixes group on the first strip',
    g.map(x => x.name).sort(), ['shArf1', 'shCtrl'])
}

console.log('\nSAYING WHAT THE NAME MEANS')
{
  const s = ['HFpEF__1', 'HFpEF_13', 'CL_1', 'CL_20']
  check('before the first separator', groupsByRule(s, 'before-first'),
    ['HFpEF', 'HFpEF', 'CL', 'CL'])
  check('drop trailing numbers', groupsByRule(s, 'strip-digits'),
    ['HFpEF', 'HFpEF', 'CL', 'CL'])
  check('whole name keeps every sample apart', groupsByRule(s, 'whole'), s)
  check('a custom pattern', groupsByRule(s, 'regex', '[-_]+\\d+$'),
    ['HFpEF', 'HFpEF', 'CL', 'CL'])

  // Two-part names where the two rules genuinely differ.
  const t = ['HFpEF_old_1', 'HFpEF_old_2', 'HFpEF_young_1']
  check('before-first collapses the second factor', groupsByRule(t, 'before-first'),
    ['HFpEF', 'HFpEF', 'HFpEF'])
  check('before-last keeps it', groupsByRule(t, 'before-last'),
    ['HFpEF_old', 'HFpEF_old', 'HFpEF_young'])

  // A half-typed pattern must not throw while someone is still typing it.
  check('an unfinished pattern is survivable', groupsByRule(s, 'regex', '[-_'), s)
  check('and is reported as invalid', isValidPattern('[-_'), false)
  check('a good one is valid', isValidPattern('[-_]+\\d+$'), true)
  check('an empty one is not an error', isValidPattern(''), true)
}

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll group tests passed\n')
process.exit(failed ? 1 : 0)
