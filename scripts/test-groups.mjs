// Regression tests for group detection (npm test, and in CI before deploy).
// Runs the real src/lib/groups.ts via Node's built-in TypeScript type-stripping.
import { detectGroups, isUsableDetection } from '../src/lib/groups.ts'

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

console.log(failed ? `\n${failed} test(s) failed\n` : '\nAll group tests passed\n')
process.exit(failed ? 1 : 0)
