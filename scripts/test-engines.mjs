// The R engines, run for real against R — the block loop is not otherwise covered.
//
// `scripts/test-contrasts.mjs` proves the contrast ALGEBRA without R. This
// proves the thing that algebra is handed to: that one fit per block happens,
// that a block's gene filter and dispersion come from that block alone, and
// that the numbers move in the direction blocking exists to fix.
//
// The R text is lifted out of webr.ts verbatim, with only the /work prefix
// repointed at a temp directory, because webr.ts itself cannot be imported
// outside a browser (it reads document.baseURI at module scope).
//
// Skips itself when R or the Bioconductor packages are missing, so it never
// turns a laptop without R into a red build.
import { execFileSync } from 'node:child_process'
import { mkdirSync, writeFileSync, readFileSync, existsSync, rmSync } from 'node:fs'
import { join } from 'node:path'
import { tmpdir } from 'node:os'

const WORK = join(process.env.ENGINE_WORK || tmpdir(), 'rnaseq-lab-engines')
let failed = 0
const check = (name, got, want) => {
  const ok = JSON.stringify(got) === JSON.stringify(want)
  if (!ok) failed++
  console.log(`${ok ? '  ok  ' : '  FAIL'} ${name}${ok ? '' : `\n        got  ${JSON.stringify(got)}\n        want ${JSON.stringify(want)}`}`)
}
const note = (name, v) => console.log(`  ..   ${name}: ${v}`)

function haveR() {
  try {
    const out = execFileSync('Rscript', ['-e',
      'cat(all(c("DESeq2","limma","ashr") %in% rownames(installed.packages())))'],
      { encoding: 'utf8', stdio: ['ignore', 'pipe', 'ignore'] })
    return out.trim() === 'TRUE'
  } catch { return false }
}
if (!haveR()) {
  console.log('\nENGINES\n  --   R with DESeq2 + limma + ashr not found; skipping\n')
  process.exit(0)
}

/* ---------- lift the R out of webr.ts ---------- */
const src = readFileSync(new URL('../src/lib/webr.ts', import.meta.url), 'utf8')
const lift = name => {
  const m = src.match(new RegExp(`const ${name} = String\\.raw\`([\\s\\S]*?)\`\\n`))
  if (!m) throw new Error(`could not lift ${name} out of webr.ts`)
  return m[1]
}
const RECODE = lift('RECODE_R')
const engine = name => lift(name).replace('__RECODE__', RECODE).replaceAll('/work/', `${WORK}/`)

/* ---------- a dataset whose blocks differ in within-group variance ---------- */
// The scenario blocking exists for: two blocks, the SAME true fold changes in
// each, and one of them far noisier. Under one fit the noisy block's variance
// is charged to the quiet one, because DESeq2 has a single dispersion per gene.
const NGENES = 2000, NDE = 200, REPS = 5
const BLOCKS = [{ name: 'Quiet', disp: 0.02 }, { name: 'Noisy', disp: 0.60 }]
const AGES = ['young', 'old']

let seed = 12345
const rnd = () => {
  seed = (seed * 1103515245 + 12345) & 0x7fffffff
  return seed / 0x7fffffff
}
const gauss = () => {
  let u = 0, v = 0
  while (!u) u = rnd()
  while (!v) v = rnd()
  return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * v)
}
// Gamma (Marsaglia-Tsang) then Poisson-ish rounding: overdispersed counts whose
// variance is set by `disp`, which is all this needs.
const gamma = k => {
  if (k < 1) return gamma(k + 1) * Math.pow(rnd(), 1 / k)
  const d = k - 1 / 3, c = 1 / Math.sqrt(9 * d)
  for (;;) {
    const x = gauss(), v = Math.pow(1 + c * x, 3)
    if (v <= 0) continue
    const u = rnd()
    if (Math.log(u) < 0.5 * x * x + d - d * v + d * Math.log(v)) return d * v
  }
}
const nb = (mu, disp) => Math.max(0, Math.round(gamma(1 / disp) * mu * disp))

const samples = []
for (const b of BLOCKS) for (const a of AGES) for (let r = 1; r <= REPS; r++)
  samples.push({ sample: `${b.name}_${a}_${r}`, group: `${b.name}_${a}`, block: b.name, disp: b.disp, age: a })

const groupLevels = BLOCKS.flatMap(b => AGES.map(a => `${b.name}_${a}`))
const rows = [`gene_id,${samples.map(s => `"${s.sample}"`).join(',')}`]
for (let g = 0; g < NGENES; g++) {
  const base = 50 + Math.floor(rnd() * 400)
  const de = g < NDE                       // the same genes move in BOTH blocks
  const cells = samples.map(s => nb(base * (de && s.age === 'old' ? 2 : 1), s.disp))
  rows.push(`"GENE${String(g).padStart(4, '0')}",${cells.join(',')}`)
}
const countsCsv = rows.join('\n') + '\n'

/* ---------- write what runAnalysis writes ---------- */
rmSync(WORK, { recursive: true, force: true })
mkdirSync(WORK, { recursive: true })
const idxOf = g => groupLevels.indexOf(g) + 1
function writeInputs(blocked) {
  writeFileSync(join(WORK, 'counts.csv'), countsCsv)
  writeFileSync(join(WORK, 'coldata.csv'), 'sample,group,block\n' +
    samples.map(s => `"${s.sample}","${s.group}","${blocked ? s.block : ''}"`).join('\n') + '\n')
  writeFileSync(join(WORK, 'levels.txt'), groupLevels.join('\n') + '\n')
  const con = BLOCKS.map(b => ({
    id: `${b.name}_old_vs_${b.name}_young`,
    plus: String(idxOf(`${b.name}_old`)), minus: String(idxOf(`${b.name}_young`)),
    block: blocked ? b.name : '',
  }))
  writeFileSync(join(WORK, 'contrasts.csv'), 'id,plus,minus,block\n' +
    con.map(c => `"${c.id}","${c.plus}","${c.minus}","${c.block}"`).join('\n') + '\n')
  return con
}
// `cat(...)` rather than letting R auto-print, so stdout is the engine's return
// value verbatim — the same string webR's evalRString hands back — and not R's
// `[1] "..."` rendering of it.
const runR = code => {
  const f = join(WORK, 'run.R')
  writeFileSync(f, `cat(${code})\n`)
  return execFileSync('Rscript', ['--vanilla', f], { encoding: 'utf8', maxBuffer: 1 << 28 }).trim()
}
const sigCount = i => {
  const lines = readFileSync(join(WORK, `deg_${i}.csv`), 'utf8').trim().split('\n')
  const head = lines[0].split(',').map(h => h.replace(/"/g, ''))
  const padj = head.indexOf('padj')
  return lines.slice(1).filter(l => {
    const v = l.split(',')[padj]
    return v && v !== 'NA' && Number(v) < 0.05
  }).length
}

console.log('\nENGINES — one fit per block, run against real R')

/* ---------- DESeq2 ---------- */
const DESEQ = engine('DESEQ_R')
writeInputs(true)
const outBlocked = runR(DESEQ)
const fitBlocked = readFileSync(join(WORK, 'fit.txt'), 'utf8').trim().split('\n')
const degBlocked = [sigCount(1), sigCount(2)]

writeInputs(false)
const outPooled = runR(DESEQ)
const fitPooled = readFileSync(join(WORK, 'fit.txt'), 'utf8').trim().split('\n')
const degPooled = [sigCount(1), sigCount(2)]

check('blocked: one fit per block is reported', fitBlocked[0].startsWith('2 fits, one per block'), true)
check('blocked: each block reports its own gene filter',
  fitBlocked.filter(l => /^(Quiet|Noisy):/.test(l)).length, 2)
check('unblocked: still one fit', fitPooled[0].startsWith('one fit'), true)
check('unblocked: one filter line for everything',
  fitPooled.filter(l => /^all:/.test(l)).length, 1)
check('every contrast produced a table', [existsSync(join(WORK, 'deg_1.csv')), existsSync(join(WORK, 'deg_2.csv'))], [true, true])

// The summary line is how runAnalysis learns each contrast's DEG count, and it
// is keyed by contrast id — a block loop that wrote the right files under the
// wrong ids would still leave every count at 0 in the bundle.
const idsIn = out => out.split('|').filter(Boolean).map(kv => kv.slice(0, kv.lastIndexOf('='))).sort()
check('blocked: the summary names both contrasts by id',
  idsIn(outBlocked), ['Noisy_old_vs_Noisy_young', 'Quiet_old_vs_Quiet_young'])
check('unblocked: the same two ids come back',
  idsIn(outPooled), ['Noisy_old_vs_Noisy_young', 'Quiet_old_vs_Quiet_young'])

note('DESeq2 blocked   [Quiet, Noisy]', JSON.stringify(degBlocked))
note('DESeq2 one fit   [Quiet, Noisy]', JSON.stringify(degPooled))
note(`truth is ${NDE} DE genes in each`, '')

// The point of the change: pooling costs the quiet block its power.
check('pooling costs the quiet block DE genes', degPooled[0] < degBlocked[0], true)
check('and the quiet block loses more than the noisy one gains or loses',
  (degBlocked[0] - degPooled[0]) > Math.abs(degBlocked[1] - degPooled[1]), true)
check('blocked recovers most of the truth in the quiet block', degBlocked[0] > NDE * 0.8, true)

/* ---------- the MLE column ---------- */
const head1 = readFileSync(join(WORK, 'deg_1.csv'), 'utf8').split('\n')[0].replace(/"/g, '')
check('the unshrunk MLE ships beside the shrunk estimate',
  ['log2FoldChange', 'lfcSE', 'log2FoldChange_MLE', 'lfcSE_MLE'].every(c => head1.split(',').includes(c)), true)
const r2 = readFileSync(join(WORK, 'deg_1.csv'), 'utf8').trim().split('\n')
const cols = r2[0].replace(/"/g, '').split(',')
const iSh = cols.indexOf('log2FoldChange'), iMle = cols.indexOf('log2FoldChange_MLE')
const shrinkPairs = r2.slice(1).map(l => l.split(',')).filter(c => c[iSh] !== 'NA' && c[iMle] !== 'NA')
check('shrinkage actually moved the estimates',
  shrinkPairs.some(c => Math.abs(Number(c[iSh]) - Number(c[iMle])) > 1e-6), true)
check('and shrinkage pulls toward zero, never away',
  shrinkPairs.every(c => Math.abs(Number(c[iSh])) <= Math.abs(Number(c[iMle])) + 1e-3), true)

/* ---------- normalized_counts depends on blocking ---------- */
writeInputs(true)
runR(DESEQ)
const normBlocked = readFileSync(join(WORK, 'norm.csv'), 'utf8').trim().split('\n')
writeInputs(false)
runR(DESEQ)
const normPooled = readFileSync(join(WORK, 'norm.csv'), 'utf8').trim().split('\n')
// CPM columns sum to 1e6; median-of-ratios columns do not.
const colSum = (lines, j) => lines.slice(1).reduce((a, l) => a + Number(l.split(',')[j]), 0)
check('a blocked run writes CPM, which means the same thing in every block',
  Math.abs(colSum(normBlocked, 2) - 1e6) < 1e3, true)
check('an unblocked run still writes median-of-ratios',
  Math.abs(colSum(normPooled, 2) - 1e6) > 1e3, true)

/* ---------- limma takes the same path ---------- */
const LIMMA = engine('LIMMA_R')
writeInputs(true)
runR(LIMMA)
const fitLimma = readFileSync(join(WORK, 'fit.txt'), 'utf8').trim().split('\n')
check('limma blocks too', fitLimma[0].startsWith('2 fits, one per block'), true)
check('limma reports a per-block filter',
  fitLimma.filter(l => /^(Quiet|Noisy):/.test(l)).length, 2)
note('limma blocked    [Quiet, Noisy]', JSON.stringify([sigCount(1), sigCount(2)]))

/* ---------- a contrast that leaves its block must not be answerable ---------- */
writeFileSync(join(WORK, 'contrasts.csv'), 'id,plus,minus,block\n' +
  `"cross","${idxOf('Quiet_old')}","${idxOf('Noisy_young')}","Quiet"\n`)
let refused = false
try { runR(DESEQ) } catch (e) { refused = /not in block/.test(String(e.stderr || e)) }
check('a contrast naming a group outside its block is refused', refused, true)

rmSync(WORK, { recursive: true, force: true })
console.log(failed ? `\n${failed} engine test(s) failed\n` : '\nAll engine tests passed\n')
process.exit(failed ? 1 : 0)
