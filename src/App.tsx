import { useMemo, useRef, useState } from 'react'
import Papa from 'papaparse'
import { runAnalysis, type AnalysisResult, type Method } from './lib/webr'
import { buildBundleFiles, zipBundle } from './lib/bundle'
import { detectGroups, isUsableDetection } from './lib/groups'

const EXPLORER_URL = 'https://jiaenlin.github.io/rnaseq-studio/'
type Step = 'upload' | 'design' | 'run' | 'result'

interface Counts { csv: string; samples: string[]; nGenes: number }

/** Assign only the two chosen groups' samples; everything else is excluded. */
function assignByGroups(
  samples: string[],
  a: { name: string; samples: string[] },
  b: { name: string; samples: string[] },
): Record<string, 'A' | 'B' | null> {
  const inA = new Set(a.samples), inB = new Set(b.samples)
  const out: Record<string, 'A' | 'B' | null> = {}
  for (const s of samples) out[s] = inA.has(s) ? 'A' : inB.has(s) ? 'B' : null
  return out
}

export default function App() {
  const [step, setStep] = useState<Step>('upload')
  const [counts, setCounts] = useState<Counts | null>(null)
  const [uploadErr, setUploadErr] = useState<string | null>(null)
  const fileRef = useRef<HTMLInputElement>(null)

  // design (two groups)
  const [nameA, setNameA] = useState('control')
  const [nameB, setNameB] = useState('treatment')
  // null = excluded from this comparison. A counts matrix often holds many
  // arms; forcing every column into A or B would silently pool unrelated ones.
  const [assign, setAssign] = useState<Record<string, 'A' | 'B' | null>>({})
  const [detected, setDetected] = useState<{ name: string; samples: string[] }[]>([])
  const [control, setControl] = useState<'A' | 'B'>('A')

  // run params
  const [method, setMethod] = useState<Method>('limma')
  const [project, setProject] = useState('My RNA-seq analysis')
  const [species, setSpecies] = useState('human')
  const [log, setLog] = useState<string[]>([])
  const [running, setRunning] = useState(false)
  const [result, setResult] = useState<AnalysisResult | null>(null)
  const [zipUrl, setZipUrl] = useState<string | null>(null)
  const [runErr, setRunErr] = useState<string | null>(null)

  const onFile = (f: File | undefined) => {
    if (!f) return
    setUploadErr(null)
    f.text().then(text => {
      const parsed = Papa.parse<string[]>(text.trim(), { skipEmptyLines: true })
      const rows = parsed.data as string[][]
      if (!rows.length || rows[0].length < 3) {
        setUploadErr('Expected a gene × sample matrix (first column = gene, ≥2 sample columns).'); return
      }
      const samples = rows[0].slice(1).map(s => String(s).trim())
      const csv = Papa.unparse(rows)                         // normalize to clean CSV for R
      setCounts({ csv, samples, nGenes: rows.length - 1 })
      // Infer the design from the column names rather than splitting them in
      // half — an arbitrary split is wrong for anything but a 2-group matrix.
      const groups = detectGroups(samples)
      setDetected(groups)
      if (isUsableDetection(groups, samples)) {
        setNameA(groups[0].name)
        setNameB(groups[1].name)
        setAssign(assignByGroups(samples, groups[0], groups[1]))
      } else {
        const a: Record<string, 'A' | 'B' | null> = {}
        samples.forEach((x, i) => { a[x] = i < Math.ceil(samples.length / 2) ? 'A' : 'B' })
        setAssign(a)
      }
      setStep('design')
    }).catch(e => setUploadErr(String(e?.message || e)))
  }

  /** Choosing a detected group renames the slot and re-assigns its samples. */
  const pickDetected = (slot: 'A' | 'B', name: string) => {
    const g = detected.find(x => x.name === name)
    if (!g || !counts) return
    const other = detected.find(x => x.name === (slot === 'A' ? nameB : nameA))
    if (slot === 'A') setNameA(name); else setNameB(name)
    const a = slot === 'A' ? g : other
    const b = slot === 'A' ? other : g
    if (a && b) setAssign(assignByGroups(counts.samples, a, b))
    else setAssign(p => {
      const next = { ...p }
      for (const smp of counts.samples) if (g.samples.includes(smp)) next[smp] = slot
      return next
    })
  }

  const conditionOf = (group: 'A' | 'B') => (group === 'A' ? nameA : nameB)
  const samplesForRun = useMemo(
    () => (counts?.samples || [])
      .filter(s => assign[s] === 'A' || assign[s] === 'B')
      .map(s => ({ sample: s, condition: conditionOf(assign[s] as 'A' | 'B') })),
    [counts, assign, nameA, nameB])
  const nA = Object.values(assign).filter(g => g === 'A').length
  const nB = Object.values(assign).filter(g => g === 'B').length
  const nExcluded = (counts?.samples || []).filter(x => assign[x] !== 'A' && assign[x] !== 'B').length
  const designOk = !!nameA.trim() && !!nameB.trim() && nameA.trim() !== nameB.trim() && nA >= 2 && nB >= 2

  const doRun = async () => {
    if (!counts) return
    setRunning(true); setRunErr(null); setLog([]); setResult(null); setZipUrl(null)
    const onLog = (m: string) => setLog(prev => [...prev, m])
    try {
      const controlCond = conditionOf(control)
      const input = { countsCsv: counts.csv, samples: samplesForRun, control: controlCond, method }
      const res = await runAnalysis(input, onLog)
      const files = buildBundleFiles(input, res, { project, species, method })
      const blob = new Blob([zipBundle(files) as BlobPart], { type: 'application/zip' })
      setZipUrl(URL.createObjectURL(blob))
      setResult(res)
      setStep('result')
    } catch (e: any) {
      onLog('❌ ' + (e?.message || e))
      setRunErr(String(e?.message || e))
    } finally {
      setRunning(false)
    }
  }

  const reset = () => { setStep('upload'); setCounts(null); setResult(null); setZipUrl(null); setLog([]); setRunErr(null) }

  return (
    <div className="mx-auto flex min-h-full max-w-4xl flex-col px-4">
      <header className="flex items-center gap-2 py-4">
        <span className="grid h-9 w-9 place-items-center rounded-lg bg-indigo-500 font-bold text-white">L</span>
        <div>
          <h1 className="text-lg font-semibold leading-none">RNA-seq Lab</h1>
          <p className="text-xs text-slate-400">Differential expression in your browser · nothing is uploaded</p>
        </div>
        {step !== 'upload' && <button className="btn ml-auto" onClick={reset}>Start over</button>}
      </header>

      <Steps step={step} />

      {/* key={step} remounts on each step, which re-fires the @starting-style entrance */}
      <main className="step-enter flex-1 py-4" key={step}>
        {step === 'upload' && (
          <div className="card p-6">
            <h2 className="mb-1 text-base font-semibold">1 · Upload a counts matrix</h2>
            <p className="mb-4 text-sm text-slate-500">
              A CSV/TSV with <b>genes as rows, samples as columns</b>; the first column is the gene id/symbol,
              the first row is sample names. Raw (integer) counts work best.
            </p>
            <button className="btn btn-primary" onClick={() => fileRef.current?.click()}>⭱ Choose counts file</button>
            <input ref={fileRef} type="file" accept=".csv,.tsv,.txt" className="hidden"
              onChange={e => onFile(e.target.files?.[0])} />
            {uploadErr && <p className="mt-3 text-sm text-red-500">{uploadErr}</p>}
            <p className="mt-4 text-xs text-slate-400">Best for small–moderate datasets (runs on your CPU via WebAssembly). Large datasets: use the desktop app.</p>
          </div>
        )}

        {step === 'design' && counts && (
          <div className="space-y-4">
            <div className="card p-6">
              <h2 className="mb-1 text-base font-semibold">2 · Define two groups</h2>
              <p className="mb-4 text-sm text-slate-500">{counts.samples.length} samples · {counts.nGenes.toLocaleString()} genes detected. Assign each sample to a group and pick the control.</p>
              {detected.length > 1 && detected.length < counts.samples.length && (
                <div className="mb-4 rounded-lg border border-slate-200 bg-slate-50 p-3 dark:border-slate-700 dark:bg-slate-800/50">
                  <p className="mb-2 text-xs text-slate-500">
                    <b>{detected.length} groups</b> detected from your sample names. Pick the two to
                    compare — samples in any other group are left out of this run.
                  </p>
                  <div className="grid gap-3 sm:grid-cols-2">
                    <label className="text-sm">Group A
                      <select
                        className="input mt-1 w-full py-1"
                        value={detected.some(g => g.name === nameA) ? nameA : ''}
                        onChange={e => pickDetected('A', e.target.value)}
                      >
                        <option value="">— choose —</option>
                        {detected.map(g => (
                          <option key={g.name} value={g.name}>{g.name} (n={g.samples.length})</option>
                        ))}
                      </select>
                    </label>
                    <label className="text-sm">Group B
                      <select
                        className="input mt-1 w-full py-1"
                        value={detected.some(g => g.name === nameB) ? nameB : ''}
                        onChange={e => pickDetected('B', e.target.value)}
                      >
                        <option value="">— choose —</option>
                        {detected.map(g => (
                          <option key={g.name} value={g.name}>{g.name} (n={g.samples.length})</option>
                        ))}
                      </select>
                    </label>
                  </div>
                </div>
              )}
              <div className="mb-4 grid gap-3 sm:grid-cols-2">
                <label className="text-sm">Group A name
                  <input className="input mt-1 w-full" value={nameA} onChange={e => setNameA(e.target.value)} /></label>
                <label className="text-sm">Group B name
                  <input className="input mt-1 w-full" value={nameB} onChange={e => setNameB(e.target.value)} /></label>
              </div>
              <label className="mb-3 flex flex-wrap items-center gap-2 text-sm text-slate-600 dark:text-slate-300">
                Control / reference:
                <select className="input py-1" value={control} onChange={e => setControl(e.target.value as 'A' | 'B')}>
                  <option value="A">{nameA || 'Group A'}</option>
                  <option value="B">{nameB || 'Group B'}</option>
                </select>
                <span className="text-xs text-slate-400">results read as “other vs control”</span>
              </label>
              <div className="max-h-72 overflow-auto rounded-lg border border-slate-100 dark:border-slate-800">
                <table className="w-full text-sm">
                  <tbody>
                    {counts.samples.map(s => (
                      <tr key={s} className="border-t border-slate-100 first:border-0 dark:border-slate-800">
                        <td className="px-3 py-1.5 font-mono">{s}</td>
                        <td className="px-3 py-1.5 text-right">
                          <div className="inline-flex overflow-hidden rounded-lg border border-slate-300 dark:border-slate-600">
                            {([['A', nameA || 'A'], ['B', nameB || 'B'], [null, 'exclude']] as const).map(([g, lbl]) => (
                              <button key={String(g)} onClick={() => setAssign(p => ({ ...p, [s]: g }))}
                                className={`px-3 py-1 text-xs ${assign[s] === g
                                  ? (g === null ? 'bg-slate-400 text-white' : 'bg-indigo-500 text-white')
                                  : 'bg-white text-slate-600 dark:bg-slate-800 dark:text-slate-300'}`}>
                                {lbl}
                              </button>
                            ))}
                          </div>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <p className="mt-2 text-xs text-slate-400">
                {nameA || 'A'}: {nA} · {nameB || 'B'}: {nB}
                {nExcluded > 0 && ` · ${nExcluded} excluded`}
                {designOk ? '' : ' · each group needs ≥ 2 samples and distinct names'}
              </p>
            </div>
            <div className="flex justify-between">
              <button className="btn" onClick={() => setStep('upload')}>← Back</button>
              <button className="btn btn-primary" disabled={!designOk} onClick={() => setStep('run')}>Next →</button>
            </div>
          </div>
        )}

        {(step === 'run' || (step === 'result' && running)) && (
          <div className="space-y-4">
            <div className="card p-6">
              <h2 className="mb-1 text-base font-semibold">3 · Run</h2>
              <div className="grid gap-3 sm:grid-cols-3">
                <label className="text-sm">Method
                  <select className="input mt-1 w-full" value={method} onChange={e => setMethod(e.target.value as Method)}>
                    <option value="limma">limma-voom (fast)</option>
                    <option value="deseq2">DESeq2 (gold standard)</option>
                  </select></label>
                <label className="text-sm">Project name
                  <input className="input mt-1 w-full" value={project} onChange={e => setProject(e.target.value)} /></label>
                <label className="text-sm">Species
                  <input className="input mt-1 w-full" value={species} onChange={e => setSpecies(e.target.value)} /></label>
              </div>
              <button className="btn btn-primary mt-4" disabled={running} onClick={doRun}>
                {running ? 'Running…' : `Run ${method === 'limma' ? 'limma-voom' : 'DESeq2'}`}
              </button>
              {!running && <button className="btn ml-2" onClick={() => setStep('design')}>← Back</button>}
              {runErr && <p className="mt-3 text-sm text-red-500">Failed: {runErr}</p>}
            </div>
            {log.length > 0 && (
              <pre className="card max-h-72 overflow-auto p-4 text-xs text-slate-600 dark:text-slate-300">{log.join('\n')}</pre>
            )}
          </div>
        )}

        {step === 'result' && result && !running && (
          <div className="space-y-4">
            <div className="card p-6 text-center">
              <div className="text-sm uppercase tracking-wide text-slate-400">Analysis complete</div>
              <div className="mt-1 text-3xl font-bold text-indigo-600">{result.nDeg.toLocaleString()}</div>
              <div className="text-sm text-slate-500">DEGs at padj &lt; 0.05 · {result.numerator} vs {result.denominator} · {method === 'limma' ? 'limma-voom' : 'DESeq2'}</div>
              <div className="mt-5 flex flex-wrap justify-center gap-2">
                {zipUrl && <a className="btn btn-primary" href={zipUrl} download={`${sanitize(project)}_bundle.zip`}>⭳ Download bundle (.zip)</a>}
                <a className="btn" href={EXPLORER_URL} target="_blank" rel="noopener noreferrer">Open RNA-seq Studio ↗</a>
              </div>
              <p className="mt-3 text-xs text-slate-400">Explore it: open RNA-seq Studio and drop the downloaded .zip onto the page.</p>
            </div>
            <div className="flex justify-center">
              <button className="btn" onClick={reset}>Analyze another dataset</button>
            </div>
          </div>
        )}
      </main>

      <footer className="border-t border-slate-200 py-3 text-center text-xs text-slate-400 dark:border-slate-700">
        Runs locally in your browser via webR (R 4.6.0) · your data never leaves this device
      </footer>
    </div>
  )
}

function Steps({ step }: { step: Step }) {
  const items: { id: Step; label: string }[] = [
    { id: 'upload', label: 'Upload' }, { id: 'design', label: 'Design' },
    { id: 'run', label: 'Run' }, { id: 'result', label: 'Results' },
  ]
  const idx = items.findIndex(i => i.id === step)
  return (
    <div className="flex items-center gap-2 text-xs">
      {items.map((it, i) => (
        <div key={it.id} className="flex items-center gap-2">
          <span className={`step-dot ${i <= idx ? 'bg-indigo-500 text-white' : 'bg-slate-200 text-slate-500 dark:bg-slate-700'}`}>{i + 1}</span>
          <span className={i <= idx ? 'font-medium text-slate-700 dark:text-slate-200' : 'text-slate-400'}>{it.label}</span>
          {i < items.length - 1 && <span className="mx-1 h-px w-6 bg-slate-200 dark:bg-slate-700" />}
        </div>
      ))}
    </div>
  )
}

const sanitize = (s: string) => s.trim().replace(/[^A-Za-z0-9._-]+/g, '_').replace(/^_+|_+$/g, '') || 'rnaseq'
