import { useMemo, useRef, useState } from 'react'
import Papa from 'papaparse'
import {
  runAnalysis, getWebR, ensureObjectPackages, installDtu,
  type AnalysisResult, type Method, type ContrastRequest,
} from './lib/webr'
import {
  buildBundleFiles, referenceGroup, referenceGroupFor, zipBundle, type IsoformOutput,
} from './lib/bundle'
import { parseMatrix, probeFromCsv, type Probe } from './lib/matrix'
import {
  isTranscriptMatrix, readLongRead, applySqanti, describe as describeLongRead,
  type LongReadInput,
} from './lib/longread'
import { runDte, dtuFromPipeline, FILTER_NOTE, type DexseqTables } from './lib/dtu'
import { readRObject, isRObjectFile } from './lib/robj'
import {
  detectFactors, pairwiseContrasts, withinFactorContrasts, interactionContrast,
  blockedContrasts,
  type Design, type ContrastSpec, type ContrastScheme,
} from './lib/design'
import { suggestBlockFactor, separationByFactor, SEPARATION_THRESHOLD } from './lib/blocking'

const EXPLORER_URL = 'https://jiaenlin.github.io/rnaseq-studio/'

function Flow({ at }: { at: 'convert' | 'done' }) {
  const steps: [string, string][] = [
    ['counts / DESeq2 object', 'genes × samples'],
    ['RNA-seq Lab', at === 'done' ? 'converted it' : 'converts it — you are here'],
    ['bundle .zip', at === 'done' ? 'ready to download' : 'the studio’s input format'],
    ['RNA-seq Studio', at === 'done' ? 'open it there next' : 'where you explore it'],
  ]
  const lit = at === 'done' ? 3 : 1
  return (
    <ol className="flex flex-wrap items-stretch gap-1.5">
      {steps.map(([name, what], i) => (
        <li key={name} className="flex items-stretch gap-1.5">
          <div className={`rounded-lg px-2.5 py-1.5 text-left ${i === lit
            ? 'bg-indigo-50 ring-1 ring-indigo-400 dark:bg-indigo-500/15'
            : 'bg-slate-100 dark:bg-slate-800/60'}`}>
            <div className={`text-[11.5px] font-semibold leading-tight ${i === lit
              ? 'text-indigo-700 dark:text-indigo-300'
              : 'text-slate-700 dark:text-slate-200'}`}>{name}</div>
            <div className="text-[10.5px] leading-tight text-slate-400">{what}</div>
          </div>
          {i < steps.length - 1 && <span className="self-center text-[11px] text-slate-300">&rarr;</span>}
        </li>
      ))}
    </ol>
  )
}

type Step = 'upload' | 'design' | 'run' | 'result'

interface Counts {
  countsCsv: string
  /** thinned numeric counts, for lib/blocking.ts */
  probe: Probe
  samples: string[]
  nGenes: number
  geneNames: Map<string, string> | null
  /** where it came from, for the note under the sample count */
  origin: string
  /** sample table that arrived with an R object, if any */
  colData: Record<string, Record<string, string>>
  colDataColumns: string[]
}

const EXCLUDED = '—'

export default function App() {
  const [step, setStep] = useState<Step>('upload')
  const [counts, setCounts] = useState<Counts | null>(null)
  const [uploadErr, setUploadErr] = useState<string | null>(null)
  const [uploadBusy, setUploadBusy] = useState(false)
  const fileRef = useRef<HTMLInputElement>(null)
  /** Transcript counts exactly as uploaded. The usage test models these, not the gene sums. */
  const txCountsRef = useRef<string | null>(null)

  // design
  const [design, setDesign] = useState<Design | null>(null)
  const [factorNames, setFactorNames] = useState<string[]>([])
  const [refs, setRefs] = useState<string[]>([])
  const [groupOf, setGroupOf] = useState<Record<string, string>>({})
  /**
   * Which comparisons to export — `null` until somebody says.
   *
   * `null` and `new Set()` are DIFFERENT states, and the difference is the
   * whole reason "Clear" can exist. This used to be one `Set`, with empty
   * meaning "the reader has not chosen, so offer all of them", which made
   * "the reader wants none of them" unrepresentable — `toggle` had to refuse
   * the last deselection (`return next.size ? next : new Set([id])`) or
   * unticking the final box would silently retick every box on the page.
   *
   * So: `null` means undecided and the default below fills in; a Set means
   * exactly these, empty included. Anything that INVALIDATES the choice — a new
   * reference level, a regrouped sample, a new file — sets it back to `null`
   * rather than to empty, because those change which contrasts exist at all.
   */
  const [chosen, setChosen] = useState<Set<string> | null>(null)

  /**
   * The factor fitted separately, if any. -1 is "one fit over everything" —
   * the default, and what every design smaller than an atlas wants.
   *
   * Held as an INDEX, because that is what survives the reader renaming a
   * factor. Keeping the name here instead meant typing "tissue" over "factor1"
   * silently switched blocking off: the stored name matched nothing, `blocking`
   * fell back to '', and the page went from 110 within-tissue comparisons to 95
   * that include Kidney-vs-Liver, with no message and nothing to click to
   * explain it. A factor's identity is its position; its name is a label on it.
   */
  const [blockFactorIdx, setBlockFactorIdx] = useState(-1)
  const [scheme, setScheme] = useState<ContrastScheme>('all-pairs')

  // run params
  // DESeq2 by default. It is what the studio runs for any comparison the
  // reader asks for later, so defaulting to it means the bundle's own tables
  // and anything computed on top of them come from the same engine.
  const [method, setMethod] = useState<Method>('deseq2')
  const [shrink, setShrink] = useState<'none' | 'apeglm'>('none')
  const [project, setProject] = useState('My RNA-seq analysis')
  const [species, setSpecies] = useState('human')
  const [log, setLog] = useState<string[]>([])
  const [running, setRunning] = useState(false)
  const [result, setResult] = useState<AnalysisResult | null>(null)
  const [zipUrl, setZipUrl] = useState<string | null>(null)
  const [runErr, setRunErr] = useState<string | null>(null)
  const [saved, setSaved] = useState(false)
  /** The isoform layer, when the uploaded matrix was transcript-level. */
  const [longRead, setLongRead] = useState<LongReadInput | null>(null)
  const [wantIsoform, setWantIsoform] = useState(true)
  const [sqantiNote, setSqantiNote] = useState<string | null>(null)
  const sqantiRef = useRef<HTMLInputElement>(null)
  const dexseqRef = useRef<HTMLInputElement>(null)
  /** The pipeline's own DEXSeq tables, when supplied. */
  const [dexseq, setDexseq] = useState<DexseqTables | null>(null)
  const [dexseqNote, setDexseqNote] = useState<string | null>(null)

  const onLog = (m: string) => setLog(prev => [...prev, m])

  /** Seed the design state from a set of sample names (+ any colData columns). */
  const seedDesign = (
    samples: string[],
    colData: Record<string, Record<string, string>>,
    colDataColumns: string[],
    probe: Probe,
  ) => {
    // A sample table that came with the object beats anything guessed from
    // names: it is what the pipeline actually modelled.
    const fromColData = colDataColumns.length > 0
    let d: Design
    if (fromColData) {
      const factors = colDataColumns.map((c, i) => {
        const values = samples.map(s => colData[s]?.[c] ?? '')
        return { name: c || `factor${i + 1}`, levels: [...new Set(values)], values }
      }).filter(f => f.levels.length > 1)
      const groups = samples.map((_, si) => factors.map(f => f.values[si]).join('_'))
      const groupLevels = [...new Set(groups)]
      const cells = new Set(groups)
      const expected = factors.reduce((a, f) => a * f.levels.length, 1)
      d = {
        factors, groups, groupLevels,
        factorial: factors.length > 1 && cells.size === expected,
        balanced: false,
      }
    } else {
      d = detectFactors(samples)
    }
    setDesign(d)
    setFactorNames(d.factors.map(f => f.name))
    setRefs(d.factors.map(f => f.levels[0]))
    // Suggested from the DATA, never from the shape of the names — see
    // lib/blocking.ts. With no probe this is null, i.e. one fit for everything.
    const suggested = suggestBlockFactor(d, probe)
    setBlockFactorIdx(suggested ? suggested.index : -1)
    const g: Record<string, string> = {}
    samples.forEach((s, i) => { g[s] = d.groups[i] })
    setGroupOf(g)
    setChosen(null)               // filled by the effect-free default below
    return d
  }

  const onFile = async (f: File | undefined) => {
    if (!f) return
    setUploadErr(null); setUploadBusy(true); setLog([])
    try {
      // EVERY upload path clears the isoform layer first. The R-object branch
      // used to leave a previous transcript upload's layer in place: uploading
      // transcript_counts.tsv and then an .rds shipped a bundle whose
      // transcripts.csv and dtu_*.csv came from the first file and whose
      // raw_counts.csv came from the second, with platform: 'long-read'
      // asserting they belonged together. Nothing errored.
      setLongRead(null); setSqantiNote(null); txCountsRef.current = null
    setDexseq(null); setDexseqNote(null)
      if (isRObjectFile(f.name)) {
        // nf-core's own object: counts, sample table and design in one file.
        const webR = await getWebR(onLog)
        await ensureObjectPackages(webR, /\.rdata|\.rda$/i.test(f.name), onLog)
        const obj = await readRObject(f, webR, onLog)
        setCounts({
          countsCsv: obj.countsCsv, samples: obj.samples, nGenes: obj.nGenes,
          geneNames: null, origin: obj.source,
          colData: obj.colData, colDataColumns: obj.colDataColumns,
          probe: probeFromCsv(obj.countsCsv),
        })
        seedDesign(obj.samples, obj.colData, obj.colDataColumns, probeFromCsv(obj.countsCsv))
      } else {
        const text = await f.text()
        const parsed = Papa.parse<string[]>(text.trim(), { skipEmptyLines: true })
        const rows = parsed.data as string[][]
        const m = parseMatrix(rows)
        // A transcript matrix is quantified per ISOFORM. The gene layer is
        // summed from it here so both levels come from one file — see
        // lib/longread.ts. Everything downstream then behaves as before.
        const lr = isTranscriptMatrix(rows[0].map(h => String(h ?? '').trim()))
          ? readLongRead(m, rows) : null
        setLongRead(lr)
        if (lr) onLog(`Isoform layer: ${describeLongRead(lr)}`)
        setCounts({
          countsCsv: lr ? lr.geneCountsCsv : m.countsCsv,
          samples: m.samples, nGenes: lr ? lr.nGenes : m.nGenes,
          geneNames: lr ? lr.geneNames
            : m.geneNames ? new Map(m.geneIds.map((id, i) => [id, m.geneNames![i]])) : null,
          origin: lr
            ? `transcript matrix (${lr.transcripts.length.toLocaleString()} isoforms, summed to ${lr.nGenes.toLocaleString()} genes)`
            : m.annotationColumns.length > 1
              ? `matrix (${m.annotationColumns.join(' + ')} read as annotation)`
              : 'matrix',
          colData: {}, colDataColumns: [],
          probe: lr ? probeFromCsv(lr.geneCountsCsv) : m.probe,
        })
        seedDesign(m.samples, {}, [], lr ? probeFromCsv(lr.geneCountsCsv) : m.probe)
        // The transcript counts are kept as read; the DTU engine needs them raw.
        // The ROUNDED, transcript-keyed matrix — not parsed.countsCsv, which is
        // headed `gene_id` and still fractional.
        txCountsRef.current = lr ? lr.txCountsCsv : null
      }
      setStep('design')
    } catch (e: any) {
      setUploadErr(String(e?.message || e))
    } finally {
      setUploadBusy(false)
    }
  }

  /**
   * The optional SQANTI3 classification, joined onto transcripts already read.
   *
   * Separate from the counts upload because it is a separate file from the
   * pipeline: out/cohort/sqanti/cohort_classification.txt. Without it every
   * isoform is "annotated or not"; with it they carry FSM / ISM / NIC / NNC,
   * which is the vocabulary the studio colours by.
   */
  async function onSqanti(f?: File) {
    if (!f || !longRead) return
    try {
      const text = await f.text()
      const parsed = Papa.parse<string[]>(text.trim(), { skipEmptyLines: true, delimiter: '\t' })
      const { matched, total, vocabulary } = applySqanti(
        longRead.transcripts, parsed.data as string[][])
      setLongRead({ ...longRead, vocabulary })
      setSqantiNote(matched === 0
        ? `No transcript in that file matches this matrix — it is probably from a different run.`
        : `${matched.toLocaleString()} of ${total.toLocaleString()} transcripts classified.`)
    } catch (e: any) {
      setSqantiNote(String(e?.message || e))
    }
  }

  /**
   * The pipeline's own DEXSeq result.
   *
   * DTU is not computed in this app — DEXSeq cannot load in webR (see
   * src/lib/dtu.ts) — so the usage numbers come from the run that already
   * produced them, unaltered. Both files are accepted at once: the transcript
   * table carries the test, the gene table carries DEXSeq's perGeneQValue.
   */
  async function onDexseq(files?: FileList | null) {
    if (!files?.length) return
    try {
      const next: DexseqTables = { transcript: [] }
      for (const f of Array.from(files)) {
        const rows = Papa.parse<string[]>((await f.text()).trim(),
          { skipEmptyLines: true, delimiter: '\t' }).data as string[][]
        const head = (rows[0] ?? []).map(h => String(h ?? '').trim().toLowerCase())
        if (head.includes('featureid')) next.transcript = rows
        else if (head.includes('geneid') && head.some(h => /qval/.test(h))) next.gene = rows
      }
      if (!next.transcript.length) {
        setDexseqNote('No results_dtu_transcript.tsv among those files — that is the one '
          + 'carrying the test. Look in out/de_analysis/<contrast>/.')
        return
      }
      setDexseq(next)
      setDexseqNote(`DEXSeq table read: ${(next.transcript.length - 1).toLocaleString()} transcripts`
        + (next.gene ? `, and ${(next.gene.length - 1).toLocaleString()} per-gene q-values.`
          : '. No results_dtu_gene.tsv, so the bundle will carry no per-gene q-value.'))
    } catch (e: any) {
      setDexseqNote(String(e?.message || e))
    }
  }

  /* ---------- derived design ---------- */

  const named: Design | null = useMemo(() => {
    if (!design) return null
    return { ...design, factors: design.factors.map((f, i) => ({ ...f, name: factorNames[i] || f.name })) }
  }, [design, factorNames])

  const activeSamples = useMemo(
    () => (counts?.samples ?? []).filter(s => groupOf[s] && groupOf[s] !== EXCLUDED),
    [counts, groupOf])

  const groupLevels = useMemo(() => {
    const seen = new Set<string>()
    const out: string[] = []
    for (const s of activeSamples) { const g = groupOf[s]; if (!seen.has(g)) { seen.add(g); out.push(g) } }
    return out
  }, [activeSamples, groupOf])

  const groupSizes = useMemo(() => {
    const m = new Map<string, number>()
    for (const s of activeSamples) m.set(groupOf[s], (m.get(groupOf[s]) ?? 0) + 1)
    return m
  }, [activeSamples, groupOf])

  /**
   * The reference as a GROUP label, not a factor level. Passing `refs[0]`
   * straight to `pairwiseContrasts` names a denominator no sample carries
   * whenever a name position is constant — see `referenceGroupFor`.
   */
  const refGroup = useMemo(
    () => (named ? referenceGroupFor(named, refs, groupLevels) : ''),
    [named, refs, groupLevels])

  /** Blocking is only meaningful once there are two factors to separate. */
  const blockable = useMemo(
    () => (named && named.factors.length > 1 ? named.factors : []),
    [named])
  const blockIdx = blockFactorIdx >= 0 && blockFactorIdx < blockable.length ? blockFactorIdx : -1

  /**
   * How far apart each factor's levels sit — the number behind the suggestion.
   *
   * Shown rather than kept private, because "fit these separately" is a claim
   * about the data and the reader is entitled to the evidence. It is also the
   * only thing on the page that distinguishes eleven tissues from six
   * timepoints, which look identical from the sample names.
   */
  const seps = useMemo(
    () => (named && counts?.probe ? separationByFactor(named, counts.probe) : []),
    [named, counts])
  const sepOf = (i: number) => seps.find(x => x.index === i)?.separation
  const blocking = blockIdx >= 0 ? blockable[blockIdx].name : ''

  /**
   * The plan when a blocking factor is chosen: one fit per level, and only
   * contrasts that live inside one. Kept separate from `available` so the
   * budget downgrade has somewhere to be reported from.
   */
  const plan = useMemo(
    () => (named && blocking
      // The reference level belongs to the factor being compared INSIDE a block,
      // which is the one that is not the block.
      ? blockedContrasts(named, blocking, { scheme, reference: refs[blockIdx === 0 ? 1 : 0] })
      : null),
    [named, blocking, blockIdx, scheme, refs])

  /** Every contrast worth offering, given the factors and their references. */
  const available: ContrastSpec[] = useMemo(() => {
    if (!named) return []
    if (plan) {
      // No interaction and no cross-block pair: neither is answerable by any
      // one of these fits, and offering a comparison nothing can compute is
      // how a design page starts lying.
      return plan.contrasts.filter(c =>
        [c.numerator, c.denominator].every(g => (groupSizes.get(g) ?? 0) >= 2))
    }
    const within = named.factors.length > 1 ? withinFactorContrasts(named, refs) : []
    const base = within.length
      ? within
      : pairwiseContrasts(groupLevels, refGroup || groupLevels[0])
    const ix = named.factors.length > 1 ? interactionContrast(named, refs) : null
    const all = ix ? [...base, ix] : base
    // Only offer contrasts whose groups actually survive the exclusions, and
    // that have enough samples on both sides to estimate anything.
    return all.filter(c => {
      const gs = c.kind === 'interaction'
        ? named.groupLevels
        : [c.numerator, c.denominator]
      return gs.every(g => (groupSizes.get(g) ?? 0) >= 2)
    })
  }, [named, refs, groupLevels, groupSizes, refGroup, plan])

  // Default selection: everything pairwise, plus the interaction if present.
  const effectiveChosen = useMemo(() => {
    if (chosen) return chosen
    // Everything by default, interaction included: the whole point of
    // detecting it is that someone would not have thought to ask for it.
    return new Set(available.map(c => c.id))
  }, [chosen, available])

  const selected = available.filter(c => effectiveChosen.has(c.id))
  const designOk = groupLevels.length >= 2 && selected.length >= 1 &&
    groupLevels.every(g => (groupSizes.get(g) ?? 0) >= 2)

  const toggle = (id: string) => setChosen(() => {
    const next = new Set(effectiveChosen)
    if (next.has(id)) next.delete(id); else next.add(id)
    return next
  })

  /**
   * All of them, or none.
   *
   * Untangling twelve checkboxes one click at a time to export the one
   * comparison you came for is the reason this is here — on a 2x3 design
   * `available` is nine contrasts plus an interaction.
   *
   * "All" goes back to `null` rather than writing every id out. It means the
   * same thing today and keeps meaning it when the set of available contrasts
   * changes underneath: pick a different reference level after selecting all
   * explicitly, and a written-out list would leave you with a page of unticked
   * boxes and no idea why.
   */
  const nSelected = selected.length
  const selectAll = () => setChosen(null)
  const clearAll = () => setChosen(new Set())

  /* ---------- run ---------- */

  const doRun = async () => {
    if (!counts || !named) return
    setRunning(true); setRunErr(null); setResult(null); setZipUrl(null)
    try {
      const covariates = named.factors.map(f => f.name)
      const bi = blockIdx
      const samples = activeSamples.map(s => {
        const si = counts.samples.indexOf(s)
        const rec: Record<string, string> = { sample: s, group: groupOf[s] }
        named.factors.forEach(f => { rec[f.name] = f.values[si] ?? '' })
        // The block travels as its own column, not as one of the covariates:
        // R partitions on it before fitting anything.
        rec.block = bi >= 0 ? (named.factors[bi].values[si] ?? '') : ''
        return rec as { sample: string; group: string; block: string }
      })
      const requests: ContrastRequest[] = selected.map(c => {
        if (c.kind !== 'interaction') {
          return { ...c, plus: [c.numerator], minus: [c.denominator], block: c.block ?? '' }
        }
        // (A1B1 - A0B1) - (A1B0 - A0B0) written over group means.
        const [fa, fb] = named.factors
        const aAlt = fa.levels.find(l => l !== refs[0])!
        const bAlt = fb.levels.find(l => l !== refs[1])!
        const g = (a: string, b: string) => {
          const i = fa.values.findIndex((v, k) => v === a && fb.values[k] === b)
          if (i < 0) throw new Error(
            `The interaction needs a ${a}/${b} sample and there is none — ` +
            `that cell of the design is empty.`)
          return named.groups[i]
        }
        return {
          ...c,
          plus: [g(aAlt, bAlt), g(refs[0], refs[1])],
          minus: [g(refs[0], bAlt), g(aAlt, refs[1])],
        }
      })

      const input = {
        countsCsv: counts.countsCsv, samples, groupLevels, contrasts: requests, method,
        shrink: method === 'deseq2' ? shrink : 'none',
      }
      const res = await runAnalysis(input, onLog)

      // The isoform layer, after the gene fit and never instead of it. One
      // DTE+DTU run per PAIRWISE contrast: DEXSeq's design is a two-level
      // usage test, so an interaction coefficient has no DTU counterpart and
      // is skipped rather than approximated.
      // Built from `wantIsoform`, not merely from `longRead`. Otherwise
      // unticking the box still shipped transcripts.csv and a meta.json naming
      // a DTU engine and filter that were never applied.
      // Built from `wantIsoform`, not merely from `longRead`. Otherwise
      // unticking the box still shipped transcripts.csv and a meta.json naming
      // a filter that was never applied.
      const runIsoform = !!(longRead && wantIsoform && txCountsRef.current)
      const isoform: IsoformOutput = {
        transcripts: runIsoform ? longRead!.transcripts : undefined,
        vocabulary: longRead?.vocabulary,
        byContrast: {},
        dtuByContrast: {},
      }
      if (runIsoform && longRead && txCountsRef.current) {
        const pairs = requests.filter(c => c.kind === 'pairwise')
        const webR = await getWebR(onLog)
        await installDtu(webR, onLog)
        for (const c of pairs) {
          if (c.plus.length !== 1 || c.minus.length !== 1) {
            onLog(`Isoform layer: skipping "${c.label}" — it is not a two-group comparison.`)
            continue
          }
          onLog(`Isoform layer: ${c.label}`)
          const pairSamples = samples
            .filter(x => x.group === c.plus[0] || x.group === c.minus[0])
            .map(x => ({ sample: x.sample, group: x.group }))
          const r = await runDte(webR, {
            txCountsCsv: txCountsRef.current,
            samples: pairSamples,
            numerator: c.plus[0], denominator: c.minus[0],
            shrink: method === 'deseq2' ? shrink : 'none',
          }, onLog)
          isoform.byContrast[c.id] = r
          onLog(`  ${r.nTested.toLocaleString()} transcripts tested · `
            + `${r.nDte} differentially expressed (FDR<0.05)`)

          // DTU is DEXSeq's, computed by the pipeline. Only the FIRST pairwise
          // comparison gets it: one uploaded table is one contrast, and
          // attaching it to a second pair would label another comparison's
          // numbers with this one's name.
          if (dexseq && c.id === pairs[0]?.id) {
            const d = dtuFromPipeline(dexseq, txCountsRef.current, pairSamples,
              c.plus[0], c.minus[0])
            isoform.dtuByContrast![c.id] = d.dtuCsv
            if (d.flipped) (isoform.dtuFlipped ??= []).push(c.id)
            d.notes.forEach(onLog)
          }
        }
      }

      const files = buildBundleFiles(input, res, {
        project, species, method, covariates,
        shrink: method === 'deseq2' ? shrink : 'none',
        blockFactor: blocking || undefined,
        // The reference the reader actually picked. Not derivable from the
        // contrasts — see referenceGroup.
        control: referenceGroup({ factors: named.factors, groupLevels }, refs),
        geneNames: counts.geneNames ?? undefined,
        isoform: isoform.transcripts ? isoform : undefined,
        transcriptCountsCsv: txCountsRef.current ?? undefined,
      })
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

  const reset = () => {
    setStep('upload'); setCounts(null); setDesign(null); setResult(null)
    setZipUrl(null); setLog([]); setRunErr(null); setSaved(false); setChosen(null)
    setLongRead(null); setSqantiNote(null); txCountsRef.current = null
    setDexseq(null); setDexseqNote(null)
  }

  return (
    <div className="mx-auto flex min-h-full max-w-4xl flex-col px-4">
      <header className="flex items-center gap-2 py-4">
        <span className="grid h-9 w-9 place-items-center rounded-lg bg-indigo-500 font-bold text-white">L</span>
        <div>
          <h1 className="text-lg font-semibold leading-none">RNA-seq Lab</h1>
          <p className="text-xs text-slate-400">Turns a counts matrix or a DESeq2 object into RNA-seq Studio&rsquo;s input file · nothing is uploaded</p>
        </div>
        {step !== 'upload' && <button className="btn ml-auto" onClick={reset}>Start over</button>}
      </header>

      <Steps step={step} />

      <main className="step-enter flex-1 py-4" key={step}>
        {step === 'upload' && (
          <div className="card p-6">
            <p className="mb-4 text-sm">
              This page does one thing: it converts your counts into the{' '}
              <code className="rounded bg-slate-100 px-1 py-0.5 text-[12px] dark:bg-slate-800">bundle.zip</code>{' '}
              that <a className="underline" href={EXPLORER_URL} target="_blank" rel="noopener noreferrer">
                RNA-seq Studio</a>{' '}opens. You read and plot your results there, not here.
            </p>
            <div className="mb-5"><Flow at="convert" /></div>
            <h2 className="mb-1 text-base font-semibold">1 · Upload counts</h2>
            <p className="mb-3 text-sm text-slate-500">
              Either a <b>CSV/TSV matrix</b> (genes as rows, samples as columns) or an{' '}
              <b>R object from nf-core/rnaseq</b>.
            </p>
            <ul className="mb-4 space-y-1 text-[13px] text-slate-500">
              <li>· <code className="font-mono text-[12px]">salmon.merged.gene_counts.tsv</code> — the
                nf-core matrix. Its <code className="font-mono text-[12px]">gene_id</code> and{' '}
                <code className="font-mono text-[12px]">gene_name</code> columns are recognised as
                annotation, not as two extra samples.</li>
              <li>· <code className="font-mono text-[12px]">deseq2.dds.RData</code> — nf-core&rsquo;s
                DESeq2 object. Carries the counts <em>and</em> the sample table, so the design does not
                have to be guessed from sample names.</li>
              <li>· <code className="font-mono text-[12px]">*.SummarizedExperiment.rds</code></li>
              <li>· <code className="font-mono text-[12px]">cohort/transcript_counts.tsv</code> — Oxford
                Nanopore <b>wf-transcriptomes</b>. Quantified per isoform, so the lab sums it to genes
                <em>and</em> keeps the transcript layer: it adds isoform-level DESeq2 to the bundle. A PacBio or StringTie transcript matrix works the same way.</li>
            </ul>
            <button className="btn btn-primary" disabled={uploadBusy} onClick={() => fileRef.current?.click()}>
              {uploadBusy ? 'Reading…' : '⭱ Choose counts file or R object'}
            </button>
            <input ref={fileRef} type="file" accept=".csv,.tsv,.txt,.rds,.RData,.rda" className="hidden"
              onChange={e => onFile(e.target.files?.[0])} />
            {uploadErr && <p className="mt-3 text-sm text-red-500">{uploadErr}</p>}
            {log.length > 0 && (
              <pre className="mt-3 max-h-40 overflow-auto rounded-lg bg-slate-50 p-3 text-xs text-slate-600 dark:bg-slate-800/60 dark:text-slate-300">{log.join('\n')}</pre>
            )}
            <p className="mt-4 text-xs text-slate-400">
              Reading an <code className="font-mono">.RData</code> loads R in the browser first (a one-time
              download, then cached). Best for small–moderate datasets.
            </p>
          </div>
        )}

        {step === 'design' && counts && named && (
          <div className="space-y-4">
            <div className="card p-6">
              <h2 className="mb-1 text-base font-semibold">2 · Design</h2>
              <p className="mb-4 text-sm text-slate-500">
                {counts.samples.length} samples · {counts.nGenes.toLocaleString()} genes · read from{' '}
                {counts.origin}.
                {counts.colDataColumns.length > 0 && ' The sample table came with the object.'}
              </p>

              {named.factors.length > 1 ? (
                <div className="mb-4 rounded-lg border border-indigo-200 bg-indigo-50/60 p-3 dark:border-indigo-500/30 dark:bg-indigo-500/10">
                  <p className="mb-2 text-xs text-slate-600 dark:text-slate-300">
                    <b>{named.factors.map(f => f.levels.length).join(' × ')} factorial design</b> detected
                    {named.factorial ? '' : ' (some combinations are missing)'}. Name each factor and pick
                    its reference level — every comparison is then read as “other vs reference”.
                  </p>
                  <div className="grid gap-3 sm:grid-cols-2">
                    {named.factors.map((f, i) => (
                      <div key={i} className="rounded-lg bg-white/70 p-2.5 dark:bg-slate-800/60">
                        <input
                          className="input mb-1.5 w-full py-1 text-sm font-medium"
                          value={factorNames[i] ?? ''}
                          onChange={e => setFactorNames(p => p.map((v, k) => (k === i ? e.target.value : v)))}
                        />
                        <label className="flex items-center gap-2 text-xs text-slate-500">
                          reference
                          <select className="input flex-1 py-0.5 text-xs" value={refs[i] ?? ''}
                            onChange={e => { setRefs(p => p.map((v, k) => (k === i ? e.target.value : v))); setChosen(null) }}>
                            {f.levels.map(l => <option key={l} value={l}>{l}</option>)}
                          </select>
                        </label>
                        <p className="mt-1 text-[11px] text-slate-400">{f.levels.join(' · ')}</p>
                      </div>
                    ))}
                  </div>

                  {/* SEPARATE FITS. Below the factors because it is a
                      consequence of them: you decide what the factors are, then
                      whether one of them is too big a jump to model across. */}
                  <div className="mt-3 rounded-lg bg-white/70 p-2.5 dark:bg-slate-800/60">
                    <label className="flex flex-wrap items-center gap-2 text-xs text-slate-600 dark:text-slate-300">
                      <b>Fit separately by</b>
                      <select className="input py-0.5 text-xs" value={String(blockIdx)}
                        onChange={e => { setBlockFactorIdx(Number(e.target.value)); setChosen(null) }}>
                        <option value="-1">nothing — one fit over every group</option>
                        {blockable.map((f, i) => (
                          <option key={i} value={String(i)}>
                            {f.name} ({f.levels.length} fits){
                              sepOf(i) != null ? ` · separation ${sepOf(i)!.toFixed(1)}` : ''}
                          </option>
                        ))}
                      </select>
                      {plan && (
                        <select className="input py-0.5 text-xs" value={scheme}
                          onChange={e => { setScheme(e.target.value as ContrastScheme); setChosen(null) }}>
                          <option value="all-pairs">every pair, later vs earlier</option>
                          <option value="vs-reference">each level vs the reference</option>
                          <option value="consecutive">consecutive levels only</option>
                        </select>
                      )}
                    </label>

                    {plan ? (
                      <>
                        <p className="mt-1.5 text-[11px] leading-relaxed text-slate-500 dark:text-slate-400">
                          <b className="tabular-nums">{plan.contrasts.length} comparisons</b>
                          {plan.perBlock > 0 && <> — {plan.perBlock} in each of {plan.blocks.length} {blocking} levels</>}.
                          Each level is fitted on its own, so its dispersion comes from samples
                          like it. Nothing compares one {blocking} to another: no single fit
                          answers that, and the fold changes are comparable between them anyway.
                          {sepOf(blockIdx) != null && (
                            <> Its levels explain <b>{sepOf(blockIdx)!.toFixed(0)}×</b> more variance
                              between them than within, against a threshold of {SEPARATION_THRESHOLD};
                              that is why this is offered and the other factor is not.</>
                          )}
                        </p>
                        {plan.scheme !== plan.requested && (
                          <p className="mt-1 text-[11px] text-amber-600 dark:text-amber-400">
                            Every pair would be too many tables to put in one bundle, so this is
                            each level against the reference instead. Pick fewer levels, or export
                            in more than one run, if you need all of them.
                          </p>
                        )}
                      </>
                    ) : (
                      <p className="mt-1.5 text-[11px] leading-relaxed text-slate-500 dark:text-slate-400">
                        One fit over every group, which is right when the groups are variants of
                        one experiment. Separate fits are for a factor whose levels are not
                        comparable — different tissues, different cell lines — because DESeq2
                        estimates one dispersion per gene across whatever is in the fit, so the
                        noisiest level sets the variance for the quietest.
                        {seps.length > 0 && (
                          <> Nothing here separates enough to need it: {seps.map(x =>
                            `${x.name} ${x.separation.toFixed(1)}`).join(', ')} — all under{' '}
                            {SEPARATION_THRESHOLD}.</>
                        )}
                      </p>
                    )}
                  </div>
                </div>
              ) : (
                <label className="mb-4 flex flex-wrap items-center gap-2 text-sm text-slate-600 dark:text-slate-300">
                  Control / reference group:
                  <select className="input py-1" value={refs[0] ?? groupLevels[0] ?? ''}
                    onChange={e => { setRefs([e.target.value]); setChosen(null) }}>
                    {groupLevels.map(g => <option key={g} value={g}>{g} (n={groupSizes.get(g) ?? 0})</option>)}
                  </select>
                  <span className="text-xs text-slate-400">results read as “other vs reference”</span>
                </label>
              )}

              <div className="mb-1.5 flex flex-wrap items-center gap-x-3 gap-y-1">
                <h3 className="text-sm font-semibold">Comparisons to export</h3>
                {/* A 2x3 design offers nine contrasts plus an interaction, and
                    getting to the one you came for meant nine clicks. */}
                {available.length > 1 && (
                  <span className="ml-auto flex items-center gap-2 text-xs">
                    <span className="tabular-nums text-slate-400">
                      {nSelected} of {available.length} selected
                    </span>
                    <button type="button" className="btn px-2 py-0.5 text-xs"
                      disabled={nSelected === available.length} onClick={selectAll}>Select all</button>
                    <button type="button" className="btn px-2 py-0.5 text-xs"
                      disabled={nSelected === 0} onClick={clearAll}>Clear</button>
                  </span>
                )}
              </div>
              <p className="mb-2 text-xs text-slate-500">
                Each one becomes a <code className="font-mono">deg_*.csv</code> in the bundle; the studio
                lets you switch between them.
              </p>
              <div className="mb-4 space-y-1.5">
                {available.length === 0 && (
                  <p className="text-xs text-amber-600">No comparison has ≥2 samples on both sides yet.</p>
                )}
                {/* Reachable for the first time now that Clear exists, and the
                    reason the Next button is dead — said where the boxes are,
                    not only in the small print under the sample table. */}
                {available.length > 0 && nSelected === 0 && (
                  <p className="text-xs text-amber-600">
                    Nothing selected, so the bundle would carry no results. Tick one, or Select all.
                  </p>
                )}
                {available.map(c => (
                  <label key={c.id}
                    className={`flex cursor-pointer items-start gap-2.5 rounded-lg border p-2.5 text-sm ${
                      effectiveChosen.has(c.id)
                        ? 'border-indigo-300 bg-indigo-50/60 dark:border-indigo-500/40 dark:bg-indigo-500/10'
                        : 'border-slate-200 dark:border-slate-700'}`}>
                    <input type="checkbox" className="mt-0.5" checked={effectiveChosen.has(c.id)}
                      onChange={() => toggle(c.id)} />
                    <span className="min-w-0 flex-1">
                      <span className="font-medium">{c.label}</span>
                      {c.kind === 'interaction' && (
                        <span className="ml-2 rounded bg-amber-100 px-1.5 py-0.5 text-[10px] font-semibold uppercase tracking-wide text-amber-700 dark:bg-amber-500/20 dark:text-amber-300">
                          interaction
                        </span>
                      )}
                      <span className="block font-mono text-[11px] text-slate-400">{c.id}</span>
                    </span>
                  </label>
                ))}
              </div>
              {available.some(c => c.kind === 'interaction') && (
                <p className="mb-4 rounded-lg bg-amber-50 p-2.5 text-[12px] text-amber-800 dark:bg-amber-500/10 dark:text-amber-200">
                  <b>The interaction</b> asks whether one factor&rsquo;s effect <em>depends on</em> the
                  other. No pairwise comparison answers it, and in a 2×2 study it is often the real
                  question. A gene significant here has a different response in one arm than the other —
                  it is not a fold change and should not be read as one.
                </p>
              )}

              <details className="mb-3">
                <summary className="cursor-pointer text-sm font-medium text-slate-600 dark:text-slate-300">
                  Sample assignment ({activeSamples.length} in · {counts.samples.length - activeSamples.length} excluded)
                </summary>
                <div className="mt-2 max-h-72 overflow-auto rounded-lg border border-slate-100 dark:border-slate-800">
                  <table className="w-full text-sm">
                    <tbody>
                      {counts.samples.map(s => (
                        <tr key={s} className="border-t border-slate-100 first:border-0 dark:border-slate-800">
                          <td className="px-3 py-1.5 font-mono text-[13px]">{s}</td>
                          <td className="px-3 py-1.5 text-right">
                            <select className="input py-0.5 text-xs" value={groupOf[s] ?? EXCLUDED}
                              onChange={e => { setGroupOf(p => ({ ...p, [s]: e.target.value })); setChosen(null) }}>
                              {[...new Set([...named.groupLevels, EXCLUDED])].map(g => (
                                <option key={g} value={g}>{g === EXCLUDED ? 'exclude' : g}</option>
                              ))}
                            </select>
                          </td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
              </details>

              <p className="text-xs text-slate-400">
                {groupLevels.map(g => `${g}: ${groupSizes.get(g) ?? 0}`).join(' · ')}
                {designOk ? ''
                  : available.length > 0 && nSelected === 0
                    ? ' · tick at least one comparison to export'
                    : ' · every group needs ≥ 2 samples, and at least one comparison'}
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
              <p className="mb-3 text-sm text-slate-500">
                {selected.length} comparison{selected.length === 1 ? '' : 's'} from one model fit.
              </p>
              {longRead && (
                <div className="mb-4 rounded-lg border border-emerald-300/60 bg-emerald-50/60 p-3 text-[13px] dark:border-emerald-700/50 dark:bg-emerald-950/30">
                  <label className="flex items-start gap-2 font-medium">
                    <input type="checkbox" className="mt-1" checked={wantIsoform}
                      onChange={e => setWantIsoform(e.target.checked)} />
                    <span>Also run the isoform layer</span>
                  </label>
                  <p className="mt-1 pl-6 text-slate-600 dark:text-slate-300">
                    {describeLongRead(longRead)}. Adds transcript-level <b>DESeq2</b> per
                    comparison, keeping {FILTER_NOTE}.
                  </p>
                  <div className="mt-2 pl-6">
                    <button className="btn btn-ghost text-xs" onClick={() => sqantiRef.current?.click()}>
                      ⭱ Add SQANTI3 categories (optional)
                    </button>
                    <input ref={sqantiRef} type="file" accept=".txt,.tsv,.csv" className="hidden"
                      onChange={e => onSqanti(e.target.files?.[0])} />
                    <p className="mt-1 text-xs text-slate-500">
                      <code className="font-mono text-[11px]">cohort/sqanti/cohort_classification.txt</code>{' '}
                      — labels each isoform FSM / ISM / NIC / NNC. Without it the studio can still show
                      novel versus annotated, but not which kind of novel.
                    </p>
                    {sqantiNote && <p className="mt-1 text-xs text-emerald-700 dark:text-emerald-400">{sqantiNote}</p>}
                  </div>
                  <div className="mt-3 pl-6">
                    <button className="btn btn-ghost text-xs" onClick={() => dexseqRef.current?.click()}>
                      ⭱ Add your pipeline’s DEXSeq result (optional)
                    </button>
                    <input ref={dexseqRef} type="file" accept=".tsv,.txt,.csv" multiple className="hidden"
                      onChange={e => onDexseq(e.target.files)} />
                    <p className="mt-1 text-xs text-slate-500">
                      <code className="font-mono text-[11px]">out/de_analysis/&lt;contrast&gt;/results_dtu_transcript.tsv</code>
                      {' '}and{' '}
                      <code className="font-mono text-[11px]">results_dtu_gene.tsv</code> — pick both.
                      Differential transcript <b>usage</b> is not computed here: DEXSeq cannot load
                      in the browser (it needs Rsamtools, which has no WebAssembly build), and a
                      different engine’s numbers under DEXSeq’s name would be worse than none.
                      Supply these and the bundle carries the real thing.
                    </p>
                    {dexseqNote && <p className="mt-1 text-xs text-emerald-700 dark:text-emerald-400">{dexseqNote}</p>}
                  </div>
                </div>
              )}
              <div className="grid gap-3 sm:grid-cols-3">
                <label className="text-sm">Method
                  <select className="input mt-1 w-full" value={method} onChange={e => setMethod(e.target.value as Method)}>
                    <option value="deseq2">DESeq2 (gold standard)</option>
                    <option value="limma">limma-voom (fast)</option>
                  </select></label>
                <label className="text-sm">Project name
                  <input className="input mt-1 w-full" value={project} onChange={e => setProject(e.target.value)} /></label>
                <label className="text-sm">Species
                  <input className="input mt-1 w-full" value={species} onChange={e => setSpecies(e.target.value)} /></label>
              </div>

              {/* Shrinkage is apeglm or nothing. There is no third option and
                  there will not be one: ashr was here and failed in both
                  directions on real data — see lib/webr.ts. */}
              {method === 'deseq2' && (
                <div className="mt-3">
                  <label className="text-sm">Fold-change shrinkage
                    <select className="input mt-1 w-full" value={shrink}
                      onChange={e => setShrink(e.target.value as 'none' | 'apeglm')}>
                      <option value="none">none — report the maximum likelihood estimate</option>
                      <option value="apeglm">apeglm — shrink low-information estimates</option>
                    </select></label>
                  <p className="mt-1 text-[11px] leading-relaxed text-slate-500 dark:text-slate-400">
                    {shrink === 'none'
                      ? <>Fold changes are reported as fitted. A gene with few counts can show a
                          very large one on very little evidence — check its lfcSE before
                          believing it.</>
                      : <>apeglm pulls estimates toward zero in proportion to how little the data
                          says, so a 20-fold change measured on 30 counts stops outranking a
                          two-fold change measured on thousands. It needs a coefficient rather
                          than a contrast, so each comparison costs an extra Wald re-test —
                          a few seconds per comparison, and the model is not refitted. The
                          unshrunk estimate is exported beside it either way.</>}
                  </p>
                </div>
              )}
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
            <div className="card p-6">
              <div className="text-sm uppercase tracking-wide text-slate-400">Converted &mdash; the lab&rsquo;s job is done</div>
              <h2 className="mt-1 text-base font-semibold">
                {sanitize(project)}_bundle.zip is ready for RNA-seq Studio
              </h2>
              <div className="mb-4 mt-2 overflow-hidden rounded-lg border border-slate-200 dark:border-slate-700">
                <table className="w-full text-sm">
                  <tbody>
                    {result.contrasts.map(c => (
                      <tr key={c.id} className="border-t border-slate-100 first:border-0 dark:border-slate-800">
                        <td className="px-3 py-1.5">
                          {c.label}
                          {c.kind === 'interaction' && (
                            <span className="ml-2 rounded bg-amber-100 px-1.5 py-0.5 text-[10px] font-semibold uppercase text-amber-700 dark:bg-amber-500/20 dark:text-amber-300">interaction</span>
                          )}
                        </td>
                        <td className="px-3 py-1.5 text-right font-mono text-[13px] tabular-nums">
                          {c.nDeg.toLocaleString()} <span className="text-slate-400">DEG</span>
                        </td>
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
              <p className="mb-4 text-xs text-slate-400">
                padj &lt; 0.05 · {method === 'limma' ? 'limma-voom' : 'DESeq2'} — a sanity check on the run.
                The tables, volcano and enrichment are in the studio.
              </p>
              <div className="mb-5"><Flow at="done" /></div>

              <div className="grid gap-2.5">
                <Handoff n={1} title="Save the bundle" done={saved}>
                  {zipUrl && (
                    <a className={saved ? 'btn' : 'btn btn-primary'} href={zipUrl}
                      download={`${sanitize(project)}_bundle.zip`} onClick={() => setSaved(true)}>
                      ⭳ Download {sanitize(project)}_bundle.zip{saved ? ' again' : ''}
                    </a>
                  )}
                </Handoff>
                <Handoff n={2} title="Open it in RNA-seq Studio">
                  <a className={saved ? 'btn btn-primary' : 'btn'} href={EXPLORER_URL}
                    target="_blank" rel="noopener noreferrer">
                    Open RNA-seq Studio &rarr;
                  </a>
                  <p className="mt-2 text-xs text-slate-400">
                    Drop the downloaded .zip onto its page. Your data never left this device, and
                    it does not leave it there either.
                  </p>
                </Handoff>
              </div>
            </div>
            <div className="flex justify-center">
              <button className="btn" onClick={reset}>Convert another dataset</button>
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

function Handoff({ n, title, done, children }: {
  n: number; title: string; done?: boolean; children: React.ReactNode
}) {
  return (
    <div className="flex gap-3 rounded-xl bg-slate-50 p-3.5 dark:bg-slate-800/50">
      <span className={`grid h-[22px] w-[22px] flex-none place-items-center rounded-md text-[11px] font-bold text-white ${done ? 'bg-emerald-500' : 'bg-indigo-500'}`}>
        {done ? '✓' : n}
      </span>
      <div className="min-w-0 flex-1">
        <div className="mb-2 text-sm font-semibold">{title}</div>
        {children}
      </div>
    </div>
  )
}

function Steps({ step }: { step: Step }) {
  const items: { id: Step; label: string }[] = [
    { id: 'upload', label: 'Upload' }, { id: 'design', label: 'Design' },
    { id: 'run', label: 'Run' }, { id: 'result', label: 'Bundle' },
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
