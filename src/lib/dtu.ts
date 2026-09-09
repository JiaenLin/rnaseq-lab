// Transcript-level differential expression (DTE), and differential transcript
// usage (DTU) read from the pipeline that computed it.
//
// A SEPARATE R path on purpose. The gene-level DESEQ_R in webr.ts is verified
// number-for-number against the studio and against the cluster; nothing here
// edits it.
//
// ── WHY DTU IS NOT COMPUTED HERE ────────────────────────────────────────────
//
// DTU is the question long reads exist for — a gene can be perfectly flat while
// its dominant isoform swaps — and DEXSeq is the established tool for it, the
// one Oxford Nanopore's own wf-transcriptomes runs.
//
// DEXSeq cannot load in webR, and the blocker is upstream of this project:
//
//   1. DEXSeq declares `import(Rsamtools)` in its NAMESPACE, so `library(DEXSeq)`
//      needs Rsamtools present — even though DEXSeq calls not one Rsamtools
//      symbol (grepped across all eight of its R files: zero hits).
//   2. Rsamtools has no WebAssembly build anywhere: repo.r-wasm.org and three
//      r-universes, at R 4.4, 4.5 and 4.6.
//   3. Building it fails, and not for a reason this repo can fix. Rsamtools
//      links against Rhtslib, and the Rhtslib wasm build ships a `libhts.a`
//      whose members are HOST objects — `wasm-ld: archive member 'hts.o' is
//      neither Wasm object file nor LLVM bitcode`, for every member. htslib is
//      not actually cross-compiled for Emscripten. (Build log: the
//      "Build the WASM binaries" step of build-wasm.yml, 2026-09-09.)
//   4. DEXSeq also reaches XML through geneplotter -> annotate, whose NAMESPACE
//      has `importFrom(XML, ...)`. XML's configure fails its libxml2 link test
//      under emscripten as well.
//
// So rather than substitute a different engine and report its numbers under
// DTU's name, this app READS the DEXSeq result the pipeline already produced.
// wf-transcriptomes writes it whenever `de_analysis` is on:
//
//   out/de_analysis/<contrast>/results_dtu_transcript.tsv   per-transcript
//   out/de_analysis/<contrast>/results_dtu_gene.tsv         perGeneQValue
//
// Those are real, unmodified DEXSeq numbers. With no such file the bundle
// simply carries no DTU — the studio already degrades to the gene and DTE
// layers — rather than carrying a different test wearing DEXSeq's label.

import type { Shrink } from './webr.ts'

export interface DteInput {
  /** `transcript_id,<sample>… ` raw counts, rounded, as read from the pipeline. */
  txCountsCsv: string
  samples: { sample: string; group: string }[]
  numerator: string
  denominator: string
  shrink?: Shrink
}

export interface DteResult {
  /** transcript_id,baseMean,log2FoldChange,lfcSE,pvalue,padj */
  dteCsv: string
  nTested: number
  nDte: number
  notes: string[]
}

/**
 * The prefilter, stated once.
 *
 * Taken from the cluster's stage 03 so the two are comparable: at least ten
 * counts in total and at least three counts in at least two samples. A
 * transcript with four reads across six libraries cannot support a test.
 */
export const FILTER_NOTE =
  'at least 10 counts in total and at least 3 counts in at least 2 samples'

const R = String.raw`webr::eval_js('0')
suppressMessages(library(DESeq2))

cnt <- read.csv("/work/tx_counts.csv", row.names = 1, check.names = FALSE)
cd  <- read.csv("/work/tx_coldata.csv", stringsAsFactors = FALSE)

cnt <- cnt[, cd$sample, drop = FALSE]
# Already integers: the caller rounds bambu's fractional expected counts once,
# before the gene matrix is summed and before this matrix is written, so the
# gene layer and the transcript layer are sums of the SAME numbers. Coerced
# rather than rounded here, so a non-integer arriving is a loud failure and not
# a second, different rounding.
cnt <- as.matrix(cnt)
storage.mode(cnt) <- "integer"
if (anyNA(cnt)) stop("transcript counts contain NA after coercion; expected pre-rounded integers")

grp <- factor(cd$group, levels = c("__DEN__", "__NUM__"))

keep <- rowSums(cnt) >= 10 & rowSums(cnt >= 3) >= 2
cnt  <- cnt[keep, , drop = FALSE]
notes <- sprintf("filter: %d of %d transcripts kept (%s)",
                 nrow(cnt), length(keep), "__FILTERNOTE__")

dds <- DESeqDataSetFromMatrix(cnt, data.frame(grp = grp, row.names = colnames(cnt)), ~ grp)
dds <- tryCatch(DESeq(dds, quiet = TRUE),
                error = function(e) suppressWarnings(DESeq(dds, fitType = "mean", quiet = TRUE)))
res <- results(dds, contrast = c("grp", "__NUM__", "__DEN__"))
if (identical("__SHRINK__", "apeglm")) {
  res <- tryCatch(lfcShrink(dds, coef = resultsNames(dds)[length(resultsNames(dds))],
                            type = "apeglm", res = res, quiet = TRUE),
                  error = function(e) { notes <<- c(notes, "apeglm failed; reporting the MLE"); res })
}
dte <- data.frame(transcript_id = rownames(res), baseMean = res$baseMean,
                  log2FoldChange = res$log2FoldChange, lfcSE = res$lfcSE,
                  pvalue = res$pvalue, padj = res$padj)
write.csv(dte, "/work/dte.csv", row.names = FALSE, na = "NA")
writeLines(notes, "/work/dte_notes.txt")
sprintf("%d|%d", nrow(cnt), sum(!is.na(dte$padj) & dte$padj < 0.05))`

const csvEsc = (s: string) => JSON.stringify(String(s))

/** Transcript-level DESeq2 for one pair of groups. Same engine as the gene layer. */
export async function runDte(
  webR: any, input: DteInput, onLog: (m: string) => void,
): Promise<DteResult> {
  const enc = new TextEncoder(); const dec = new TextDecoder()
  try { await webR.FS.mkdir('/work') } catch { /* exists */ }

  await webR.FS.writeFile('/work/tx_counts.csv', enc.encode(input.txCountsCsv))
  await webR.FS.writeFile('/work/tx_coldata.csv', enc.encode('sample,group\n' +
    input.samples.map(s => `${csvEsc(s.sample)},${csvEsc(s.group)}`).join('\n') + '\n'))

  const summary: string = await webR.evalRString(
    R.replaceAll('__NUM__', input.numerator)
     .replaceAll('__DEN__', input.denominator)
     .replaceAll('__SHRINK__', input.shrink ?? 'none')
     .replaceAll('__FILTERNOTE__', FILTER_NOTE))

  const [nTested, nDte] = summary.split('|').map(n => parseInt(n, 10) || 0)
  let notes: string[] = []
  try {
    notes = dec.decode(await webR.FS.readFile('/work/dte_notes.txt')).split('\n').filter(Boolean)
  } catch { /* none */ }
  notes.forEach(onLog)

  return { dteCsv: dec.decode(await webR.FS.readFile('/work/dte.csv')), nTested, nDte, notes }
}

// ── Reading the pipeline's DEXSeq result ────────────────────────────────────

export interface DexseqTables {
  /** `results_dtu_transcript.tsv`: featureID, groupID, log2FoldChange, pvalue, padj */
  transcript: string[][]
  /** `results_dtu_gene.tsv`: GENEID, qval — DEXSeq's own perGeneQValue. */
  gene?: string[][]
}

const norm = (s: string) => (s ?? '').trim().toLowerCase().replace(/[_.]/g, '')
const col = (h: string[], ...names: string[]) => {
  const want = new Set(names.map(norm))
  return h.findIndex(x => want.has(norm(x)))
}
const clean = (v: string | undefined) => {
  const t = (v ?? '').trim()
  return t === 'NA' || t === 'NaN' || t === 'null' ? '' : t
}

export interface DtuFromPipeline {
  dtuCsv: string
  nTested: number
  nDtu: number
  nGeneQ: number
  notes: string[]
}

/**
 * Turn the pipeline's DEXSeq tables into the bundle's `dtu_<contrast>.csv`.
 *
 * Every statistic here is DEXSeq's, carried through unaltered: the effect, the
 * p-value, the FDR and the per-gene q-value. The only numbers this function
 * computes are the two observed usage SHARES, which are not a test — they are
 * each transcript's fraction of its gene's counts, and they exist so the studio
 * can draw the mix without recomputing anything the test depended on.
 */
export function dtuFromPipeline(
  tables: DexseqTables,
  txCountsCsv: string,
  samples: { sample: string; group: string }[],
  numerator: string,
  denominator: string,
): DtuFromPipeline {
  const th = tables.transcript[0].map(x => String(x ?? '').trim())
  const iTx = col(th, 'featureID', 'transcript_id')
  const iGene = col(th, 'groupID', 'gene_id')
  const iLfc = col(th, 'log2FoldChange', 'log2fold')
  const iP = col(th, 'pvalue')
  const iQ = col(th, 'padj')
  if (iTx < 0 || iGene < 0 || iP < 0) {
    throw new Error(
      'That does not look like a DEXSeq transcript table. Expected featureID, groupID and ' +
      'pvalue columns — wf-transcriptomes writes it to ' +
      'out/de_analysis/<contrast>/results_dtu_transcript.tsv.')
  }

  // DEXSeq's own perGeneQValue, when the gene table came with it.
  const geneQ = new Map<string, string>()
  if (tables.gene?.length) {
    const gh = tables.gene[0].map(x => String(x ?? '').trim())
    const gi = col(gh, 'GENEID', 'gene_id', 'groupID')
    const qi = col(gh, 'qval', 'padj', 'qvalue')
    if (gi >= 0 && qi >= 0) {
      for (const r of tables.gene.slice(1)) {
        const g = clean(r[gi]); if (g) geneQ.set(g, clean(r[qi]))
      }
    }
  }

  // Observed shares, from the same matrix the bundle ships.
  const lines = txCountsCsv.trim().split(/\r?\n/)
  const head = lines[0].split(',').map(x => x.replace(/^"|"$/g, ''))
  const sampleCol = new Map(head.map((h, i) => [h, i] as const))
  const numIdx = samples.filter(s => s.group === numerator).map(s => sampleCol.get(s.sample)!)
    .filter(i => i != null)
  const denIdx = samples.filter(s => s.group === denominator).map(s => sampleCol.get(s.sample)!)
    .filter(i => i != null)
  const counts = new Map<string, number[]>()
  for (let r = 1; r < lines.length; r++) {
    const c = lines[r].split(',')
    counts.set(c[0].replace(/^"|"$/g, ''), c.map(Number))
  }
  const geneTotal = new Map<string, number[]>()
  const geneOf = new Map<string, string>()
  for (const row of tables.transcript.slice(1)) {
    const t = clean(row[iTx]), g = clean(row[iGene])
    if (!t || !g) continue
    geneOf.set(t, g)
    const v = counts.get(t); if (!v) continue
    let acc = geneTotal.get(g)
    if (!acc) { acc = new Array(v.length).fill(0); geneTotal.set(g, acc) }
    for (let i = 1; i < v.length; i++) acc[i] += Number.isFinite(v[i]) ? v[i] : 0
  }
  const share = (t: string, idx: number[]) => {
    const v = counts.get(t), tot = geneTotal.get(geneOf.get(t) ?? '')
    if (!v || !tot) return ''
    const xs = idx.map(i => (tot[i] > 0 ? v[i] / tot[i] : NaN)).filter(Number.isFinite)
    return xs.length ? (xs.reduce((a, b) => a + b, 0) / xs.length).toFixed(6) : ''
  }

  const out = ['transcript_id,gene_id,usage_effect,pvalue,padj,gene_padj,mean_usage_num,mean_usage_den']
  let nDtu = 0, nTested = 0
  for (const row of tables.transcript.slice(1)) {
    const t = clean(row[iTx]), g = clean(row[iGene])
    if (!t) continue
    nTested++
    const padj = iQ >= 0 ? clean(row[iQ]) : ''
    if (padj && Number(padj) < 0.05) nDtu++
    out.push([t, g, iLfc >= 0 ? clean(row[iLfc]) : '', clean(row[iP]), padj,
      geneQ.get(g) ?? '', share(t, numIdx), share(t, denIdx)].map(csvEsc).join(','))
  }

  const notes = [
    `DEXSeq from the pipeline: ${nTested.toLocaleString()} transcripts, ` +
    `${nDtu.toLocaleString()} with changed usage (FDR<0.05)`,
    geneQ.size
      ? `per-gene q-values: ${geneQ.size.toLocaleString()} from DEXSeq's perGeneQValue`
      : 'no results_dtu_gene.tsv supplied, so the bundle carries no per-gene q-value',
  ]

  /**
   * The shares are computed from the UPLOADED matrix; DEXSeq ran on the
   * pipeline's. If those two disagree about which isoforms a gene has, the
   * statistics stay right and the shares silently do not — a gene truncated to
   * one isoform reads as 100% in both groups beside a significant p-value,
   * which looks like a result and is an artifact of the upload.
   */
  const missing = [...geneOf.keys()].filter(t => !counts.has(t)).length
  if (missing) {
    notes.push(
      `${missing.toLocaleString()} of ${nTested.toLocaleString()} transcripts in the DEXSeq ` +
      `table are not in the counts matrix, so they carry no usage share. The statistics are ` +
      `unaffected; the shares are drawn from the matrix you uploaded.`)
  }
  const singleton = [...geneTotal.keys()].filter(g =>
    [...geneOf.values()].filter(x => x === g).length === 1).length
  if (singleton) {
    notes.push(
      `${singleton.toLocaleString()} gene(s) have a single isoform in the uploaded matrix, so ` +
      `their shares are 100% in both groups. If DEXSeq tested them, it saw more isoforms than ` +
      `this matrix has — check that the counts and the DEXSeq table come from the same run.`)
  }
  return { dtuCsv: out.join('\n') + '\n', nTested, nDtu, nGeneQ: geneQ.size, notes }
}
