// webR loader + differential-expression engines (limma-voom, DESeq2), all in-browser.

export type Method = 'limma' | 'deseq2'

export interface AnalysisInput {
  countsCsv: string                                   // raw counts CSV (gene rows, sample cols; col 1 = gene id)
  samples: { sample: string; condition: string }[]
  control: string                                     // reference condition
  method: Method
}

export interface AnalysisResult {
  degCsv: string
  normCsv: string
  nDeg: number
  numerator: string
  denominator: string
}

const WEBR_URL = 'https://webr.r-wasm.org/v0.6.0/webr.mjs'
// Our own WASM repo (hosts locfit for DESeq2), served alongside the app on Pages.
const LOCFIT_REPO = new URL('wasm/', document.baseURI).href.replace(/\/$/, '')

let webRPromise: Promise<any> | null = null
const installed = new Set<Method>()

async function getWebR(onLog: (m: string) => void): Promise<any> {
  if (!webRPromise) {
    webRPromise = (async () => {
      if (!self.crossOriginIsolated) onLog('⚠ not cross-origin isolated yet — the page may reload once.')
      onLog('Loading webR (R 4.6.0)…')
      const mod: any = await import(/* @vite-ignore */ WEBR_URL)
      const webR = new mod.WebR()
      await webR.init()
      onLog('webR ready.')
      return webR
    })()
  }
  return webRPromise
}

async function ensurePackages(webR: any, method: Method, onLog: (m: string) => void) {
  if (installed.has(method)) return
  if (method === 'limma') {
    onLog('Installing limma… (first run downloads a few MB, then cached)')
    await webR.installPackages(['limma'], { repos: ['https://bioc.r-universe.dev', 'https://repo.r-wasm.org'] })
  } else {
    onLog('Installing DESeq2 + apeglm… (first run downloads ~tens of MB, then cached)')
    await webR.installPackages(['DESeq2', 'apeglm'], {
      repos: [LOCFIT_REPO, 'https://bioc.r-universe.dev', 'https://repo.r-wasm.org'],
    })
  }
  installed.add(method)
}

const LIMMA_R = `local({
  suppressMessages(library(limma))
  REF <- __REF__
  counts <- as.matrix(read.csv("/work/counts.csv", row.names = 1, check.names = FALSE))
  storage.mode(counts) <- "double"
  # colClasses="character": a group named "517E2" is otherwise read as
  # scientific notation (51700) and then matches nothing.
  cd <- read.csv("/work/coldata.csv", colClasses = "character", check.names = FALSE)
  counts <- counts[, cd$sample, drop = FALSE]
  grp <- relevel(factor(cd$condition), ref = REF)
  design <- model.matrix(~ grp)
  cpm <- t(t(counts) / colSums(counts)) * 1e6
  write.csv(data.frame(gene_id = rownames(cpm), gene_name = rownames(cpm),
            round(as.data.frame(cpm), 3), check.names = FALSE), "/work/norm.csv", row.names = FALSE)
  keep <- rowSums(counts >= 10) >= max(2, min(table(grp)))
  v <- voom(counts[keep, , drop = FALSE], design)
  fit <- eBayes(lmFit(v, design))
  coefName <- colnames(design)[2]
  tt <- topTable(fit, coef = coefName, number = Inf, sort.by = "none")
  write.csv(data.frame(gene_id = rownames(tt), gene_name = rownames(tt),
            baseMean = round(2^tt$AveExpr, 3), log2FoldChange = round(tt$logFC, 4),
            lfcSE = NA, pvalue = tt$P.Value, padj = tt$adj.P.Val), "/work/deg.csv", row.names = FALSE)
  sprintf("%s|%s|%d", sub("^grp", "", coefName), REF, sum(tt$adj.P.Val < 0.05, na.rm = TRUE))
})`

const DESEQ_R = `local({
  suppressMessages(library(DESeq2))
  REF <- __REF__
  counts <- round(as.matrix(read.csv("/work/counts.csv", row.names = 1, check.names = FALSE)))
  storage.mode(counts) <- "integer"
  cd <- read.csv("/work/coldata.csv", colClasses = "character", check.names = FALSE)
  rownames(cd) <- cd$sample
  counts <- counts[, cd$sample, drop = FALSE]
  cd$condition <- relevel(factor(cd$condition), ref = REF)
  dds <- DESeqDataSetFromMatrix(counts, cd, ~condition)
  # Default parametric dispersion fit never calls locfit; on the rare fallback
  # (locfit is a stub in this build) retry with the locfit-free "mean" fit.
  dds <- tryCatch(DESeq(dds, quiet = TRUE),
                  error = function(e) suppressWarnings(DESeq(dds, fitType = "mean", quiet = TRUE)))
  res <- as.data.frame(results(dds))
  nc <- counts(dds, normalized = TRUE)
  write.csv(data.frame(gene_id = rownames(nc), gene_name = rownames(nc),
            round(as.data.frame(nc), 3), check.names = FALSE), "/work/norm.csv", row.names = FALSE)
  write.csv(data.frame(gene_id = rownames(res), gene_name = rownames(res),
            baseMean = round(res$baseMean, 3), log2FoldChange = round(res$log2FoldChange, 4),
            lfcSE = round(res$lfcSE, 4), pvalue = res$pvalue, padj = res$padj), "/work/deg.csv", row.names = FALSE)
  sprintf("%s|%s|%d", tail(levels(cd$condition), 1), REF, sum(res$padj < 0.05, na.rm = TRUE))
})`

export async function runAnalysis(input: AnalysisInput, onLog: (m: string) => void): Promise<AnalysisResult> {
  const webR = await getWebR(onLog)
  await ensurePackages(webR, input.method, onLog)

  try { await webR.FS.mkdir('/work') } catch { /* exists */ }
  const enc = new TextEncoder()
  await webR.FS.writeFile('/work/counts.csv', enc.encode(input.countsCsv))
  const coldata = 'sample,condition\n' +
    input.samples.map(s => `${JSON.stringify(s.sample)},${JSON.stringify(s.condition)}`).join('\n') + '\n'
  await webR.FS.writeFile('/work/coldata.csv', enc.encode(coldata))

  onLog(`Running ${input.method === 'limma' ? 'limma-voom' : 'DESeq2'}…`)
  const script = (input.method === 'limma' ? LIMMA_R : DESEQ_R).replace('__REF__', JSON.stringify(input.control))
  const summary = await webR.evalRString(script)

  const dec = new TextDecoder()
  const degCsv = dec.decode(await webR.FS.readFile('/work/deg.csv'))
  const normCsv = dec.decode(await webR.FS.readFile('/work/norm.csv'))
  const [numerator, denominator, nDeg] = summary.split('|')
  onLog(`Done — ${nDeg} DEGs at padj < 0.05.`)
  return { degCsv, normCsv, nDeg: parseInt(nDeg, 10) || 0, numerator, denominator }
}
