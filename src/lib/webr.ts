// webR loader + differential-expression engines (limma-voom, DESeq2), all in-browser.

export type Method = 'limma' | 'deseq2'

/**
 * One question to ask of the fitted model.
 *
 * `plus` and `minus` are group labels. A simple comparison is one group on each
 * side; an INTERACTION is two on each:
 *
 *   (KO_Cold - Ctrl_Cold) - (KO_Thermo - Ctrl_Thermo)
 *     plus  = [KO_Cold,   Ctrl_Thermo]
 *     minus = [Ctrl_Cold, KO_Thermo]
 *
 * Weights are UNIT (+1 per plus group, -1 per minus group), never averaged.
 * For a one-vs-one comparison the two are identical; for an interaction they are
 * not, and dividing by the group count would report exactly half the true
 * effect size with correct-looking p-values — the difference between "the cold
 * response is blunted by 1.5 log2" and "by 0.75". Averaging weights would only
 * be right for a "mean of several groups vs control" contrast, which this app
 * never builds.
 *
 * Writing both forms as one object is why this needs only ONE model fit: under
 * a cell-means design (`~ 0 + group`) every one of them is a linear combination
 * of group means. Fitting `~ genotype * temperature` separately just to reach
 * the interaction would answer the same question from a second, differently
 * parameterised fit — and two fits invite two sets of numbers that disagree at
 * the last decimal for reasons nobody wants to explain to a reviewer.
 */
export interface ContrastRequest {
  id: string
  label: string
  numerator: string
  denominator: string
  kind: 'pairwise' | 'interaction'
  plus: string[]
  minus: string[]
}

export interface AnalysisInput {
  countsCsv: string                                   // gene_id,<sample>… raw counts
  samples: { sample: string; group: string; [covariate: string]: string }[]
  groupLevels: string[]                               // every group, in display order
  contrasts: ContrastRequest[]
  method: Method
}

export interface ContrastResult {
  id: string
  label: string
  numerator: string
  denominator: string
  kind: 'pairwise' | 'interaction'
  degCsv: string
  nDeg: number
}

export interface AnalysisResult {
  contrasts: ContrastResult[]
  normCsv: string
}

const WEBR_URL = 'https://webr.r-wasm.org/v0.6.0/webr.mjs'
// Our own WASM repo (hosts locfit for DESeq2), served alongside the app on Pages.
const LOCFIT_REPO = new URL('wasm/', document.baseURI).href.replace(/\/$/, '')

let webRPromise: Promise<any> | null = null
const installed = new Set<string>()

export async function getWebR(onLog: (m: string) => void): Promise<any> {
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

async function install(webR: any, key: string, pkgs: string[], onLog: (m: string) => void, note: string) {
  if (installed.has(key)) return
  onLog(note)
  await webR.installPackages(pkgs, {
    repos: [LOCFIT_REPO, 'https://bioc.r-universe.dev', 'https://repo.r-wasm.org'],
  })
  installed.add(key)
}

/**
 * Packages needed to OPEN an nf-core R object. A DESeqDataSet cannot be
 * deserialised without the classes that define it, so reading
 * `deseq2.dds.RData` needs DESeq2 present even though nothing is fitted yet.
 */
export function ensureObjectPackages(webR: any, needsDESeq2: boolean, onLog: (m: string) => void) {
  return needsDESeq2
    ? install(webR, 'deseq2', ['DESeq2', 'apeglm'], onLog,
        'Installing DESeq2 to open the object… (first run downloads tens of MB, then cached)')
    : install(webR, 'se', ['SummarizedExperiment'], onLog,
        'Installing SummarizedExperiment to open the object…')
}

const ensurePackages = (webR: any, method: Method, onLog: (m: string) => void) =>
  method === 'limma'
    ? install(webR, 'limma', ['limma'], onLog,
        'Installing limma… (first run downloads a few MB, then cached)')
    : install(webR, 'deseq2', ['DESeq2', 'apeglm'], onLog,
        'Installing DESeq2 + apeglm… (first run downloads ~tens of MB, then cached)')

// Group labels are recoded to g1..gN before they reach R. Real labels contain
// "+", "-" and spaces ("517E2+RSL3"), which make.names() mangles into something
// that no longer matches the contrast we asked for — silently, and only for the
// groups with awkward names.
const RECODE_R = String.raw`
  cd <- read.csv("/work/coldata.csv", colClasses = "character", check.names = FALSE)
  lv <- readLines("/work/levels.txt")
  idx <- match(cd$group, lv)
  if (anyNA(idx)) stop("a sample carries a group that is not in the level list")
  grp <- factor(paste0("g", idx), levels = paste0("g", seq_along(lv)))
  counts <- as.matrix(read.csv("/work/counts.csv", row.names = 1, check.names = FALSE))
  counts <- counts[, cd$sample, drop = FALSE]
  con <- read.csv("/work/contrasts.csv", colClasses = "character", check.names = FALSE)
  ids <- function(s) if (nzchar(s)) paste0("g", strsplit(s, ";", fixed = TRUE)[[1]]) else character(0)
`

const LIMMA_R = String.raw`local({
  suppressMessages(library(limma))
  __RECODE__
  storage.mode(counts) <- "double"
  design <- model.matrix(~ 0 + grp)
  colnames(design) <- levels(grp)

  cpm <- t(t(counts) / colSums(counts)) * 1e6
  write.csv(data.frame(gene_id = rownames(cpm), gene_name = rownames(cpm),
            round(as.data.frame(cpm), 3), check.names = FALSE), "/work/norm.csv", row.names = FALSE)

  keep <- rowSums(counts >= 10) >= max(2, min(table(grp)))
  v <- voom(counts[keep, , drop = FALSE], design)
  fit <- lmFit(v, design)

  out <- character(0)
  for (i in seq_len(nrow(con))) {
    p <- ids(con$plus[i]); m <- ids(con$minus[i])
    cv <- setNames(rep(0, ncol(design)), colnames(design))
    cv[p] <- cv[p] + 1
    cv[m] <- cv[m] - 1
    f2 <- eBayes(contrasts.fit(fit, cv))
    tt <- topTable(f2, number = Inf, sort.by = "none")
    write.csv(data.frame(gene_id = rownames(tt), gene_name = rownames(tt),
              baseMean = round(2^tt$AveExpr, 3), log2FoldChange = round(tt$logFC, 4),
              lfcSE = NA, pvalue = tt$P.Value, padj = tt$adj.P.Val),
              sprintf("/work/deg_%d.csv", i), row.names = FALSE)
    out <- c(out, sprintf("%s=%d", con$id[i], sum(tt$adj.P.Val < 0.05, na.rm = TRUE)))
  }
  paste(out, collapse = "|")
})`

const DESEQ_R = String.raw`local({
  suppressMessages(library(DESeq2))
  __RECODE__
  counts <- round(counts); storage.mode(counts) <- "integer"
  cd2 <- data.frame(grp = grp); rownames(cd2) <- cd$sample
  # ~ 0 + grp is a cell-means model: one coefficient per group, so every
  # comparison, interaction included, is a linear combination of them.
  dds <- DESeqDataSetFromMatrix(counts, cd2, ~ 0 + grp)
  dds <- tryCatch(DESeq(dds, quiet = TRUE),
                  error = function(e) suppressWarnings(DESeq(dds, fitType = "mean", quiet = TRUE)))

  nc <- counts(dds, normalized = TRUE)
  write.csv(data.frame(gene_id = rownames(nc), gene_name = rownames(nc),
            round(as.data.frame(nc), 3), check.names = FALSE), "/work/norm.csv", row.names = FALSE)

  rn <- resultsNames(dds)
  out <- character(0)
  for (i in seq_len(nrow(con))) {
    p <- paste0("grp", ids(con$plus[i])); m <- paste0("grp", ids(con$minus[i]))
    cv <- setNames(rep(0, length(rn)), rn)
    cv[p] <- cv[p] + 1
    cv[m] <- cv[m] - 1
    res <- as.data.frame(results(dds, contrast = cv))
    write.csv(data.frame(gene_id = rownames(res), gene_name = rownames(res),
              baseMean = round(res$baseMean, 3), log2FoldChange = round(res$log2FoldChange, 4),
              lfcSE = round(res$lfcSE, 4), pvalue = res$pvalue, padj = res$padj),
              sprintf("/work/deg_%d.csv", i), row.names = FALSE)
    out <- c(out, sprintf("%s=%d", con$id[i], sum(res$padj < 0.05, na.rm = TRUE)))
  }
  paste(out, collapse = "|")
})`

const csvEsc = (s: string) => JSON.stringify(String(s))

export async function runAnalysis(
  input: AnalysisInput, onLog: (m: string) => void, webRIn?: any,
): Promise<AnalysisResult> {
  const webR = webRIn ?? await getWebR(onLog)
  await ensurePackages(webR, input.method, onLog)

  const idxOf = (g: string) => input.groupLevels.indexOf(g) + 1
  const bad = input.contrasts.find(c => [...c.plus, ...c.minus].some(g => idxOf(g) === 0))
  if (bad) throw new Error(`Contrast "${bad.label}" names a group that no sample has.`)

  try { await webR.FS.mkdir('/work') } catch { /* exists */ }
  const enc = new TextEncoder()
  await webR.FS.writeFile('/work/counts.csv', enc.encode(input.countsCsv))

  const coldata = 'sample,group\n' +
    input.samples.map(s => `${csvEsc(s.sample)},${csvEsc(s.group)}`).join('\n') + '\n'
  await webR.FS.writeFile('/work/coldata.csv', enc.encode(coldata))
  await webR.FS.writeFile('/work/levels.txt', enc.encode(input.groupLevels.join('\n') + '\n'))

  // Contrasts travel as 1-based level indices, so no group label is ever parsed by R.
  const rows = input.contrasts.map(c =>
    [csvEsc(c.id), csvEsc(c.plus.map(idxOf).join(';')), csvEsc(c.minus.map(idxOf).join(';'))].join(','))
  await webR.FS.writeFile('/work/contrasts.csv', enc.encode('id,plus,minus\n' + rows.join('\n') + '\n'))

  onLog(`Running ${input.method === 'limma' ? 'limma-voom' : 'DESeq2'} — ${input.contrasts.length} contrast(s)…`)
  const script = (input.method === 'limma' ? LIMMA_R : DESEQ_R).replace('__RECODE__', RECODE_R)
  const summary: string = await webR.evalRString(script)

  const nDegOf = new Map(summary.split('|').filter(Boolean).map(kv => {
    const i = kv.lastIndexOf('=')
    return [kv.slice(0, i), parseInt(kv.slice(i + 1), 10) || 0] as const
  }))

  const dec = new TextDecoder()
  const contrasts: ContrastResult[] = []
  for (let i = 0; i < input.contrasts.length; i++) {
    const c = input.contrasts[i]
    contrasts.push({
      id: c.id, label: c.label, numerator: c.numerator, denominator: c.denominator, kind: c.kind,
      degCsv: dec.decode(await webR.FS.readFile(`/work/deg_${i + 1}.csv`)),
      nDeg: nDegOf.get(c.id) ?? 0,
    })
  }
  const normCsv = dec.decode(await webR.FS.readFile('/work/norm.csv'))
  onLog(`Done — ${contrasts.map(c => `${c.label}: ${c.nDeg}`).join(' · ')} DEGs at padj < 0.05.`)
  return { contrasts, normCsv }
}
