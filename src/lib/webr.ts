import { blocksOf, cpmCsv, mapLimit, poolSize, splitByBlock } from './parallel'
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
  /**
   * Which fit answers this contrast. Every request in a run either carries one
   * or none carries one; a mixture would mean two contrasts sharing a label and
   * not a model. Absent is the whole dataset as a single fit — the default, and
   * the only thing that existed before blocking.
   */
  block?: string
}

export interface AnalysisInput {
  countsCsv: string                                   // gene_id,<sample>… raw counts
  /** `block` names the fit a sample belongs to; absent or '' means one fit for all. */
  samples: { sample: string; group: string; block?: string; [covariate: string]: string | undefined }[]
  groupLevels: string[]                               // every group, in display order
  contrasts: ContrastRequest[]
  method: Method
  /**
   * Fold-change shrinkage. apeglm, or nothing — never ashr.
   *
   * ashr was here and is gone. It failed in both directions on real data:
   * measured over the 110 tables of an 11-tissue ageing atlas, it shrank
   * well-measured significant effects to 0.67 of their MLE, and left 287 genes
   * carrying |log2FC| > 5 — up to 25.4, a 44-million-fold change — because
   * those are genes with zero counts in one group where the estimate is
   * unbounded and its normal-mixture prior cannot pull it back. On the same
   * contrast apeglm gave 0.86 on the well-measured genes and 0.002 on the
   * unbounded one. None of ashr's settings changed either number.
   */
  shrink?: Shrink
}

/** Shrinkage estimators this app will use. There is deliberately no third. */
export type Shrink = 'none' | 'apeglm'

export interface ContrastResult {
  id: string
  label: string
  numerator: string
  denominator: string
  kind: 'pairwise' | 'interaction'
  degCsv: string
  nDeg: number
  /** The fit that answered it, carried through to meta.json. */
  block?: string
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

/**
 * An ADDITIONAL webR instance — its own worker, its own R, its own heap.
 *
 * Deliberately not the shared singleton above: that one is reused across
 * uploads and holds any state a previous run left behind, which is exactly what
 * a pool must not share. These are created for one run and closed after it.
 */
async function spawnWebR(mod: any): Promise<any> {
  const webR = new mod.WebR({ interactive: false })
  await webR.init()
  return webR
}

/** The module handle, so extra instances do not re-import it. */
let webRModule: Promise<any> | null = null
const getModule = () => (webRModule ??= import(/* @vite-ignore */ WEBR_URL))

export async function getWebR(onLog: (m: string) => void): Promise<any> {
  if (!webRPromise) {
    webRPromise = (async () => {
      if (!self.crossOriginIsolated) onLog('⚠ not cross-origin isolated yet — the page may reload once.')
      onLog('Loading webR (R 4.6.0)…')
      const mod: any = await getModule()
      /**
       * There is no memory or CPU dial to turn here, and it is worth saying so
       * rather than leaving someone to look for one. webR is a single wasm32
       * thread: no BLAS, no parallelism, a 4 GB address ceiling, and no heap
       * option in its API. The browser will be several times slower than native
       * R on the same data and the only real lever is doing less work.
       *
       * What CAN be chosen is the transport. With cross-origin isolation the
       * channel is backed by SharedArrayBuffer; without it webR falls back to a
       * service-worker channel where every eval and every FS write is an
       * order of magnitude dearer — which on a run that writes a 29 MB matrix
       * and reads back 110 tables is the difference between slow and unusable.
       * GitHub Pages cannot send COOP/COEP, so coi-serviceworker.js installs a
       * worker and reloads once; until that lands, isolation is off.
       */
      const webR = new mod.WebR({ interactive: false })
      await webR.init()
      onLog(self.crossOriginIsolated
        ? 'webR ready (cross-origin isolated — fast SharedArrayBuffer channel).'
        : 'webR ready — NOT cross-origin isolated, so the slow channel is in use. '
          + 'Reload the page once to let the service worker enable it.')
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
    ? install(webR, 'deseq2', ['DESeq2', 'ashr'], onLog,
        'Installing DESeq2 to open the object… (first run downloads tens of MB, then cached)')
    : install(webR, 'se', ['SummarizedExperiment'], onLog,
        'Installing SummarizedExperiment to open the object…')
}

/**
 * Packages for one POOL WORKER, quietly and without the shared `installed` set.
 *
 * That set is a per-module memo of what the singleton has installed; a fresh
 * instance has installed nothing, so consulting it would skip the install and
 * leave the worker without DESeq2. Each instance has its own virtual
 * filesystem and has to be furnished separately.
 */
/** Packages for an engine, plus apeglm only when shrinkage is actually asked for. */
export const packagesFor = (method: Method, shrink: Shrink): string[] =>
  method === 'limma' ? ['limma']
    : shrink === 'apeglm' ? ['DESeq2', 'apeglm'] : ['DESeq2']

/**
 * DESeq2 + satuRn, for the isoform layer.
 *
 * satuRn comes from bioc.r-universe.dev, which is already in the repo list, and
 * every one of its dependencies has a WebAssembly build. DEXSeq does not: it
 * imports Rsamtools, which wraps htslib and is not built for wasm anywhere, so
 * it fails at `library()` after the download. Checked against both repo indexes,
 * not assumed.
 */
export const installDtu = (webR: any, onLog: (m: string) => void) =>
  install(webR, 'dtu', ['DESeq2', 'satuRn'], onLog,
    'Installing DESeq2 and satuRn for the isoform layer… (first run downloads tens of MB, then cached)')

export const installFor = (webR: any, method: Method, shrink: Shrink) =>
  webR.installPackages(packagesFor(method, shrink), {
    repos: [LOCFIT_REPO, 'https://bioc.r-universe.dev', 'https://repo.r-wasm.org'],
  })

const ensurePackages = (webR: any, method: Method, shrink: Shrink, onLog: (m: string) => void) =>
  method === 'limma'
    ? install(webR, 'limma', ['limma'], onLog,
        'Installing limma… (first run downloads a few MB, then cached)')
    : install(webR, `deseq2:${shrink}`, packagesFor(method, shrink), onLog,
        `Installing ${packagesFor(method, shrink).join(' + ')}… `
        + '(first run downloads ~tens of MB, then cached)')

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

  # BLOCKS. One fit per level, contrasts only within it. An unblocked run is
  # the single block "", which takes every sample and every group — the same
  # object the old code built, so its numbers are unchanged.
  if (is.null(cd$block)) cd$block <- ""
  blocks <- unique(cd$block)
  if (is.null(con$block)) con$block <- ""

  # ONE FIT, THEN PER-CONTRAST EXTRACTION. This is the architecture, and it is
  # not negotiable per contrast.
  #
  # ~ 0 + grp is a cell-means model: one coefficient per group over EVERY
  # sample, so every comparison the app can build — pairwise, pooled-side,
  # interaction — is a linear combination of those group means, extracted from
  # the single fit. Dispersion is estimated once across all groups, which is
  # what DESeq2's own FAQ recommends for multi-group data ("the final dispersion
  # value will incorporate the within-group variability across all groups"), and
  # it means every contrast in a run is answered by the same model rather than
  # by a different one each time.
  #
  # Do NOT narrow the fit to the groups a contrast happens to name. Two
  # contrasts from one run would then come from two different models, with two
  # dispersion estimates, and the same comparison would change its p-value
  # depending on what else was ticked beside it. The contrasts are the cheap
  # part; keep them that way and pay the fit once.
  #
  # A BLOCK is not that. It is a declared partition of the SAMPLES, fixed
  # before any contrast is chosen, so every comparison inside a block still
  # comes from that block's one fit and still cannot be moved by what is ticked
  # beside it. What blocking buys is the thing one global fit cannot give: a
  # dispersion estimated among comparable samples. DESeq2 has one dispersion
  # per gene and no way to vary it by group, so eleven tissues in one fit means
  # brown fat setting the variance for hypothalamus.
  fitnote <- if (length(blocks) > 1)
    sprintf("%d fits, one per block: %d samples, %d groups", length(blocks), ncol(counts), nlevels(grp))
  else sprintf("one fit: %d samples, %d groups", ncol(counts), nlevels(grp))
`

const LIMMA_R = String.raw`local({
  suppressMessages(library(limma))
  __RECODE__
  storage.mode(counts) <- "double"

  # CPM over EVERY sample, before any block is taken apart. A share of the
  # library is comparable between blocks in a way no per-block model quantity
  # is; see the DESeq2 path below for why that matters.
  writeNorm <- __WRITENORM__
  if (writeNorm) {
    cpm <- t(t(counts) / colSums(counts)) * 1e6
    write.csv(data.frame(gene_id = rownames(cpm), gene_name = rownames(cpm),
              round(as.data.frame(cpm), 3), check.names = FALSE), "/work/norm.csv", row.names = FALSE)
    rm(cpm)
  }

  out <- character(0)
  for (b in blocks) {
    sel <- which(cd$block == b)
    cb <- counts[, sel, drop = FALSE]
    gb <- droplevels(grp[sel])
    rows <- which(con$block == b)
    if (!length(rows)) next
    if (nlevels(gb) < 2)
      stop(sprintf("block %s has %d group(s); its contrasts cannot be fitted", b, nlevels(gb)))

    design <- model.matrix(~ 0 + gb)
    colnames(design) <- levels(gb)
    # The filter is per block: a gene silent in this tissue should not be
    # fitted here, nor sit in this block's multiple-testing correction.
    keep <- rowSums(cb >= 10) >= max(2, min(table(gb)))
    fitnote <- c(fitnote, sprintf("%s: %d samples, %d groups, %d of %d genes",
                                  if (nzchar(b)) b else "all", length(sel),
                                  nlevels(gb), sum(keep), length(keep)))
    if (!writeNorm) writeLines(rownames(cb)[keep], "/work/kept.txt")
    v <- voom(cb[keep, , drop = FALSE], design)
    fit <- lmFit(v, design)

    for (i in rows) {
      p <- ids(con$plus[i]); m <- ids(con$minus[i])
      if (!all(c(p, m) %in% colnames(design)))
        stop(sprintf("contrast %s names a group that is not in block %s", con$id[i], b))
      cv <- setNames(rep(0, ncol(design)), colnames(design))
      cv[p] <- cv[p] + 1
      cv[m] <- cv[m] - 1
      f2 <- eBayes(contrasts.fit(fit, cv))
      tt <- topTable(f2, number = Inf, sort.by = "none")
      write.csv(data.frame(gene_id = rownames(tt), gene_name = rownames(tt),
                baseMean = round(2^tt$AveExpr, 3), log2FoldChange = round(tt$logFC, 4),
                lfcSE = NA, pvalue = signif(tt$P.Value, 4), padj = signif(tt$adj.P.Val, 4)),
                sprintf("/work/deg_%d.csv", i), row.names = FALSE)
      out <- c(out, sprintf("%s=%d", con$id[i], sum(tt$adj.P.Val < 0.05, na.rm = TRUE)))
    }
  }
  writeLines(fitnote, "/work/fit.txt")
  paste(out, collapse = "|")
})`

const DESEQ_R = String.raw`local({
  suppressMessages(library(DESeq2))
  __RECODE__
  counts <- round(counts); storage.mode(counts) <- "integer"
  SHRINK <- "__SHRINK__"

  # Where the time goes, reported back so a slow browser run can be diagnosed
  # instead of guessed at. Native R on 275 samples: fits 38 s, shrinkage 46 s,
  # results() 17 s, IO 5 s.
  .t0 <- Sys.time(); .acc <- c()
  tick <- function(lbl) {
    now <- Sys.time()
    .acc[[lbl]] <<- (if (is.null(.acc[[lbl]])) 0 else .acc[[lbl]]) +
      as.numeric(difftime(now, .t0, units = "secs"))
    .t0 <<- now
  }
  tick("read+prep")

  # DROP GENES NO BLOCK WILL FIT, ONCE, BEFORE ANYTHING ELSE.
  #
  # Each block filters again below, so this changes no result — it only removes
  # rows that every block would have removed anyway. On the 11-tissue atlas that
  # is 10,966 of 34,514 (32%), including 1,742 that are zero in all 275 samples.
  # They were being carried through the whole run: held in the matrix, copied
  # into each block's subset, and written into normalized_counts.csv. In wasm,
  # where the address space is 4 GB and there is no swap, peak memory is the
  # thing most likely to end a large run outright.
  #
  # It also makes normalized_counts.csv mean one thing — genes detectable in at
  # least one block — instead of depending on whether the run was blocked.
  # LIBRARY SIZE IS WHAT WAS SEQUENCED, so the CPM denominator is taken over
  # every gene, before anything is dropped. Only the totals are kept, not a
  # second copy of the matrix — the prefilter below exists to hold peak memory
  # down and computing CPM first would undo it. limma already normalises this
  # way, and lib/parallel.ts computes the same denominator when a pool assembles
  # the matrix, so all three agree to the last digit.
  libTotals <- colSums(counts)
  libTotals[libTotals <= 0] <- 1

  anyKeep <- rep(FALSE, nrow(counts))
  for (b in blocks) {
    sel <- which(cd$block == b)
    if (!length(sel)) next
    gb <- droplevels(grp[sel])
    anyKeep <- anyKeep | (rowSums(counts[, sel, drop = FALSE] >= 10) >= max(2, min(table(gb))))
  }
  fitnote <- c(fitnote, sprintf("%d of %d genes are detectable in at least one block",
                                sum(anyKeep), length(anyKeep)))
  counts <- counts[anyKeep, , drop = FALSE]
  tick("prefilter")

  # WHAT normalized_counts.csv HOLDS DEPENDS ON WHETHER THE RUN IS BLOCKED.
  #
  # DESeq2's normalized counts are median-of-ratios against a reference built
  # from the samples in THAT object. Two blocks are two objects, so liver's
  # normalized counts sit on liver's reference and heart's on heart's, with no
  # conversion between them — putting them in one matrix would produce a file
  # that looks cross-comparable and is not. So a blocked run writes CPM, which
  # is each gene's share of its library and means the same thing in every
  # block. An unblocked run keeps median-of-ratios exactly as before.
  # meta.json records which, in counts_unit.
  # writeNorm is FALSE when this session is ONE WORKER OF SEVERAL. A worker
  # holds only its own block's columns, so it cannot write a matrix spanning
  # every sample and must not try; lib/parallel.ts assembles CPM from the whole
  # file instead. It writes the genes it fitted, and the caller unions those
  # across workers to decide which rows that matrix carries.
  writeNorm <- __WRITENORM__
  if (writeNorm && length(blocks) > 1) {
    cpm <- t(t(counts) / libTotals) * 1e6
    write.csv(data.frame(gene_id = rownames(cpm), gene_name = rownames(cpm),
              round(as.data.frame(cpm), 3), check.names = FALSE), "/work/norm.csv", row.names = FALSE)
    rm(cpm)
  }

  out <- character(0)
  shrunk <- 0L
  ncon <- 0L
  for (b in blocks) {
    sel <- which(cd$block == b)
    cb <- counts[, sel, drop = FALSE]
    gb <- droplevels(grp[sel])
    rows <- which(con$block == b)
    if (!length(rows)) next
    if (nlevels(gb) < 2)
      stop(sprintf("block %s has %d group(s); its contrasts cannot be fitted", b, nlevels(gb)))

    # Drop genes nothing can be said about, WITHIN THIS BLOCK. A 34,514-row
    # mouse annotation carried 1,742 rows that are zero in all 275 samples and
    # ~10,800 that never reach 10 counts in a group's worth of them; every one
    # of those was still fitted, and each only made the correction harsher.
    # Applied per block the rule also stops a gene expressed in one tissue from
    # entering another tissue's correction at all.
    keep <- rowSums(cb >= 10) >= max(2, min(table(gb)))
    cb <- cb[keep, , drop = FALSE]
    fitnote <- c(fitnote, sprintf("%s: %d samples, %d groups, %d of %d genes",
                                  if (nzchar(b)) b else "all", length(sel),
                                  nlevels(gb), sum(keep), length(keep)))
    cd2 <- data.frame(grp = gb); rownames(cd2) <- cd$sample[sel]
    # ~ 0 + grp is a cell-means model: one coefficient per group, so every
    # comparison, interaction included, is a linear combination of them.
    dds <- DESeqDataSetFromMatrix(cb, cd2, ~ 0 + grp)
    dds <- tryCatch(DESeq(dds, quiet = TRUE),
                    error = function(e) suppressWarnings(DESeq(dds, fitType = "mean", quiet = TRUE)))
    tick("DESeq() fits")

    if (!writeNorm) writeLines(rownames(cb), "/work/kept.txt")
    if (writeNorm && length(blocks) == 1) {
      # Median-of-ratios over EVERY gene, not just the fitted ones.
      #
      # This used to write counts(dds, normalized = TRUE), which is the matrix
      # AFTER the expression filter — so an unblocked bundle carried ~15k genes
      # in normalized_counts.csv while a blocked one carried all 34k, and the
      # same engine shipped a different gene set depending on a setting that has
      # nothing to do with which genes exist. Applying the size factors to the
      # unfiltered matrix is the same normalisation over the full annotation, so
      # a gene the model declined to test can still be plotted.
      nc <- sweep(counts[, sel, drop = FALSE], 2, sizeFactors(dds), "/")
      write.csv(data.frame(gene_id = rownames(nc), gene_name = rownames(nc),
                round(as.data.frame(nc), 3), check.names = FALSE), "/work/norm.csv", row.names = FALSE)
      rm(nc)
    }

    sf <- sizeFactors(dds)

    # THE FILTER STATISTIC MUST BE PER CONTRAST.
    #
    # results() screens out genes with no chance of significance using the mean
    # of normalised counts, and sets their padj to NA. Under one fit over every
    # group that mean spans EVERY sample, which is the wrong denominator for a
    # comparison between two of them: on an 11-tissue matrix a brain-only gene
    # keeps a healthy mean, survives the screen, and is handed a p-value from
    # samples where it reads zero — while a liver-only gene has its mean diluted
    # elevenfold and can come back NA despite being a strong liver signal.
    # Measured on a synthetic 4-group set with two tissue-specific blocks: the
    # global filter let all 731 foreign genes through, the per-contrast filter
    # excluded all 731 and lost none of the 759 genes that were real for the
    # contrast being asked.
    #
    # Blocking shrinks the problem — a block is already comparable samples —
    # but does not retire it: five ages of one tissue still differ, so the mean
    # over all five is still the wrong denominator for two of them.
    cmean <- function(cv) {
      lev <- sub("^grp", "", names(cv)[abs(cv) > 1e-12])
      inC <- as.character(gb) %in% lev
      if (!any(inC)) stop("a contrast names no sample")
      if (is.null(sf)) rowMeans(counts(dds, normalized = TRUE)[, inC, drop = FALSE])
      else rowMeans(sweep(counts(dds)[, inC, drop = FALSE], 2, sf[inC], "/"))
    }

    rn <- resultsNames(dds)

    # SHRINKAGE IS apeglm OR NOTHING.
    #
    # apeglm needs a COEFFICIENT and a cell-means design has none — every
    # comparison is a contrast between two of its coefficients. The way across
    # is to relevel so the denominator is the reference and re-run ONLY the Wald
    # test: nbinomWaldTest reuses the dispersions already estimated, so this is
    # a reparameterisation of the same model, not a second fit.
    #
    # One relevel serves every contrast sharing that denominator, so five age
    # levels need four relevels for all ten pairs rather than ten. Measured on
    # one tissue: 2.2 s per relevel, 2.8 s per apeglm call.
    #
    # p-values, padj and baseMean still come from the cell-means results()
    # above, with its per-contrast filter. Only the effect size and its standard
    # error are replaced — exactly the swap lfcShrink is for.
    shrinkFits <- new.env(parent = emptyenv())
    apeglmFor <- function(numLvl, denLvl) {
      if (SHRINK != "apeglm") return(NULL)
      if (is.null(shrinkFits[[denLvl]])) {
        d2 <- dds
        design(d2) <- ~ grp
        d2$grp <- relevel(d2$grp, ref = denLvl)
        assign(denLvl, nbinomWaldTest(d2, quiet = TRUE), envir = shrinkFits)
      }
      d2 <- shrinkFits[[denLvl]]
      cf <- paste0("grp_", numLvl, "_vs_", denLvl)
      if (!(cf %in% resultsNames(d2))) return(NULL)
      tryCatch(suppressMessages(lfcShrink(d2, coef = cf, type = "apeglm", quiet = TRUE)),
               error = function(e) NULL)
    }

    for (i in rows) {
      p <- paste0("grp", ids(con$plus[i])); m <- paste0("grp", ids(con$minus[i]))
      # A contrast is answered by ONE fit, so every group it names must be in
      # this block. Without this an unknown name indexes cv with NA and the
      # comparison comes back as noise rather than as an error.
      if (!all(c(p, m) %in% rn))
        stop(sprintf("contrast %s names a group that is not in block %s", con$id[i], b))
      cv <- setNames(rep(0, length(rn)), rn)
      cv[p] <- cv[p] + 1
      cv[m] <- cv[m] - 1
      bm <- cmean(cv)
      res <- results(dds, contrast = cv, filter = bm)
      tick("results()")

      # THE MLE IS ALWAYS EXPORTED, shrunk or not. Comparing effect sizes
      # BETWEEN fits — which is the whole of the cross-block view in the studio
      # — has to be done on estimates that were not each pulled toward their
      # own fit's prior by a different amount.
      mle <- res$log2FoldChange; mleSE <- res$lfcSE

      # Only a one-group-per-side comparison is a coefficient. A pooled side or
      # an interaction is not, so it keeps the MLE and the note says how many.
      simple <- length(ids(con$plus[i])) == 1 && length(ids(con$minus[i])) == 1
      sh <- if (simple) apeglmFor(sub("^grp", "", p), sub("^grp", "", m)) else NULL
      if (!is.null(sh)) {
        j <- match(rownames(res), rownames(sh))
        res$log2FoldChange <- sh$log2FoldChange[j]
        res$lfcSE <- sh$lfcSE[j]
        shrunk <- shrunk + 1L
      }
      tick("lfcShrink (apeglm)")

      # p-values at 4 significant figures. R writes them at full double
      # precision — 0.0435007582036718 is eighteen characters where five would
      # do — and across 110 tables of 17k rows those two columns are a third of
      # the bundle. Four figures rather than three keeps ranked lists from
      # gaining ties, and no threshold anyone applies can tell the difference.
      write.csv(data.frame(gene_id = rownames(res), gene_name = rownames(res),
                baseMean = round(bm, 3), log2FoldChange = round(res$log2FoldChange, 4),
                lfcSE = round(res$lfcSE, 4),
                pvalue = signif(res$pvalue, 4), padj = signif(res$padj, 4),
                log2FoldChange_MLE = round(mle, 4), lfcSE_MLE = round(mleSE, 4)),
                sprintf("/work/deg_%d.csv", i), row.names = FALSE)
      out <- c(out, sprintf("%s=%d", con$id[i], sum(res$padj < 0.05, na.rm = TRUE)))
      ncon <- ncon + 1L
    }
  }
  fitnote <- c(fitnote, if (SHRINK == "apeglm")
    sprintf("per-contrast filter; apeglm shrinkage on %d of %d contrasts", shrunk, ncon)
    else sprintf("per-contrast filter; no shrinkage — fold changes are the MLE (%d contrasts)", ncon))
  tick("write tables")
  fitnote <- c(fitnote, paste0("time: ", paste(sprintf("%s %.1fs", names(.acc), unlist(.acc)),
                                               collapse = " | ")))
  writeLines(fitnote, "/work/fit.txt")
  paste(out, collapse = "|")
})`

const csvEsc = (s: string) => JSON.stringify(String(s))

/**
 * One R session, over whatever samples and contrasts it is handed.
 *
 * Both paths below go through here: the single-instance run passes everything,
 * and each pool worker passes one block. `writeNorm` says which of them owns
 * normalized_counts.csv — a worker holding 25 of 275 columns cannot write it.
 */
async function runSession(
  webR: any,
  input: AnalysisInput,
  writeNorm: boolean,
): Promise<{ degCsv: string[]; nDeg: Map<string, number>; notes: string[]; kept: string[] }> {
  const idxOf = (g: string) => input.groupLevels.indexOf(g) + 1
  const enc = new TextEncoder()
  const dec = new TextDecoder()
  try { await webR.FS.mkdir('/work') } catch { /* exists */ }
  await webR.FS.writeFile('/work/counts.csv', enc.encode(input.countsCsv))
  await webR.FS.writeFile('/work/coldata.csv', enc.encode('sample,group,block\n' +
    input.samples.map(s =>
      `${csvEsc(s.sample)},${csvEsc(s.group)},${csvEsc(s.block ?? '')}`).join('\n') + '\n'))
  await webR.FS.writeFile('/work/levels.txt', enc.encode(input.groupLevels.join('\n') + '\n'))
  const rows = input.contrasts.map(c =>
    [csvEsc(c.id), csvEsc(c.plus.map(idxOf).join(';')), csvEsc(c.minus.map(idxOf).join(';')),
     csvEsc(c.block ?? '')].join(','))
  await webR.FS.writeFile('/work/contrasts.csv',
    enc.encode('id,plus,minus,block\n' + rows.join('\n') + '\n'))

  const script = (input.method === 'limma' ? LIMMA_R : DESEQ_R)
    .replaceAll('__RECODE__', RECODE_R)
    .replaceAll('__WRITENORM__', writeNorm ? 'TRUE' : 'FALSE')
    .replaceAll('__SHRINK__', input.shrink ?? 'none')
  const summary: string = await webR.evalRString(script)

  const nDeg = new Map(summary.split('|').filter(Boolean).map(kv => {
    const i = kv.lastIndexOf('=')
    return [kv.slice(0, i), parseInt(kv.slice(i + 1), 10) || 0] as const
  }))
  let notes: string[] = []
  try {
    notes = dec.decode(await webR.FS.readFile('/work/fit.txt')).split('\n').filter(Boolean)
  } catch { /* engine wrote no note */ }
  let kept: string[] = []
  if (!writeNorm) {
    try {
      kept = dec.decode(await webR.FS.readFile('/work/kept.txt')).split('\n').filter(Boolean)
    } catch { /* nothing fitted */ }
  }
  const degCsv: string[] = []
  for (let i = 0; i < input.contrasts.length; i++) {
    degCsv.push(dec.decode(await webR.FS.readFile(`/work/deg_${i + 1}.csv`)))
  }
  return { degCsv, nDeg, notes, kept }
}

export async function runAnalysis(
  input: AnalysisInput, onLog: (m: string) => void, webRIn?: any,
): Promise<AnalysisResult> {
  const idxOf = (g: string) => input.groupLevels.indexOf(g) + 1
  const bad = input.contrasts.find(c => [...c.plus, ...c.minus].some(g => idxOf(g) === 0))
  if (bad) throw new Error(`Contrast "${bad.label}" names a group that no sample has.`)

  // A contrast is answered by ONE fit, so every group it names has to live in
  // the block it was filed under. Caught here rather than in R, where it
  // surfaces as a missing coefficient halfway through an 11-fit run.
  const blockOfGroup = new Map<string, string>()
  for (const s of input.samples) blockOfGroup.set(s.group, s.block ?? '')
  const strayed = input.contrasts.find(c =>
    [...c.plus, ...c.minus].some(g => blockOfGroup.get(g) !== (c.block ?? '')))
  if (strayed) throw new Error(
    `Contrast "${strayed.label}" spans more than one block, so no single fit answers it.`)

  const blocks = blocksOf(input.samples)
  const shrink: Shrink = input.shrink ?? 'none'
  const engine = input.method === 'limma' ? 'limma-voom' : 'DESeq2'
  const workers = webRIn ? 1 : poolSize(blocks.length)

  /* ---------- one instance: unblocked runs, and small blocked ones ---------- */
  if (workers <= 1) {
    const webR = webRIn ?? await getWebR(onLog)
    await ensurePackages(webR, input.method, shrink, onLog)
    onLog(`Running ${engine} — ${input.contrasts.length} contrast(s)`
      + (blocks.length > 1 ? ` across ${blocks.length} blocks (${blocks.length} fits)…` : '…'))
    const r = await runSession(webR, input, true)
    r.notes.forEach(onLog)
    const contrasts = input.contrasts.map((c, i) => ({
      id: c.id, label: c.label, numerator: c.numerator, denominator: c.denominator,
      kind: c.kind, block: c.block, degCsv: r.degCsv[i], nDeg: r.nDeg.get(c.id) ?? 0,
    }))
    const normCsv = new TextDecoder().decode(await webR.FS.readFile('/work/norm.csv'))
    onLog(`Done — ${contrasts.length} contrast(s) at padj < 0.05.`)
    return { contrasts, normCsv }
  }

  /* ---------- a pool: one worker per block, several at a time ---------- */
  //
  // INSTALL ONCE, THEN FAN OUT. The first instance fetches DESeq2 over the
  // network; the rest find those tarballs in the browser's HTTP cache. Doing it
  // the other way costs N times the download on a cold cache and puts N
  // concurrent installs on the critical path — the one part of this that is
  // reasoned rather than measured, and the cheapest thing to keep off it.
  const jobs = splitByBlock(input.countsCsv, input.samples, input.contrasts)
  onLog(`Running ${engine} — ${input.contrasts.length} contrasts over ${jobs.length} blocks `
    + `on ${workers} workers (${blocks.length} independent fits).`)

  const mod = await getModule()
  const lead = webRIn ?? await getWebR(onLog)
  await ensurePackages(lead, input.method, shrink, onLog)

  const pool: any[] = [lead]
  if (workers > 1) {
    onLog(`Starting ${workers - 1} more R worker(s)…`)
    const extra = await Promise.all(
      Array.from({ length: workers - 1 }, async () => {
        const w = await spawnWebR(mod)
        await installFor(w, input.method, shrink)
        return w
      }))
    pool.push(...extra)
  }

  const degCsv = new Array<string>(input.contrasts.length)
  const nDeg = new Map<string, number>()
  const keptAll = new Set<string>()
  let done = 0
  try {
    const results = await mapLimit(jobs, pool.length, async (job, slot) => {
      const sub: AnalysisInput = {
        countsCsv: job.countsCsv,
        samples: input.samples.filter(s => (s.block ?? '') === job.block),
        groupLevels: input.groupLevels,
        contrasts: job.contrasts.map(c => c.request),
        method: input.method,
      }
      const r = await runSession(pool[slot], sub, false)
      onLog(`  [${++done}/${jobs.length}] ${job.block}: `
        + (r.notes.find(n => n.startsWith(job.block)) ?? `${job.contrasts.length} contrast(s)`))
      return { job, r }
    })
    // Assembled by the ORIGINAL index, so the order never depends on which
    // worker finished first.
    for (const { job, r } of results) {
      job.contrasts.forEach((c, k) => { degCsv[c.index] = r.degCsv[k] })
      r.nDeg.forEach((v, k) => nDeg.set(k, v))
      r.kept.forEach(g => keptAll.add(g))
    }
  } finally {
    // Close only what this run created; the shared singleton outlives it.
    await Promise.all(pool.slice(1).map(w => w.close?.()))
  }

  const missing = degCsv.findIndex(x => x === undefined)
  if (missing >= 0) {
    throw new Error(
      `No table came back for "${input.contrasts[missing].label}". `
      + `The run is incomplete, so no bundle was built.`)
  }

  const contrasts: ContrastResult[] = input.contrasts.map((c, i) => ({
    id: c.id, label: c.label, numerator: c.numerator, denominator: c.denominator,
    kind: c.kind, block: c.block, degCsv: degCsv[i], nDeg: nDeg.get(c.id) ?? 0,
  }))
  onLog(`Assembling normalized counts over ${keptAll.size.toLocaleString()} genes…`)
  const normCsv = cpmCsv(input.countsCsv, keptAll)
  onLog(`Done — ${contrasts.length} contrasts from ${jobs.length} fits on ${pool.length} workers.`)
  return { contrasts, normCsv }
}
