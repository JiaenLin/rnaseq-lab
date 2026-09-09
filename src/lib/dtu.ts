// Transcript-level differential expression (DTE) and differential transcript
// usage (DTU), in webR.
//
// A SEPARATE R path on purpose. The gene-level DESEQ_R in webr.ts is verified
// number-for-number against the studio and against the cluster; nothing here
// edits it. This module opens its own session, on its own matrix, and writes its
// own files.
//
// DTE answers "is this isoform present at a different level?" — DESeq2 on
// transcript rows, same engine, same filter, same shrinkage.
//
// DTU answers the question long reads exist for: "did the MIX change?" A gene
// can be flat while its isoforms swap.
//
// WHY satuRn AND NOT DEXSeq. DEXSeq is what Oxford Nanopore's own
// wf-transcriptomes runs, so it was the first choice and a bundle built with it
// would have been directly comparable to the cluster's own
// `results_dtu_transcript.tsv`. It cannot run here: DEXSeq `Imports: Rsamtools`,
// Rsamtools wraps htslib, and **no WebAssembly build of it exists** — not on
// repo.r-wasm.org, not on bioc.r-universe.dev. `library(DEXSeq)` fails at
// namespace load with "there is no package called 'Rsamtools'", which is what a
// live run on the deployed app reported after the gene-level fit had already
// finished.
//
// Of the two alternatives, DRIMSeq needs only locfit (which this app already
// builds and hosts) and satuRn needs nothing at all. satuRn is chosen: it is
// built for exactly this scale, and its quasi-binomial model of each isoform's
// share is the same question DEXSeq asks.
//
// The cost is stated rather than hidden: **satuRn's numbers are not DEXSeq's**,
// so a bundle built here is no longer arithmetically comparable to the
// cluster's table. meta.json names the engine for this reason.

import type { Shrink } from './webr.ts'

export interface DtuInput {
  /** `transcript_id,<sample>… ` raw counts, as read from the pipeline. */
  txCountsCsv: string
  /** transcript_id -> gene_id. Every row of the matrix must be here. */
  txToGene: Map<string, string>
  samples: { sample: string; group: string }[]
  /** Exactly two groups. DEXSeq's `~ sample + exon + condition:exon` is a two-level test. */
  numerator: string
  denominator: string
  shrink?: Shrink
}

export interface DtuResult {
  /** DESeq2 on transcript rows: transcript_id,baseMean,log2FoldChange,lfcSE,pvalue,padj */
  dteCsv: string
  /** transcript_id,gene_id,usage_effect,pvalue,padj,gene_padj,mean_usage_num,mean_usage_den */
  dtuCsv: string
  nTested: number
  nDte: number
  nDtu: number
  notes: string[]
}

/**
 * The prefilter, stated once and applied to BOTH tests.
 *
 * Taken from the cluster's stage 03 so the two are comparable: at least ten
 * counts in total and at least three counts in at least two samples. On this
 * cohort it took 55,852 transcripts to 22,275. It is not a nicety — DEXSeq's
 * cost is roughly linear in features, and a transcript with four reads across
 * six libraries cannot support a usage test whatever the engine.
 */
export const FILTER_NOTE =
  'at least 10 counts in total and at least 3 counts in at least 2 samples'

const R = String.raw`webr::eval_js('0')
suppressMessages({ library(DESeq2); library(satuRn); library(SummarizedExperiment) })

cnt <- read.csv("/work/tx_counts.csv", row.names = 1, check.names = FALSE)
cd  <- read.csv("/work/tx_coldata.csv", stringsAsFactors = FALSE)
t2g <- read.csv("/work/tx2gene.csv", stringsAsFactors = FALSE)

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

gene <- t2g$gene_id[match(rownames(cnt), t2g$transcript_id)]
gene[is.na(gene)] <- rownames(cnt)[is.na(gene)]

# ---- DTE: DESeq2 on transcript rows -------------------------------------
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

# ---- DTU: satuRn, a quasi-binomial model of each isoform's SHARE ---------
# A gene with ONE surviving transcript carries no usage information at all --
# its single isoform is 100% of the gene in every sample by construction -- and
# no usage model can fit it. Dropping them here is what keeps the run from
# failing on a third of the matrix.
multi <- gene %in% names(which(table(gene) > 1))
notes <- c(notes, sprintf("DTU: %d transcripts in %d multi-isoform genes",
                          sum(multi), length(unique(gene[multi]))))
dtu <- data.frame(transcript_id = character(0), gene_id = character(0),
                  usage_effect = numeric(0), pvalue = numeric(0), padj = numeric(0),
                  gene_padj = numeric(0), mean_usage_num = numeric(0),
                  mean_usage_den = numeric(0))
if (sum(multi) > 1) {
  cm <- cnt[multi, , drop = FALSE]
  gm <- gene[multi]

  tx  <- data.frame(isoform_id = rownames(cm), gene_id = gm, row.names = rownames(cm))
  cd2 <- data.frame(group = grp, row.names = colnames(cm))
  se  <- SummarizedExperiment(assays = list(counts = cm), colData = cd2, rowData = tx)
  se  <- satuRn::fitDTU(object = se, formula = ~ 0 + group, parallel = FALSE, verbose = FALSE)

  design <- model.matrix(~ 0 + group, data = cd2)
  colnames(design) <- levels(grp)
  L <- matrix(0, nrow = ncol(design), ncol = 1,
              dimnames = list(colnames(design), "contrast"))
  L["__DEN__", 1] <- -1
  L["__NUM__", 1] <- 1
  se <- satuRn::testDTU(object = se, contrasts = L,
                        diagplot1 = FALSE, diagplot2 = FALSE, sort = FALSE)
  sr <- rowData(se)[["fitDTUResult_contrast"]]

  # satuRn reports a regular FDR and an EMPIRICAL one, and recommends the
  # empirical for FDR control. On a small design the empirical null can fail to
  # fit and come back all-NA, so it is used when it exists and the fallback is
  # recorded rather than silently taken.
  padj <- sr$empirical_FDR
  pval <- sr$empirical_pval
  which_p <- "satuRn empirical"
  if (all(is.na(padj))) {
    padj <- sr$regular_FDR; pval <- sr$pval
    which_p <- "satuRn regular (the empirical null did not fit)"
  }
  notes <- c(notes, paste("DTU p-values:", which_p))

  # A per-gene q-value, built rather than borrowed: Simes across a gene's own
  # transcript p-values, then BH across genes. The minimum transcript padj is
  # NOT a gene-level q and must not be labelled as one.
  gid <- as.character(rowData(se)$gene_id)
  simes <- tapply(pval, gid, function(p) {
    p <- sort(p[!is.na(p)])
    if (!length(p)) return(NA_real_)
    min(length(p) * p / seq_along(p), 1)
  })
  gq <- p.adjust(simes, method = "BH")

  # Observed usage per group: each transcript's share of its gene's total.
  # Written out rather than left to the reader, because a proportion the app
  # recomputed could disagree with the test that was actually run.
  #
  # NAMES MATTER HERE AND ARE NOT AUTOMATIC. ifelse copies the attributes of its
  # TEST argument, so ifelse(tot > 0, ...) returns a matrix carrying tot's
  # rownames -- which are GENE ids, duplicated, because tot was built by rowsum.
  # rowMeans then inherits those, and a transcript-id lookup into a gene-named
  # vector is NA on every row, with no error anywhere.
  gtot <- rowsum(cm, gm)
  share <- function(which) {
    s <- cm[, which, drop = FALSE]
    tot <- gtot[match(gm, rownames(gtot)), which, drop = FALSE]
    stats::setNames(rowMeans(ifelse(tot > 0, s / tot, NA_real_), na.rm = TRUE),
                    rownames(cm))
  }
  un <- share(grp == "__NUM__"); ud <- share(grp == "__DEN__")
  ids <- as.character(rowData(se)$isoform_id)
  dtu <- data.frame(
    transcript_id  = ids,
    gene_id        = gid,
    # satuRn's effect is on the quasi-binomial LOGIT scale -- a change in log
    # odds of usage, not a log2 fold change. Named for what it is; meta.json
    # records the scale so nothing downstream reads it as a fold change.
    usage_effect   = as.numeric(sr$estimates),
    pvalue         = as.numeric(pval),
    padj           = as.numeric(padj),
    gene_padj      = as.numeric(gq[gid]),
    mean_usage_num = as.numeric(un[ids]),
    mean_usage_den = as.numeric(ud[ids]))
}
write.csv(dtu, "/work/dtu.csv", row.names = FALSE, na = "NA")

writeLines(notes, "/work/dtu_notes.txt")
sprintf("%d|%d|%d", nrow(cnt),
        sum(!is.na(dte$padj) & dte$padj < 0.05),
        sum(!is.na(dtu$padj) & dtu$padj < 0.05))`

const csvEsc = (s: string) => JSON.stringify(String(s))

/** Run DTE and DTU for one pair of groups. `webR` must already have DTU_PACKAGES. */
export async function runDtu(
  webR: any, input: DtuInput, onLog: (m: string) => void,
): Promise<DtuResult> {
  const enc = new TextEncoder(); const dec = new TextDecoder()
  try { await webR.FS.mkdir('/work') } catch { /* exists */ }

  await webR.FS.writeFile('/work/tx_counts.csv', enc.encode(input.txCountsCsv))
  await webR.FS.writeFile('/work/tx_coldata.csv', enc.encode('sample,group\n' +
    input.samples.map(s => `${csvEsc(s.sample)},${csvEsc(s.group)}`).join('\n') + '\n'))
  await webR.FS.writeFile('/work/tx2gene.csv', enc.encode('transcript_id,gene_id\n' +
    [...input.txToGene].map(([t, g]) => `${csvEsc(t)},${csvEsc(g)}`).join('\n') + '\n'))

  onLog('Running DESeq2 on transcripts, then satuRn for usage.')

  const summary: string = await webR.evalRString(
    R.replaceAll('__NUM__', input.numerator)
     .replaceAll('__DEN__', input.denominator)
     .replaceAll('__SHRINK__', input.shrink ?? 'none')
     .replaceAll('__FILTERNOTE__', FILTER_NOTE))

  const [nTested, nDte, nDtu] = summary.split('|').map(n => parseInt(n, 10) || 0)
  let notes: string[] = []
  try {
    notes = dec.decode(await webR.FS.readFile('/work/dtu_notes.txt')).split('\n').filter(Boolean)
  } catch { /* none */ }
  notes.forEach(onLog)

  return {
    dteCsv: dec.decode(await webR.FS.readFile('/work/dte.csv')),
    dtuCsv: dec.decode(await webR.FS.readFile('/work/dtu.csv')),
    nTested, nDte, nDtu, notes,
  }
}
