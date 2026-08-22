// Reading nf-core/rnaseq's R objects directly, instead of asking for a CSV.
//
// nf-core/rnaseq already builds everything this app rebuilds. Its output holds
//
//   star_salmon/deseq2_qc/deseq2.dds.RData          a DESeqDataSet
//   star_salmon/salmon.merged.gene.SummarizedExperiment.rds
//   star_salmon/salmon.merged.gene_counts.tsv       the matrix
//
// The RData is strictly better than the TSV: it carries the counts AND the
// sample table AND the design formula that produced them, so the design does
// not have to be re-guessed from sample names. It is also the object the
// pipeline actually analysed, which removes a whole class of "the CSV I
// exported is not quite what was modelled" mismatch.
//
// We can read it because this app already runs R: webR loads the object with
// R's own `load()`/`readRDS()`, which is the only thing that reliably parses
// R serialization. No JS RData parser is involved.

export interface RObjectImport {
  countsCsv: string                      // gene_id,<sample>… — raw integer counts
  samples: string[]
  /** colData from the object: sample -> {column: value}. Empty when it carried none. */
  colData: Record<string, Record<string, string>>
  colDataColumns: string[]
  /** the object's own design, e.g. "~ condition", when it had one */
  designFormula: string | null
  nGenes: number
  source: 'DESeqDataSet' | 'SummarizedExperiment' | 'matrix' | 'data.frame'
}

export const isRObjectFile = (name: string) => /\.(rdata|rda|rds)$/i.test(name)

/**
 * R that finds the interesting object in whatever was loaded and writes it out
 * as plain CSV. Written to be defensive: an .RData can hold several objects,
 * and nf-core's is not guaranteed to keep the same variable name forever, so we
 * pick by CLASS rather than by name.
 */
const READ_R = String.raw`local({
  path <- "/work/robj.bin"
  isRds <- __IS_RDS__

  obj <- NULL
  if (isRds) {
    obj <- readRDS(path)
  } else {
    e <- new.env()
    nms <- load(path, envir = e)
    # Prefer a DESeqDataSet, then any SummarizedExperiment, then a matrix.
    rank <- function(x) {
      cl <- class(x)
      if (any(cl == "DESeqDataSet")) 3L
      else if (any(cl %in% c("SummarizedExperiment", "RangedSummarizedExperiment"))) 2L
      else if (is.matrix(x) || is.data.frame(x)) 1L else 0L
    }
    best <- 0L
    for (n in nms) {
      v <- get(n, envir = e)
      r <- rank(v)
      if (r > best) { best <- r; obj <- v }
    }
    if (is.null(obj)) stop("No counts-like object found in that .RData")
  }

  cls <- class(obj)
  src <- if (any(cls == "DESeqDataSet")) "DESeqDataSet"
         else if (any(cls %in% c("SummarizedExperiment", "RangedSummarizedExperiment"))) "SummarizedExperiment"
         else if (is.data.frame(obj)) "data.frame" else "matrix"

  design <- ""
  cd <- NULL
  if (src %in% c("DESeqDataSet", "SummarizedExperiment")) {
    suppressMessages(library(SummarizedExperiment))
    # 'counts' is what DESeq2 names it; fall back to the first assay.
    an <- assayNames(obj)
    which <- if ("counts" %in% an) "counts" else an[1]
    m <- as.matrix(assay(obj, which))
    cd <- as.data.frame(colData(obj), stringsAsFactors = FALSE)
    if (src == "DESeqDataSet") {
      d <- tryCatch(DESeq2::design(obj), error = function(e) NULL)
      if (!is.null(d)) design <- paste(deparse(d), collapse = " ")
    }
  } else {
    m <- as.matrix(obj)
  }

  storage.mode(m) <- "double"
  if (is.null(rownames(m))) rownames(m) <- paste0("gene_", seq_len(nrow(m)))
  if (is.null(colnames(m))) colnames(m) <- paste0("sample_", seq_len(ncol(m)))

  write.csv(data.frame(gene_id = rownames(m), as.data.frame(m), check.names = FALSE),
            "/work/robj_counts.csv", row.names = FALSE)

  # colData: drop columns DESeq2 adds for its own bookkeeping, and any column
  # that is constant or unique-per-sample — neither can define a group.
  if (!is.null(cd) && ncol(cd) > 0) {
    drop <- grepl("^(sizeFactor|replaceable)$", colnames(cd))
    cd <- cd[, !drop, drop = FALSE]
    keep <- vapply(cd, function(x) {
      u <- length(unique(as.character(x)))
      u > 1 && u < length(x)
    }, logical(1))
    cd <- cd[, keep, drop = FALSE]
  }
  if (is.null(cd) || ncol(cd) == 0) {
    write.csv(data.frame(sample = colnames(m)), "/work/robj_coldata.csv", row.names = FALSE)
  } else {
    out <- data.frame(sample = rownames(cd), lapply(cd, as.character),
                      stringsAsFactors = FALSE, check.names = FALSE)
    write.csv(out, "/work/robj_coldata.csv", row.names = FALSE)
  }

  sprintf("%s|%d|%s", src, nrow(m), design)
})`

/**
 * Load an .RData / .rds through webR and return counts + sample table.
 * `evalR` is injected so this module does not own the webR lifecycle.
 */
export async function readRObject(
  file: File,
  webR: any,
  onLog: (m: string) => void,
): Promise<RObjectImport> {
  const isRds = /\.rds$/i.test(file.name)
  onLog(`Reading ${file.name} with R…`)

  try { await webR.FS.mkdir('/work') } catch { /* exists */ }
  await webR.FS.writeFile('/work/robj.bin', new Uint8Array(await file.arrayBuffer()))

  const summary: string = await webR.evalRString(READ_R.replace('__IS_RDS__', isRds ? 'TRUE' : 'FALSE'))
  const [source, nGenesStr, designFormula] = summary.split('|')

  const dec = new TextDecoder()
  const countsCsv = dec.decode(await webR.FS.readFile('/work/robj_counts.csv'))
  const coldataCsv = dec.decode(await webR.FS.readFile('/work/robj_coldata.csv'))

  const { samples, colData, columns } = parseColData(coldataCsv)
  const headerSamples = countsCsv.split('\n', 1)[0].split(',').slice(1).map(unquote)

  onLog(`Read a ${source}: ${headerSamples.length} samples, ${nGenesStr} genes.`)
  if (columns.length) onLog(`Sample table columns: ${columns.join(', ')}`)
  if (designFormula) onLog(`Design carried by the object: ${designFormula}`)

  return {
    countsCsv,
    samples: headerSamples.length ? headerSamples : samples,
    colData,
    colDataColumns: columns,
    designFormula: designFormula || null,
    nGenes: parseInt(nGenesStr, 10) || 0,
    source: source as RObjectImport['source'],
  }
}

const unquote = (s: string) => s.trim().replace(/^"|"$/g, '')

/** Minimal CSV split that respects quotes — colData values can contain commas. */
export function splitCsvLine(line: string): string[] {
  const out: string[] = []
  let cur = '', inQ = false
  for (let i = 0; i < line.length; i++) {
    const ch = line[i]
    if (inQ) {
      if (ch === '"') { if (line[i + 1] === '"') { cur += '"'; i++ } else inQ = false }
      else cur += ch
    } else if (ch === '"') inQ = true
    else if (ch === ',') { out.push(cur); cur = '' }
    else cur += ch
  }
  out.push(cur)
  return out.map(s => s.trim())
}

export function parseColData(csv: string) {
  const lines = csv.trim().split(/\r?\n/)
  const header = splitCsvLine(lines[0]).map(unquote)
  const columns = header.slice(1)
  const samples: string[] = []
  const colData: Record<string, Record<string, string>> = {}
  for (const line of lines.slice(1)) {
    const cells = splitCsvLine(line).map(unquote)
    const s = cells[0]
    if (!s) continue
    samples.push(s)
    const rec: Record<string, string> = {}
    columns.forEach((c, i) => { rec[c] = cells[i + 1] ?? '' })
    colData[s] = rec
  }
  return { samples, colData, columns }
}
