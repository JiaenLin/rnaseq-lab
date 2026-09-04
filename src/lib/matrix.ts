// Reading a counts matrix, and telling gene annotation apart from samples.
//
// The naive rule — "column 0 is the gene, everything else is a sample" — breaks
// on the most common file anyone will bring here. nf-core/rnaseq writes
// `salmon.merged.gene_counts.tsv` with TWO leading annotation columns:
//
//     gene_id            gene_name   Ctrl_Cold_1   Ctrl_Cold_2   …
//     ENSMUSG00000000001 Gnai3       1234.5        1102.0        …
//
// Under that rule `gene_name` becomes a sample, so a 19-sample delivery reads as
// 20 samples with one nonsense column full of gene symbols. It is visible in the
// UI only as an off-by-one in a count nobody double-checks.
//
// So a column is a sample if its VALUES look numeric — not because of where it
// sits or what it is called. Names are used only to break ties.

export interface ParsedMatrix {
  geneIds: string[]
  geneNames: string[] | null      // null when the file carries no symbol column
  samples: string[]
  /** CSV in the canonical `gene_id,<sample>…` shape R reads. */
  countsCsv: string
  /** Annotation columns that were recognised and set aside. */
  annotationColumns: string[]
  nGenes: number
  /**
   * A thinned copy of the counts, row-major genes x samples.
   *
   * Every row is walked here anyway, so keeping every Nth one costs almost
   * nothing and is the only numeric view the rest of the app gets — the counts
   * otherwise exist solely as a CSV string on its way to R. `lib/blocking.ts`
   * reads it to decide whether a factor's levels are far enough apart to
   * deserve separate fits, which is a question about the DATA and cannot be
   * answered from sample names.
   */
  probe: Probe
}

/** Thinned counts, row-major: `values[g * nSamples + s]`. */
export interface Probe {
  values: Float64Array
  nGenes: number
  nSamples: number
}

/** Rows kept for the probe. 2,000 x 275 is 4 MB and plenty to compare levels. */
const PROBE_ROWS = 2000

/** Column names that are annotation even if they somehow parse as numeric. */
const ANNOTATION_RX =
  /^(gene[_.]?id|gene[_.]?name|gene[_.]?symbol|symbol|name|transcript[_.]?id|tx[_.]?id|id|entrez([_.]?(gene|id))?|ensembl([_.]?id)?|refseq|description|biotype|gene[_.]?biotype|chr|chromosome|start|end|strand|length|gene[_.]?length|width|locus)$/i

/** A gene symbol column is the one we want to keep as `gene_name`. */
const SYMBOL_RX = /^(gene[_.]?name|gene[_.]?symbol|symbol|name)$/i

const isNumericish = (v: string) => {
  const t = (v ?? '').trim()
  if (t === '' || t === 'NA' || t === 'NaN' || t === 'null') return true   // missing, not text
  return Number.isFinite(Number(t))
}

/**
 * Decide which columns hold counts.
 *
 * Sampling 200 rows rather than all of them keeps a 57k-gene matrix responsive;
 * a column of gene symbols announces itself in the first handful of rows.
 */
function classify(header: string[], rows: string[][]): { sampleIdx: number[]; annIdx: number[] } {
  const probe = rows.slice(0, 200)
  const sampleIdx: number[] = []
  const annIdx: number[] = []
  for (let c = 0; c < header.length; c++) {
    const name = (header[c] ?? '').trim()
    const values = probe.map(r => r[c] ?? '')
    const numeric = values.length > 0 && values.every(isNumericish)
    // Name wins over content: an all-numeric `entrez_id` column is still an id.
    if (ANNOTATION_RX.test(name) || !numeric) annIdx.push(c)
    else sampleIdx.push(c)
  }
  return { sampleIdx, annIdx }
}

const csvCell = (v: string) => (/[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v)

/**
 * Parse an already-split table (rows of cells) into the canonical shape.
 * Throws with a readable message rather than producing a subtly wrong matrix.
 */
export function parseMatrix(rows: string[][]): ParsedMatrix {
  if (!rows.length) throw new Error('The file is empty.')
  const header = rows[0].map(h => String(h ?? '').trim())
  const body = rows.slice(1).filter(r => r.some(c => String(c ?? '').trim() !== ''))
  if (!body.length) throw new Error('The file has a header but no data rows.')

  const { sampleIdx, annIdx } = classify(header, body)

  /**
   * A sample name carrying a comma, quote or newline breaks the CSV contract.
   *
   * `csvCell` below quotes it correctly and R reads it back correctly, so this
   * looks safe and is not: every OTHER reader of the bundle splits on commas
   * — this module's own `reshapeCounts`, and the studio's counts parser — and
   * a quoted comma silently shifts every column right of it by one. A header of
   * `gene_id,"ctrl,rep1",ctrl2` is seven fields over rows of six, and nothing
   * downstream notices.
   *
   * Refused at the door rather than escaped, because the alternative is making
   * two parsers in two repositories quote-aware and staying that way.
   */
  const badName = sampleIdx.map(i => header[i] ?? '').filter(n => /[",\n\r]/.test(n))
  if (badName.length) {
    throw new Error(
      `Sample name(s) contain a comma, quote or line break, which a CSV column header cannot ` +
      `carry safely: ${badName.slice(0, 3).map(n => JSON.stringify(n)).join(', ')}` +
      `${badName.length > 3 ? ` and ${badName.length - 3} more` : ''}. Rename the columns and re-upload.`)
  }

  /**
   * Duplicate sample names silently DELETE a library.
   *
   * R selects columns by name — `counts[, cd$sample]` — so two columns called
   * `s1` both resolve to the first one, and the second library's counts are
   * replaced by the first's with no error at any layer. Verified in R: a matrix
   * whose second column reads 999/888 comes back reading 10/30.
   *
   * Duplicate GENE ids fail loudly on their own (`read.csv(row.names = 1)`
   * refuses them), but the message names neither the file nor the offender, so
   * they are caught here too.
   */
  const dupOf = (xs: string[]) => {
    const seen = new Set<string>(), dup = new Set<string>()
    for (const x of xs) (seen.has(x) ? dup : seen).add(x)
    return [...dup]
  }
  const dupSample = dupOf(sampleIdx.map(i => header[i] ?? ''))
  if (dupSample.length) {
    throw new Error(
      `Duplicate sample column name(s): ${dupSample.slice(0, 5).join(', ')}` +
      `${dupSample.length > 5 ? ` and ${dupSample.length - 5} more` : ''}. ` +
      `Every column would be read as the first one of its name, so a library would be ` +
      `silently replaced. Give each sample a unique name.`)
  }

  if (sampleIdx.length < 2) {
    throw new Error(
      `Found ${sampleIdx.length} numeric column(s); a counts matrix needs at least 2 samples. ` +
      `Columns read as annotation: ${annIdx.map(i => header[i] || `#${i + 1}`).join(', ') || 'none'}.`)
  }

  // The gene id is the first annotation column; the symbol column, if any, is
  // whichever annotation column is named like a symbol.
  const idIdx = annIdx.length ? annIdx[0] : -1
  const symIdx = annIdx.find(i => SYMBOL_RX.test(header[i] ?? '')) ?? -1
  // `gene_id` + `gene_name` both present: keep the id as the key, not the symbol.
  const keyIdx = idIdx === symIdx && annIdx.length > 1 ? annIdx[1] : idIdx

  const geneIds = body.map((r, i) => (keyIdx >= 0 ? String(r[keyIdx] ?? '').trim() : '') || `row_${i + 1}`)
  const geneNames = symIdx >= 0 && symIdx !== keyIdx
    ? body.map(r => String(r[symIdx] ?? '').trim())
    : null

  const samples = sampleIdx.map(i => header[i] || `sample_${i + 1}`)
  const out: string[] = ['gene_id,' + samples.map(csvCell).join(',')]
  for (let r = 0; r < body.length; r++) {
    const row = body[r]
    out.push(csvCell(geneIds[r]) + ',' + sampleIdx.map(i => {
      const v = String(row[i] ?? '').trim()
      return v === '' || v === 'NA' || v === 'NaN' ? '0' : v
    }).join(','))
  }

  const dupGene = dupOf(geneIds)
  if (dupGene.length) {
    throw new Error(
      `Duplicate gene id(s): ${dupGene.slice(0, 5).join(', ')}` +
      `${dupGene.length > 5 ? ` and ${dupGene.length - 5} more` : ''}. ` +
      `R refuses a matrix with repeated row names, so the run would fail after the upload.`)
  }

  // Thinned evenly across the file rather than from the top: the first rows of
  // a sorted annotation are not a sample of it.
  const stride = Math.max(1, Math.ceil(body.length / PROBE_ROWS))
  const keptRows: number[] = []
  for (let r = 0; r < body.length; r += stride) keptRows.push(r)
  const probe: Probe = {
    values: new Float64Array(keptRows.length * sampleIdx.length),
    nGenes: keptRows.length,
    nSamples: sampleIdx.length,
  }
  for (let k = 0; k < keptRows.length; k++) {
    const row = body[keptRows[k]]
    for (let c = 0; c < sampleIdx.length; c++) {
      const v = Number(String(row[sampleIdx[c]] ?? '').trim())
      probe.values[k * sampleIdx.length + c] = Number.isFinite(v) ? v : 0
    }
  }

  return {
    geneIds,
    geneNames,
    samples,
    countsCsv: out.join('\n') + '\n',
    annotationColumns: annIdx.map(i => header[i] || `#${i + 1}`),
    nGenes: body.length,
    probe,
  }
}

/**
 * A probe from an already-canonical `gene_id,<sample>…` CSV.
 *
 * The nf-core path never touches `parseMatrix` — R opens the object and hands
 * back a CSV directly — so without this, blocking would simply never be
 * suggested for the one input format that most often carries many tissues.
 * Values are unquoted by position, which is safe here because this CSV was
 * written by us or by R's write.csv over numeric columns.
 */
export function probeFromCsv(csv: string): Probe {
  const lines = csv.trim().split(/\r?\n/)
  const nSamples = Math.max(0, (lines[0]?.split(',').length ?? 1) - 1)
  const body = lines.length - 1
  if (nSamples < 1 || body < 1) return { values: new Float64Array(0), nGenes: 0, nSamples: 0 }
  const stride = Math.max(1, Math.ceil(body / PROBE_ROWS))
  const rows: number[] = []
  for (let r = 0; r < body; r += stride) rows.push(r)
  const values = new Float64Array(rows.length * nSamples)
  for (let k = 0; k < rows.length; k++) {
    const cells = lines[rows[k] + 1].split(',')
    for (let c = 0; c < nSamples; c++) {
      const v = Number(cells[c + 1])
      values[k * nSamples + c] = Number.isFinite(v) ? v : 0
    }
  }
  return { values, nGenes: rows.length, nSamples }
}
