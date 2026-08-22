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
}

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

  return {
    geneIds,
    geneNames,
    samples,
    countsCsv: out.join('\n') + '\n',
    annotationColumns: annIdx.map(i => header[i] || `#${i + 1}`),
    nGenes: body.length,
  }
}
