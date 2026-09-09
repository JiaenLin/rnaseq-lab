// The long-read (isoform) layer.
//
// WHICH FILES. This module is written against Oxford Nanopore
// `epi2me-labs/wf-transcriptomes` output, which is what the cluster produces.
// Three files, all under the run's `out/` directory:
//
//   cohort/transcript_counts.tsv       REQUIRED. Counts per transcript per sample,
//                                      with 13 leading annotation columns.
//   cohort/transcript_metadata.tsv     optional. The same 13 annotation columns
//                                      without the counts. Only needed if you
//                                      bring counts from somewhere else.
//   cohort/sqanti/cohort_classification.txt
//                                      optional. SQANTI3's structural category
//                                      per transcript, plus ORF/NMD when the run
//                                      had `sqanti_skip_orf: false`.
//
// `cohort/gene_counts.tsv` is NOT needed: the gene matrix is summed from the
// transcript matrix here, so both levels come from one fit on one file and
// cannot disagree about what a gene contains.
//
// A PacBio/IsoSeq or StringTie matrix works too as long as the first column is a
// transcript id and some column names the gene; everything below degrades to
// "annotated transcript, no category" rather than refusing.

import type { ParsedMatrix } from './matrix.ts'

/** One transcript's annotation, as the bundle's `transcripts.csv` row. */
export interface TranscriptAnn {
  transcript_id: string
  gene_id: string
  gene_name: string
  transcript_name: string
  /** What a reader should SEE. See `displayName`. */
  display_name: string
  /** SQANTI3 structural category, '' when unclassified. */
  structural_category: string
  novel: boolean
}

/**
 * Which vocabulary `structural_category` is written in.
 *
 * Two incompatible ones reach that column. bambu's own `txClassDescription`
 * says `newWithin`, `newLastJunction:newJunction:newLastExon`, `newGene-spliced`;
 * SQANTI3 says `full-splice_match`, `novel_in_catalog`, `novel_not_in_catalog`.
 * A reader colouring by FSM/ISM/NIC/NNC maps the first set to nothing at all,
 * and — worse — `applySqanti` leaves unmatched rows alone, so a partial join
 * produces a column that is half each. The bundle records which, and 'mixed'
 * is a real answer rather than an embarrassment.
 */
export type CategoryVocabulary = 'bambu' | 'sqanti' | 'mixed' | 'none'

export interface LongReadInput {
  /** transcript_id -> annotation, in matrix row order. */
  transcripts: TranscriptAnn[]
  /** `gene_id,<sample>…` raw counts, summed within gene from the transcript matrix. */
  geneCountsCsv: string
  /**
   * `transcript_id,<sample>…` raw counts, ROUNDED, keyed by transcript id.
   *
   * Not `parsed.countsCsv`: that one is headed `gene_id` (matrix.ts writes one
   * header for every matrix) and still fractional. Shipping it would put a file
   * called transcript_counts.csv into a bundle claiming a `gene_id` column, and
   * numbers no engine actually used.
   */
  txCountsCsv: string
  /** gene_id -> symbol, for the gene layer's own annotation. */
  geneNames: Map<string, string>
  /** transcript_id -> gene_id, for the DTU offset. */
  txToGene: Map<string, string>
  nGenes: number
  /** See `CategoryVocabulary`. Set to 'bambu' or 'none' on read; `applySqanti` updates it. */
  vocabulary: CategoryVocabulary
}

const norm = (s: string) => s.trim().toLowerCase().replace(/[_.]/g, '')

const csvCell = (v: string) => (/[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v)

/**
 * Find a column by any of several accepted names. -1 when absent.
 *
 * `skip` is the key column, and it is not optional politeness. In
 * wf-transcriptomes the key is `TXNAME`, which normalises to `txname` — the
 * same string as the alias `tx_name` for the transcript's readable NAME. Without
 * the skip, `findIndex` returns column 0 and every transcript's display name
 * becomes its own accession: 3,999 isoforms, none of them readable, and no error
 * anywhere. Measured on this cohort before the skip existed.
 */
const col = (header: string[], skip: number, ...names: string[]) => {
  const want = new Set(names.map(norm))
  return header.findIndex((h, i) => i !== skip && want.has(norm(h)))
}

/**
 * Is this a transcript-level matrix rather than a gene-level one?
 *
 * Decided on the KEY COLUMN'S NAME, not on the ids in it. A bambu cohort matrix
 * is keyed by `TXNAME` and holds a mixture of `ENSMUST…` and `BambuTx…`; a
 * StringTie one is keyed by `transcript_id`. Sniffing the ids instead would call
 * a gene matrix of `ENSMUSG…` transcript-level the moment someone renamed a
 * column, and the whole isoform layer would be built on genes.
 */
export function isTranscriptMatrix(header: string[]): boolean {
  const key = norm(header[0] ?? '')
  return key === 'txname' || key === 'transcriptid' || key === 'txid'
    || key === 'transcript' || key === 'isoform'
}

/**
 * The name a reader sees, in three tiers.
 *
 *   1. Annotated            -> the reference transcript name, `Nppb-201`.
 *   2. Novel, known gene    -> `<symbol>-novel-<n>`, n from the caller's counter.
 *   3. Novel gene           -> the transcript id itself; there is no better name,
 *                              and inventing a locus string from data we do not
 *                              have here would be worse than the accession.
 *
 * The accession is never lost — it stays the key of every table — exactly as
 * src/lib/symbols.ts does one level up. NOTHING IS MERGED: two models may
 * legitimately land on the same display name, and both keep their own row.
 */
export function displayName(
  transcript_id: string, gene_name: string, transcript_name: string, novelIndex: number,
): string {
  const tn = transcript_name.trim()
  if (tn && tn !== 'NA') return tn
  const gn = gene_name.trim()
  if (gn && gn !== 'NA') return `${gn}-novel-${novelIndex}`
  return transcript_id
}

const clean = (v: string | undefined) => {
  const t = (v ?? '').trim()
  return t === 'NA' || t === 'NaN' || t === 'null' ? '' : t
}

/**
 * Read the isoform layer out of a transcript matrix.
 *
 * `parsed` has already separated annotation columns from samples, but it keeps
 * only the key and one symbol column — every other annotation column is
 * discarded there and is exactly what this layer needs. So the raw rows are
 * passed in again rather than re-parsed.
 */
export function readLongRead(parsed: ParsedMatrix, rows: string[][]): LongReadInput {
  const header = rows[0].map(h => String(h ?? '').trim())
  const body = rows.slice(1).filter(r => r.some(c => String(c ?? '').trim() !== ''))

  const iTx = 0
  const iGene = col(header, iTx, 'GENEID', 'gene_id', 'associated_gene')
  const iGeneName = col(header, iTx, 'gene_name', 'gene_symbol', 'symbol')
  const iTxName = col(header, iTx, 'transcript_name')
  const iNovelTx = col(header, iTx, 'novelTranscript')
  const iCategory = col(header, iTx, 'structural_category', 'txClassDescription')

  if (iGene < 0) {
    throw new Error(
      'A transcript matrix needs a gene column so transcripts can be grouped into genes. ' +
      'Expected GENEID (wf-transcriptomes) or gene_id. Columns read as annotation: ' +
      `${parsed.annotationColumns.join(', ') || 'none'}.`)
  }

  // Novel counters are per gene and follow FILE ORDER, which bambu writes in
  // genomic order within a gene. Stable across rebuilds of the same run, which
  // is what makes `Nppb-novel-1` safe to quote in a figure.
  const novelSeen = new Map<string, number>()
  const transcripts: TranscriptAnn[] = []
  const txToGene = new Map<string, string>()
  const geneNames = new Map<string, string>()

  for (const r of body) {
    const transcript_id = clean(r[iTx]) || 'unknown'
    const gene_id = clean(r[iGene]) || transcript_id
    const gene_name = iGeneName >= 0 ? clean(r[iGeneName]) : ''
    const transcript_name = iTxName >= 0 ? clean(r[iTxName]) : ''
    const category = iCategory >= 0 ? clean(r[iCategory]) : ''
    // `novelTranscript` is bambu's own flag and is authoritative when present.
    // Otherwise a missing reference transcript name is the tell.
    const novel = iNovelTx >= 0
      ? /^(true|t|yes|1)$/i.test(clean(r[iNovelTx]))
      : !transcript_name

    let n = 0
    if (novel) {
      n = (novelSeen.get(gene_id) ?? 0) + 1
      novelSeen.set(gene_id, n)
    }
    transcripts.push({
      transcript_id, gene_id, gene_name, transcript_name,
      display_name: displayName(transcript_id, gene_name, transcript_name, n),
      structural_category: category,
      novel,
    })
    txToGene.set(transcript_id, gene_id)
    if (gene_name && !geneNames.has(gene_id)) geneNames.set(gene_id, gene_name)
  }

  return {
    transcripts,
    ...sumToGenes(parsed, transcripts),
    geneNames,
    txToGene,
    vocabulary: iCategory >= 0 && transcripts.some(t => t.structural_category) ? 'bambu' : 'none',
  }
}

/**
 * Sum a transcript matrix into a gene matrix.
 *
 * Deliberately NOT read from `cohort/gene_counts.tsv`, even though that file
 * exists beside it. bambu writes the two independently; if the gene layer came
 * from one file and the transcript layer from another, a DTU result could
 * disagree with the gene-level fold change it sits beside and there would be no
 * way to tell which was right. One matrix, one derivation, one answer.
 *
 * Counts are rounded: bambu distributes reads across transcripts by expectation
 * and writes fractional values (36.69048), and DESeq2 requires integers. The
 * rounding happens ONCE, here, on the gene total — not per transcript and then
 * summed, which would accumulate the error.
 */
function sumToGenes(
  parsed: ParsedMatrix, transcripts: TranscriptAnn[],
): { geneCountsCsv: string; txCountsCsv: string; nGenes: number } {
  const lines = parsed.countsCsv.trim().split(/\r?\n/)
  const nSamples = parsed.samples.length
  if (lines.length - 1 !== transcripts.length) {
    // The two layers are tied together HERE and nowhere else: row r of the
    // counts is transcript r. If those ever disagree, every transcript's counts
    // land on the wrong gene and nothing downstream can notice.
    throw new Error(
      `Internal: ${lines.length - 1} count rows but ${transcripts.length} transcripts. ` +
      `Please report this — the isoform layer would have been misaligned.`)
  }

  const totals = new Map<string, Float64Array>()
  const order: string[] = []
  const txOut = ['transcript_id,' + parsed.samples.map(csvCell).join(',')]

  for (let r = 1; r < lines.length; r++) {
    // parsed.countsCsv is written by matrix.ts, which QUOTES a key containing a
    // comma. A naive split then shifts every sample one column to the left, on
    // that row only, silently. splitCsv honours the quoting its writer applied.
    const cells = splitCsv(lines[r])
    const t = transcripts[r - 1]
    let acc = totals.get(t.gene_id)
    if (!acc) { acc = new Float64Array(nSamples); totals.set(t.gene_id, acc); order.push(t.gene_id) }
    // ROUNDED ONCE, HERE. bambu writes expected counts (36.69048) and both
    // engines model integers. Rounding per transcript and then summing, and
    // rounding the gene sum, give different gene totals; doing it once at the
    // transcript level means the gene matrix, the transcript matrix shipped in
    // the bundle, and the matrix the usage test is handed are all sums of the SAME
    // integers, and a reader can reproduce one from the other.
    const rounded = new Array<number>(nSamples)
    for (let sIdx = 0; sIdx < nSamples; sIdx++) {
      const v = Number(cells[sIdx + 1])
      const n = Number.isFinite(v) ? Math.round(v) : 0
      rounded[sIdx] = n
      acc[sIdx] += n
    }
    txOut.push(csvCell(t.transcript_id) + ',' + rounded.join(','))
  }

  const out = ['gene_id,' + parsed.samples.map(csvCell).join(',')]
  for (const gene of order) {
    out.push(csvCell(gene) + ',' + Array.from(totals.get(gene)!, v => String(v)).join(','))
  }
  return {
    geneCountsCsv: out.join('\n') + '\n',
    txCountsCsv: txOut.join('\n') + '\n',
    nGenes: order.length,
  }
}

/** Split one CSV line, honouring the double-quoting matrix.ts's writer applies. */
function splitCsv(line: string): string[] {
  if (!line.includes('"')) return line.split(',')
  const out: string[] = []
  let cur = '', q = false
  for (let i = 0; i < line.length; i++) {
    const c = line[i]
    if (q) {
      if (c === '"') { if (line[i + 1] === '"') { cur += '"'; i++ } else q = false }
      else cur += c
    } else if (c === '"') q = true
    else if (c === ',') { out.push(cur); cur = '' }
    else cur += c
  }
  out.push(cur)
  return out
}

/**
 * Join SQANTI3's classification onto the transcripts already read.
 *
 * `cohort_classification.txt` is tab-separated with `isoform` as its first
 * column and `structural_category` among the rest; with ORF prediction on it
 * also carries `coding`, `ORF_length` and `predicted_NMD`. Rows it does not
 * mention are left alone rather than dropped — SQANTI classifies the models the
 * run gave it, and a mismatch means the wrong run's file, which is worth saying.
 */
export function applySqanti(
  transcripts: TranscriptAnn[], rows: string[][],
): { matched: number; total: number; vocabulary: CategoryVocabulary } {
  const header = rows[0].map(h => String(h ?? '').trim())
  const iId = col(header, -1, 'isoform', 'transcript_id')
  const iCat = col(header, -1, 'structural_category')
  if (iId < 0 || iCat < 0) {
    throw new Error(
      'That does not look like a SQANTI3 classification table: expected an `isoform` ' +
      'column and a `structural_category` column. wf-transcriptomes writes it to ' +
      'out/cohort/sqanti/cohort_classification.txt.')
  }
  const cat = new Map<string, string>()
  for (const r of rows.slice(1)) {
    const id = clean(r[iId])
    if (id) cat.set(id, clean(r[iCat]))
  }
  let matched = 0
  for (const t of transcripts) {
    const c = cat.get(t.transcript_id)
    if (c) { t.structural_category = c; matched++ }
  }
  // A join that reached some rows and not others leaves the column speaking two
  // languages. Say so rather than let the bundle claim it is all SQANTI.
  const vocabulary: CategoryVocabulary =
    matched === 0 ? 'none'
      : transcripts.some(t => t.structural_category && !cat.has(t.transcript_id)) ? 'mixed'
        : 'sqanti'
  return { matched, total: transcripts.length, vocabulary }
}

/** The bundle's `transcripts.csv`. */
export function transcriptsCsv(transcripts: readonly TranscriptAnn[]): string {
  const head = 'transcript_id,gene_id,gene_name,transcript_name,display_name,structural_category,novel'
  const rows = transcripts.map(t => [
    t.transcript_id, t.gene_id, t.gene_name, t.transcript_name,
    t.display_name, t.structural_category, t.novel ? 'TRUE' : 'FALSE',
  ].map(csvCell).join(','))
  return [head, ...rows].join('\n') + '\n'
}

/** What the UI states before a run, so nobody discovers the shape afterwards. */
const VOCAB_LABEL: Record<CategoryVocabulary, string> = {
  bambu: 'bambu classes', sqanti: 'SQANTI categories',
  mixed: 'mixed bambu/SQANTI classes', none: '',
}

export function describe(input: LongReadInput): string {
  const n = input.transcripts.length
  const novel = input.transcripts.filter(t => t.novel).length
  const named = input.transcripts.filter(t => t.display_name !== t.transcript_id).length
  const cats = new Set(input.transcripts.map(t => t.structural_category).filter(Boolean))
  return `${n.toLocaleString()} transcripts over ${input.nGenes.toLocaleString()} genes`
    + ` · ${novel.toLocaleString()} novel`
    + ` · ${named.toLocaleString()} with a readable name`
    + (cats.size ? ` · ${cats.size} ${VOCAB_LABEL[input.vocabulary] || 'categories'}`
      : ' · no structural categories')
}
