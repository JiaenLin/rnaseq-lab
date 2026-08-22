import { zipSync } from 'fflate'
import type { AnalysisInput, AnalysisResult, Method } from './webr'

export interface BundleParams {
  project: string
  species: string
  method: Method
  /** extra per-sample columns (the recovered factors) to carry into samples.csv */
  covariates?: string[]
  /** gene_id -> symbol, when the source carried both */
  geneNames?: Map<string, string>
  countsUnitNote?: string
}

const csvCell = (v: string) => (/[",\n]/.test(v) ? `"${v.replace(/"/g, '""')}"` : v)

/**
 * Normalize a matrix to the bundle's `gene_id,gene_name,…` shape.
 *
 * When the source carried a symbol column we put the real symbols back here;
 * otherwise the id is duplicated, which is what the studio expects from a file
 * with no separate symbol. Without the map an nf-core matrix keyed by Ensembl
 * id puts `ENSMUSG00000000001` on every plot axis instead of a gene name.
 */
function reshapeCounts(csv: string, geneNames?: Map<string, string>): string {
  const lines = csv.trim().split(/\r?\n/)
  if (!lines.length) return csv
  const head = lines[0].split(',')
  const second = (head[1] ?? '').replace(/^"|"$/g, '')
  const hasSymbol = /^(gene_name|symbol|name)$/i.test(second)
  if (hasSymbol && !geneNames) return csv

  const header = ['gene_id', 'gene_name', ...head.slice(hasSymbol ? 2 : 1)].join(',')
  const rows = lines.slice(1).map(l => {
    const c = l.split(',')
    const id = (c[0] ?? '').replace(/^"|"$/g, '')
    const sym = geneNames?.get(id) || (hasSymbol ? (c[1] ?? '').replace(/^"|"$/g, '') : id)
    return [c[0], csvCell(sym), ...c.slice(hasSymbol ? 2 : 1)].join(',')
  })
  return [header, ...rows].join('\n') + '\n'
}

/** Put real gene symbols into a DEG table whose first two columns are id,id. */
function reshapeDeg(csv: string, geneNames: Map<string, string>): string {
  const lines = csv.trim().split(/\r?\n/)
  if (lines.length < 2) return csv
  const rows = lines.slice(1).map(l => {
    const c = l.split(',')
    const id = (c[0] ?? '').replace(/^"|"$/g, '')
    const sym = geneNames.get(id)
    return sym ? [c[0], csvCell(sym), ...c.slice(2)].join(',') : l
  })
  return [lines[0], ...rows].join('\n') + '\n'
}

/** Assemble the RNA-seq Studio bundle from an analysis result. */
export function buildBundleFiles(
  input: AnalysisInput, result: AnalysisResult, params: BundleParams,
): Record<string, Uint8Array> {
  const enc = new TextEncoder()
  const covariates = params.covariates ?? []

  // The studio lands on the first contrast and offers the rest in a selector.
  // A pairwise one belongs first: an interaction table read without the
  // pairwise ones beside it is easy to mistake for a plain fold change.
  const ordered = [...result.contrasts].sort(
    (a, b) => (a.kind === b.kind ? 0 : a.kind === 'pairwise' ? -1 : 1))
  const control = ordered.find(c => c.kind === 'pairwise')?.denominator ?? input.groupLevels[0]

  const meta = {
    schema: 1,
    project: params.project || 'RNA-seq analysis',
    species: params.species || 'unknown',
    created: new Date().toISOString().slice(0, 10),
    engine: params.method === 'limma' ? 'webr-limma-voom' : 'webr-deseq2',
    control,
    conditions: input.groupLevels,
    /**
     * What `gene_id` actually holds, which is not always a symbol.
     *
     * This was the constant 'symbol'. It is wrong exactly when `geneNames` is
     * present, because that map EXISTS to carry symbols separately — an
     * nf-core matrix is keyed by Ensembl accession and the accession stays in
     * gene_id while the symbol goes to gene_name. Nothing in the studio reads
     * this field today, so nothing broke; meta.json is the documented contract
     * between the two apps, and a contract that misdescribes its own first
     * column is worth more than nothing being broken yet.
     */
    gene_id_type: params.geneNames?.size ? 'ensembl' : 'symbol',
    counts_unit: params.countsUnitNote ?? (params.method === 'limma'
      ? 'CPM (library-size normalized)'
      : 'DESeq2 normalized (median-of-ratios)'),
    n_samples: input.samples.length,
    contrasts: ordered.map(c => ({
      id: c.id,
      numerator: c.numerator,
      denominator: c.denominator,
      label: c.label,
      deg_file: `deg_${c.id}.csv`,
      n_deg: c.nDeg,
      padj_threshold: 0.05,
      lfc_threshold: 1,
      // Not in schema v1; ignored by a reader that does not know it, and the
      // only way to tell an interaction table from a pairwise one after export.
      kind: c.kind,
    })),
  }

  // samples.csv carries the recovered factors as covariate columns. The studio's
  // SampleRow already allows them, and without them a 2x2 design arrives as four
  // opaque labels with no way to colour a plot by genotype or by temperature.
  const sampleHeader = ['sample', 'condition', ...covariates].join(',')
  const sampleRows = input.samples.map(s =>
    [csvCell(s.sample), csvCell(s.group), ...covariates.map(c => csvCell(s[c] ?? ''))].join(','))

  const files: Record<string, Uint8Array> = {
    'meta.json': enc.encode(JSON.stringify(meta, null, 2) + '\n'),
    'samples.csv': enc.encode([sampleHeader, ...sampleRows].join('\n') + '\n'),
    'normalized_counts.csv': enc.encode(reshapeCounts(result.normCsv, params.geneNames)),
    // The raw matrix, so the studio can run DESeq2 on pairs this run did not
    // export — DESeq2 models raw counts, not normalized ones.
    'raw_counts.csv': enc.encode(reshapeCounts(input.countsCsv, params.geneNames)),
  }
  for (const c of ordered) {
    files[`deg_${c.id}.csv`] = enc.encode(
      params.geneNames ? reshapeDeg(c.degCsv, params.geneNames) : c.degCsv)
  }
  return files
}

export function zipBundle(files: Record<string, Uint8Array>): Uint8Array {
  return zipSync(files, { level: 6 })
}
