import { zipSync } from 'fflate'
import type { AnalysisInput, AnalysisResult, Method } from './webr'

export interface BundleParams {
  project: string
  species: string
  method: Method
  /** extra per-sample columns (the recovered factors) to carry into samples.csv */
  covariates?: string[]
  /**
   * The factor whose levels were fitted separately, when the run was blocked.
   *
   * Recorded because it changes what the bundle means and nothing else in it
   * says so: the DEG tables come from one fit per level rather than one fit
   * overall, and `normalized_counts.csv` is CPM rather than median-of-ratios
   * because per-block normalized counts are not comparable between blocks.
   */
  blockFactor?: string
  /** gene_id -> symbol, when the source carried both */
  geneNames?: Map<string, string>
  countsUnitNote?: string
  /**
   * The group every reference level points at — see `referenceGroup`.
   *
   * Passed in rather than guessed from the contrasts, because the contrasts
   * cannot say it. Each one holds the reference of the factor it varies and
   * carries the OTHER factors at whatever level that cell happens to be, so
   * "the denominator of the first pairwise contrast" is the first level of the
   * second factor — a fact about the order the sample names were written in,
   * not about anything the reader chose.
   */
  control?: string
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

/**
 * The group all the chosen references point at.
 *
 * On a 2x2 with genotype = Ctrl and temperature = Thermo that is Ctrl_Thermo,
 * and it is the one thing about the design that no contrast records. Each
 * within-factor contrast pins the reference of the factor it varies and holds
 * the others at whatever level its cell is, so `Ctrl_Cold` and `Ctrl_Thermo`
 * are both denominators of perfectly correct contrasts. Reading `control` off
 * the first of them made it the first level of the second factor — the order
 * the sample names happened to be written in — so picking Thermo as the
 * reference produced a bundle that recorded Ctrl_Cold.
 *
 * Joined with '_' because that is what `withinFactorContrasts` joins with, and
 * its names are the ones matched against the real groups. Checked against
 * `groupLevels` rather than trusted: a label no sample carries would leave the
 * studio's comparison bar seeded with nothing, so an existing group is a better
 * answer than an invented name.
 */
export function referenceGroup(
  design: { factors: readonly { levels: string[] }[]; groupLevels: readonly string[] },
  refs: readonly string[],
): string {
  const groups = design.groupLevels
  if (design.factors.length > 1) {
    const combined = refs.slice(0, design.factors.length).join('_')
    if (groups.includes(combined)) return combined
  }
  if (refs[0] && groups.includes(refs[0])) return refs[0]
  return groups[0] ?? ''
}

/**
 * The denominator pairwise contrasts are built against — a GROUP label.
 *
 * `refs` holds factor LEVELS; groups hold whole labels. They are the same
 * string only when the name has exactly one varying position. Upload one
 * tissue of a `Tissue_age_rep` study and the tissue position is constant, so
 * `detectFactors` rightly drops it and leaves one factor (`008w`…`104w`) while
 * the groups stay `Liver_008w`…`Liver_104w`. Handing `refs[0]` to
 * `pairwiseContrasts` then names a denominator no sample carries: every
 * contrast fails the caller's group-size filter, and the page offers ZERO
 * comparisons with Run disabled and nothing on screen saying why.
 *
 * So map the level back to the group that carries it — which also keeps the
 * reference dropdown meaningful, since picking `104w` must give `Liver_104w`
 * and not merely the first group in the file.
 */
export function referenceGroupFor(
  design: {
    factors: readonly { levels: string[]; values: string[] }[]
    groups: readonly string[]
    groupLevels: readonly string[]
  },
  refs: readonly string[],
  groupLevels: readonly string[] = design.groupLevels,
): string {
  if (design.factors.length === 1) {
    const i = design.factors[0].values.findIndex(v => v === refs[0])
    if (i >= 0 && groupLevels.includes(design.groups[i])) return design.groups[i]
  }
  return referenceGroup({ factors: design.factors, groupLevels }, refs)
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
  const control = params.control
    || ordered.find(c => c.kind === 'pairwise')?.denominator
    || input.groupLevels[0]

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
    counts_unit: params.countsUnitNote ?? (params.method === 'limma' || params.blockFactor
      ? 'CPM (library-size normalized)'
      : 'DESeq2 normalized (median-of-ratios)'),
    /**
     * Not in schema v1; ignored by a reader that does not know it. A studio
     * that does know it can say that a comparison BETWEEN blocks was never
     * fitted, rather than leaving the reader to infer it from the contrast list.
     */
    block_factor: params.blockFactor ?? null,
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
      block: c.block ?? null,
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
