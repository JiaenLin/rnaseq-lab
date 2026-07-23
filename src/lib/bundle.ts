import { zipSync } from 'fflate'
import type { AnalysisInput, AnalysisResult, Method } from './webr'

export interface BundleParams {
  project: string
  species: string
  method: Method
}

// Assemble the RNA-seq Studio bundle files from an analysis result.
export function buildBundleFiles(
  input: AnalysisInput, result: AnalysisResult, params: BundleParams,
): Record<string, Uint8Array> {
  const enc = new TextEncoder()
  const contrastId = `${result.numerator}_vs_${result.denominator}`
  const conditions = Array.from(new Set([input.control, ...input.samples.map(s => s.condition)]))
  const meta = {
    schema: 1,
    project: params.project || 'RNA-seq analysis',
    species: params.species || 'unknown',
    created: new Date().toISOString().slice(0, 10),
    engine: params.method === 'limma' ? 'webr-limma-voom' : 'webr-deseq2',
    control: input.control,
    conditions,
    gene_id_type: 'symbol',
    counts_unit: params.method === 'limma'
      ? 'CPM (library-size normalized)'
      : 'DESeq2 normalized (median-of-ratios)',
    n_samples: input.samples.length,
    contrasts: [{
      id: contrastId,
      numerator: result.numerator,
      denominator: result.denominator,
      label: `${result.numerator} vs ${result.denominator}`,
      deg_file: `deg_${contrastId}.csv`,
      n_deg: result.nDeg,
      padj_threshold: 0.05,
      lfc_threshold: 1,
    }],
  }
  const samplesCsv = 'sample,condition\n'
    + input.samples.map(s => `${s.sample},${s.condition}`).join('\n') + '\n'
  return {
    'meta.json': enc.encode(JSON.stringify(meta, null, 2) + '\n'),
    'samples.csv': enc.encode(samplesCsv),
    'normalized_counts.csv': enc.encode(result.normCsv),
    [`deg_${contrastId}.csv`]: enc.encode(result.degCsv),
  }
}

export function zipBundle(files: Record<string, Uint8Array>): Uint8Array {
  return zipSync(files, { level: 6 })
}
