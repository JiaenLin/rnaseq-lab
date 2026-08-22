# RNA-seq Lab

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**One job: turn bulk RNA-seq counts into the file
[RNA-seq Studio](https://jiaenlin.github.io/rnaseq-studio/) opens — nothing is uploaded.**

Upload a counts matrix **or an nf-core/rnaseq DESeq2 object**, confirm the design, and the lab
writes a Studio **bundle** (`.zip`) to download and drop on the Studio. That is the whole
product — you read and plot your results there, not here.

Building that bundle means running the differential expression, because the Studio's format
carries `deg_<contrast>.csv`. So **limma-voom** or **DESeq2** runs entirely client-side via
[webR](https://docs.r-wasm.org/webr/latest/) (R 4.6.0 compiled to WebAssembly) to produce it.
The DE run is the mechanism, not a second product.

### The family — one job each

| App | Takes | Produces |
|---|---|---|
| **rnaseq-service** | raw FASTQ | an analysis request + nf-core sample sheet |
| **rnaseq-lab** (here) | a bulk counts matrix | `bundle.zip` for rnaseq-studio |
| **rnaseq-studio** | `bundle.zip` | the figures you read |
| **scrnaseq-lab** | an annotated `.h5ad` / `.rds` | `bundle.zip` for scrnaseq-studio |
| **scrnaseq-studio** | `bundle.zip` | the figures you read |

👉 **[Open the app](https://jiaenlin.github.io/rnaseq-lab/)**

## How it works

1. **Upload** a CSV/TSV counts matrix, or an R object straight from nf-core/rnaseq.
2. **Design** — the lab reads the factors out of your sample names (or out of the object's own
   sample table) and offers the comparisons that design actually asks.
3. **Run** — webR loads R + the DE package on first use (cached after), runs the analysis on
   your CPU, and produces a bundle. Your data never leaves the browser.
4. **Explore** — download the `.zip` and open it in RNA-seq Studio.

## What it reads

| File | Notes |
|---|---|
| any CSV/TSV counts matrix | genes as rows, samples as columns |
| `salmon.merged.gene_counts.tsv` | nf-core's matrix. Its `gene_id` **and** `gene_name` columns are recognised as annotation — the naive "column 0 is the gene, the rest are samples" rule turns `gene_name` into a 20th sample full of gene symbols |
| `deseq2.dds.RData` | nf-core's DESeq2 object: counts **and** the sample table **and** the design, so nothing has to be guessed from sample names |
| `*.SummarizedExperiment.rds` | gene- or transcript-level |

R objects are opened by R itself, in webR — the same engine that runs the DE. No JavaScript
RData parser is involved.

## Designs it understands

A one-factor study (`WT` vs `KO`) is the easy case. The one that needs help is factorial:

```
KO_Cold_1 …  KO_Thermo_1 …  Ctrl_Cold_1 …  Ctrl_Thermo_1 …
```

That is 2×2, and it asks **five** questions, not one:

| | |
|---|---|
| `KO_Thermo vs Ctrl_Thermo` | knockout effect at thermoneutrality |
| `KO_Cold vs Ctrl_Cold` | knockout effect under cold |
| `Ctrl_Cold vs Ctrl_Thermo` | cold response in controls |
| `KO_Cold vs KO_Thermo` | cold response in knockouts |
| **interaction** | does the knockout effect *depend on* temperature? |

The lab detects the factors, names them, lets you set each one's reference level, and exports
every comparison you tick as its own `deg_*.csv`. The recovered factors also travel into
`samples.csv` as covariate columns, so the Studio can colour by genotype or by temperature
rather than by four opaque group labels.

**Two things it will not do**, both deliberate:

- It never offers a comparison that moves two factors at once (`KO+Cold` vs `Ctrl+Thermo`).
  A hit there cannot be attributed to either factor.
- It withholds the interaction when any cell of the grid is empty. The coefficient still
  exists, but it is not estimated from the cells the label implies.

All contrasts come from **one model fit**. Under a cell-means design (`~ 0 + group`) every
comparison — interaction included — is a linear combination of group means, so there is no
second, differently-parameterised fit to disagree with the first.

## Methods

- **limma-voom** — no extra setup; works out of the box in webR.
- **DESeq2** (+ apeglm) — needs one CRAN dependency (`locfit`) that has no public WASM binary,
  so we build it ourselves. The binary is served from `/<repo>/wasm/` and added to webR's
  install repos.

## Develop

```bash
npm install
npm run dev      # dev server sets COOP/COEP so SharedArrayBuffer works locally
npm run build && npm run preview
```

Stack: React + TypeScript + Vite + Tailwind; webR from CDN; zip via `fflate`. No backend.
Deploys to GitHub Pages on push to `main`.

## Related

- **[RNA-seq Studio](https://jiaenlin.github.io/rnaseq-studio/)**
  ([source](https://github.com/JiaenLin/rnaseq-studio)) — explore the bundle this app produces.
- **[RNA-seq Service](https://jiaenlin.github.io/rnaseq-service/)**
  ([source](https://github.com/JiaenLin/rnaseq-service)) — start from raw FASTQ files instead of
  a count matrix: scan your sequencing folder and generate an analysis request.
