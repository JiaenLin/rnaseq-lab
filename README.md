# RNA-seq Lab

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**One job: turn a gene counts matrix into the file
[RNA-seq Studio](https://jiaenlin.github.io/rnaseq-studio/) opens — nothing is uploaded.**

Upload a counts matrix, define two groups, and the lab writes a Studio **bundle** (`.zip`) to
download and drop on the Studio. That is the whole product — you read and plot your results
there, not here.

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

1. **Upload** a CSV/TSV counts matrix (genes × samples; first column = gene, first row = samples).
2. **Design** — assign samples to two groups and pick the control.
3. **Run** — webR loads R + the DE package on first use (cached after), runs the analysis on
   your CPU, and produces a bundle. Your data never leaves the browser.
4. **Explore** — download the `.zip` and open it in RNA-seq Studio.

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
