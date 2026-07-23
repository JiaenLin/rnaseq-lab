# RNA-seq Lab

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Differential expression in your browser — nothing is uploaded.**

Upload a counts matrix, define two groups, and run **limma-voom** or **DESeq2** entirely
client-side via [webR](https://docs.r-wasm.org/webr/latest/) (R 4.6.0 compiled to
WebAssembly). The result is an RNA-seq Studio **bundle** you can download and explore in the
[RNA-seq Studio](https://jiaenlin.github.io/rnaseq-studio/) viewer.

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
