# TGVR Overview

TGVR stands for `Thermodynamics-Guided Variant Reclassification`.

It is the active modular development workflow under [workflow/tgvr](/home/paul/cdkl5-variants/workflow/tgvr), intended to gradually replace repeated notebook logic with reusable, stage-oriented code while keeping the original repository root as the scientific baseline.

## Goal

TGVR is designed to:

- build reusable variant-processing stages instead of relying on large monolithic notebooks
- preserve the scientific logic of the original CDKL5 study workflow
- keep legacy root-level files untouched as the baseline
- move forward in a way that can later support proteins beyond `CDKL5`

## Current Maturity

Current status by stage:

- `01_variant_curation`
  The most mature and reproducible TGVR stage today.
- `02_folding`
  Reset to a clean DDGun-first preparation stage. Input preparation is implemented; backend execution is not yet fully integrated.
- `03_binding`
  Conceptually defined and documented, but not currently exposed as a runnable public TGVR script in this tree.

## Working Model

TGVR is organized around per-gene workspaces and explicit stages.

The high-level flow is:

1. initialize a gene-specific TGVR workspace
2. build a master variant dataset
3. prepare downstream folding inputs
4. later integrate folding and binding evidence into broader reclassification analysis

## Layout

```text
workflow/tgvr/
  00_data/
  config/
  outputs/
  scripts/
  src/
```

Key parts:

- `00_data/`
  Stable per-gene workflow inputs such as manual curated variants and canonical reference FASTA files.
- `config/`
  TGVR-native gene configuration.
- `scripts/`
  Main user-facing command entrypoints.
- `outputs/`
  Generated workflow outputs, downloaded stage caches, fetched Palmetto runs, and logs.
- `src/`
  Reusable stage logic.

## CDKL5 As The First Real Target

`CDKL5` is currently the most developed TGVR target.
That work includes:

- a mature `01_variant_curation` path
- a cached `1kgp` backend for stable local reuse
- a Palmetto-backed raw `1kgp` path for the heavier legacy-style reproduction
- a current working kinase-domain default of `1-302`

## Next Read

- Open `Installation` for the environment setup.
- Open `Usage` for common commands.
- Open the stage pages for operational details.
