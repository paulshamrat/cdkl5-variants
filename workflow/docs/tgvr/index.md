# What The Workflow Does

TGVR stands for `Thermodynamics-Guided Variant Reclassification`.

It is the active workflow under `workflow/tgvr`.

## Main Job

The workflow is meant to take a gene such as `CDKL5` and:

1. gather variant evidence from sources such as `ClinVar`, `1kgp`, and `gnomAD`
2. merge those inputs into one filtered master dataset
3. prepare downstream inputs for folding analysis
4. later support binding-analysis workflows

## Current Stages

### 1. Variant Curation

This is the main working stage today.

It:

- builds the master variant dataset
- supports manual curated variants
- writes summary tables for the full protein and selected residue ranges

### 2. Folding

This stage currently prepares `DDGun` inputs from the curated variant set.

It does not yet provide the full end-to-end folding backend.

### 3. Binding

This stage is planned, but it is not yet exposed as a public runnable script in the current tree.

## Where Things Live

```text
workflow/tgvr/
  00_data/    stable inputs
  config/     gene configuration
  outputs/    generated results
  scripts/    commands you run
  src/        implementation code
```

## Practical Takeaway

If you only need the part that works best right now, focus on:

- `Setup`
- `Run The Workflow`
- `Current Results`
