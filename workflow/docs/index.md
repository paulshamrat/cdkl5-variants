# CDKL5 Workflow Docs

This site is the simple local guide to the active workflow under `workflow/tgvr`.

## What This Workflow Does

The workflow is being built to do three practical things:

1. collect and merge variant evidence for a target gene
2. prepare folding-analysis inputs for selected variants
3. support later binding-analysis steps as that stage matures

Today, the most complete part is variant curation.

## Start Here

- `What The Workflow Does`
  Short overview of the workflow and the current stages.
- `Setup`
  How to create the environment and run the tools.
- `Run The Workflow`
  The main commands you actually use.
- `Current Results`
  The current verified `CDKL5` output summary.

## Current Reality

Right now:

- `01_variant_curation` is the main working stage
- `02_folding` prepares DDGun inputs
- `03_binding` is planned but not yet exposed as a public runnable script

## Local Preview

To serve the docs locally:

```bash
mkdocs serve -f workflow/docs/mkdocs.yml
```

Then open `http://127.0.0.1:8000/`.
