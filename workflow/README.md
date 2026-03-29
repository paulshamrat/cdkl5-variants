# Workflow

This folder now has two clearly separated active/archive parts:

- [tgvr](/home/paul/cdkl5-variants/workflow/tgvr): active forward development for Thermodynamics-Guided Variant Reclassification
- [legacy_modularization](/home/paul/cdkl5-variants/workflow/legacy_modularization): older exact-reproduction modularization work kept as archive/reference
- [docs](/home/paul/cdkl5-variants/workflow/docs): shared memory, stage docs, and figures

## Active Path

Current work should move forward under:

- [workflow/tgvr](/home/paul/cdkl5-variants/workflow/tgvr)
- visual progress map: [workflow/docs/figures/tgvr_progress_map.drawio](/home/paul/cdkl5-variants/workflow/docs/figures/tgvr_progress_map.drawio)
- main TGVR environment: [workflow/tgvr/environment.tgvr.yml](/home/paul/cdkl5-variants/workflow/tgvr/environment.tgvr.yml)

The current active stage index is:

- [workflow/tgvr/README.md](/home/paul/cdkl5-variants/workflow/tgvr/README.md)
- [workflow/docs/stages/01_variant_curation.md](/home/paul/cdkl5-variants/workflow/docs/stages/01_variant_curation.md)
- [workflow/docs/stages/02_folding.md](/home/paul/cdkl5-variants/workflow/docs/stages/02_folding.md)

Current TGVR `01_variant_curation` supports two explicit master-build paths:

- `--stage master` for a master dataset without manual curated variants
- `--stage master --use-manual` for a master dataset that includes `workflow/tgvr/00_data/<gene>/01_variant_curation/manual/curated_variants.csv`
- important: `--stage all` currently means the non-manual path; for the manual-inclusive `CDKL5` path, run `--stage all`, then rerun `--stage master --use-manual`, then the summary commands

For a brand-new gene, TGVR now auto-initializes its own per-gene config and canonical UniProt FASTA from the supplied `GENE UNIPROT_ID` on first run.

Current TGVR `02_folding` is now package-usable as a stage under `workflow/tgvr`:

Current TGVR `02_folding` has been deliberately cleaned up to a fresh-start stage:

- TGVR no longer presents `02_folding` as a broad multi-method implemented stage
- current active direction is `DDGun` first
- current commands are for stage status and DDGun-first input preparation from the `01_variant_curation` master dataset
- DDGun source is present for inspection under `workflow/tgvr/third_party/ddgun/`, but the backend is intentionally not installed yet

Current TGVR uses one main environment for the package itself:

- create it with `conda env create -f workflow/tgvr/environment.tgvr.yml`
- activate it with `conda activate tgvr`
- keep heavier method runtimes such as `DDGun` as optional backends instead of forcing them into the base environment immediately

## Archived Path

The older modularization branch was moved into:

- [workflow/legacy_modularization](/home/paul/cdkl5-variants/workflow/legacy_modularization)

That archive contains the earlier `config/`, `main.py`, `src/`, `tests/`, and `outputs/` layout that was used for exact reproduction of the legacy notebook-driven workflow.
