# TGVR Installation

TGVR currently uses one main Conda environment for local orchestration.

## Main Environment

Environment file:

- `workflow/tgvr/environment.tgvr.yml`

Create and activate it from the repository root:

```bash
conda env create -f workflow/tgvr/environment.tgvr.yml
conda activate tgvr
```

## What This Environment Covers

The main `tgvr` environment is intended for:

- workflow orchestration
- data ingestion and staging
- master-table building
- summary generation
- preparation for heavier downstream tools

It is intentionally not trying to bundle every scientific backend immediately.

## Optional Backends

Some heavier tools are still expected to remain optional or external:

- `DDGun`
- Palmetto-backed raw `1kgp` processing
- external binding predictors and partner-model workflows

This keeps the base environment lighter and avoids forcing large, harder-to-reproduce scientific installations into the core TGVR setup too early.

## Public Variant-Curation Entry Point

Variant curation now uses one public script:

- `workflow/tgvr/scripts/run_variant_curation.py`

That single entrypoint covers:

- live `ClinVar`
- `1kgp` cached mode
- `1kgp` public VCF-backed mode
- `1kgp` Palmetto-backed mode
- live `gnomAD`
- master table generation
- summary tables

Use the public VCF-backed `1kgp` path when you want a Palmetto-free run. Use the Palmetto mode only when you explicitly want the legacy raw backend.

## Local Browser Docs

To browse these docs locally in a browser:

```bash
pip install mkdocs
mkdocs serve -f workflow/docs/mkdocs.yml
```

Then open:

```text
http://127.0.0.1:8000
```
