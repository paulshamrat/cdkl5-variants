# TGVR Installation

TGVR currently uses one main Conda environment for local orchestration.

## Main Environment

Environment file:

- [environment.tgvr.yml](/home/paul/cdkl5-variants/workflow/tgvr/environment.tgvr.yml)

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
