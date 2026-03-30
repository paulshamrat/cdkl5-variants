# 02 Folding

`02_folding` is currently a fresh-start TGVR folding stage with a deliberately narrow scope.

## Current Status

This stage is not yet a full end-to-end folding workflow.

What is currently reproducible:

- loading the TGVR master mutation dataset
- selecting an optional residue range
- preparing DDGun-ready mutation inputs

What is not yet complete:

- a fully integrated DDGun runtime backend
- result ingestion and downstream scientific summaries

## Main Runner

- `workflow/tgvr/scripts/run_folding.py`

## Reusable Implementation

- `workflow/tgvr/src/folding/core.py`

## Current Methods

The active method family right now is `DDGun`.

Supported preparation modes:

- `ddgun_seq`
- `ddgun_str`

## Commands

Show current stage status:

```bash
python workflow/tgvr/scripts/run_folding.py GENE UNIPROT_ID --stage status
```

Prepare DDGun input:

```bash
python workflow/tgvr/scripts/run_folding.py GENE UNIPROT_ID --stage prepare-ddgun --method ddgun_seq --range START-END --label "Feature name"
```

## Why The Scope Is Narrow Right Now

A real DDGun backend still requires heavier supporting assets such as `hh-suite` and the `uniclust30_2018_08` database.

That makes full local installation relatively heavy, so the current TGVR state stops at preparation rather than pretending the whole backend is already reproducible.

## Immediate Next Step

The next major improvement for this stage is to define a reproducible DDGun execution backend:

- local install
- container
- or Palmetto-backed execution
