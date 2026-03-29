# TGVR Architecture

TGVR is being built as a modular workflow rather than a notebook-only pipeline.

## Design Principles

- keep the repository root as the untouched legacy baseline
- place active development under `workflow/`
- make stage boundaries explicit
- separate reusable implementation from user-facing entrypoints
- keep heavy outputs and downloaded runtime artifacts local by default

## Layering

There are currently three important layers in the repository:

- `legacy study workflow`
  The original root-level notebook-driven workflow.
- `legacy modularization`
  The archived exact-reproduction modularization stored under `workflow/legacy_modularization/`.
- `active development workflow`
  The forward-moving TGVR area under `workflow/tgvr/`.

## TGVR Structure

```text
workflow/tgvr/
  config/
  internal/
  scripts/
  src/
```

Responsibilities:

- `scripts/`
  Small public command entrypoints for users.
- `src/`
  Reusable implementation modules.
- `internal/`
  Helper code that supports public scripts but is not intended as the main public interface.
- `config/`
  Gene-specific configuration and reusable workflow settings.

## Runtime Data Model

TGVR organizes work by gene and by stage.

Typical paths:

```text
workflow/tgvr/00_data/<gene>/<stage>/
workflow/tgvr/outputs/<gene>/<stage>/
```

This allows:

- repeatable per-gene execution
- stage-local inputs and outputs
- cleaner progression from curation to downstream modeling

## Current Public Surface

The currently exposed runnable surface is centered on:

- variant curation
- folding preparation
- Palmetto management for the raw `1kgp` path

The binding stage is documented conceptually, but it is not currently exposed through a public script in `workflow/tgvr/scripts/`.

## Current State Summary

- `01_variant_curation` is the strongest stage today
- `02_folding` is prepared for DDGun-first work but not yet end-to-end
- `03_binding` remains a planned modular stage whose public implementation still needs to be completed
