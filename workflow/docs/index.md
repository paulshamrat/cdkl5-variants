# Workflow Docs

This is the local documentation home for the workflow development area under `workflow/`.

## What Lives Here

- `TGVR`
  The active development workflow for Thermodynamics-Guided Variant Reclassification.
- `Project Status`
  The public-safe summary of the current workflow state and verified results.
- `Reference Notes`
  Internal notes and legacy maps that remain in the repo but are not surfaced in the public browser docs.

## Recommended Reading Order

1. Open `TGVR > Overview` to understand the active workflow.
2. Open `TGVR > Usage` for the main runnable commands.
3. Open the individual TGVR stage pages for operational details.
4. Use `PROJECT_STATUS.md` for the public-safe current-state summary.

## Current State

The repository now has three important layers:

- `legacy study workflow`
  The original notebook-driven baseline at the repository root.
- `legacy modularization`
  The archived exact-reproduction modularization under `workflow/legacy_modularization/`.
- `active development workflow`
  The current TGVR development area under `workflow/tgvr/`.

The active stage with the strongest current support is `01_variant_curation`.
`02_folding` is a DDGun-first preparation stage, and `03_binding` is currently documented as the planned modular binding stage rather than a fully implemented runnable stage.

## Browse In A Local Browser

Install MkDocs if needed:

```bash
pip install mkdocs
```

Then from the repository root run:

```bash
mkdocs serve -f workflow/docs/mkdocs.yml -a 127.0.0.1:8001
```

Open:

```text
http://127.0.0.1:8001
```

## Internal Files

Some internal working files may remain locally in `workflow/docs/` for private continuity, but they should not be tracked in the public branch:

- `MEMORY.md`
- `notebook.md`
- `legacy_maps/`
- `publications/`
