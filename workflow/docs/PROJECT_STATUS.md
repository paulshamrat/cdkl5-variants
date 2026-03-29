# Project Status

This file is the public-safe status summary for the workflow development branch.

## Current Structure

The repository currently uses three layers:

- `legacy study workflow`
  The original notebook-driven baseline at the repository root.
- `legacy modularization`
  The archived exact-reproduction modularization under `workflow/legacy_modularization/`.
- `active development workflow`
  The current TGVR development area under `workflow/tgvr/`.

## Current TGVR State

Current stage maturity:

- `01_variant_curation`
  Most mature and actively verified stage.
- `02_folding`
  DDGun-first preparation stage.
- `03_binding`
  Planned/documented stage, not currently exposed as a public runnable script.

## Current Verified CDKL5 Variant-Curation Result

Current verified end-to-end `CDKL5` rerun results:

- full-length total: `162`
- kinase-domain total (`1-302`): `114`

Full-length class counts:

- `Benign` 22
- `Benign/Likely benign` 15
- `Conflicting classifications of pathogenicity` 7
- `Likely benign` 10
- `Likely pathogenic` 26
- `Pathogenic` 10
- `Pathogenic/Likely pathogenic` 24
- `Uncertain significance` 48

## Public Documentation

Public browser-safe docs live under:

- `workflow/docs/tgvr/`

Run them locally with:

```bash
mkdocs serve -f workflow/docs/mkdocs.yml -a 127.0.0.1:8001
```

Then open:

```text
http://127.0.0.1:8001/
```

## Private Notes

Personal notes, private infrastructure details, and local reference material should stay outside the public branch history.
