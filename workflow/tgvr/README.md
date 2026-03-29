# TGVR

This folder contains the active development workflow for `Thermodynamics-Guided Variant Reclassification (TGVR)`.

The public-safe browser documentation lives under:

- [workflow/docs/tgvr](/home/paul/cdkl5-variants/workflow/docs/tgvr)

Use that documentation for:

- workflow overview
- environment setup
- canonical `CDKL5` variant-curation commands
- current stage status

## Main Entry Points

Current public entry scripts:

- [run_variant_curation.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_variant_curation.py)
- [run_folding.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_folding.py)
- [manage_1kgp_palmetto.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/manage_1kgp_palmetto.py)

## Browser Docs

Run locally from the repository root:

```bash
mkdocs serve -f workflow/docs/mkdocs.yml -a 127.0.0.1:8001
```

Then open:

```text
http://127.0.0.1:8001/
```

Expected current `CDKL5` result after the manual-inclusive path:

- full-length variants: `162`
- kinase-domain variants `1-302`: `114`
- note: the earlier UniProt-style kinase annotation had been discussed as `13-297`, but TGVR should use `1-302` for now
- retained 1KGP benign rows in the final master: `8`

</details>

<details open>
<summary><strong>02 Folding</strong></summary>

This stage has been intentionally reset to a clean fresh start.

Current meaning of `02_folding`:

- TGVR folding is **not** being treated as broadly implemented yet
- the first active method family is `DDGun`
- the immediate goal is reproducible `DDGun` setup plus DDGun-first input preparation from TGVR `01_variant_curation`
- DDGun source is present for inspection under `workflow/tgvr/third_party/ddgun/`, but the backend is intentionally not installed yet because full setup requires `hh-suite` and a large `uniclust30` database
- the main `tgvr` environment does **not** install DDGun

Main stage guide:

- [workflow/docs/stages/02_folding.md](/home/paul/cdkl5-variants/workflow/docs/stages/02_folding.md)

Shared runner:

- [scripts/run_folding.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_folding.py)

Main commands:

```bash
python workflow/tgvr/scripts/run_folding.py CDKL5 O76039 --stage status
python workflow/tgvr/scripts/run_folding.py CDKL5 O76039 --stage prepare-ddgun --method ddgun_seq --range 1-302 --label "Kinase domain"
```

Current `CDKL5` locations:

- data: `workflow/tgvr/00_data/cdkl5/02_folding/`
- outputs: `workflow/tgvr/outputs/cdkl5/02_folding/`

</details>

<details>
<summary><strong>03 Binding</strong></summary>

This stage is planned for later return.

The earlier binding scaffold is not part of the current minimal TGVR tree, so the active focus remains:

- `01_variant_curation`
- `02_folding` fresh start with `DDGun`

</details>

## Notes

- `assess_local_gene_bootstrap.py` is older exploratory code and should not be treated as the main stage runner.
- Update this file whenever a new stage is added or when a stage changes materially.
- For now, treat `workflow/docs/stages/03_binding.md` as archived planning context, not as an active TGVR stage in the current minimal tree.
