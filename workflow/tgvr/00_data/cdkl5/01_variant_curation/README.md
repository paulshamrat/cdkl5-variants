Seeded workflow-local inputs for `CDKL5` variant curation.

These files were restored from the legacy top-level [`00_data`](/home/paul/cdkl5-variants/00_data) area so the TGVR workflow can be tested from within [`workflow/tgvr`](/home/paul/cdkl5-variants/workflow/tgvr) without depending on external copies.

Included seeds:

- `manual/curated_variants.xlsx`
  Source: legacy `00_data/hector2017.xlsx`
- `1kgp/1kgp_cdkl5_grch38.xlsx`
  Source: legacy `00_data/1kgp_cdkl5_grch38.xlsx`
- `gnomad/gnomad.csv`
  Source: legacy `00_data/gnomAD.csv`

These are workflow-scoped seed inputs and benchmarks. Generated outputs should continue to go under [`workflow/tgvr/outputs`](/home/paul/cdkl5-variants/workflow/tgvr/outputs).
