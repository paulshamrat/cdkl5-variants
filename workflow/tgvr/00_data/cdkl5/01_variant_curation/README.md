Stable workflow-local inputs for `CDKL5` variant curation.

These files are the preserved inputs that should remain in [`workflow/tgvr/00_data`](/home/paul/cdkl5-variants/workflow/tgvr/00_data) between runs.

Current stable inputs:

- `manual/curated_variants.xlsx`
  Source: legacy `00_data/hector2017.xlsx`
- `../reference/canonical_uniprot_O76039.fasta`
  Source: UniProt canonical sequence for `O76039`

Everything generated at run time now belongs under [`workflow/tgvr/outputs`](/home/paul/cdkl5-variants/workflow/tgvr/outputs), including:

- downloaded ClinVar and gnomAD caches
- cached `1kgp` workbooks
- public `phase3` VCF source downloads
- fetched Palmetto run directories
- logs, summaries, and final master outputs
