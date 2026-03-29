# TGVR Usage

TGVR is currently driven through explicit stage scripts under [workflow/tgvr/scripts](/home/paul/cdkl5-variants/workflow/tgvr/scripts).

## Current Public Entry Points

- [run_variant_curation.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_variant_curation.py)
- [run_folding.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/run_folding.py)
- [manage_1kgp_palmetto.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/manage_1kgp_palmetto.py)

At the moment, there is no public `run_binding.py` script in the tree, so the binding stage should be treated as planned/documented work rather than a runnable public CLI.

Most stage commands take `GENE UNIPROT_ID`, for example:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage all
```

## 01 Variant Curation

This section is the main begin-to-finish workflow for `CDKL5` variant curation.

Typical `CDKL5` variant-curation layout:

```text
workflow/tgvr
├── 00_data
│   └── cdkl5
│       ├── 01_variant_curation
│       │   ├── 1kgp
│       │   ├── clinvar
│       │   ├── gnomad
│       │   └── manual
│       │       └── curated_variants.csv
│       └── reference
└── outputs
    └── cdkl5
        └── 01_variant_curation
            ├── logs
            └── master
```

### Optional Prerequisite: Manual Curated File

If you want `master --use-manual` to include manually curated variants, install the manual file before running the main staged workflow.

Install a CSV:

```bash
cd ~/cdkl5-variants
conda activate tgvr

python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --manual-file /path/to/curated_variants.csv
```

Or install an Excel file:

```bash
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --manual-file /path/to/curated_variants.xlsx
```

This copies the file into:

```text
workflow/tgvr/00_data/cdkl5/01_variant_curation/manual/
```

Expected minimal columns:

```text
Mutation,Consequence,Source
```

Example preview:

```csv
Mutation,Consequence,Source
R178W,Likely pathogenic,Hector2017
G20R,Uncertain significance,Manual review
```

### Canonical End-To-End Run

Use this command block when you want the full reproducible path from Palmetto-backed raw `1kgp` regeneration through the final TGVR summaries:

```bash
cd ~/cdkl5-variants
conda activate tgvr

rm -f ~/.ssh/palmetto.sock
ssh -M -S ~/.ssh/palmetto.sock -o ControlPersist=6h -fN $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock -O check $USER@slogin.palmetto.clemson.edu
ssh -S ~/.ssh/palmetto.sock $USER@slogin.palmetto.clemson.edu "hostname && whoami"

python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 setup
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 submit
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 status
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 log
python workflow/tgvr/scripts/manage_1kgp_palmetto.py CDKL5 fetch

python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage clinvar
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode cached
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage gnomad
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --use-manual
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-full
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-range --range 1-302 --label "Kinase domain"
```

### Cached-Only Local Rerun

If the cached `1kgp` workbook is already trusted and does not need regeneration, you can skip the Palmetto block and run:

```bash
cd ~/cdkl5-variants
conda activate tgvr

python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage clinvar
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage 1kgp --1kgp-mode cached
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage gnomad
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage master --use-manual
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-full
python workflow/tgvr/scripts/run_variant_curation.py CDKL5 O76039 --stage summary-range --range 1-302 --label "Kinase domain"
```

The `all` shortcut is still available, but the explicit staged form is preferred because the master depends on a correct `1kgp` stage.

Important note:

- `--stage all` currently means the non-manual path
- if you want manual curated variants included, rerun the `master` step explicitly with `--use-manual`
- the manual curated file is not part of the `1kgp` stage itself; it is read only by the `master --use-manual` step from `workflow/tgvr/00_data/<gene>/01_variant_curation/manual/curated_variants.csv` (or an installed `.xlsx` equivalent)

### Variant-Curation Outputs

Current verified outputs from the staged rerun:

| Output | Value |
| --- | ---: |
| Full-length total | `162` |
| Kinase-domain total (`1-302`) | `114` |
| Benign | `22` |
| Benign/Likely benign | `15` |
| Conflicting classifications of pathogenicity | `7` |
| Likely benign | `10` |
| Likely pathogenic | `26` |
| Pathogenic | `10` |
| Pathogenic/Likely pathogenic | `24` |
| Uncertain significance | `48` |

### Palmetto 1KGP Details

The heavier raw `1kgp` reproduction path is a Palmetto-backed subworkflow used only when you need to regenerate the cached `1kgp` workbook before the local TGVR master build.

Public manager:

- [manage_1kgp_palmetto.py](/home/paul/cdkl5-variants/workflow/tgvr/scripts/manage_1kgp_palmetto.py)

Detailed stage notes:

- [01 Variant Curation: 1KGP](/home/paul/cdkl5-variants/workflow/docs/tgvr/stages/01_variant_curation_1kgp.md)

## 02 Folding

Show stage status:

```bash
python workflow/tgvr/scripts/run_folding.py CDKL5 O76039 --stage status
```

Prepare DDGun input:

```bash
python workflow/tgvr/scripts/run_folding.py CDKL5 O76039 --stage prepare-ddgun --method ddgun_seq --range 1-302 --label "Kinase domain"
```
