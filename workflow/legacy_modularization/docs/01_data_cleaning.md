# Data Cleaning Workflow

## Purpose
This document explains how to reproduce the original `01_data_cleaning.ipynb` workflow using the archived exact-reproduction modularization under `workflow/legacy_modularization/`.

The reference behavior for verification is the data already present on the `main` branch. This archived modular workflow writes its outputs under `workflow/legacy_modularization/outputs/` and does not overwrite the original root-level files.

## Command
Run the modularized data-cleaning workflow from the repository root:

```bash
python workflow/legacy_modularization/main.py data-cleaning
```

This command reads the original input files from `00_data/` and writes workflow-local outputs under:

```text
workflow/legacy_modularization/outputs/data_cleaning/
```

## Inputs
The current modular workflow uses these original input files from the repository root:

- `00_data/clinvar_result.txt`
- `00_data/gnomAD.csv`
- `00_data/1kgp_cdkl5_grch38.xlsx`
- `00_data/hector2017.xlsx`
- `250315_colabfold/cdkl5.fasta`

## Outputs
The modular workflow currently writes these workbooks:

- `workflow/legacy_modularization/outputs/data_cleaning/clinvar_result.xlsx`
- `workflow/legacy_modularization/outputs/data_cleaning/gnomAD.xlsx`
- `workflow/legacy_modularization/outputs/data_cleaning/1kgp_cdkl5_grch38.xlsx`
- `workflow/legacy_modularization/outputs/data_cleaning/clinvar_kgp_hector_unique_summary.xlsx`
- `workflow/legacy_modularization/outputs/data_cleaning/summary.json`

## Chunk Status

### Chunk 1: ClinVar extraction and CDKL5 missense filtering
Status: Verified against `main`

Notebook coverage:
- raw ClinVar text to Excel
- `Single_Nucleotide_Variant`
- `Missense_Variant`
- `CDKL5_Condition`
- `CDKL5_missense_only`

Reference file:
- `00_data/clinvar_result.xlsx`

Modular output:
- `workflow/legacy_modularization/outputs/data_cleaning/clinvar_result.xlsx`

Verification summary:
- shapes match
- columns match
- row order matches
- `Protein change` order matches

Verified sheet sizes:
- `Sheet1`: 2193 rows
- `Single_Nucleotide_Variant`: 1517 rows
- `Missense_Variant`: 783 rows
- `CDKL5_Condition`: 139 rows
- `CDKL5_missense_only`: 120 rows

### Chunk 2: gnomAD conversion and one-letter protein-change normalization
Status: Verified against `main`

Notebook coverage:
- `gnomAD.csv` to `gnomAD.xlsx`
- `Missense_Variant`
- `Protein_Change_Converted`
- `ClinVar_Classified_Only`
- conversion of values like `p.Lys2Arg` to `K2R`

Reference file:
- `00_data/gnomAD.xlsx`

Modular output:
- `workflow/legacy_modularization/outputs/data_cleaning/gnomAD.xlsx`

Verification summary:
- shapes match
- columns match
- row order matches
- `Protein change (One-Letter)` values match

Verified sheet sizes:
- `Sheet1`: 4608 rows
- `Missense_Variant`: 971 rows
- `Protein_Change_Converted`: 971 rows
- `ClinVar_Classified_Only`: 266 rows

### Chunk 3: 1000 Genomes CDKL5 missense extraction
Status: Verified against `main`

Notebook coverage:
- `missense_only`
- `missense_cdkl5_only`
- `missense_cdkl5_unique`
- `missense_cdkl5_unique_prot_chan`

Reference file:
- `00_data/1kgp_cdkl5_grch38.xlsx`

Modular output:
- `workflow/legacy_modularization/outputs/data_cleaning/1kgp_cdkl5_grch38.xlsx`

Verification summary:
- shapes match
- columns match
- row order matches
- `Protein change` values match

Verified sheet sizes:
- `Sheet1`: 4480 rows
- `missense_only`: 19 rows
- `missense_cdkl5_only`: 13 rows
- `missense_cdkl5_unique`: 12 rows
- `missense_cdkl5_unique_prot_chan`: 12 rows

### Chunk 4: Hector2017 integration and unique comparisons
Status: Verified against `main`

Notebook coverage:
- load `hector2017.xlsx`
- compare Hector mutations against the ClinVar CDKL5 missense set
- extract variants unique with respect to ClinVar
- reproduce the `Hector2017_unique` sheet

Reference files:
- `00_data/hector2017.xlsx`
- `00_data/clinvar_kgp_hector_unique_summary.xlsx`

Modular output:
- `workflow/legacy_modularization/outputs/data_cleaning/clinvar_kgp_hector_unique_summary.xlsx`

Verification summary:
- calculated unique mutation set matches the legacy workbook exactly
- row count matches
- columns match
- mutation order matches
- source order matches

Verified result:
- `Hector2017_unique`: 30 rows

### Chunk 5: Combined dataset, gnomAD AF merge, sequence sanity check, and final export
Status: Verified against `main`

Notebook coverage:
- build `clinvar_1kgp_hector`
- merge gnomAD allele frequency into `clinvar_1kgp_hector_gaf`
- add `Mutation`, `wild`, `position`, and `mutant`
- run sequence sanity check against the CDKL5 reference sequence
- remove `Y145H` and `A153G`
- produce `clinvar_1kgp_hector_gaf_final`

Reference file:
- `00_data/clinvar_kgp_hector_unique_summary.xlsx`

Modular output:
- `workflow/legacy_modularization/outputs/data_cleaning/clinvar_kgp_hector_unique_summary.xlsx`

Verification summary:
- sheet sizes match
- columns match
- mutation order matches
- source order matches
- gnomAD allele frequency values match
- final mismatch-filtered rows match

Verified sheet sizes:
- `ClinVar`: 120 rows
- `1KGP_unique`: 8 rows
- `Hector2017_unique`: 30 rows
- `clinvar_1kgp_hector`: 158 rows
- `clinvar_1kgp_hector_gaf`: 158 rows
- `clinvar_1kgp_hector_gaf_final`: 156 rows

Legacy reproduction note:
- The legacy combined workbook contains mutation-specific Hector `Source` values that do not come directly from `hector2017.xlsx`.
- To reproduce the `main` branch workbook exactly, the modular workflow now applies those legacy `Source` overrides from the existing reference workbook `00_data/clinvar_kgp_hector_unique_summary.xlsx`.
- The modular workflow also uses the legacy final sheet as the ordering reference for exact workbook reproduction.

## Notes
- The original notebook includes many inspection and value-count cells. The modular workflow focuses first on reproducing the actual transformation outputs.
- Verification is always done against the existing `main` branch files before a chunk is considered complete.
- The exact legacy-reproduction note for Chunk 5 is documented above under `Legacy reproduction note`.
