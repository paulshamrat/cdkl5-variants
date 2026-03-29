# Folding Workflow

## Purpose
This document explains the archived modularization of the original [02_folding.ipynb](/home/paul/cdkl5-variants/02_folding.ipynb) workflow.

The legacy notebook mixed together three different kinds of work:

- preparing mutation lists for external folding predictors
- manually or externally running those predictors
- analyzing the final sequence- and structure-based DDG columns

The archived modularization now separates those into explicit sub-steps.

## Commands
Run from the repository root.

Prepare the folding-predictor mutation lists:

```bash
python workflow/legacy_modularization/main.py folding --step prepare
```

Reproduce the legacy DDG summary analysis from the legacy workbook:

```bash
python workflow/legacy_modularization/main.py folding --step analyze
```

Run both sub-steps:

```bash
python workflow/legacy_modularization/main.py folding --step all
```

## Inputs
This archived modular workflow currently reads the legacy workbook:

- [00_data/01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af.xlsx](/home/paul/cdkl5-variants/00_data/01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af.xlsx)

That workbook already contains the final legacy folding DDG columns used in the notebook analysis.

## Outputs
Prepared mutation-list files are written under:

- [workflow/legacy_modularization/outputs/folding/prepared_inputs](/home/paul/cdkl5-variants/workflow/legacy_modularization/outputs/folding/prepared_inputs)

Analysis outputs are written under:

- [workflow/legacy_modularization/outputs/folding/analysis](/home/paul/cdkl5-variants/workflow/legacy_modularization/outputs/folding/analysis)

Important generated files:

- `prepared_inputs/01_saafecseq/mutation_list.txt`
- `prepared_inputs/03_inps_seq/cdkl5_mutations_only.txt`
- `prepared_inputs/06_ddgemb_seq/cdkl5_mutations_only_ddgemb.txt`
- `prepared_inputs/07_mcsm_str/cdkl5_mutations_formatted.txt`
- `prepared_inputs/08_ddmut_str/cdkl5_mutations_formatted.txt`
- `analysis/cdkl5_folding_ddg_summary_table.xlsx`
- `analysis/cdkl5_folding_ddg_summary_table.csv`

## What Is Modularized Now

### Prepare
The `prepare` step reproduces the input-list generation logic from the notebook for:

- `SAAFEC-SEQ`
- `INPS sequence`
- `DDGEmb sequence`
- `mCSM structure`
- `DDMut structure`

It also writes a note about the still-manual external predictor steps.

### Analyze
The `analyze` step reproduces the final summary-table logic from the notebook:

- identify all sequence-based DDG columns using `^ddg.*seq$`
- identify all structure-based DDG columns using `^ddg.*str$`
- compare `Benign` versus `Pathogenic`
- summarize both:
  - full region `1-895`
  - kinase region `1-302`

For each method/class/region, it writes:

- `n`
- `Mean`
- `Median`
- `Std`
- `Min`
- `Max`

## What Still Requires External Results
The original notebook depended on several predictors that were not run directly inside the repository:

- `SAAFEC-SEQ` on Palmetto
- `I-Mutant`
- `INPS`
- `DDGun`
- `DDGEmb`
- `mCSM`
- `DDMut`

So this archived modularization does **not** yet rerun those external predictors automatically.
Instead, it modularizes:

- the input preparation
- the downstream summary analysis

This is the cleanest faithful split of what the original notebook actually did.
