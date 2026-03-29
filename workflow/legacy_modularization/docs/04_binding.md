# 04 Binding

This archived module captures the reproducible parts of the legacy [04_binding.ipynb](/home/paul/cdkl5-variants/04_binding.ipynb) workflow.

## What It Reproduces

- preparation of mutation-list inputs for the legacy binding predictors
- summary analysis of binding DDG columns already present in the legacy final binding workbook

## Commands

Run from the repository root:

```bash
python workflow/legacy_modularization/main.py binding --step prepare
python workflow/legacy_modularization/main.py binding --step analyze
python workflow/legacy_modularization/main.py binding --step all
```

## Inputs

Preparation uses:

- [workflow/legacy_modularization/config/config.yaml](/home/paul/cdkl5-variants/workflow/legacy_modularization/config/config.yaml)
- [04_binding/clinvar_1kgp_hector_gaf_final.xlsx](/home/paul/cdkl5-variants/04_binding/clinvar_1kgp_hector_gaf_final.xlsx)

Analysis uses:

- [04_binding/analysis/clinvar_1kgp_hector_gaf_final_binding.xlsx](/home/paul/cdkl5-variants/04_binding/analysis/clinvar_1kgp_hector_gaf_final_binding.xlsx)

## Outputs

- [workflow/legacy_modularization/outputs/binding/prepared_inputs](/home/paul/cdkl5-variants/workflow/legacy_modularization/outputs/binding/prepared_inputs)
- [workflow/legacy_modularization/outputs/binding/analysis](/home/paul/cdkl5-variants/workflow/legacy_modularization/outputs/binding/analysis)

## Current Scope

The archived modularization currently automates:

- `SAAMBE-3D` mutation list preparation
- `FoldX` individual list template preparation
- `mCSM-PPI2` mutation list preparation
- `DDMutPPI` mutation list preparation
- binding DDG summary table generation from the legacy final workbook

It does **not** yet automate the original notebook's external predictor submissions or HPC execution.

