# Folding Manual Inputs

This archived modularization reproduces the input-preparation and analysis logic from `02_folding.ipynb`.

What is automated now:
- creation of mutation-list files for the legacy folding predictors
- summary analysis of sequence/structure DDG columns already present in the legacy workbook

What remains external or manual from the original notebook:
- SAAFEC-SEQ execution on Palmetto or another suitable environment
- I-Mutant sequence/structure submissions
- INPS sequence and structure submissions
- DDGun sequence and structure submissions
- DDGEmb sequence submissions
- mCSM structure submissions
- DDMut structure submissions

If those predictor outputs are collected later, they should be stored in workflow-local folders under:
- /home/paul/cdkl5-variants/workflow/legacy_modularization/outputs/folding/prepared_inputs
