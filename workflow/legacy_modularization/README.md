# Workflow

This folder contains the modularized workflow layer for analyzing CDKL5 missense variants. It is intentionally separated from the legacy notebook-first study files at the repository root.

The goal is to keep the original workflow intact while developing a reusable controller layer here in `workflow/`.

## Layout

- `config/`: workflow-specific configuration
- `docs/`: memory and supporting documentation
- `src/`: reusable modular code
- `tests/`: workflow-specific tests
- `outputs/`: generated files created by the modular workflow

## 0. Setup Environment

Before starting, ensure you have the `cdkl5-activation` environment active:
```bash
conda activate cdkl5-activation
```

## 1. Project Configuration

The entire workflow is driven by `config/config.yaml`. Open this file to configure:
- `uniprot_id`: The target protein (default: `O76039`).
- `residue_range`: The range of residues considered (default: `[1, 302]`).
- `paths`: Key directories for raw data and predictor results.

> [!TIP]
> The workflow may read legacy data from the repository root using `../...` paths, but generated files should go into `workflow/outputs/` by default. This keeps the original baseline untouched.

## 2. Automated Initial Processing

Run the orchestration script to clean the raw data and generate input files for external tools:

```bash
python workflow/main.py
```

You can also `cd workflow` and run:

```bash
python main.py
```

This script performs the following:
1.  **Data Cleaning**: Reads `00_data/clinvar_result.txt` and writes a workflow-local cleaned table under `workflow/outputs/processed/`.
2.  **PolyPhen-2 Batch**: Generates a workflow-local batch file under `workflow/outputs/predictor_inputs/`.
3.  **SAAMBE-3D Batch**: Generates a workflow-local mutation list under `workflow/outputs/predictor_inputs/`.

Archived exact-reproduction commands now available from the repo root:

```bash
python workflow/legacy_modularization/main.py data-cleaning
python workflow/legacy_modularization/main.py folding --step all
python workflow/legacy_modularization/main.py binding --step all
```

## 3. Running External Predictors

### PolyPhen-2
1.  Go to the [PolyPhen-2 Webserver](http://genetics.bwh.harvard.edu/pph2/).
2.  Upload the generated `batch_submission.txt`.
3.  Download the results as a TSV file and save it to the path specified in `config.yaml` as `polyphen2_results`.

### SAAMBE-3D (on Palmetto HPC)
1.  Transfer the generated SAAMBE mutation list from `workflow/outputs/predictor_inputs/` to Palmetto.
2.  Run the SAAMBE-3D script (refer to `04_binding.ipynb` for the exact SSH/salloc commands).
3.  Transfer the resulting `.out` files into `workflow/outputs/predictor_results/saambe_3d/`.

## 4. Reclassification Analysis

Once you have the external results, use the refactored Jupyter notebooks to complete the analysis.

- **04_binding.ipynb**: Merges SAAMBE-3D results back into the variant database.
- **06_variant_reclass_pathogenicity.ipynb**: Merges PolyPhen-2 results.
- **05_variant_reclass_ddG_FoldingBinding.ipynb**: Performs the final thermodynamic reclassification using the collected data.

Archived stage docs:

- [01_data_cleaning.md](/home/paul/cdkl5-variants/workflow/legacy_modularization/docs/01_data_cleaning.md)
- [02_folding.md](/home/paul/cdkl5-variants/workflow/legacy_modularization/docs/02_folding.md)
- [04_binding.md](/home/paul/cdkl5-variants/workflow/legacy_modularization/docs/04_binding.md)

The notebooks are now modular:
```python
from src.analysis.energy_analysis import calculate_max_ddg
# Reclassification logic is now handled by src/
```

## 5. Verification & Testing

You can verify that the core utilities are working correctly by running the built-in tests:
```bash
python workflow/tests/test_utils.py
```
