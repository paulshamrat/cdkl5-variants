# Tutorial: High-Throughput CDKL5 Variant Analysis

This tutorial guides you through using the modularized workflow for analyzing CDKL5 missense variants. The workflow integrates data cleaning, external predictor submission (PolyPhen-2, SAAMBE-3D), and thermodynamic reclassification.

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
> Use repository-relative paths in the config. The tools will automatically resolve them to absolute paths on your system.

## 2. Automated Initial Processing

Run the orchestration script to clean the raw data and generate input files for external tools:

```bash
python main.py
```

This script performs the following:
1.  **Data Cleaning**: Reads `00_data/clinvar_result.txt` and formats it into `01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af.xlsx`.
2.  **PolyPhen-2 Batch**: Generates `05_pathogenicity/01_polyphen2/batch_submission.txt`.
3.  **SAAMBE-3D Batch**: Generates `04_binding/mutations_list.txt`.

## 3. Running External Predictors

### PolyPhen-2
1.  Go to the [PolyPhen-2 Webserver](http://genetics.bwh.harvard.edu/pph2/).
2.  Upload the generated `batch_submission.txt`.
3.  Download the results as a TSV file and save it to the path specified in `config.yaml` as `polyphen2_results`.

### SAAMBE-3D (on Palmetto HPC)
1.  Transfer `04_binding/mutations_list.txt` to Palmetto.
2.  Run the SAAMBE-3D script (refer to `04_binding.ipynb` for the exact SSH/salloc commands).
3.  Transfer the resulting `.out` files back to the `04_binding/outputs/saambe_3d/` directory.

## 4. Reclassification Analysis

Once you have the external results, use the refactored Jupyter notebooks to complete the analysis.

- **04_binding.ipynb**: Merges SAAMBE-3D results back into the variant database.
- **06_variant_reclass_pathogenicity.ipynb**: Merges PolyPhen-2 results.
- **05_variant_reclass_ddG_FoldingBinding.ipynb**: Performs the final thermodynamic reclassification using the collected data.

The notebooks are now modular:
```python
from src.analysis.energy_analysis import calculate_max_ddg
# Reclassification logic is now handled by src/
```

## 5. Verification & Testing

You can verify that the core utilities are working correctly by running the built-in tests:
```bash
python tests/test_utils.py
```
