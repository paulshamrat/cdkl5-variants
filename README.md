
# CDKL5 Variants Analysis

This repository provides an end-to-end workflow for the analysis of CDKL5 gene variants, including variant annotation, stability prediction, CDKL5 and it's substrate structural modeling, and reclassification. The workflow integrates state-of-the-art bioinformatics, machine learning, and deep learning tools for comprehensive variant effect assessment, substrate interaction modeling, and protein-protein docking.

## Features

- **Variant Annotation & Classification:**  
  - Integrates ClinVar, gnomAD, and other population/clinical databases.
  - Supports germline classification, reclassification, and value count summaries.
  - Ready for integration with AI/ML-based variant effect predictors.

- **Stability Prediction:**  
  - Sequence-based and structure-based ΔΔG (ddG) predictions using multiple methods (SAAFEC-SEQ, I-Mutant, INPS, DDGun, mCSM, DDMut, DynaMut, DDGEmb, etc.).
  - Compatible with AI/ML and deep learning-based stability predictors (e.g., DDGEmb, AlphaFold2-based, and future models).
  - Automated merging and comparison of results across methods.
  - Publication-ready violin plots and statistical summaries.

- **Structural Modeling & AI Integration:**  
  - Substrate and variant structure prediction using ColabFold (AlphaFold2, deep learning-based modeling).
  - Protein-protein docking with HADDOCK3 (physics-based docking tool).
  - Analysis of interaction interfaces and energetic changes.
  - Designed for easy extension with generative AI, ML-based scoring, and structure prediction tools.

- **Visualization & Reporting:**  
  - Multi-panel, customizable plots (PDF export, font/legend control).
  - Jupyter notebooks for reproducible analysis and figure generation.
  - Data cleaning, inspection, and summary tables.
  - Ready for integration with AI/ML-based visualization and analytics.

## Directory Structure

```
cdkl5-variants/
├── 250315_colabfold/                # ColabFold/AlphaFold2 structure prediction scripts and results
├── 250519_energy/                   # Main energy, stability, and docking analysis
│   ├── 00_data/
│   ├── 01_data_cleaning.ipynb
│   ├── 02_folding.ipynb
│   ├── 03_haddock/
│   │   ├── 250401_01_docking_data_processing.ipynb
│   │   ├── 250401_02_docking_automation_test.ipynb
│   │   ├── 250728_haddock_energies.py
│   │   └── ...
│   ├── 04_binding.ipynb
│   ├── 05_pathogenicity
│   ├── 250630_reclassification_of_mutation_final2_last_minute_corr.ipynb
│   ├── 250703_pathogenicity_2.ipynb
│   └── ...
├── README.md
```

## Getting Started

1. **Clone the repository:**
   ```bash
   git clone https://github.com/paulshamrat/ColabMDA.git
   cd cdkl5-variants
   ```

2. **Install dependencies:**  
   - Python 3.8+  
   - pandas, matplotlib, seaborn, openpyxl, scikit-learn, tensorflow/torch (for AI/ML), etc.  
   - (See notebooks for full requirements.)

3. **Run the analysis:**  
   - Open the relevant notebook(s) in Jupyter or VS Code for your analysis goal:
     - `250519_energy/01_data_cleaning.ipynb`: Data cleaning and inspection.
     - `250519_energy/02_folding.ipynb`: Folding and stability analysis.
     - `250519_energy/03_haddock/250401_01_docking_data_processing.ipynb`: Docking data processing (HADDOCK3).
     - `250519_energy/03_haddock/250401_02_docking_automation_test.ipynb`: Docking automation and batch processing.
     - `250519_energy/03_haddock/250728_haddock_energies.py`: Energetic analysis of docking results.
     - `250519_energy/04_binding.ipynb`: Binding analysis.
     - `250519_energy/250630_reclassification_of_mutation_final2_last_minute_corr.ipynb`: Final variant reclassification.
     - `250519_energy/250703_pathogenicity_2.ipynb`: Pathogenicity analysis.
     - `250315_colabfold/`: Structure prediction scripts and batch jobs (see shell scripts and FASTA files).
     - Additional notebooks/scripts for focused tasks (e.g., plotting, data exploration, ML integration).
   - Follow the notebook cells for data loading, cleaning, prediction, and visualization.
   - Integrate or extend with your own AI/ML models as needed.

4. **Modeling & Docking:**  
   - Use ColabFold for structure prediction (deep learning-based, see provided scripts in `250315_colabfold/`).
   - Run HADDOCK3 for docking (physics-based docking tool; see `250519_energy/03_haddock/`).
   - Ready for integration with generative AI and ML-based structure/docking tools.

## Notebooks

- `250519_energy/01_data_cleaning.ipynb`: Data cleaning and inspection.
- `250519_energy/02_folding.ipynb`: Folding and stability analysis.
- `250519_energy/03_haddock/250401_01_docking_data_processing.ipynb`: Docking data processing (HADDOCK3).
- `250519_energy/03_haddock/250401_02_docking_automation_test.ipynb`: Docking automation and batch processing.
- `250519_energy/03_haddock/250728_haddock_energies.py`: Energetic analysis of docking results.
- `250519_energy/04_binding.ipynb`: Binding analysis.
- `250519_energy/250630_reclassification_of_mutation_final2_last_minute_corr.ipynb`: Final variant reclassification.
- `250519_energy/250703_pathogenicity_2.ipynb`: Pathogenicity analysis.
- `250315_colabfold/`: Structure prediction scripts and batch jobs (see shell scripts and FASTA files).
- Additional notebooks/scripts for modeling, docking, ML integration, and downstream analysis.


## Data

- **Main variant data:**  
  `250519_energy/00_data/01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af_noddg.xlsx`
- **Folding ΔΔG dataset:**  
  `250519_energy/00_data/01_cdkl5_clinvar_gaf_1kgp_hctr_comb_unq_af.xlsx`
- **Binding ΔΔG datasets:**  
  - `250519_energy/04_binding/clinvar_1kgp_hector_gaf_final_binding.xlsx`  
  - `250519_energy/04_binding/clinvar_1kgp_hector_gaf_final_binding_znf219_111_only.xlsx`  
  - Method-specific subfolders:  
    `01_saambeseq/`, `02_saambe3d/`, `03_foldx/`, `04_mutabind2/`, `05_mcsmppi/`, `06_beatmusic/`, `07_ddmutppi/`, `08_bindprofx/`, `10_isee/`  
  - Additional results and analysis in `analysis/` and `outputs/` subfolders.
- **Prediction results and intermediate files:**  
  Located throughout `250519_energy/` and its subfolders.

## Citation

If you use this workflow or code, please cite the relevant tools and this repository.

---

**Contact:**  
For questions or contributions, please open an issue or contact the repository maintainer.
