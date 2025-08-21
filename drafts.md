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
   - pandas, matplotlib, seaborn, openpyxl, scikit-learn, tensorflow/torch, etc.  
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


4. **Modeling & Docking:**  
   - Use ColabFold for structure prediction (deep learning-based, see provided scripts in `250315_colabfold/`).
   - Run HADDOCK3 for docking (physics-based docking tool; see `250519_energy/03_haddock/`).


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
