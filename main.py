import os
import sys
import pandas as pd
from src.data_utils import load_config, get_absolute_path, add_mutation_components
from src.predictors.saambe_handler import generate_saambe_mutation_list
from src.predictors.polyphen2_handler import generate_polyphen2_batch
from src.analysis.energy_analysis import calculate_max_ddg

def main():
    # 1. Load config
    config = load_config()
    uniprot_id = config['project']['uniprot_id']
    
    # 2. Paths
    data_dir = get_absolute_path(config['paths']['data_dir'])
    raw_txt = os.path.join(data_dir, config['files']['raw_clinvar'])
    cleaned_xlsx = get_absolute_path(config['files']['processed_variants'])
    
    print(f"Starting modularized workflow for {config['project']['name']}...")
    
    # 3. Data Cleaning
    if os.path.exists(raw_txt):
        print(f"Loading raw data from {raw_txt}...")
        df = pd.read_csv(raw_txt, sep='\t', low_memory=False)
        # Basic filtering for missense/variants if needed
        df = add_mutation_components(df)
        df.to_excel(cleaned_xlsx, index=False)
        print(f"Cleaned data saved to {cleaned_xlsx}")
    else:
        print(f"Raw data not found at {raw_txt}, skipping cleaning step.")
        if os.path.exists(cleaned_xlsx):
            df = pd.read_excel(cleaned_xlsx)
        else:
            print("No data available to process.")
            return

    # 4. Predictor Input Generation
    saambe_muts_path = get_absolute_path(os.path.join(config['paths']['binding_dir'], "mutations_list.txt"))
    generate_saambe_mutation_list(df, saambe_muts_path)
    
    pp2_batch_path = get_absolute_path(os.path.join(config['paths']['polyphen2_results'], "batch_submission.txt"))
    generate_polyphen2_batch(df, uniprot_id, pp2_batch_path)
    
    print("\n--- Next Steps ---")
    print(f"1. Run SAAMBE-3D on Palmetto using the mutation list at: {saambe_muts_path}")
    print(f"2. Submit the PolyPhen-2 batch file at: {pp2_batch_path}")
    print("3. Once results are available, run reclassification analysis using the refactored notebooks.")

if __name__ == "__main__":
    main()
