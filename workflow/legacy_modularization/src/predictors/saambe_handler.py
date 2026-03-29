import pandas as pd
import glob
import os
from src.data_utils import get_absolute_path

def generate_saambe_mutation_list(df, output_path, chain="A"):
    """
    Generates a mutation list file for SAAMBE-3D.
    Format: CHAIN POSITION WILD MUTANT
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        for _, row in df.iterrows():
            wt, pos, mt = row['wild'], row['position'], row['mutant']
            if pd.isna(wt) or pd.isna(pos) or pd.isna(mt):
                continue
            f.write(f"{chain} {int(pos)} {wt} {mt}\n")
    print(f"Wrote mutations to {output_path}")

def merge_saambe_results(df, outputs_dir):
    """
    Parses SAAMBE-3D .out files and merges DDG values into the DataFrame.
    """
    out_files = glob.glob(os.path.join(outputs_dir, "*.out"))
    
    for out_file in out_files:
        # derive column name from filename, e.g., "ddg_1a22_str"
        base = os.path.basename(out_file).rsplit('.out', 1)[0]
        col_name = f'ddg_{base}_str'

        # read output, expected format: Wild Position Mutant ddG(kcal/mol)
        try:
            out_df = pd.read_csv(out_file, sep=r'\s+', header=0).rename(columns={'ddG(kcal/mol)': 'ddG'})
            
            for _, row in out_df.iterrows():
                mask = (
                    (df['wild'] == row['Wild']) &
                    (df['position'] == row['Position']) &
                    (df['mutant'] == row['Mutant'])
                )
                df.loc[mask, col_name] = row['ddG']
        except Exception as e:
            print(f"Error parsing {out_file}: {e}")

    return df
