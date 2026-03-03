import pandas as pd
import os
from src.data_utils import get_absolute_path

def generate_polyphen2_batch(df, uniprot_id, output_path):
    """
    Generates a batch submission file for the PolyPhen-2 webserver.
    Format: UNIPROT_ID POSITION WILD MUTANT
    """
    lines = []
    for _, row in df.iterrows():
        wt, pos, mt = row["wild"], row["position"], row["mutant"]
        if pd.isna(wt) or pd.isna(pos) or pd.isna(mt):
            continue
        lines.append(f"{uniprot_id} {int(pos)} {wt} {mt}")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as fo:
        fo.write("\n".join(lines))
    print(f"Wrote {len(lines)} variants to {output_path}")

def merge_polyphen2_results(variants_df, polyphen_results_path):
    """
    Merges PolyPhen-2 results into the variants DataFrame.
    Assumes polyphen_results_path points to an Excel file with 'pos', 'prediction', etc.
    """
    pp2_df = pd.read_excel(polyphen_results_path)
    
    # Normalize key columns
    variants_df['position'] = variants_df['position'].astype(int)
    pp2_df['pos'] = pp2_df['pos'].astype(int)

    # Select relevant columns and rename for consistency
    cols_to_pull = ['pos', 'prediction', 'pph2_prob', 'pph2_FPR', 'pph2_TPR']
    pp2_small = pp2_df[cols_to_pull].rename(columns={'pos': 'position'})

    # Drop duplicates to ensure 1:1 mapping if multiple entries exist for same position
    pp2_unique = pp2_small.drop_duplicates(subset='position', keep='first')

    # Map back to original dataframe
    lookup = pp2_unique.set_index('position')
    for col in ['prediction', 'pph2_prob', 'pph2_FPR', 'pph2_TPR']:
        variants_df[col] = variants_df['position'].map(lookup[col])

    return variants_df
