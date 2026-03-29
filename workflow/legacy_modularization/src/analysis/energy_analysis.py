import pandas as pd
import numpy as np

def calculate_max_ddg(df, str_cols):
    """Calculates the maximum absolute DDG across specified structural predictors."""
    df['ddG_Fmax'] = df[str_cols].abs().max(axis=1)
    return df

def get_variant_data_dicts(df, str_cols):
    """ Splits dataframe into benign and pathogenic dictionaries for plotting. """
    methods = [c.replace('_str', '') for c in str_cols]
    
    ben_df = df[df['Germline classification'] == 'Benign']
    path_df = df[df['Germline classification'] == 'Pathogenic']
    
    def to_dict(data_df):
        res = {}
        for _, r in data_df.iterrows():
            res[r['mutation']] = {
                'pos': r['position'],
                **{c.replace('_str',''): abs(r[c]) for c in str_cols},
                'max': r['ddG_Fmax']
            }
        return res
    
    return to_dict(ben_df), to_dict(path_df)

def calculate_reclassification_threshold(df):
    """Calculates the threshold for reclassification based on max benign and min pathogenic max DDG."""
    benign_max = df[df['Germline classification'] == 'Benign']['ddG_Fmax'].max()
    pathogenic_min = df[df['Germline classification'] == 'Pathogenic']['ddG_Fmax'].min()
    
    if pd.isna(benign_max) or pd.isna(pathogenic_min):
        return None
        
    return (benign_max + pathogenic_min) / 2
