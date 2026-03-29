import os
import yaml
import pandas as pd

def get_workflow_root():
    """Returns the absolute path to the workflow/ directory."""
    return os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def load_config(config_path=None):
    """Loads the project configuration from a YAML file."""
    if config_path is None:
        config_path = os.path.join(get_workflow_root(), "config", "config.yaml")
    
    with open(config_path, 'r') as f:
        return yaml.safe_load(f)

def get_absolute_path(relative_path, base_dir=None):
    """Converts a workflow-relative path to an absolute path."""
    if base_dir is None:
        base_dir = get_workflow_root()
    return os.path.abspath(os.path.join(base_dir, relative_path))

def parse_mutation(mutation_str):
    """
    Parses a mutation string like 'F13S' into (wild, position, mutant).
    Example: 'F13S' -> ('F', 13, 'S')
    """
    if not isinstance(mutation_str, str) or len(mutation_str) < 3:
        return None, None, None
    
    wild = mutation_str[0]
    mutant = mutation_str[-1]
    try:
        position = int(mutation_str[1:-1])
        return wild, position, mutant
    except ValueError:
        return None, None, None

def add_mutation_components(df, mutation_col=None):
    """Adds 'wild', 'position', and 'mutant' columns to a DataFrame."""
    if mutation_col is None:
        # Try common column names
        for col in ['Mutation', 'mutation', 'Protein change', 'Protein Change']:
            if col in df.columns:
                mutation_col = col
                break
    
    if mutation_col is None or mutation_col not in df.columns:
        # Add empty columns if no source found to avoid KeyErrors later
        for col in ['wild', 'position', 'mutant']:
            if col not in df.columns:
                df[col] = None
        return df

    parsed = df[mutation_col].apply(parse_mutation)
    df['wild'] = parsed.apply(lambda x: x[0] if x else None)
    df['position'] = parsed.apply(lambda x: x[1] if x else None)
    df['mutant'] = parsed.apply(lambda x: x[2] if x else None)
    
    # Also ensure 'Mutation' column exists for downstream compatibility if we found it under another name
    if 'Mutation' not in df.columns and mutation_col != 'Mutation':
         df['Mutation'] = df[mutation_col]

    return df

def save_to_excel(df, file_path, sheet_name='Sheet1'):
    """Utility to save a DataFrame to Excel, creating parent directories if needed."""
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    df.to_excel(file_path, index=False)
