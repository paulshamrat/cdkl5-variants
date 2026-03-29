import json
import os
import re
from pathlib import Path

import pandas as pd

from src.data_utils import get_absolute_path, load_config


def ensure_directory(path):
    os.makedirs(path, exist_ok=True)
    return path


def write_lines(path, lines):
    ensure_directory(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
        if lines:
            handle.write("\n")


def load_binding_input_dataframe(config):
    binding_cfg = config["binding"]
    input_workbook = get_absolute_path(binding_cfg["input_workbook"])
    return pd.read_excel(input_workbook), input_workbook


def load_binding_analysis_dataframe(config):
    binding_cfg = config["binding"]
    analysis_workbook = get_absolute_path(binding_cfg["analysis_workbook"])
    return pd.read_excel(analysis_workbook), analysis_workbook


def prepare_binding_inputs(config):
    df, input_workbook = load_binding_input_dataframe(config)
    binding_cfg = config["binding"]
    output_root = get_absolute_path(binding_cfg["prepared_inputs_dir"])
    chain_id = binding_cfg.get("structure_chain", "A")

    valid_mutation_df = df.dropna(subset=["wild", "position", "mutant"]).copy()
    valid_mutation_df["wild"] = valid_mutation_df["wild"].astype(str).str.upper()
    valid_mutation_df["mutant"] = valid_mutation_df["mutant"].astype(str).str.upper()
    valid_mutation_df["position"] = valid_mutation_df["position"].astype(int)
    valid_mutation_df["mutation"] = valid_mutation_df["mutation"].astype(str)

    mutation_strings = (
        valid_mutation_df["mutation"].dropna().astype(str).drop_duplicates().tolist()
    )
    saambe3d_lines = [
        f"{chain_id} {row.position} {row.wild} {row.mutant}"
        for row in valid_mutation_df.itertuples(index=False)
    ]
    chain_mutations = [f"{chain_id} {mutation}" for mutation in mutation_strings]
    foldx_lines = [
        f"{row.wild}{chain_id}{row.position}{row.mutant};"
        for row in valid_mutation_df.itertuples(index=False)
    ]

    files = {
        "saambe3d_mutation_list": os.path.join(output_root, "02_saambe3d", "mutations_list.txt"),
        "mcsmppi_mutation_list": os.path.join(output_root, "05_mcsmppi", "mutations_for_mcsmppi2.txt"),
        "ddmutppi_mutation_list": os.path.join(output_root, "07_ddmutppi", "clinvar_1kgp_hector_gaf_final_mutlist.txt"),
        "foldx_individual_list": os.path.join(output_root, "03_foldx", "individual_list_template.txt"),
    }

    write_lines(files["saambe3d_mutation_list"], saambe3d_lines)
    write_lines(files["mcsmppi_mutation_list"], chain_mutations)
    write_lines(files["ddmutppi_mutation_list"], chain_mutations)
    write_lines(files["foldx_individual_list"], foldx_lines)

    manual_note_path = os.path.join(output_root, "MANUAL_INPUTS_REQUIRED.md")
    manual_note = "\n".join(
        [
            "# Binding Manual Inputs",
            "",
            "This archived modularization reproduces the input-preparation and analysis logic from `04_binding.ipynb`.",
            "",
            "What is automated now:",
            "- creation of mutation-list files for SAAMBE-3D, FoldX, mCSM-PPI2, and DDMutPPI",
            "- summary analysis of binding DDG columns already present in the legacy final binding workbook",
            "",
            "What remains external or manual from the original notebook:",
            "- SAAMBE-3D execution on Palmetto or another suitable environment",
            "- mCSM-PPI2 server submissions",
            "- FoldX complex modeling and AnalyseComplex execution",
            "- DDMutPPI server submissions",
            "- any additional partner-specific result collection",
            "",
            "Legacy final binding workbook used for analysis:",
            f"- {Path(load_binding_analysis_dataframe(config)[1]).as_posix()}",
        ]
    )
    write_lines(manual_note_path, manual_note.splitlines())

    summary = {
        "input_workbook": input_workbook,
        "prepared_inputs_dir": output_root,
        "mutation_count": int(len(mutation_strings)),
        "files": files,
        "manual_note": manual_note_path,
    }
    summary_path = os.path.join(output_root, "summary.json")
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return summary


def parse_binding_column(column_name):
    match = re.match(r"^ddg_(?P<partner>.+)_str_(?P<method>[^_]+)$", str(column_name))
    if match:
        return match.group("partner"), match.group("method")
    return str(column_name), "unknown"


def collect_binding_summary(df, cols, region_lo, region_hi, region_label, classes):
    rows = []
    sub = df[df["position"].between(region_lo, region_hi) & df["Germline classification"].isin(classes)].copy()
    for col in cols:
        partner, method = parse_binding_column(col)
        for cls in classes:
            vals = pd.to_numeric(
                sub.loc[sub["Germline classification"] == cls, col],
                errors="coerce",
            ).dropna()
            if len(vals) == 0:
                continue
            rows.append(
                {
                    "Partner": partner,
                    "Method": method,
                    "Column": col,
                    "Region": region_label,
                    "Class": cls,
                    "n": int(len(vals)),
                    "Mean": round(float(vals.mean()), 3),
                    "Median": round(float(vals.median()), 3),
                    "Std": round(float(vals.std()), 3) if len(vals) > 1 else 0.0,
                    "Min": round(float(vals.min()), 3),
                    "Max": round(float(vals.max()), 3),
                }
            )
    return rows


def analyze_binding_results(config):
    df, analysis_workbook = load_binding_analysis_dataframe(config)
    binding_cfg = config["binding"]
    output_dir = get_absolute_path(binding_cfg["analysis_dir"])
    ensure_directory(output_dir)

    ddg_str_cols = [col for col in df.columns if re.match(r"^ddg_.*_str(?:_.+)?$", str(col))]
    classes = binding_cfg.get("analysis_classes", ["Benign", "Pathogenic"])
    full_range = binding_cfg["analysis_ranges"]["full"]
    kinase_range = binding_cfg["analysis_ranges"]["kinase"]

    all_rows = []
    all_rows += collect_binding_summary(df, ddg_str_cols, full_range[0], full_range[1], "Full", classes)
    all_rows += collect_binding_summary(df, ddg_str_cols, kinase_range[0], kinase_range[1], "Kinase", classes)
    summary_df = pd.DataFrame(all_rows)

    summary_excel = os.path.join(output_dir, "cdkl5_binding_ddg_summary_table.xlsx")
    summary_csv = os.path.join(output_dir, "cdkl5_binding_ddg_summary_table.csv")
    summary_df.to_excel(summary_excel, index=False)
    summary_df.to_csv(summary_csv, index=False)

    metadata = {
        "analysis_workbook": analysis_workbook,
        "analysis_dir": output_dir,
        "ddg_str_columns": ddg_str_cols,
        "classes": classes,
        "ranges": {
            "full": full_range,
            "kinase": kinase_range,
        },
        "rows": int(len(summary_df)),
        "files": {
            "summary_excel": summary_excel,
            "summary_csv": summary_csv,
        },
    }
    metadata_path = os.path.join(output_dir, "summary.json")
    with open(metadata_path, "w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return metadata


def run_binding_workflow(step="all"):
    config = load_config()
    results = {}
    if step in {"prepare", "all"}:
        results["prepare"] = prepare_binding_inputs(config)
    if step in {"analyze", "all"}:
        results["analyze"] = analyze_binding_results(config)
    return results
