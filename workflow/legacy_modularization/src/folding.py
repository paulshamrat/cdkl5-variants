import json
import os
import re
from pathlib import Path

import pandas as pd

from src.data_utils import get_absolute_path, load_config


def load_folding_input_dataframe(config):
    folding_cfg = config["folding"]
    input_workbook = get_absolute_path(folding_cfg["input_workbook"])
    return pd.read_excel(input_workbook), input_workbook


def ensure_directory(path):
    os.makedirs(path, exist_ok=True)
    return path


def write_lines(path, lines):
    ensure_directory(os.path.dirname(path))
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))
        if lines:
            handle.write("\n")


def prepare_folding_inputs(config):
    df, input_workbook = load_folding_input_dataframe(config)
    folding_cfg = config["folding"]
    output_root = get_absolute_path(folding_cfg["prepared_inputs_dir"])
    chain_id = folding_cfg.get("structure_chain", "A")

    valid_mutation_df = df.dropna(subset=["wild", "position", "mutant"]).copy()
    valid_mutation_df["wild"] = valid_mutation_df["wild"].astype(str).str.upper()
    valid_mutation_df["mutant"] = valid_mutation_df["mutant"].astype(str).str.upper()
    valid_mutation_df["position"] = valid_mutation_df["position"].astype(int)

    mutation_strings = (
        valid_mutation_df["mutation"].dropna().astype(str).drop_duplicates().tolist()
    )
    mutation_triplets = [
        f"{row.wild} {row.position} {row.mutant}"
        for row in valid_mutation_df.itertuples(index=False)
    ]
    chain_mutations = [f"{chain_id} {mutation}" for mutation in mutation_strings]

    files = {
        "saafecseq_mutation_list": os.path.join(output_root, "01_saafecseq", "mutation_list.txt"),
        "inps_seq_mutation_list": os.path.join(output_root, "03_inps_seq", "cdkl5_mutations_only.txt"),
        "ddgemb_seq_mutation_list": os.path.join(output_root, "06_ddgemb_seq", "cdkl5_mutations_only_ddgemb.txt"),
        "mcsm_str_mutation_list": os.path.join(output_root, "07_mcsm_str", "cdkl5_mutations_formatted.txt"),
        "ddmut_str_mutation_list": os.path.join(output_root, "08_ddmut_str", "cdkl5_mutations_formatted.txt"),
    }

    write_lines(files["saafecseq_mutation_list"], mutation_triplets)
    write_lines(files["inps_seq_mutation_list"], mutation_strings)
    write_lines(files["ddgemb_seq_mutation_list"], mutation_strings)
    write_lines(files["mcsm_str_mutation_list"], chain_mutations)
    write_lines(files["ddmut_str_mutation_list"], chain_mutations)

    manual_note_path = os.path.join(output_root, "MANUAL_INPUTS_REQUIRED.md")
    manual_note = "\n".join(
        [
            "# Folding Manual Inputs",
            "",
            "This archived modularization reproduces the input-preparation and analysis logic from `02_folding.ipynb`.",
            "",
            "What is automated now:",
            "- creation of mutation-list files for the legacy folding predictors",
            "- summary analysis of sequence/structure DDG columns already present in the legacy workbook",
            "",
            "What remains external or manual from the original notebook:",
            "- SAAFEC-SEQ execution on Palmetto or another suitable environment",
            "- I-Mutant sequence/structure submissions",
            "- INPS sequence and structure submissions",
            "- DDGun sequence and structure submissions",
            "- DDGEmb sequence submissions",
            "- mCSM structure submissions",
            "- DDMut structure submissions",
            "",
            "If those predictor outputs are collected later, they should be stored in workflow-local folders under:",
            f"- {Path(output_root).as_posix()}",
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


def collect_ddg_summary(df, cols, region_lo, region_hi, region_label, method_type, classes):
    rows = []
    sub = df[df["position"].between(region_lo, region_hi) & df["Germline classification"].isin(classes)].copy()
    for col in cols:
        for cls in classes:
            vals = pd.to_numeric(
                sub.loc[sub["Germline classification"] == cls, col],
                errors="coerce",
            ).dropna()
            if len(vals) == 0:
                continue
            rows.append(
                {
                    "Method": col,
                    "Type": method_type,
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


def analyze_folding_results(config):
    df, input_workbook = load_folding_input_dataframe(config)
    folding_cfg = config["folding"]
    output_dir = get_absolute_path(folding_cfg["analysis_dir"])
    ensure_directory(output_dir)

    ddg_seq_cols = [col for col in df.columns if re.match(r"(?i)^ddg.*seq$", str(col))]
    ddg_str_cols = [col for col in df.columns if re.match(r"(?i)^ddg.*str$", str(col))]
    classes = folding_cfg.get("analysis_classes", ["Benign", "Pathogenic"])
    full_range = folding_cfg["analysis_ranges"]["full"]
    kinase_range = folding_cfg["analysis_ranges"]["kinase"]

    all_rows = []
    all_rows += collect_ddg_summary(df, ddg_seq_cols, full_range[0], full_range[1], "Full", "Sequence", classes)
    all_rows += collect_ddg_summary(df, ddg_seq_cols, kinase_range[0], kinase_range[1], "Kinase", "Sequence", classes)
    all_rows += collect_ddg_summary(df, ddg_str_cols, full_range[0], full_range[1], "Full", "Structure", classes)
    all_rows += collect_ddg_summary(df, ddg_str_cols, kinase_range[0], kinase_range[1], "Kinase", "Structure", classes)
    summary_df = pd.DataFrame(all_rows)

    summary_excel = os.path.join(output_dir, "cdkl5_folding_ddg_summary_table.xlsx")
    summary_csv = os.path.join(output_dir, "cdkl5_folding_ddg_summary_table.csv")
    summary_df.to_excel(summary_excel, index=False)
    summary_df.to_csv(summary_csv, index=False)

    metadata = {
        "input_workbook": input_workbook,
        "analysis_dir": output_dir,
        "ddg_seq_columns": ddg_seq_cols,
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


def run_folding_workflow(step="all"):
    config = load_config()
    results = {}
    if step in {"prepare", "all"}:
        results["prepare"] = prepare_folding_inputs(config)
    if step in {"analyze", "all"}:
        results["analyze"] = analyze_folding_results(config)
    return results
