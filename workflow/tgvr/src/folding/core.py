import json
import re
import traceback
from datetime import datetime
from pathlib import Path

import pandas as pd
import yaml


TGVR_ROOT = Path(__file__).resolve().parents[2]
REPO_ROOT = TGVR_ROOT.parents[1]

DDGUN_METHODS = {
    "ddgun_seq": {
        "description": "DDGun sequence mode",
        "input_filename": "mutations.txt",
    },
    "ddgun_str": {
        "description": "DDGun structure mode",
        "input_filename": "mutations.txt",
    },
}


class RunLogger:
    def __init__(self, gene_name, stage_name):
        logs_dir = folding_output_dir(gene_name) / "logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_log_path = logs_dir / f"{timestamp}_{stage_name}.log"
        self.latest_log_path = logs_dir / "latest.log"
        self._run_handle = self.run_log_path.open("w", encoding="utf-8")
        self._latest_handle = self.latest_log_path.open("w", encoding="utf-8")

    def log(self, message):
        print(message, flush=True)
        self._run_handle.write(message + "\n")
        self._latest_handle.write(message + "\n")
        self._run_handle.flush()
        self._latest_handle.flush()

    def close(self):
        self._run_handle.close()
        self._latest_handle.close()


def display_path(path):
    path = Path(path).resolve()
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def gene_config_path(gene_name):
    return TGVR_ROOT / "config" / "genes" / f"{gene_name.upper()}.yaml"


def load_gene_config(gene_name):
    path = gene_config_path(gene_name)
    if not path.exists():
        raise FileNotFoundError(
            f"No TGVR gene config found for {gene_name.upper()}: {display_path(path)}"
        )
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def folding_data_dir(gene_name):
    return TGVR_ROOT / "00_data" / gene_name.lower() / "02_folding"


def folding_output_dir(gene_name):
    return TGVR_ROOT / "outputs" / gene_name.lower() / "02_folding"


def variant_master_path(gene_name):
    return (
        TGVR_ROOT
        / "outputs"
        / gene_name.lower()
        / "01_variant_curation"
        / "master"
        / "10_master_dataset_final.csv"
    )


def ensure_folding_dirs(gene_name):
    data_root = folding_data_dir(gene_name)
    output_root = folding_output_dir(gene_name)
    for subdir in ["ddgun", "logs"]:
        (output_root / subdir).mkdir(parents=True, exist_ok=True)
    for subdir in ["ddgun", "reference", "results", "env"]:
        (data_root / subdir).mkdir(parents=True, exist_ok=True)
    return data_root, output_root


def load_master_dataframe(gene_name):
    path = variant_master_path(gene_name)
    if not path.exists():
        raise FileNotFoundError(
            f"Variant-curation master dataset not found for {gene_name.upper()}: {display_path(path)}. "
            "Run 01_variant_curation first."
        )
    df = pd.read_csv(path)
    required = {"Mutation", "wild", "position", "mutant", "Germline classification"}
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Master dataset is missing required columns: {', '.join(missing)}")
    return df, path


def parse_range_text(range_text):
    match = re.fullmatch(r"(\d+)-(\d+)", str(range_text).strip())
    if not match:
        raise ValueError("Range must be formatted as START-END, for example 13-297.")
    start, end = int(match.group(1)), int(match.group(2))
    if start < 1 or end < start:
        raise ValueError("Range must satisfy 1 <= START <= END.")
    return start, end


def write_lines(path, lines):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    content = "\n".join(lines)
    if content:
        content += "\n"
    path.write_text(content, encoding="utf-8")


def ddgun_prepare(gene_name, uniprot_id, logger, method, range_text=None, label=None):
    if method not in DDGUN_METHODS:
        raise ValueError(f"Unsupported folding method for fresh-start TGVR 02_folding: {method}")

    gene_upper = gene_name.upper()
    gene_config = load_gene_config(gene_upper)
    df, master_path = load_master_dataframe(gene_name)
    ensure_folding_dirs(gene_name)

    scope_label = "Full length"
    if range_text is not None:
        start, end = parse_range_text(range_text)
        scope_label = label or f"Residues {start}-{end}"
        df = df[df["position"].between(start, end)].copy()

    valid = df.dropna(subset=["Mutation", "wild", "position", "mutant"]).copy()
    valid["Mutation"] = valid["Mutation"].astype(str)
    valid["wild"] = valid["wild"].astype(str).str.upper()
    valid["mutant"] = valid["mutant"].astype(str).str.upper()
    valid["position"] = valid["position"].astype(int)
    valid = valid.drop_duplicates(subset=["Mutation"]).copy()

    data_root = folding_data_dir(gene_name)
    output_root = folding_output_dir(gene_name)
    method_root = output_root / "ddgun"

    mutation_lines = valid["Mutation"].tolist()
    prepared_input = method_root / DDGUN_METHODS[method]["input_filename"]
    write_lines(prepared_input, mutation_lines)

    reference_fasta = TGVR_ROOT / gene_config["reference"]["fasta"]
    note_path = method_root / "README.md"
    note = "\n".join(
        [
            "# TGVR 02_folding DDGun Start",
            "",
            f"Gene: `{gene_upper}`",
            f"UniProt: `{uniprot_id.upper()}`",
            f"Method: `{method}`",
            f"Scope: `{scope_label}`",
            "",
            "This fresh-start TGVR folding stage currently begins with DDGun only.",
            "",
            "Prepared from:",
            f"- {display_path(master_path)}",
            f"- reference FASTA: {display_path(reference_fasta)}",
            "",
            "Prepared mutation list:",
            f"- {display_path(prepared_input)}",
            "",
            "Place future DDGun outputs under:",
            f"- {display_path(data_root / 'results')}",
            "",
            "Notes:",
            "- `ddgun_seq` will need a sequence FASTA and mutation list.",
            "- `ddgun_str` will need a structure/PDB, chain, and mutation list.",
            "- installation/container/backend setup is not finalized yet; this stage is intentionally minimal.",
        ]
    )
    note_path.write_text(note + "\n", encoding="utf-8")

    summary = {
        "gene": gene_upper,
        "uniprot_id": uniprot_id.upper(),
        "method": method,
        "scope": scope_label,
        "source_master_dataset": str(master_path),
        "reference_fasta": str(reference_fasta),
        "mutation_count": int(len(mutation_lines)),
        "prepared_input": str(prepared_input),
        "note": str(note_path),
    }
    summary_path = method_root / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    logger.log(f"Prepared DDGun input for {len(mutation_lines)} unique mutations in scope: {scope_label}.")
    logger.log(f"Prepared input: {display_path(prepared_input)}")
    logger.log(f"Summary file: {display_path(summary_path)}")
    return summary


def status_summary(gene_name, logger):
    summary = {
        "stage": "02_folding",
        "status": "fresh_start",
        "active_method_family": "DDGun",
        "supported_methods_now": sorted(DDGUN_METHODS.keys()),
        "next_step": "set up a reproducible DDGun backend before scientific runs",
    }
    logger.log("TGVR 02_folding status:")
    logger.log("  - fresh-start stage: yes")
    logger.log("  - active method family: DDGun")
    logger.log("  - supported methods now: ddgun_seq, ddgun_str")
    logger.log("  - next step: reproducible DDGun setup/backend")
    return summary


def run_folding_stage(gene_name, uniprot_id, stage, range_text=None, label=None, method=None):
    logger = RunLogger(gene_name.upper(), stage.replace("-", "_"))
    ensure_folding_dirs(gene_name)
    try:
        logger.log(f"Gene workspace: {display_path(folding_output_dir(gene_name))}")
        logger.log(f"Run log: {display_path(logger.run_log_path)}")

        config = load_gene_config(gene_name.upper())
        expected_uniprot = str(config["reference"].get("uniprot_accession", "")).upper()
        if expected_uniprot and expected_uniprot != uniprot_id.upper():
            raise ValueError(
                f"UniProt ID mismatch for {gene_name.upper()}: expected {expected_uniprot}, received {uniprot_id.upper()}"
            )

        result = {}
        if stage == "status":
            result["status"] = status_summary(gene_name, logger)
        elif stage == "prepare-ddgun":
            if not method:
                raise ValueError("Use --method with prepare-ddgun. Supported: ddgun_seq, ddgun_str")
            logger.log("Preparing DDGun folding input...")
            logger.log("  - fresh-start TGVR 02_folding currently supports DDGun-first preparation only")
            result["prepare_ddgun"] = ddgun_prepare(gene_name, uniprot_id, logger, method=method, range_text=range_text, label=label)
        else:
            raise ValueError(f"Unsupported TGVR 02_folding stage after cleanup: {stage}")
        return result
    except Exception as exc:
        logger.log(f"ERROR: {exc}")
        logger.log(traceback.format_exc().rstrip())
        raise
    finally:
        logger.close()
