import argparse
import os
import subprocess
from pathlib import Path

import pandas as pd


SCRIPT_PATH = Path(__file__).resolve()
TGVR_ROOT = SCRIPT_PATH.parents[2]

DEFAULT_REMOTE_ROOT = os.environ.get("TGVR_PALMETTO_REMOTE_ROOT", "/home/YOUR_USERNAME/tgvr")
DEFAULT_PALMETTO_TARGET = os.environ.get(
    "TGVR_PALMETTO_TARGET", "YOUR_USERNAME@slogin.palmetto.clemson.edu"
)

AA_3_TO_1 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Glu": "E",
    "Gln": "Q",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
    "Ter": "*",
}


def stage_data_dir(gene_name):
    return TGVR_ROOT / "00_data" / gene_name.lower() / "01_variant_curation"


def run_ssh(sock_path, target, remote_command):
    return subprocess.run(
        ["ssh", "-S", sock_path, target, remote_command],
        text=True,
        capture_output=True,
        check=True,
    )


def convert_three_letter_change(change):
    if not isinstance(change, str):
        return None
    import re
    match = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)", change)
    if not match:
        return None
    wild, position, mutant = match.groups()
    wild_one = AA_3_TO_1.get(wild)
    mutant_one = AA_3_TO_1.get(mutant)
    if not wild_one or not mutant_one:
        return None
    return f"{wild_one}{position}{mutant_one}"


def build_parser():
    parser = argparse.ArgumentParser(
        description="Fetch a completed Palmetto raw 1KGP run back into local TGVR and build the cached workbook used by the variant-curation pipeline."
    )
    parser.add_argument("gene", help="Currently only CDKL5 is supported")
    parser.add_argument(
        "--remote-root",
        default=DEFAULT_REMOTE_ROOT,
        help=f"Remote TGVR root on Palmetto (default: {DEFAULT_REMOTE_ROOT})",
    )
    parser.add_argument(
        "--remote-run-dir",
        default=None,
        help="Explicit remote run directory to fetch. If omitted, the latest grch38_allvar_noid_* run under the remote TGVR root is used.",
    )
    parser.add_argument(
        "--palmetto-target",
        default=DEFAULT_PALMETTO_TARGET,
        help=f"SSH target used with the Palmetto bridge (default: {DEFAULT_PALMETTO_TARGET})",
    )
    parser.add_argument(
        "--socket",
        default=os.path.expanduser("~/.ssh/palmetto.sock"),
        help="Path to the authenticated Palmetto bridge control socket",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    gene_upper = args.gene.upper()
    if gene_upper != "CDKL5":
        raise SystemExit("fetch_1kgp_palmetto.py currently supports only CDKL5")

    if args.remote_run_dir:
        remote_run_dir = args.remote_run_dir
    else:
        latest_cmd = (
            f"ls -1dt {args.remote_root}/outputs/{gene_upper.lower()}/01_variant_curation/1kgp/grch38_allvar_noid_* "
            "2>/dev/null | head -n 1"
        )
        result = run_ssh(args.socket, args.palmetto_target, latest_cmd)
        remote_run_dir = result.stdout.strip()
        if not remote_run_dir:
            raise SystemExit("No remote 1KGP Palmetto run directory was found.")

    remote_tsv = f"{remote_run_dir}/cdkl5.GRCh38.all_variants_noid.tsv"
    remote_xlsx = f"{remote_run_dir}/cdkl5.GRCh38.all_variants_noid.xlsx"
    remote_log_dir = f"{remote_run_dir}/logs"

    local_stage_dir = stage_data_dir(gene_upper)
    local_palmetto_dir = local_stage_dir / "1kgp" / "palmetto_runs" / Path(remote_run_dir).name
    local_palmetto_dir.mkdir(parents=True, exist_ok=True)

    tsv_text = run_ssh(args.socket, args.palmetto_target, f"cat {remote_tsv}").stdout
    (local_palmetto_dir / "cdkl5.GRCh38.all_variants_noid.tsv").write_text(tsv_text, encoding="utf-8")

    xlsx_bytes = subprocess.run(
        ["ssh", "-S", args.socket, args.palmetto_target, f"cat {remote_xlsx}"],
        capture_output=True,
        check=True,
    ).stdout
    (local_palmetto_dir / "cdkl5.GRCh38.all_variants_noid.xlsx").write_bytes(xlsx_bytes)

    logs_listing = run_ssh(args.socket, args.palmetto_target, f"ls -1 {remote_log_dir} 2>/dev/null || true").stdout
    if logs_listing.strip():
        for log_name in [line.strip() for line in logs_listing.splitlines() if line.strip()]:
            log_text = run_ssh(args.socket, args.palmetto_target, f"cat {remote_log_dir}/{log_name}").stdout
            (local_palmetto_dir / log_name).write_text(log_text, encoding="utf-8")

    raw_df = pd.read_csv(local_palmetto_dir / "cdkl5.GRCh38.all_variants_noid.tsv", sep="\t")
    missense_only_df = raw_df[raw_df["Consequence"] == "missense_variant"].copy()
    missense_gene_only_df = missense_only_df[missense_only_df["Gene"].astype(str).str.upper() == gene_upper].copy()
    missense_unique_df = missense_gene_only_df.drop_duplicates().copy()
    missense_unique_prot_df = missense_unique_df.copy()
    missense_unique_prot_df["Protein change 3L"] = missense_unique_prot_df["HGVSp"].astype(str).str.extract(r":p\.(.+)")
    missense_unique_prot_df["Protein change"] = missense_unique_prot_df["Protein change 3L"].apply(convert_three_letter_change)
    missense_unique_prot_df = missense_unique_prot_df[missense_unique_prot_df["Protein change"].notna()].copy()

    workbook_path = local_stage_dir / "1kgp" / f"1kgp_{gene_upper.lower()}_grch38.xlsx"
    workbook_path.parent.mkdir(parents=True, exist_ok=True)
    with pd.ExcelWriter(workbook_path) as writer:
        raw_df.to_excel(writer, sheet_name="Sheet1", index=False)
        missense_only_df.to_excel(writer, sheet_name="missense_only", index=False)
        missense_gene_only_df.to_excel(writer, sheet_name="missense_cdkl5_only", index=False)
        missense_unique_df.to_excel(writer, sheet_name="missense_cdkl5_unique", index=False)
        missense_unique_prot_df.to_excel(writer, sheet_name="missense_cdkl5_unique_prot_chan", index=False)

    print(f"Remote run dir: {remote_run_dir}")
    print(f"Local fetched run dir: {local_palmetto_dir}")
    print(f"Updated cached workbook: {workbook_path}")
    print(f"Counts: raw={len(raw_df)} missense={len(missense_only_df)} gene_only={len(missense_gene_only_df)} unique={len(missense_unique_df)} prot={len(missense_unique_prot_df)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
