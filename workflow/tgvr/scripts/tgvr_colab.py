import argparse
import json
import re
import shutil
import sys
from pathlib import Path

import pandas as pd
import requests


SCRIPT_PATH = Path(__file__).resolve()


DEFAULT_DRIVE_ROOT = "/content/drive/MyDrive/tgvr_colab"
ENSEMBL_REST_BASE = "https://rest.ensembl.org"
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


def build_drive_workspace(drive_root):
    drive_root = Path(drive_root).expanduser()
    paths = {
        "root": drive_root,
        "data": drive_root / "data",
        "outputs": drive_root / "outputs",
        "logs": drive_root / "logs",
        "cache": drive_root / "cache",
        "scripts": drive_root / "scripts",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def colab_stage_output_dir(drive_root, gene_name):
    return (
        Path(drive_root).expanduser()
        / "outputs"
        / gene_name.lower()
        / "01_variant_curation"
        / "1kgp"
    )


def convert_three_letter_change(change):
    if not isinstance(change, str):
        return None
    match = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)", change)
    if not match:
        return None
    wild, position, mutant = match.groups()
    if wild not in AA_3_TO_1 or mutant not in AA_3_TO_1:
        return None
    return f"{AA_3_TO_1[wild]}{position}{AA_3_TO_1[mutant]}"


def choose_transcript_consequence(consequences, gene_symbol, canonical_transcript):
    canonical_root = canonical_transcript.split(".", 1)[0] if canonical_transcript else None
    candidates = []
    for consequence in consequences or []:
        if (consequence.get("gene_symbol") or "").upper() != gene_symbol.upper():
            continue
        if not str(consequence.get("biotype", "")).startswith("protein_coding"):
            continue
        transcript_id = consequence.get("transcript_id") or ""
        transcript_root = transcript_id.split(".", 1)[0] if transcript_id else ""
        rank = 1
        if canonical_root and transcript_root == canonical_root:
            rank = 0
        has_hgvsp_penalty = 0 if consequence.get("hgvsp") else 1
        flags_penalty = len(consequence.get("flags") or [])
        candidates.append((rank, has_hgvsp_penalty, flags_penalty, transcript_id, consequence))

    if not candidates:
        return None

    candidates.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    return candidates[0][4]


def extract_1kgp_af(vep_entry, alt_allele):
    for colocated in vep_entry.get("colocated_variants", []):
        frequencies = colocated.get("frequencies", {})
        alt_data = frequencies.get(alt_allele)
        if isinstance(alt_data, dict) and alt_data.get("af") is not None:
            return alt_data.get("af")
    return None


def build_colab_1kgp_workbook(gene_name, drive_root, logger=None):
    gene_upper = gene_name.upper()
    stage_dir = colab_stage_output_dir(drive_root, gene_upper)
    stage_dir.mkdir(parents=True, exist_ok=True)

    session = requests.Session()
    session.headers.update({"Content-Type": "application/json", "Accept": "application/json"})

    def log(message):
        if logger:
            logger(message)

    def get_json(url):
        response = session.get(url, timeout=120)
        response.raise_for_status()
        return response.json()

    def post_json(url, payload):
        response = session.post(url, json=payload, timeout=120)
        response.raise_for_status()
        return response.json()

    log(f"Looking up {gene_upper} in Ensembl...")
    lookup_url = f"{ENSEMBL_REST_BASE}/lookup/symbol/homo_sapiens/{gene_upper}?content-type=application/json"
    lookup_payload = get_json(lookup_url)

    chromosome = lookup_payload["seq_region_name"]
    start = lookup_payload["start"]
    end = lookup_payload["end"]
    canonical_transcript = lookup_payload.get("canonical_transcript")

    log(f"Fetching 1KGP overlap variants for {gene_upper} in {chromosome}:{start}-{end}...")
    region_url = (
        f"{ENSEMBL_REST_BASE}/overlap/region/human/"
        f"{chromosome}:{start}-{end}?feature=variation;variant_set=1kg_3;content-type=application/json"
    )
    region_payload = get_json(region_url)

    gene_rows = []
    for row in region_payload:
        alleles = row.get("alleles") or []
        gene_rows.append(
            {
                "CHROM": f"chr{row.get('seq_region_name')}",
                "POS": row.get("start"),
                "REF": alleles[0] if len(alleles) >= 1 else None,
                "ALT": alleles[1] if len(alleles) >= 2 else None,
                "RefGenome": row.get("assembly_name"),
                "AF": pd.NA,
                "Gene": gene_upper,
                "Consequence": row.get("consequence_type"),
                "HGVSc": pd.NA,
                "HGVSp": pd.NA,
                "dbSNP ID": row.get("id"),
            }
        )
    gene_rows_df = pd.DataFrame(gene_rows)

    missense_region_rows = [
        row
        for row in region_payload
        if row.get("consequence_type") == "missense_variant" and row.get("id")
    ]
    missense_ids = [row["id"] for row in missense_region_rows]

    vep_payload = []
    if missense_ids:
        log(f"Annotating {len(missense_ids)} missense rsIDs through Ensembl VEP...")
        for index in range(0, len(missense_ids), 200):
            chunk = missense_ids[index : index + 200]
            vep_payload.extend(
                post_json(f"{ENSEMBL_REST_BASE}/vep/human/id?hgvs=1", {"ids": chunk})
            )

    vep_by_id = {entry["input"]: entry for entry in vep_payload if "input" in entry}

    missense_rows = []
    for region_row in missense_region_rows:
        rsid = region_row.get("id")
        vep_entry = vep_by_id.get(rsid, {})
        transcript = choose_transcript_consequence(
            vep_entry.get("transcript_consequences", []),
            gene_upper,
            canonical_transcript,
        )
        if transcript is None:
            continue
        allele_string = vep_entry.get("allele_string") or "/".join(region_row.get("alleles", []))
        ref = allele_string.split("/")[0] if allele_string else None
        alt = transcript.get("variant_allele")
        missense_rows.append(
            {
                "CHROM": f"chr{region_row.get('seq_region_name')}",
                "POS": region_row.get("start"),
                "REF": ref,
                "ALT": alt,
                "RefGenome": region_row.get("assembly_name"),
                "AF": extract_1kgp_af(vep_entry, alt),
                "Gene": transcript.get("gene_symbol"),
                "Consequence": region_row.get("consequence_type"),
                "HGVSc": transcript.get("hgvsc"),
                "HGVSp": transcript.get("hgvsp"),
                "dbSNP ID": rsid,
            }
        )

    columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "RefGenome",
        "AF",
        "Gene",
        "Consequence",
        "HGVSc",
        "HGVSp",
        "dbSNP ID",
    ]
    missense_only_df = pd.DataFrame(missense_rows)
    if missense_only_df.empty:
        missense_only_df = pd.DataFrame(columns=columns)

    missense_gene_only_df = missense_only_df[
        missense_only_df["Gene"].astype(str).str.upper() == gene_upper
    ].copy()
    missense_unique_df = missense_gene_only_df.drop_duplicates().copy()
    missense_unique_prot_df = missense_unique_df.copy()
    missense_unique_prot_df["Protein change 3L"] = missense_unique_prot_df["HGVSp"].astype(str).str.extract(
        r":p\.(.+)"
    )
    missense_unique_prot_df["Protein change"] = missense_unique_prot_df["Protein change 3L"].apply(
        convert_three_letter_change
    )
    missense_unique_prot_df = missense_unique_prot_df[
        missense_unique_prot_df["Protein change"].notna()
    ].copy()

    workbook_path = stage_dir / f"1kgp_{gene_upper.lower()}_grch38.xlsx"
    summary_path = stage_dir / f"1kgp_{gene_upper.lower()}_grch38.summary.json"

    with pd.ExcelWriter(workbook_path) as writer:
        gene_rows_df.to_excel(writer, sheet_name="Sheet1", index=False)
        missense_only_df.to_excel(writer, sheet_name="missense_only", index=False)
        missense_gene_only_df.to_excel(writer, sheet_name="missense_cdkl5_only", index=False)
        missense_unique_df.to_excel(writer, sheet_name="missense_cdkl5_unique", index=False)
        missense_unique_prot_df.to_excel(
            writer, sheet_name="missense_cdkl5_unique_prot_chan", index=False
        )

    summary = {
        "gene": gene_upper,
        "assembly": "GRCh38",
        "ensembl_gene_id": lookup_payload.get("id"),
        "canonical_transcript": canonical_transcript,
        "region": f"{chromosome}:{start}-{end}",
        "counts": {
            "gene_rows": int(len(gene_rows_df)),
            "missense_only": int(len(missense_only_df)),
            "missense_gene_only": int(len(missense_gene_only_df)),
            "missense_gene_unique": int(len(missense_unique_df)),
            "missense_gene_unique_prot_chan": int(len(missense_unique_prot_df)),
        },
        "files": {
            "xlsx": str(workbook_path),
            "summary_json": str(summary_path),
        },
    }
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def ensure_drive_is_mounted(drive_root):
    drive_root = Path(drive_root).expanduser()
    mydrive = drive_root.parent
    if not mydrive.exists():
        raise SystemExit(
            "Google Drive does not appear to be mounted. In a Colab notebook cell, run:\n"
            "from google.colab import drive\n"
            "drive.mount('/content/drive')"
        )
    return drive_root


def copy_colab_files(drive_root):
    drive_root = ensure_drive_is_mounted(drive_root)
    build_drive_workspace(drive_root)

    destination = drive_root / "scripts" / "tgvr_colab.py"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SCRIPT_PATH, destination)

    return {
        "root": str(drive_root),
        "scripts_dir": str((drive_root / "scripts").resolve()),
        "runner": str(destination.resolve()),
    }


def build_parser():
    parser = argparse.ArgumentParser(
        description="Drive-backed Colab helper for TGVR."
    )
    parser.add_argument(
        "--drive-root",
        default=DEFAULT_DRIVE_ROOT,
        help=f"Drive-backed TGVR workspace root (default: {DEFAULT_DRIVE_ROOT})",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser(
        "init",
        help="Create the Drive workspace and copy the Colab TGVR scripts into it.",
    )

    run_parser = subparsers.add_parser(
        "run-1kgp",
        help="Build a Colab-safe Drive-backed 1KGP workbook without Palmetto.",
    )
    run_parser.add_argument("gene", help="Gene symbol, for example CDKL5")

    return parser


def log(message):
    print(message, flush=True)


def main():
    parser = build_parser()
    args = parser.parse_args()

    drive_root = ensure_drive_is_mounted(args.drive_root)

    if args.command == "init":
        copied = copy_colab_files(drive_root)
        print("Created Drive-backed TGVR workspace and copied Colab helpers:")
        print(json.dumps(copied, indent=2))
        print("Next command:")
        print(f"  cd {drive_root}")
        print(f"  python3 scripts/tgvr_colab.py --drive-root {drive_root} run-1kgp CDKL5")
        return 0

    if args.command == "run-1kgp":
        build_drive_workspace(drive_root)
        summary = build_colab_1kgp_workbook(args.gene, drive_root, logger=log)
        print("Done.")
        print(json.dumps(summary, indent=2))
        print("Workbook:", Path(summary["files"]["xlsx"]).resolve())
        print("Summary:", Path(summary["files"]["summary_json"]).resolve())
        return 0

    raise SystemExit(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
