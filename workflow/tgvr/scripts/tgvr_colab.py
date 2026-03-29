import argparse
import json
import re
import shutil
import subprocess
import sys
import urllib.request
from pathlib import Path

import pandas as pd
import requests


SCRIPT_PATH = Path(__file__).resolve()


DEFAULT_WORKSPACE_ROOT = "/content/tgvr_colab"
ENSEMBL_REST_BASE = "https://rest.ensembl.org"
CLINVAR_VARIANT_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
DEFAULT_LEGACY_1KGP_COUNTS = {
    "CDKL5": {
        "gene_rows": 4480,
        "missense_only": 19,
        "missense_gene_only": 13,
        "missense_gene_unique": 12,
        "missense_gene_unique_prot_chan": 12,
    }
}
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


def build_workspace(workspace_root):
    workspace_root = Path(workspace_root).expanduser()
    paths = {
        "root": workspace_root,
        "data": workspace_root / "data",
        "outputs": workspace_root / "outputs",
        "logs": workspace_root / "logs",
        "cache": workspace_root / "cache",
        "scripts": workspace_root / "scripts",
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def gene_data_root(workspace_root, gene_name):
    return (
        Path(workspace_root).expanduser()
        / "data"
        / gene_name.lower()
        / "01_variant_curation"
    )


def gene_output_root(workspace_root, gene_name):
    return (
        Path(workspace_root).expanduser()
        / "outputs"
        / gene_name.lower()
        / "01_variant_curation"
    )


def colab_stage_output_dir(workspace_root, gene_name):
    return gene_output_root(workspace_root, gene_name) / "1kgp"


def clinvar_output_dir(workspace_root, gene_name):
    return gene_output_root(workspace_root, gene_name) / "clinvar"


def manual_data_dir(workspace_root, gene_name):
    return gene_data_root(workspace_root, gene_name) / "manual"


def clinvar_data_dir(workspace_root, gene_name):
    return gene_data_root(workspace_root, gene_name) / "clinvar"


def ensure_gene_workspace_layout(workspace_root, gene_name):
    gene_upper = gene_name.upper()
    paths = {
        "data_root": gene_data_root(workspace_root, gene_upper),
        "manual_dir": manual_data_dir(workspace_root, gene_upper),
        "clinvar_data_dir": clinvar_data_dir(workspace_root, gene_upper),
        "clinvar_output_dir": clinvar_output_dir(workspace_root, gene_upper),
        "kgp_output_dir": colab_stage_output_dir(workspace_root, gene_upper),
    }
    for path in paths.values():
        path.mkdir(parents=True, exist_ok=True)
    return paths


def download_file(url, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "codex-cdkl5-variants/1.0"})
    with urllib.request.urlopen(request) as response, open(destination, "wb") as handle:
        shutil.copyfileobj(response, handle)
    return destination


def run_checked(command):
    result = subprocess.run(
        command,
        check=True,
        text=True,
        capture_output=True,
    )
    return result.stdout


def detect_vcf_contig_prefix(vcf_path):
    header = run_checked(["bcftools", "view", "-h", str(vcf_path)])
    if "##contig=<ID=chrX>" in header:
        return "chr"
    if "##contig=<ID=X>" in header:
        return ""
    if "##contig=<ID=chr1>" in header:
        return "chr"
    return ""


def workbook_output_paths(stage_dir, gene_upper, suffix=""):
    stem = f"1kgp_{gene_upper.lower()}_grch38"
    if suffix:
        stem = f"{stem}_{suffix}"
    return {
        "xlsx": stage_dir / f"{stem}.xlsx",
        "summary_json": stage_dir / f"{stem}.summary.json",
    }


def write_1kgp_workbook(
    gene_upper,
    lookup_payload,
    gene_rows_df,
    missense_only_df,
    stage_dir,
    suffix="",
):
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

    outputs = workbook_output_paths(stage_dir, gene_upper, suffix=suffix)
    with pd.ExcelWriter(outputs["xlsx"]) as writer:
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
        "canonical_transcript": lookup_payload.get("canonical_transcript"),
        "region": f"{lookup_payload['seq_region_name']}:{lookup_payload['start']}-{lookup_payload['end']}",
        "counts": {
            "gene_rows": int(len(gene_rows_df)),
            "missense_only": int(len(missense_only_df)),
            "missense_gene_only": int(len(missense_gene_only_df)),
            "missense_gene_unique": int(len(missense_unique_df)),
            "missense_gene_unique_prot_chan": int(len(missense_unique_prot_df)),
        },
        "files": {
            "xlsx": str(outputs["xlsx"]),
            "summary_json": str(outputs["summary_json"]),
        },
    }
    outputs["summary_json"].write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def install_manual_file(gene_name, manual_source, workspace_root):
    source_path = Path(manual_source).expanduser().resolve()
    if not source_path.exists():
        raise FileNotFoundError(f"Manual curated file not found: {source_path}")

    suffix = source_path.suffix.lower()
    if suffix not in {".csv", ".xlsx"}:
        raise ValueError("Manual curated file must be .csv or .xlsx")

    target_dir = manual_data_dir(workspace_root, gene_name)
    target_dir.mkdir(parents=True, exist_ok=True)
    target_name = "curated_variants.csv" if suffix == ".csv" else "curated_variants.xlsx"
    target_path = target_dir / target_name
    shutil.copy2(source_path, target_path)
    return target_path


def ensure_clinvar_bulk(gene_name, workspace_root, requested_path=None, logger=None):
    if requested_path:
        path = Path(requested_path).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Required ClinVar bulk file not found: {path}")
        return path

    bulk_path = clinvar_data_dir(workspace_root, gene_name) / "clinvar_variant_summary.txt.gz"
    if bulk_path.exists():
        return bulk_path

    if logger:
        logger(f"ClinVar bulk file not found. Downloading latest variant_summary.txt.gz to {bulk_path}")
    return download_file(CLINVAR_VARIANT_SUMMARY_URL, bulk_path).resolve()


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


def build_colab_1kgp_workbook(gene_name, workspace_root, logger=None):
    gene_upper = gene_name.upper()
    ensure_gene_workspace_layout(workspace_root, gene_upper)
    stage_dir = colab_stage_output_dir(workspace_root, gene_upper)
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

    missense_only_df = pd.DataFrame(missense_rows)
    return write_1kgp_workbook(
        gene_upper,
        lookup_payload,
        gene_rows_df,
        missense_only_df,
        stage_dir,
        suffix="",
    )


def build_colab_clinvar_outputs(gene_name, workspace_root, clinvar_bulk_path, logger=None):
    gene_upper = gene_name.upper()
    ensure_gene_workspace_layout(workspace_root, gene_upper)
    output_dir = clinvar_output_dir(workspace_root, gene_upper)
    output_dir.mkdir(parents=True, exist_ok=True)

    if logger:
        logger(f"Reading ClinVar bulk file for {gene_upper} from {clinvar_bulk_path}...")

    usecols = ["GeneSymbol", "Type", "Name", "PhenotypeList", "ClinicalSignificance"]
    chunks = []
    for chunk in pd.read_csv(
        clinvar_bulk_path,
        sep="\t",
        compression="gzip",
        usecols=usecols,
        low_memory=False,
        chunksize=200000,
    ):
        sub = chunk[chunk["GeneSymbol"].astype(str).str.upper() == gene_upper]
        if not sub.empty:
            chunks.append(sub)

    gene_df = pd.concat(chunks, ignore_index=True) if chunks else pd.DataFrame(columns=usecols)
    snv_df = gene_df[
        gene_df["Type"].astype(str).str.contains("single nucleotide variant", case=False, na=False)
    ].copy()
    missense_df = snv_df[
        snv_df["Name"].astype(str).str.contains(r"\(p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)", regex=True, na=False)
    ].copy()
    condition_df = missense_df[
        missense_df["PhenotypeList"].astype(str).str.contains(gene_upper, case=False, na=False)
    ].copy()
    final_df = condition_df.drop_duplicates().copy()

    stage_files = {
        "gene_rows": output_dir / "01_gene_rows.csv",
        "snv_rows": output_dir / "02_snv_rows.csv",
        "missense_rows": output_dir / "03_missense_rows.csv",
        "condition_rows": output_dir / "04_condition_rows.csv",
        "final_rows": output_dir / "05_final_deduplicated_rows.csv",
    }
    gene_df.to_csv(stage_files["gene_rows"], index=False)
    snv_df.to_csv(stage_files["snv_rows"], index=False)
    missense_df.to_csv(stage_files["missense_rows"], index=False)
    condition_df.to_csv(stage_files["condition_rows"], index=False)
    final_df.to_csv(stage_files["final_rows"], index=False)

    summary = {
        "gene": gene_upper,
        "source": str(Path(clinvar_bulk_path).resolve()),
        "condition_filter": {
            "mode": "contains",
            "value": gene_upper,
        },
        "counts": {
            "gene_rows": int(len(gene_df)),
            "snv_rows": int(len(snv_df)),
            "missense_rows": int(len(missense_df)),
            "condition_rows": int(len(condition_df)),
            "final_rows": int(len(final_df)),
        },
        "files": {key: str(path) for key, path in stage_files.items()},
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def build_colab_1kgp_workbook_from_vcf(gene_name, workspace_root, vcf_path, logger=None):
    gene_upper = gene_name.upper()
    ensure_gene_workspace_layout(workspace_root, gene_upper)
    stage_dir = colab_stage_output_dir(workspace_root, gene_upper)
    cache_dir = Path(workspace_root).expanduser() / "cache" / gene_upper.lower() / "1kgp_vcf"
    stage_dir.mkdir(parents=True, exist_ok=True)
    cache_dir.mkdir(parents=True, exist_ok=True)

    vcf_path = Path(vcf_path).expanduser().resolve()
    if not vcf_path.exists():
        raise FileNotFoundError(f"Required VCF file not found: {vcf_path}")

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
    canonical_transcript = lookup_payload.get("canonical_transcript")
    contig_prefix = detect_vcf_contig_prefix(vcf_path)
    region = f"{contig_prefix}{lookup_payload['seq_region_name']}:{lookup_payload['start']}-{lookup_payload['end']}"
    region_vcf = cache_dir / f"{gene_upper.lower()}.region.vcf.gz"

    log(f"Subsetting {region} from {vcf_path.name}...")
    subprocess.run(
        [
            "bcftools",
            "view",
            "-r",
            region,
            str(vcf_path),
            "-Oz",
            "-o",
            str(region_vcf),
        ],
        check=True,
    )
    subprocess.run(["tabix", "-f", "-p", "vcf", str(region_vcf)], check=True)

    raw_lines = run_checked(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%ID\n",
            str(region_vcf),
        ]
    ).splitlines()

    gene_rows = []
    rsid_to_rows = {}
    for line in raw_lines:
        chrom, pos, ref, alt, af, variant_id = line.split("\t")
        row = {
            "CHROM": chrom,
            "POS": int(pos),
            "REF": ref,
            "ALT": alt,
            "RefGenome": "GRCh38",
            "AF": None if af == "." else af,
            "Gene": gene_upper,
            "Consequence": pd.NA,
            "HGVSc": pd.NA,
            "HGVSp": pd.NA,
            "dbSNP ID": variant_id if variant_id != "." else None,
        }
        gene_rows.append(row)
        if variant_id and variant_id != ".":
            rsid_to_rows.setdefault(variant_id, []).append(row)

    gene_rows_df = pd.DataFrame(gene_rows)
    missense_rows = []
    rsids = sorted(rsid_to_rows.keys())
    vep_payload = []
    if rsids:
        log(f"Annotating {len(rsids)} VCF rsIDs through Ensembl VEP...")
        for index in range(0, len(rsids), 200):
            chunk = rsids[index : index + 200]
            vep_payload.extend(post_json(f"{ENSEMBL_REST_BASE}/vep/human/id?hgvs=1", {"ids": chunk}))

    for entry in vep_payload:
        rsid = entry.get("input")
        if not rsid or rsid not in rsid_to_rows:
            continue
        transcript = choose_transcript_consequence(
            entry.get("transcript_consequences", []),
            gene_upper,
            canonical_transcript,
        )
        if transcript is None:
            continue
        if "missense_variant" not in str(transcript.get("consequence_terms", "")):
            continue
        for source_row in rsid_to_rows[rsid]:
            missense_rows.append(
                {
                    "CHROM": source_row["CHROM"],
                    "POS": source_row["POS"],
                    "REF": source_row["REF"],
                    "ALT": source_row["ALT"],
                    "RefGenome": source_row["RefGenome"],
                    "AF": source_row["AF"],
                    "Gene": transcript.get("gene_symbol"),
                    "Consequence": "missense_variant",
                    "HGVSc": transcript.get("hgvsc"),
                    "HGVSp": transcript.get("hgvsp"),
                    "dbSNP ID": rsid,
                }
            )

    missense_only_df = pd.DataFrame(missense_rows)
    summary = write_1kgp_workbook(
        gene_upper,
        lookup_payload,
        gene_rows_df,
        missense_only_df,
        stage_dir,
        suffix="vcf",
    )
    summary["source_vcf"] = str(vcf_path)
    summary["subset_vcf"] = str(region_vcf)
    Path(summary["files"]["summary_json"]).write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def compare_1kgp_outputs(gene_name, workspace_root, summary_path=None):
    gene_upper = gene_name.upper()
    manual_csv = manual_data_dir(workspace_root, gene_upper) / "curated_variants.csv"
    manual_xlsx = manual_data_dir(workspace_root, gene_upper) / "curated_variants.xlsx"
    default_summary = colab_stage_output_dir(workspace_root, gene_upper) / f"1kgp_{gene_upper.lower()}_grch38.summary.json"
    if summary_path:
        summary_file = Path(summary_path).expanduser().resolve()
    else:
        summary_file = default_summary
    if not summary_file.exists():
        raise FileNotFoundError(f"1KGP summary file not found: {summary_file}")

    observed = json.loads(summary_file.read_text())
    expected = DEFAULT_LEGACY_1KGP_COUNTS.get(gene_upper)
    comparison = {
        "workspace_root": str(Path(workspace_root).expanduser().resolve()),
        "summary_file": str(summary_file),
        "manual": {
            "csv_exists": manual_csv.exists(),
            "xlsx_exists": manual_xlsx.exists(),
        },
        "observed_counts": observed.get("counts", {}),
        "expected_counts": expected,
        "retention_percent": {},
    }
    if expected:
        for key, expected_value in expected.items():
            observed_value = observed.get("counts", {}).get(key)
            if observed_value is None:
                continue
            comparison["retention_percent"][key] = round((observed_value / expected_value) * 100.0, 1)
    return comparison


def validate_workspace_root(workspace_root):
    workspace_root = Path(workspace_root).expanduser()
    if str(workspace_root).startswith("/content/drive/") and not Path("/content/drive/MyDrive").exists():
        raise SystemExit(
            "Google Drive does not appear to be mounted. In a Colab notebook cell, run:\n"
            "from google.colab import drive\n"
            "drive.mount('/content/drive')"
        )
    return workspace_root


def copy_colab_files(workspace_root):
    workspace_root = validate_workspace_root(workspace_root)
    build_workspace(workspace_root)

    destination = workspace_root / "scripts" / "tgvr_colab.py"
    destination.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(SCRIPT_PATH, destination)

    return {
        "root": str(workspace_root),
        "scripts_dir": str((workspace_root / "scripts").resolve()),
        "runner": str(destination.resolve()),
    }


def build_parser():
    parser = argparse.ArgumentParser(
        description="Colab-terminal helper for TGVR using a dedicated workspace."
    )
    parser.add_argument(
        "--workspace-root",
        "--drive-root",
        dest="workspace_root",
        default=DEFAULT_WORKSPACE_ROOT,
        help=(
            "TGVR workspace root. Use /content/tgvr_colab for fast ephemeral runs or "
            "/content/drive/MyDrive/tgvr_colab for persistent Google Drive runs "
            f"(default: {DEFAULT_WORKSPACE_ROOT})"
        ),
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser(
        "init",
        help="Create the workspace and copy the Colab TGVR script into it.",
    )

    manual_parser = subparsers.add_parser(
        "install-manual",
        help="Copy a manual curated variant file into the TGVR workspace.",
    )
    manual_parser.add_argument("gene", help="Gene symbol, for example CDKL5")
    manual_parser.add_argument("manual_file", help="Path to a CSV or XLSX manual curated file")

    clinvar_parser = subparsers.add_parser(
        "run-clinvar",
        help="Build ClinVar outputs for a gene inside the TGVR workspace.",
    )
    clinvar_parser.add_argument("gene", help="Gene symbol, for example CDKL5")
    clinvar_parser.add_argument(
        "--clinvar-bulk",
        default=None,
        help="Optional path to a pre-downloaded ClinVar variant_summary.txt.gz file",
    )

    run_parser = subparsers.add_parser(
        "run-1kgp",
        help="Build a Colab-safe 1KGP workbook without Palmetto.",
    )
    run_parser.add_argument("gene", help="Gene symbol, for example CDKL5")

    run_vcf_parser = subparsers.add_parser(
        "run-1kgp-vcf",
        help="Build a 1KGP workbook from a local/source VCF such as the legacy phase3 crossmap file.",
    )
    run_vcf_parser.add_argument("gene", help="Gene symbol, for example CDKL5")
    run_vcf_parser.add_argument(
        "--vcf",
        required=True,
        help="Path to a bgzipped and indexed source VCF, for example phase3.chrX.GRCh38.GT.crossmap.vcf.gz",
    )

    compare_parser = subparsers.add_parser(
        "compare-1kgp",
        help="Compare a generated 1KGP summary against the legacy/TGVR benchmark counts and report manual-data status.",
    )
    compare_parser.add_argument("gene", help="Gene symbol, for example CDKL5")
    compare_parser.add_argument(
        "--summary-path",
        default=None,
        help="Optional path to a specific 1KGP summary json to compare",
    )

    return parser


def log(message):
    print(message, flush=True)


def main():
    parser = build_parser()
    args = parser.parse_args()

    workspace_root = validate_workspace_root(args.workspace_root)

    if args.command == "init":
        copied = copy_colab_files(workspace_root)
        print("Created TGVR workspace and copied the Colab helper:")
        print(json.dumps(copied, indent=2))
        print("Suggested flow:")
        print(f"  cd {workspace_root}")
        print("  mkdir -p manual_input")
        print("  # put your manual curated file under manual_input/ if you have one")
        print(f"  python3 scripts/tgvr_colab.py --workspace-root {workspace_root} install-manual CDKL5 {workspace_root}/manual_input/curated_variants.csv")
        print(f"  python3 scripts/tgvr_colab.py --workspace-root {workspace_root} run-clinvar CDKL5")
        print(f"  python3 scripts/tgvr_colab.py --workspace-root {workspace_root} run-1kgp CDKL5")
        print(f"  # or for the legacy-style VCF-backed path")
        print(f"  python3 scripts/tgvr_colab.py --workspace-root {workspace_root} run-1kgp-vcf CDKL5 --vcf /content/1kgp_source/phase3.chrX.GRCh38.GT.crossmap.vcf.gz")
        print(f"  python3 scripts/tgvr_colab.py --workspace-root {workspace_root} compare-1kgp CDKL5")
        return 0

    if args.command == "install-manual":
        ensure_gene_workspace_layout(workspace_root, args.gene)
        installed_path = install_manual_file(args.gene, args.manual_file, workspace_root)
        print("Installed manual curated file:")
        print(installed_path.resolve())
        return 0

    if args.command == "run-clinvar":
        ensure_gene_workspace_layout(workspace_root, args.gene)
        clinvar_bulk = ensure_clinvar_bulk(args.gene, workspace_root, args.clinvar_bulk, logger=log)
        summary = build_colab_clinvar_outputs(args.gene, workspace_root, clinvar_bulk, logger=log)
        print("Done.")
        print(json.dumps(summary, indent=2))
        print("Final deduplicated rows:", Path(summary["files"]["final_rows"]).resolve())
        print("Summary:", (clinvar_output_dir(workspace_root, args.gene) / "summary.json").resolve())
        return 0

    if args.command == "run-1kgp":
        build_workspace(workspace_root)
        ensure_gene_workspace_layout(workspace_root, args.gene)
        summary = build_colab_1kgp_workbook(args.gene, workspace_root, logger=log)
        print("Done.")
        print(json.dumps(summary, indent=2))
        print("Workbook:", Path(summary["files"]["xlsx"]).resolve())
        print("Summary:", Path(summary["files"]["summary_json"]).resolve())
        return 0

    if args.command == "run-1kgp-vcf":
        build_workspace(workspace_root)
        ensure_gene_workspace_layout(workspace_root, args.gene)
        summary = build_colab_1kgp_workbook_from_vcf(args.gene, workspace_root, args.vcf, logger=log)
        print("Done.")
        print(json.dumps(summary, indent=2))
        print("Workbook:", Path(summary["files"]["xlsx"]).resolve())
        print("Summary:", Path(summary["files"]["summary_json"]).resolve())
        return 0

    if args.command == "compare-1kgp":
        comparison = compare_1kgp_outputs(args.gene, workspace_root, summary_path=args.summary_path)
        print(json.dumps(comparison, indent=2))
        return 0

    raise SystemExit(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
