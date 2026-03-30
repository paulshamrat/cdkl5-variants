import argparse
from datetime import datetime
import json
import os
import re
import shutil
import subprocess
import sys
import traceback
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd
import yaml

SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
WORKFLOW_ROOT = REPO_ROOT / "workflow"
TGVR_ROOT = WORKFLOW_ROOT / "tgvr"
PALMETTO_1KGP_TEMPLATE = """#!/bin/bash
#SBATCH --job-name=tgvr_1kgp_cdkl5_grch38_allvar_noid
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=__REMOTE_RUN_DIR__/logs/%x_%j.log

set -euo pipefail

module load anaconda3/2023.09-0
source activate cdkl51000G_env
module load biocontainers samtools/1.17 bcftools/1.17

REMOTE_ROOT=__REMOTE_ROOT__
REMOTE_RUN_DIR=__REMOTE_RUN_DIR__
ASSET_ROOT=__ASSET_ROOT__
DATA=$ASSET_ROOT/GRCh38
VCF=phase3.chrX.GRCh38.GT.crossmap.vcf.gz
CACHE=$ASSET_ROOT/vep_cache_GRCh38
SIF=$ASSET_ROOT/vep.sif
REGION=chrX:18425583-18653629

mkdir -p "$REMOTE_RUN_DIR" "$REMOTE_RUN_DIR/logs"
rm -f "$REMOTE_RUN_DIR"/*.vcf.gz "$REMOTE_RUN_DIR"/*.tbi "$REMOTE_RUN_DIR"/*.tsv "$REMOTE_RUN_DIR"/*.xlsx "$REMOTE_RUN_DIR"/*.html

echo "[1] Subsetting CDKL5 ($REGION)..."
bcftools view -r "$REGION" "$DATA/$VCF" -Oz -o "$REMOTE_RUN_DIR/cdkl5.region.GRCh38.vcf.gz"
tabix -p vcf "$REMOTE_RUN_DIR/cdkl5.region.GRCh38.vcf.gz"

echo "[2] Annotating with VEP (GRCh38)..."
apptainer exec \
  --bind "$REMOTE_RUN_DIR":/data_out \
  --bind "$CACHE":/cache \
  "$SIF" vep \
    --input_file /data_out/cdkl5.region.GRCh38.vcf.gz \
    --output_file /data_out/cdkl5.GRCh38.vep.vcf.gz \
    --format vcf \
    --vcf \
    --cache --dir_cache /cache \
    --assembly GRCh38 \
    --offline \
    --fork 4 \
    --everything \
    --compress_output bgzip \
    --force_overwrite

tabix -p vcf "$REMOTE_RUN_DIR/cdkl5.GRCh38.vep.vcf.gz"

echo "[3] Extracting AF & CSQ..."
bcftools query \
  -f '%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/CSQ\n' \
  "$REMOTE_RUN_DIR/cdkl5.GRCh38.vep.vcf.gz" \
  > "$REMOTE_RUN_DIR/cdkl5.GRCh38.raw.af.csq.tsv"

echo "[4] Parsing Gene, Consequence, HGVSc, HGVSp..."
hdr=$(zgrep '^##INFO=<ID=CSQ' "$REMOTE_RUN_DIR/cdkl5.GRCh38.vep.vcf.gz" | sed 's/.*Format: //;s/\">//')
IFS='|' read -r -a F <<< "$hdr"
for i in "${!F[@]}"; do
  [[ "${F[$i]}" == "SYMBOL" ]] && SI=$((i+1))
  [[ "${F[$i]}" == "Consequence" ]] && CI=$((i+1))
  [[ "${F[$i]}" == "HGVSc" ]] && HSC=$((i+1))
  [[ "${F[$i]}" == "HGVSp" ]] && HSP=$((i+1))
done

echo -e "CHROM\tPOS\tREF\tALT\tRefGenome\tAF\tGene\tConsequence\tHGVSc\tHGVSp" \
  > "$REMOTE_RUN_DIR/cdkl5.GRCh38.all_variants_noid.tsv"

awk -v si="$SI" -v ci="$CI" -v hsc="$HSC" -v hsp="$HSP" -F'\t' 'BEGIN{OFS="\t"}{
  af=$5; if(af==".") af="NA";
  split($6,txs,","); split(txs[1],a,"|");
  print $1,$2,$3,$4,"GRCh38",af,a[si],a[ci],a[hsc],a[hsp]
}' "$REMOTE_RUN_DIR/cdkl5.GRCh38.raw.af.csq.tsv" >> "$REMOTE_RUN_DIR/cdkl5.GRCh38.all_variants_noid.tsv"

python3 - <<'PY2'
import pandas as pd
from pathlib import Path
run_dir = Path("__REMOTE_RUN_DIR__")
tsv_path = run_dir / "cdkl5.GRCh38.all_variants_noid.tsv"
xlsx_path = run_dir / "cdkl5.GRCh38.all_variants_noid.xlsx"
pd.read_csv(tsv_path, sep="\t").to_excel(xlsx_path, index=False)
PY2

echo "DONE: $REMOTE_RUN_DIR"
ls -lh "$REMOTE_RUN_DIR"
"""

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
    "Sec": "U",
    "Pyl": "O",
    "Ter": "*",
    "Asx": "B",
    "Glx": "Z",
    "Xaa": "X",
}
AA_1_TO_3 = {
    value: key
    for key, value in AA_3_TO_1.items()
    if len(value) == 1 and key not in {"Asx", "Glx", "Xaa"}
}

ENSEMBL_REST_BASE = "https://rest.ensembl.org"
CLINVAR_VARIANT_SUMMARY_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
PHASE3_CHRX_CROSSMAP_VCF_URL = (
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/"
    "1000G_2504_high_coverage/working/phase3_liftover_nygc_dir/"
    "phase3.chrX.GRCh38.GT.crossmap.vcf.gz"
)
PHASE3_CHRX_CROSSMAP_TBI_URL = PHASE3_CHRX_CROSSMAP_VCF_URL + ".tbi"
GENE_OUTPUT_SUBDIRS = ["clinvar", "1kgp", "gnomad", "master", "logs"]
CLINVAR_COLUMNS = [
    "Name",
    "Gene(s)",
    "Protein change",
    "Condition(s)",
    "Germline classification",
    "Source",
]
DEFAULT_PALMETTO_USERNAME = os.environ.get("TGVR_PALMETTO_USERNAME", "shamrap")
DEFAULT_PALMETTO_REMOTE_ROOT = os.environ.get("TGVR_PALMETTO_REMOTE_ROOT", f"/home/{DEFAULT_PALMETTO_USERNAME}/tgvr")
DEFAULT_PALMETTO_TARGET = os.environ.get(
    "TGVR_PALMETTO_TARGET", f"{DEFAULT_PALMETTO_USERNAME}@slogin.palmetto.clemson.edu"
)
DEFAULT_PALMETTO_SOCKET = os.path.expanduser("~/.ssh/palmetto.sock")
DEFAULT_LEGACY_1KGP_ROOT = "/project/ealexov/compbio/shamrat/250419_1000Genome"
RUN_LOGGER = None


class RunLogger:
    def __init__(self, gene_name, stage_name):
        logs_dir = stage_output_dir(gene_name) / "logs"
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


def gene_config_path(gene_name):
    return TGVR_ROOT / "config" / "genes" / f"{gene_name.upper()}.yaml"


def load_gene_config(gene_name):
    path = gene_config_path(gene_name)
    if not path.exists():
        raise FileNotFoundError(
            f"No TGVR gene config found for {gene_name.upper()}: {display_path(path)}"
        )
    return yaml.safe_load(path.read_text(encoding="utf-8"))


def stage_data_dir(gene_name):
    return TGVR_ROOT / "00_data" / gene_name.lower() / "01_variant_curation"


def stage_output_dir(gene_name):
    return TGVR_ROOT / "outputs" / gene_name.lower() / "01_variant_curation"


def runtime_stage_dir(gene_name, stage_name):
    path = stage_output_dir(gene_name) / stage_name
    path.mkdir(parents=True, exist_ok=True)
    return path


def runtime_cache_dir(gene_name, stage_name):
    path = runtime_stage_dir(gene_name, stage_name) / "cache"
    path.mkdir(parents=True, exist_ok=True)
    return path


def display_path(path):
    path = Path(path).resolve()
    try:
        return str(path.relative_to(REPO_ROOT))
    except ValueError:
        return str(path)


def print_stage_counts(title, counts):
    RUN_LOGGER.log(title)
    for key, value in counts.items():
        RUN_LOGGER.log(f"  - {key.replace('_', ' ')}: {value}")


def ensure_gene_workspace(gene_name):
    gene_dir = stage_output_dir(gene_name)
    gene_dir.mkdir(parents=True, exist_ok=True)
    for subdir in GENE_OUTPUT_SUBDIRS:
        (gene_dir / subdir).mkdir(parents=True, exist_ok=True)
    return gene_dir


def ensure_gene_input_dirs(gene_name):
    input_root = stage_data_dir(gene_name)
    (input_root / "manual").mkdir(parents=True, exist_ok=True)
    return input_root


def find_manual_file(gene_name):
    gene_lower = gene_name.lower()
    manual_dir = stage_data_dir(gene_lower) / "manual"
    for name in ["curated_variants.csv", "curated_variants.xlsx"]:
        path = manual_dir / name
        if path.exists():
            return path
    return None


def maybe_install_manual_file(gene_name, manual_source):
    if manual_source is None:
        return find_manual_file(gene_name)

    source_path = Path(manual_source).expanduser().resolve()
    if not source_path.exists():
        raise FileNotFoundError(f"Manual curated file not found: {source_path}")

    suffix = source_path.suffix.lower()
    if suffix not in {".csv", ".xlsx"}:
        raise ValueError("Manual curated file must be .csv or .xlsx")

    target_name = "curated_variants.csv" if suffix == ".csv" else "curated_variants.xlsx"
    target_path = stage_data_dir(gene_name.lower()) / "manual" / target_name
    shutil.copy2(source_path, target_path)
    return target_path


def resolve_required_input(path_arg, prompt_text):
    if path_arg:
        path = Path(path_arg).expanduser().resolve()
        if not path.exists():
            raise FileNotFoundError(f"Required input file not found: {path}")
        return path

    if not sys.stdin.isatty():
        raise FileNotFoundError(prompt_text)

    response = input(prompt_text + " ").strip()
    if not response:
        raise FileNotFoundError("Required input was not provided.")
    path = Path(response).expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Required input file not found: {path}")
    return path


def download_file(url, destination):
    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(url, headers={"User-Agent": "codex-cdkl5-variants/1.0"})
    with urllib.request.urlopen(request) as response, open(destination, "wb") as handle:
        shutil.copyfileobj(response, handle)
    return destination


def api_get_json(url):
    request = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(request) as response:
        return json.loads(response.read().decode("utf-8"))


def api_get_text(url):
    request = urllib.request.Request(url, headers={"Accept": "text/plain"})
    with urllib.request.urlopen(request) as response:
        return response.read().decode("utf-8")


def api_post_json(url, payload):
    request = urllib.request.Request(
        url,
        data=json.dumps(payload).encode("utf-8"),
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        method="POST",
    )
    with urllib.request.urlopen(request) as response:
        return json.loads(response.read().decode("utf-8"))


def run_checked(command):
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    return result.stdout


def palmetto_remote_root(args):
    return args.palmetto_remote_root or DEFAULT_PALMETTO_REMOTE_ROOT


def palmetto_target(args):
    return args.palmetto_target or DEFAULT_PALMETTO_TARGET


def palmetto_socket(args):
    return os.path.expanduser(args.palmetto_socket or DEFAULT_PALMETTO_SOCKET)


def palmetto_target_username(target):
    return target.split("@", 1)[0] if "@" in target else target


def run_ssh_checked(sock_path, target, remote_command, stdin_text=None):
    return subprocess.run(
        ["ssh", "-S", sock_path, target, remote_command],
        input=stdin_text,
        text=True,
        capture_output=True,
        check=True,
    )


def render_template(template_text, replacements):
    rendered = template_text
    for key, value in replacements.items():
        rendered = rendered.replace(key, value)
    return rendered


def extract_palmetto_job_id(stdout_text):
    match = re.search(r"Submitted batch job (\d+)", stdout_text)
    return match.group(1) if match else None


def resolve_palmetto_remote_run_dir(gene_name, args):
    if args.palmetto_remote_run_dir:
        return args.palmetto_remote_run_dir
    remote_root = palmetto_remote_root(args)
    target = palmetto_target(args)
    socket_path = palmetto_socket(args)
    latest_cmd = (
        f"ls -1dt {remote_root}/outputs/{gene_name.lower()}/01_variant_curation/1kgp/grch38_allvar_noid_* "
        "2>/dev/null | head -n 1"
    )
    result = run_ssh_checked(socket_path, target, latest_cmd)
    remote_run_dir = result.stdout.strip()
    if not remote_run_dir:
        raise SystemExit("No remote 1KGP Palmetto run directory was found.")
    return remote_run_dir


def setup_palmetto_1kgp_assets(args):
    remote_root = palmetto_remote_root(args)
    target = palmetto_target(args)
    socket_path = palmetto_socket(args)
    asset_root = f"{remote_root}/resources/1kgp"
    remote_cmd = f"""
set -euo pipefail
ASSET_ROOT="{asset_root}"
LEGACY_ROOT="{DEFAULT_LEGACY_1KGP_ROOT}"
mkdir -p "$ASSET_ROOT/GRCh38"
if [ ! -f "$ASSET_ROOT/vep.sif" ]; then
  cp -a "$LEGACY_ROOT/vep.sif" "$ASSET_ROOT/vep.sif"
fi
if [ ! -d "$ASSET_ROOT/vep_cache_GRCh38" ]; then
  cp -a "$LEGACY_ROOT/vep_cache_GRCh38" "$ASSET_ROOT/vep_cache_GRCh38"
fi
if [ ! -f "$ASSET_ROOT/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz" ]; then
  cp -a "$LEGACY_ROOT/00_data/1000G_highcov/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz" "$ASSET_ROOT/GRCh38/"
fi
if [ ! -f "$ASSET_ROOT/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz.tbi" ]; then
  cp -a "$LEGACY_ROOT/00_data/1000G_highcov/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz.tbi" "$ASSET_ROOT/GRCh38/"
fi
du -sh "$ASSET_ROOT/vep.sif" "$ASSET_ROOT/vep_cache_GRCh38" "$ASSET_ROOT/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz" "$ASSET_ROOT/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz.tbi"
"""
    result = run_ssh_checked(socket_path, target, remote_cmd)
    print(f"Remote asset root: {asset_root}")
    print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip())


def submit_palmetto_1kgp(gene_name, args):
    gene_upper = gene_name.upper()
    if gene_upper != "CDKL5":
        raise SystemExit("Palmetto-backed raw 1KGP currently supports only CDKL5")

    remote_root = palmetto_remote_root(args)
    target = palmetto_target(args)
    socket_path = palmetto_socket(args)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    remote_run_dir = f"{remote_root}/outputs/cdkl5/01_variant_curation/1kgp/grch38_allvar_noid_{timestamp}"
    remote_script_path = f"{remote_root}/jobs/1kgp_cdkl5_grch38_allvar_noid_{timestamp}.sh"
    asset_root = f"{remote_root}/resources/1kgp"

    template_text = PALMETTO_1KGP_TEMPLATE
    script_text = render_template(
        template_text,
        {
            "__REMOTE_ROOT__": remote_root,
            "__REMOTE_RUN_DIR__": remote_run_dir,
            "__ASSET_ROOT__": asset_root,
        },
    )
    setup_command = (
        f"mkdir -p {remote_root}/jobs "
        f"{remote_root}/outputs/cdkl5/01_variant_curation/1kgp "
        f"{remote_run_dir}/logs"
    )
    run_ssh_checked(socket_path, target, setup_command)
    run_ssh_checked(
        socket_path,
        target,
        (
            f"test -f {asset_root}/vep.sif "
            f"-a -d {asset_root}/vep_cache_GRCh38 "
            f"-a -f {asset_root}/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz "
            f"-a -f {asset_root}/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz.tbi"
        ),
    )
    run_ssh_checked(socket_path, target, f"cat > {remote_script_path}", stdin_text=script_text)
    run_ssh_checked(socket_path, target, f"chmod +x {remote_script_path} && bash -n {remote_script_path}")
    submit_result = run_ssh_checked(socket_path, target, f"sbatch {remote_script_path}")
    job_id = extract_palmetto_job_id(submit_result.stdout)
    print(f"Remote script: {remote_script_path}")
    print(f"Remote run dir: {remote_run_dir}")
    print(f"Remote asset root: {asset_root}")
    if job_id:
        print(f"Submitted batch job {job_id}")
        print(
            "Check status with:\n"
            f"ssh -S {socket_path} {target} \"squeue -j {job_id} -o '%i %T %M %R'\""
        )
        print(
            "Check log with:\n"
            f"ssh -S {socket_path} {target} "
            f"\"sed -n '1,160p' {remote_run_dir}/logs/tgvr_1kgp_cdkl5_grch38_allvar_noid_{job_id}.log\""
        )
        print(
            "Fetch the completed run back into local TGVR with:\n"
            f"python workflow/tgvr/scripts/run_variant_curation.py {gene_upper} O76039 "
            f"--stage 1kgp --1kgp-mode palmetto --1kgp-palmetto-action fetch "
            f"--palmetto-target {target} --palmetto-remote-root {remote_root} "
            f"--palmetto-remote-run-dir {remote_run_dir}"
        )
    if submit_result.stderr.strip():
        print(submit_result.stderr.strip())


def palmetto_status(args):
    target = palmetto_target(args)
    socket_path = palmetto_socket(args)
    username = palmetto_target_username(target)
    if args.palmetto_job_id:
        cmd = f"squeue -j {args.palmetto_job_id} -o '%i %T %M %R %j'"
    else:
        cmd = f"squeue -u {username} -o '%i %T %M %R %j'"
    result = run_ssh_checked(socket_path, target, cmd)
    print(result.stdout.strip())


def palmetto_log(gene_name, args):
    target = palmetto_target(args)
    socket_path = palmetto_socket(args)
    remote_run_dir = resolve_palmetto_remote_run_dir(gene_name, args)
    if args.palmetto_job_id:
        log_path = f"{remote_run_dir}/logs/tgvr_1kgp_cdkl5_grch38_allvar_noid_{args.palmetto_job_id}.log"
    else:
        latest_log_cmd = f"ls -1t {remote_run_dir}/logs/tgvr_1kgp_cdkl5_grch38_allvar_noid_*.log 2>/dev/null | head -n 1"
        result = run_ssh_checked(socket_path, target, latest_log_cmd)
        log_path = result.stdout.strip()
        if not log_path:
            raise SystemExit("No remote 1KGP Palmetto log file was found.")
    result = run_ssh_checked(socket_path, target, f"sed -n '1,200p' {log_path}")
    print(result.stdout.rstrip())


def fetch_palmetto_1kgp(gene_name, args):
    gene_upper = gene_name.upper()
    if gene_upper != "CDKL5":
        raise SystemExit("Palmetto-backed raw 1KGP currently supports only CDKL5")

    target = palmetto_target(args)
    socket_path = palmetto_socket(args)
    remote_run_dir = resolve_palmetto_remote_run_dir(gene_name, args)
    remote_tsv = f"{remote_run_dir}/cdkl5.GRCh38.all_variants_noid.tsv"
    remote_xlsx = f"{remote_run_dir}/cdkl5.GRCh38.all_variants_noid.xlsx"
    remote_log_dir = f"{remote_run_dir}/logs"

    local_palmetto_dir = runtime_stage_dir(gene_upper, "1kgp") / "palmetto_runs" / Path(remote_run_dir).name
    local_palmetto_dir.mkdir(parents=True, exist_ok=True)

    tsv_text = run_ssh_checked(socket_path, target, f"cat {remote_tsv}").stdout
    (local_palmetto_dir / "cdkl5.GRCh38.all_variants_noid.tsv").write_text(tsv_text, encoding="utf-8")

    xlsx_bytes = subprocess.run(
        ["ssh", "-S", socket_path, target, f"cat {remote_xlsx}"],
        capture_output=True,
        check=True,
    ).stdout
    (local_palmetto_dir / "cdkl5.GRCh38.all_variants_noid.xlsx").write_bytes(xlsx_bytes)

    logs_listing = run_ssh_checked(socket_path, target, f"ls -1 {remote_log_dir} 2>/dev/null || true").stdout
    if logs_listing.strip():
        for log_name in [line.strip() for line in logs_listing.splitlines() if line.strip()]:
            log_text = run_ssh_checked(socket_path, target, f"cat {remote_log_dir}/{log_name}").stdout
            (local_palmetto_dir / log_name).write_text(log_text, encoding="utf-8")

    raw_df = pd.read_csv(local_palmetto_dir / "cdkl5.GRCh38.all_variants_noid.tsv", sep="\t")
    missense_only_df = raw_df[raw_df["Consequence"] == "missense_variant"].copy()
    missense_gene_only_df = missense_only_df[missense_only_df["Gene"].astype(str).str.upper() == gene_upper].copy()
    missense_unique_df = missense_gene_only_df.drop_duplicates().copy()
    missense_unique_prot_df = missense_unique_df.copy()
    missense_unique_prot_df["Protein change 3L"] = missense_unique_prot_df["HGVSp"].astype(str).str.extract(r":p\.(.+)")
    missense_unique_prot_df["Protein change"] = missense_unique_prot_df["Protein change 3L"].apply(convert_three_letter_change)
    missense_unique_prot_df = missense_unique_prot_df[missense_unique_prot_df["Protein change"].notna()].copy()

    workbook_path = runtime_cache_dir(gene_upper, "1kgp") / f"1kgp_{gene_upper.lower()}_grch38.xlsx"
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
    print(
        f"Counts: raw={len(raw_df)} missense={len(missense_only_df)} "
        f"gene_only={len(missense_gene_only_df)} unique={len(missense_unique_df)} "
        f"prot={len(missense_unique_prot_df)}"
    )


def run_palmetto_1kgp_action(gene_name, args):
    action = args.kgp_palmetto_action
    if action == "setup":
        setup_palmetto_1kgp_assets(args)
        return
    if action == "submit":
        submit_palmetto_1kgp(gene_name, args)
        return
    if action == "status":
        palmetto_status(args)
        return
    if action == "log":
        palmetto_log(gene_name, args)
        return
    if action == "fetch":
        fetch_palmetto_1kgp(gene_name, args)
        return
    raise SystemExit(f"Unsupported Palmetto action: {action}")


def fetch_gnomad_gene_payload(gene_name):
    query = """
    query GeneVariants($geneName: String!) {
      gene(gene_symbol: $geneName, reference_genome: GRCh38) {
        gene_id
        canonical_transcript_id
        variants(dataset: gnomad_r4) {
          consequence
          hgvsp
          variant_id
          exome {
            af
          }
          genome {
            af
          }
        }
      }
    }
    """
    payload = {"query": query, "variables": {"geneName": gene_name.upper()}}
    request = urllib.request.Request(
        GNOMAD_API_URL,
        data=json.dumps(payload).encode("utf-8"),
        headers={
            "Accept": "application/json",
            "Content-Type": "application/json",
            "User-Agent": "codex-cdkl5-variants/1.0",
        },
        method="POST",
    )
    with urllib.request.urlopen(request) as response:
        response_payload = json.loads(response.read().decode("utf-8"))
    if "errors" in response_payload:
        raise RuntimeError(f"gnomAD API error: {response_payload['errors']}")
    return response_payload


def initialize_gene_reference(gene_name, uniprot_id):
    gene_upper = gene_name.upper()
    accession = uniprot_id.upper()
    config_path = gene_config_path(gene_upper)
    reference_dir = TGVR_ROOT / "00_data" / gene_name.lower() / "reference"
    reference_dir.mkdir(parents=True, exist_ok=True)
    fasta_relative = Path("00_data") / gene_name.lower() / "reference" / f"canonical_uniprot_{accession}.fasta"
    fasta_path = TGVR_ROOT / fasta_relative

    json_url = f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(accession)}.json"
    fasta_url = f"https://rest.uniprot.org/uniprotkb/{urllib.parse.quote(accession)}.fasta"

    payload = api_get_json(json_url)
    fasta_text = api_get_text(fasta_url)

    primary_accession = str(payload.get("primaryAccession") or accession).upper()
    sequence_block = payload.get("sequence") or {}
    sequence_length = int(sequence_block.get("length") or 0)
    if sequence_length < 1:
        raise ValueError(f"UniProt record for {primary_accession} did not include a valid sequence length.")

    fasta_path.write_text(fasta_text, encoding="utf-8")
    config_path.parent.mkdir(parents=True, exist_ok=True)

    config_payload = {
        "gene": {
            "symbol": gene_upper,
        },
        "clinvar": {
            "condition_filter_mode": "contains",
            "condition_filter_value": gene_upper,
        },
        "reference": {
            "uniprot_accession": primary_accession,
            "sequence_length": sequence_length,
            "fasta": str(fasta_relative).replace("\\", "/"),
            "source": "UniProt REST API",
            "source_url": fasta_url,
        },
        "transcripts": {
            "preferred": [],
        },
        "features": {},
    }

    gene_names = payload.get("genes") or []
    if gene_names:
        gene_name_block = gene_names[0].get("geneName") or {}
        if gene_name_block.get("value"):
            config_payload["reference"]["reported_gene_name"] = gene_name_block["value"]

    uni_entry = payload.get("uniProtkbId")
    if uni_entry:
        config_payload["reference"]["uniprot_entry"] = uni_entry

    comments = payload.get("comments") or []
    for comment in comments:
        if comment.get("commentType") == "ALTERNATIVE PRODUCTS":
            isoforms = comment.get("isoforms") or []
            for isoform in isoforms:
                if isoform.get("displayed"):
                    isoform_ids = isoform.get("isoformIds") or []
                    if isoform_ids:
                        config_payload["reference"]["displayed_isoform"] = isoform_ids[0]
                        break
            if config_payload["reference"].get("displayed_isoform"):
                break

    config_path.write_text(yaml.safe_dump(config_payload, sort_keys=False), encoding="utf-8")
    return config_payload, config_path, fasta_path


def ensure_default_clinvar_bulk(gene_name, requested_path=None):
    default_path = runtime_cache_dir(gene_name, "clinvar") / "clinvar_variant_summary.txt.gz"
    target_path = Path(requested_path) if requested_path else default_path
    if target_path.exists():
        return target_path.resolve()
    if target_path.resolve() == default_path.resolve():
        legacy_path = stage_data_dir(gene_name) / "clinvar" / "clinvar_variant_summary.txt.gz"
        if legacy_path.exists():
            shutil.copy2(legacy_path, default_path)
            RUN_LOGGER.log(
                "ClinVar bulk file not found in outputs cache. Seeded it from the older workflow-local location at "
                f"{display_path(legacy_path)}"
            )
            return default_path.resolve()
        RUN_LOGGER.log(
            "ClinVar bulk file not found. Downloading latest variant_summary.txt.gz to "
            f"{display_path(default_path)}"
        )
        return download_file(CLINVAR_VARIANT_SUMMARY_URL, default_path).resolve()
    raise FileNotFoundError(f"Required input file not found: {target_path.resolve()}")


def ensure_default_gnomad_json(gene_name, requested_path=None):
    default_path = runtime_cache_dir(gene_name, "gnomad") / "raw_response.json"
    target_path = Path(requested_path) if requested_path else default_path
    if target_path.exists():
        return target_path.resolve()
    if target_path.resolve() == default_path.resolve():
        legacy_path = stage_data_dir(gene_name) / "gnomad" / "raw_response.json"
        if legacy_path.exists():
            shutil.copy2(legacy_path, default_path)
            RUN_LOGGER.log(
                "gnomAD raw response not found in outputs cache. Seeded it from the older workflow-local location at "
                f"{display_path(legacy_path)}"
            )
            return default_path.resolve()
        RUN_LOGGER.log(
            "gnomAD raw response not found. Fetching live gene payload to "
            f"{display_path(default_path)}"
        )
        payload = fetch_gnomad_gene_payload(gene_name)
        default_path.parent.mkdir(parents=True, exist_ok=True)
        default_path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return default_path.resolve()
    raise FileNotFoundError(f"Required input file not found: {target_path.resolve()}")


def ensure_default_1kgp_workbook(gene_name, requested_path=None):
    gene_lower = gene_name.lower()
    default_path = runtime_cache_dir(gene_name, "1kgp") / f"1kgp_{gene_lower}_grch38.xlsx"
    target_path = Path(requested_path) if requested_path else default_path
    if target_path.exists():
        return target_path.resolve()
    if requested_path:
        raise FileNotFoundError(f"Required input file not found: {target_path.resolve()}")

    legacy_workflow_path = stage_data_dir(gene_name) / "1kgp" / f"1kgp_{gene_lower}_grch38.xlsx"
    if legacy_workflow_path.exists():
        default_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(legacy_workflow_path, default_path)
        RUN_LOGGER.log(
            "Cached 1KGP workbook not found in outputs cache. Seeded it from the older workflow-local location at "
            f"{display_path(legacy_workflow_path)}"
        )
        return default_path.resolve()

    legacy_seed = REPO_ROOT / "00_data" / f"1kgp_{gene_lower}_grch38.xlsx"
    if legacy_seed.exists():
        default_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(legacy_seed, default_path)
        RUN_LOGGER.log(
            "Cached 1KGP workbook not found in outputs cache. Seeded it from the legacy local workbook at "
            f"{display_path(legacy_seed)}"
        )
        return default_path.resolve()

    return None


def ensure_default_1kgp_vcf(gene_name, requested_path=None):
    gene_lower = gene_name.lower()
    source_dir = runtime_stage_dir(gene_name, "1kgp") / "source"
    default_vcf_path = source_dir / f"phase3.{gene_lower}.chrX.GRCh38.GT.crossmap.vcf.gz"
    if requested_path:
        target_path = Path(requested_path).expanduser().resolve()
        if not target_path.exists():
            raise FileNotFoundError(f"Required input file not found: {target_path}")
        return target_path

    if default_vcf_path.exists():
        return default_vcf_path.resolve()

    legacy_vcf_path = stage_data_dir(gene_name) / "1kgp" / "source" / f"phase3.{gene_lower}.chrX.GRCh38.GT.crossmap.vcf.gz"
    legacy_tbi_path = Path(str(legacy_vcf_path) + ".tbi")
    if legacy_vcf_path.exists():
        source_dir.mkdir(parents=True, exist_ok=True)
        shutil.copy2(legacy_vcf_path, default_vcf_path)
        if legacy_tbi_path.exists():
            shutil.copy2(legacy_tbi_path, Path(str(default_vcf_path) + ".tbi"))
        RUN_LOGGER.log(
            "1KGP phase3 crossmap VCF not found in outputs source cache. Seeded it from the older workflow-local location at "
            f"{display_path(legacy_vcf_path)}"
        )
        return default_vcf_path.resolve()

    source_dir.mkdir(parents=True, exist_ok=True)
    RUN_LOGGER.log(
        "1KGP phase3 crossmap VCF not found in outputs source cache. Downloading public chrX crossmap source to "
        f"{display_path(default_vcf_path)}"
    )
    download_file(PHASE3_CHRX_CROSSMAP_VCF_URL, default_vcf_path)
    tbi_path = Path(str(default_vcf_path) + ".tbi")
    download_file(PHASE3_CHRX_CROSSMAP_TBI_URL, tbi_path)
    return default_vcf_path.resolve()


def detect_vcf_contig_prefix(vcf_path):
    header = run_checked(["bcftools", "view", "-h", str(vcf_path)])
    if "##contig=<ID=chrX>" in header or "##contig=<ID=chr1>" in header:
        return "chr"
    if "##contig=<ID=X>" in header:
        return ""
    return ""


def convert_three_letter_change(change):
    if not isinstance(change, str):
        return None
    match = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2}|Ter)", change)
    if not match:
        return None
    wild, position, mutant = match.groups()
    wild_one = AA_3_TO_1.get(wild)
    mutant_one = AA_3_TO_1.get(mutant)
    if not wild_one or not mutant_one:
        return None
    return f"{wild_one}{position}{mutant_one}"


def build_three_letter_hgvsp(amino_acids, protein_start):
    if not isinstance(amino_acids, str) or "/" not in amino_acids or protein_start is None:
        return None
    wild_one, mutant_one = amino_acids.split("/", 1)
    wild_three = AA_1_TO_3.get(wild_one)
    mutant_three = AA_1_TO_3.get(mutant_one)
    if not wild_three or not mutant_three:
        return None
    return f"p.{wild_three}{protein_start}{mutant_three}"


def convert_hgvsp_to_one_letter(hgvsp):
    if not isinstance(hgvsp, str):
        return None
    cleaned = hgvsp.replace("p.", "")
    return convert_three_letter_change(cleaned)


def extract_protein_change(name):
    if not isinstance(name, str):
        return None
    match = re.search(r"\(p\.([A-Z][a-z]{2}\d+[A-Z][a-z]{2}|[A-Z][a-z]{2}\d+Ter)\)", name)
    if not match:
        return None
    return convert_three_letter_change(match.group(1))


def extract_hgvsp_change(hgvsp):
    if not isinstance(hgvsp, str):
        return None
    match = re.search(r"p\.([A-Z][a-z]{2}\d+(?:[A-Z][a-z]{2}|Ter))", hgvsp)
    if not match:
        return None
    return match.group(1)


def is_simple_missense_change(change):
    if not isinstance(change, str):
        return False
    match = re.fullmatch(r"([A-Z])(\d+)([A-Z])", change)
    if not match:
        return False
    wild, _position, mutant = match.groups()
    if wild == mutant:
        return False
    return wild != "*" and mutant != "*"


def read_fasta_sequence(fasta_path):
    sequence_lines = []
    seen_header = False
    with open(fasta_path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seen_header and sequence_lines:
                    break
                seen_header = True
                continue
            sequence_lines.append(line)
    return "".join(sequence_lines)


def parse_mutation(mutation_str):
    if not isinstance(mutation_str, str) or len(mutation_str) < 3:
        return None, None, None
    wild = mutation_str[0]
    mutant = mutation_str[-1]
    try:
        position = int(mutation_str[1:-1])
    except ValueError:
        return None, None, None
    return wild, position, mutant


def add_mutation_components(df, mutation_col="Mutation"):
    parsed = df[mutation_col].apply(parse_mutation)
    df["wild"] = parsed.apply(lambda x: x[0] if x else None)
    df["position"] = parsed.apply(lambda x: x[1] if x else None)
    df["mutant"] = parsed.apply(lambda x: x[2] if x else None)
    return df


def check_reference_match(wild, position, reference_sequence):
    if wild is None or position is None:
        return "Invalid"
    if position < 1 or position > len(reference_sequence):
        return "Out of range"
    return "Match" if reference_sequence[position - 1] == wild else f"Mismatch (Seq: {reference_sequence[position - 1]})"


def build_clinvar_outputs(gene_name, clinvar_bulk_path):
    gene_upper = gene_name.upper()
    gene_config = load_gene_config(gene_upper)
    output_dir = stage_output_dir(gene_name) / "clinvar"
    output_dir.mkdir(parents=True, exist_ok=True)
    clinvar_config = gene_config.get("clinvar", {})
    condition_filter_mode = str(clinvar_config.get("condition_filter_mode", "contains")).lower()
    condition_filter_value = str(clinvar_config.get("condition_filter_value", gene_upper))

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
    snv_df = gene_df[gene_df["Type"].astype(str).str.contains("single nucleotide variant", case=False, na=False)].copy()
    missense_df = snv_df[
        snv_df["Name"].astype(str).str.contains(r"\(p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}\)", regex=True, na=False)
    ].copy()
    phenotype_series = missense_df["PhenotypeList"].astype(str)
    if condition_filter_mode == "contains":
        condition_mask = phenotype_series.str.contains(condition_filter_value, case=False, na=False)
    else:
        raise ValueError(
            f"Unsupported clinvar.condition_filter_mode for {gene_upper}: {condition_filter_mode}"
        )
    condition_df = missense_df[condition_mask].copy()
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
        "source": str(clinvar_bulk_path),
        "condition_filter": {
            "mode": condition_filter_mode,
            "value": condition_filter_value,
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
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def choose_transcript_consequence(
    transcript_consequences,
    gene_name,
    canonical_transcript,
    preferred_transcripts=None,
):
    gene_upper = gene_name.upper()
    canonical_root = canonical_transcript.split(".")[0] if canonical_transcript else None
    normalized_preferred = [str(tx).split(".")[0] for tx in (preferred_transcripts or []) if tx]
    preferred = normalized_preferred + ([canonical_root] if canonical_root and canonical_root not in normalized_preferred else [])
    preferred_index = {tx: idx for idx, tx in enumerate(preferred)}

    candidates = []
    for consequence in transcript_consequences or []:
        if consequence.get("gene_symbol") != gene_upper:
            continue
        if "missense_variant" not in consequence.get("consequence_terms", []):
            continue
        if not str(consequence.get("biotype", "")).startswith("protein_coding"):
            continue
        transcript_id = consequence.get("transcript_id")
        transcript_root = transcript_id.split(".")[0] if transcript_id else None
        rank = preferred_index.get(transcript_id, len(preferred) + 10)
        if transcript_root in preferred_index:
            rank = min(rank, preferred_index[transcript_root])
        flags = consequence.get("flags") or []
        has_hgvsp = 0 if consequence.get("hgvsp") else 1
        candidates.append((rank, has_hgvsp, len(flags), transcript_id or "", consequence))

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


def normalize_overlap_row(overlap_row, gene_name):
    alleles = overlap_row.get("alleles") or []
    ref = alleles[0] if len(alleles) >= 1 else None
    alt = alleles[1] if len(alleles) >= 2 else None
    return {
        "CHROM": f"chr{overlap_row.get('seq_region_name')}",
        "POS": overlap_row.get("start"),
        "REF": ref,
        "ALT": alt,
        "RefGenome": overlap_row.get("assembly_name"),
        "AF": pd.NA,
        "Gene": gene_name.upper(),
        "Consequence": overlap_row.get("consequence_type"),
        "HGVSc": pd.NA,
        "HGVSp": pd.NA,
        "dbSNP ID": overlap_row.get("id"),
    }


def build_1kgp_outputs(gene_name):
    gene_upper = gene_name.upper()
    output_dir = stage_output_dir(gene_name) / "1kgp"
    output_dir.mkdir(parents=True, exist_ok=True)
    gene_config = load_gene_config(gene_upper)
    preferred_transcripts = gene_config.get("transcripts", {}).get("preferred", [])

    lookup_url = (
        f"{ENSEMBL_REST_BASE}/lookup/symbol/homo_sapiens/"
        f"{urllib.parse.quote(gene_upper)}?content-type=application/json"
    )
    lookup_payload = api_get_json(lookup_url)
    chromosome = lookup_payload["seq_region_name"]
    start = lookup_payload["start"]
    end = lookup_payload["end"]
    canonical_transcript = lookup_payload.get("canonical_transcript")

    region_url = (
        f"{ENSEMBL_REST_BASE}/overlap/region/human/"
        f"{chromosome}:{start}-{end}?feature=variation;variant_set=1kg_3;content-type=application/json"
    )
    region_payload = api_get_json(region_url)

    gene_rows_df = pd.DataFrame([normalize_overlap_row(row, gene_upper) for row in region_payload])
    missense_region_rows = [row for row in region_payload if row.get("consequence_type") == "missense_variant"]
    missense_ids = [row["id"] for row in missense_region_rows if row.get("id")]

    vep_payload = []
    if missense_ids:
        vep_payload = api_post_json(f"{ENSEMBL_REST_BASE}/vep/human/id?hgvs=1", {"ids": missense_ids})
    vep_by_id = {entry["input"]: entry for entry in vep_payload}

    missense_rows = []
    for region_row in missense_region_rows:
        rsid = region_row.get("id")
        vep_entry = vep_by_id.get(rsid, {})
        transcript = choose_transcript_consequence(
            vep_entry.get("transcript_consequences", []),
            gene_upper,
            canonical_transcript,
            preferred_transcripts=preferred_transcripts,
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
    missense_gene_only_df = missense_only_df[missense_only_df["Gene"] == gene_upper].copy()
    missense_unique_df = missense_gene_only_df.drop_duplicates().copy()
    missense_unique_prot_df = missense_unique_df.copy()
    missense_unique_prot_df["Protein change 3L"] = missense_unique_prot_df["HGVSp"].apply(extract_hgvsp_change)
    missense_unique_prot_df["Protein change"] = missense_unique_prot_df["Protein change 3L"].apply(convert_three_letter_change)

    stage_files = {
        "lookup_response": output_dir / "00_gene_lookup.json",
        "region_response": output_dir / "00_region_response.json",
        "missense_vep_response": output_dir / "00_missense_vep_response.json",
        "gene_rows": output_dir / "01_gene_rows.csv",
        "missense_only": output_dir / "02_missense_only.csv",
        "missense_gene_only": output_dir / "03_missense_gene_only.csv",
        "missense_unique": output_dir / "04_missense_gene_unique.csv",
        "missense_unique_prot": output_dir / "05_missense_gene_unique_prot_chan.csv",
    }

    stage_files["lookup_response"].write_text(json.dumps(lookup_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    stage_files["region_response"].write_text(json.dumps(region_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    stage_files["missense_vep_response"].write_text(json.dumps(vep_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    gene_rows_df.to_csv(stage_files["gene_rows"], index=False)
    missense_only_df.to_csv(stage_files["missense_only"], index=False)
    missense_gene_only_df.to_csv(stage_files["missense_gene_only"], index=False)
    missense_unique_df.to_csv(stage_files["missense_unique"], index=False)
    missense_unique_prot_df.to_csv(stage_files["missense_unique_prot"], index=False)

    summary = {
        "gene": gene_upper,
        "ensembl_gene_id": lookup_payload.get("id"),
        "canonical_transcript": canonical_transcript,
        "counts": {
            "gene_rows": int(len(gene_rows_df)),
            "missense_only": int(len(missense_only_df)),
            "missense_gene_only": int(len(missense_gene_only_df)),
            "missense_gene_unique": int(len(missense_unique_df)),
            "missense_gene_unique_prot_chan": int(len(missense_unique_prot_df)),
        },
        "files": {key: str(path) for key, path in stage_files.items()},
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def build_genomic_hgvs(seq_region_name, pos, ref, alt):
    if not all(isinstance(value, str) and value for value in [seq_region_name, ref, alt]):
        return None
    if "," in alt:
        return None
    if len(ref) == 1 and len(alt) == 1:
        return f"{seq_region_name}:g.{pos}{ref}>{alt}"
    return None


def build_1kgp_outputs_vcf(gene_name, vcf_path):
    gene_upper = gene_name.upper()
    output_dir = stage_output_dir(gene_name) / "1kgp"
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = output_dir / "cache_vcf"
    cache_dir.mkdir(parents=True, exist_ok=True)
    gene_config = load_gene_config(gene_upper)
    preferred_transcripts = gene_config.get("transcripts", {}).get("preferred", [])

    lookup_url = (
        f"{ENSEMBL_REST_BASE}/lookup/symbol/homo_sapiens/"
        f"{urllib.parse.quote(gene_upper)}?content-type=application/json"
    )
    lookup_payload = api_get_json(lookup_url)
    canonical_transcript = lookup_payload.get("canonical_transcript")
    seq_region_name = lookup_payload["seq_region_name"]
    contig_prefix = detect_vcf_contig_prefix(vcf_path)
    region = f"{contig_prefix}{seq_region_name}:{lookup_payload['start']}-{lookup_payload['end']}"

    region_vcf = cache_dir / f"{gene_name.lower()}.region.vcf.gz"
    subprocess.run(
        ["bcftools", "view", "-r", region, str(vcf_path), "-Oz", "-o", str(region_vcf)],
        check=True,
    )
    subprocess.run(["tabix", "-f", "-p", "vcf", str(region_vcf)], check=True)

    query_output = run_checked(
        [
            "bcftools",
            "query",
            "-f",
            "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%ID\n",
            str(region_vcf),
        ]
    ).splitlines()

    gene_rows = []
    rsid_rows = {}
    hgvs_rows = {}
    for index, line in enumerate(query_output):
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
            "dbSNP ID": None if variant_id == "." else variant_id,
        }
        gene_rows.append(row)
        if variant_id != ".":
            rsid_rows.setdefault(variant_id, []).append(row)
        else:
            hgvs = build_genomic_hgvs(seq_region_name, pos, ref, alt)
            if hgvs:
                hgvs_rows.setdefault(hgvs, []).append(row)

    gene_rows_df = pd.DataFrame(gene_rows)

    vep_by_rsid = {}
    rsids = sorted(rsid_rows)
    if rsids:
        for start in range(0, len(rsids), 200):
            chunk = rsids[start : start + 200]
            for entry in api_post_json(f"{ENSEMBL_REST_BASE}/vep/human/id?hgvs=1", {"ids": chunk}):
                if "input" in entry:
                    vep_by_rsid[entry["input"]] = entry

    vep_by_hgvs = {}
    hgvs_ids = sorted(hgvs_rows)
    if hgvs_ids:
        for start in range(0, len(hgvs_ids), 200):
            chunk = hgvs_ids[start : start + 200]
            variants = []
            for hgvs in chunk:
                source_row = hgvs_rows[hgvs][0]
                variants.append(
                    f"{seq_region_name} {source_row['POS']} . {source_row['REF']} {source_row['ALT']} . . ."
                )
            for entry in api_post_json(f"{ENSEMBL_REST_BASE}/vep/homo_sapiens/region", {"variants": variants}):
                for consequence in entry.get("transcript_consequences") or []:
                    if not consequence.get("hgvsp"):
                        synthesized_hgvsp = build_three_letter_hgvsp(
                            consequence.get("amino_acids"),
                            consequence.get("protein_start"),
                        )
                        if synthesized_hgvsp:
                            consequence["hgvsp"] = synthesized_hgvsp
                if "input" in entry:
                    fields = str(entry["input"]).split()
                    if len(fields) >= 5:
                        input_key = build_genomic_hgvs(seq_region_name, fields[1], fields[3], fields[4])
                        if input_key:
                            vep_by_hgvs[input_key] = entry

    missense_rows = []

    def append_annotated_rows(source_rows, vep_entry, fallback_id=None):
        transcript = choose_transcript_consequence(
            vep_entry.get("transcript_consequences", []),
            gene_upper,
            canonical_transcript,
            preferred_transcripts=preferred_transcripts,
        )
        if transcript is None:
            return
        if "missense_variant" not in (transcript.get("consequence_terms") or []):
            return
        colocated = vep_entry.get("colocated_variants") or []
        colocated_id = next((item.get("id") for item in colocated if item.get("id")), None)
        for source_row in source_rows:
            variant_allele = transcript.get("variant_allele") or source_row["ALT"]
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
                    "dbSNP ID": source_row["dbSNP ID"] or colocated_id or fallback_id,
                }
            )

    for rsid, rows in rsid_rows.items():
        entry = vep_by_rsid.get(rsid)
        if entry:
            append_annotated_rows(rows, entry, fallback_id=rsid)

    for hgvs, rows in hgvs_rows.items():
        entry = vep_by_hgvs.get(hgvs)
        if entry:
            append_annotated_rows(rows, entry)

    missense_only_df = pd.DataFrame(missense_rows)
    missense_gene_only_df = missense_only_df[missense_only_df["Gene"] == gene_upper].copy()
    missense_unique_df = missense_gene_only_df.drop_duplicates().copy()
    missense_unique_prot_df = missense_unique_df.copy()
    missense_unique_prot_df["Protein change 3L"] = missense_unique_prot_df["HGVSp"].apply(extract_hgvsp_change)
    missense_unique_prot_df["Protein change"] = missense_unique_prot_df["Protein change 3L"].apply(convert_three_letter_change)
    missense_unique_prot_df = missense_unique_prot_df[missense_unique_prot_df["Protein change"].notna()].copy()

    stage_files = {
        "lookup_response": output_dir / "00_gene_lookup.json",
        "region_vcf": output_dir / "00_region_subset.vcf.gz",
        "gene_rows": output_dir / "01_gene_rows.csv",
        "missense_only": output_dir / "02_missense_only.csv",
        "missense_gene_only": output_dir / "03_missense_gene_only.csv",
        "missense_unique": output_dir / "04_missense_gene_unique.csv",
        "missense_unique_prot": output_dir / "05_missense_gene_unique_prot_chan.csv",
    }
    stage_files["lookup_response"].write_text(json.dumps(lookup_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    shutil.copy2(region_vcf, stage_files["region_vcf"])
    index_path = Path(str(region_vcf) + ".tbi")
    if index_path.exists():
        shutil.copy2(index_path, Path(str(stage_files["region_vcf"]) + ".tbi"))
    gene_rows_df.to_csv(stage_files["gene_rows"], index=False)
    missense_only_df.to_csv(stage_files["missense_only"], index=False)
    missense_gene_only_df.to_csv(stage_files["missense_gene_only"], index=False)
    missense_unique_df.to_csv(stage_files["missense_unique"], index=False)
    missense_unique_prot_df.to_csv(stage_files["missense_unique_prot"], index=False)

    summary = {
        "gene": gene_upper,
        "mode": "vcf",
        "ensembl_gene_id": lookup_payload.get("id"),
        "canonical_transcript": canonical_transcript,
        "source_vcf": str(Path(vcf_path).resolve()),
        "subset_region": region,
        "counts": {
            "gene_rows": int(len(gene_rows_df)),
            "missense_only": int(len(missense_only_df)),
            "missense_gene_only": int(len(missense_gene_only_df)),
            "missense_gene_unique": int(len(missense_unique_df)),
            "missense_gene_unique_prot_chan": int(len(missense_unique_prot_df)),
        },
        "files": {key: str(path) for key, path in stage_files.items()},
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def build_1kgp_outputs_cached(gene_name, workbook_path):
    gene_upper = gene_name.upper()
    output_dir = stage_output_dir(gene_name) / "1kgp"
    output_dir.mkdir(parents=True, exist_ok=True)

    raw_df = pd.read_excel(workbook_path, sheet_name="Sheet1")
    missense_only_df = raw_df[raw_df["Consequence"] == "missense_variant"].copy()
    missense_gene_only_df = missense_only_df[missense_only_df["Gene"].astype(str).str.upper() == gene_upper].copy()
    missense_unique_df = missense_gene_only_df.drop_duplicates().copy()
    missense_unique_prot_df = missense_unique_df.copy()
    missense_unique_prot_df["Protein change 3L"] = missense_unique_prot_df["HGVSp"].apply(extract_hgvsp_change)
    missense_unique_prot_df["Protein change"] = missense_unique_prot_df["Protein change 3L"].apply(convert_three_letter_change)
    missense_unique_prot_df = missense_unique_prot_df[missense_unique_prot_df["Protein change"].notna()].copy()

    stage_files = {
        "source_workbook": output_dir / "00_cached_source_workbook.txt",
        "gene_rows": output_dir / "01_gene_rows.csv",
        "missense_only": output_dir / "02_missense_only.csv",
        "missense_gene_only": output_dir / "03_missense_gene_only.csv",
        "missense_unique": output_dir / "04_missense_gene_unique.csv",
        "missense_unique_prot": output_dir / "05_missense_gene_unique_prot_chan.csv",
    }

    stage_files["source_workbook"].write_text(str(Path(workbook_path).resolve()) + "\n", encoding="utf-8")
    raw_df.to_csv(stage_files["gene_rows"], index=False)
    missense_only_df.to_csv(stage_files["missense_only"], index=False)
    missense_gene_only_df.to_csv(stage_files["missense_gene_only"], index=False)
    missense_unique_df.to_csv(stage_files["missense_unique"], index=False)
    missense_unique_prot_df.to_csv(stage_files["missense_unique_prot"], index=False)

    summary = {
        "gene": gene_upper,
        "mode": "cached",
        "source_workbook": str(Path(workbook_path).resolve()),
        "counts": {
            "gene_rows": int(len(raw_df)),
            "missense_only": int(len(missense_only_df)),
            "missense_gene_only": int(len(missense_gene_only_df)),
            "missense_gene_unique": int(len(missense_unique_df)),
            "missense_gene_unique_prot_chan": int(len(missense_unique_prot_df)),
        },
        "files": {key: str(path) for key, path in stage_files.items()},
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def build_gnomad_outputs(gene_name, raw_json_path):
    output_dir = stage_output_dir(gene_name) / "gnomad"
    output_dir.mkdir(parents=True, exist_ok=True)

    payload = json.loads(Path(raw_json_path).read_text(encoding="utf-8"))
    gene_data = payload["data"]["gene"]
    variants = gene_data.get("variants", [])
    gene_df = pd.DataFrame(variants)

    if gene_df.empty:
        summary = {
            "gene": gene_name.upper(),
            "counts": {"gene_variants": 0, "missense_variants": 0, "converted_variants": 0, "af_enriched_variants": 0},
            "files": {},
        }
        (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        return summary

    gene_df["exome_af"] = gene_df["exome"].apply(lambda x: None if not isinstance(x, dict) else x.get("af"))
    gene_df["genome_af"] = gene_df["genome"].apply(lambda x: None if not isinstance(x, dict) else x.get("af"))
    gene_df["max_af"] = gene_df[["exome_af", "genome_af"]].max(axis=1, skipna=True)

    missense_df = gene_df[gene_df["consequence"].astype(str) == "missense_variant"].copy()
    converted_df = missense_df.copy()
    converted_df["protein_change"] = converted_df["hgvsp"].apply(convert_hgvsp_to_one_letter)
    af_df = converted_df[converted_df["protein_change"].notna()].copy()

    stage_files = {
        "gene_variants": output_dir / "01_gene_variants.csv",
        "missense_variants": output_dir / "02_missense_variants.csv",
        "converted_variants": output_dir / "03_one_letter_converted.csv",
        "af_enriched_variants": output_dir / "04_af_enriched.csv",
    }
    gene_df.drop(columns=["exome", "genome"]).to_csv(stage_files["gene_variants"], index=False)
    missense_df.drop(columns=["exome", "genome"]).to_csv(stage_files["missense_variants"], index=False)
    converted_df.drop(columns=["exome", "genome"]).to_csv(stage_files["converted_variants"], index=False)
    af_df.drop(columns=["exome", "genome"]).to_csv(stage_files["af_enriched_variants"], index=False)

    summary = {
        "gene": gene_name.upper(),
        "canonical_transcript_id": gene_data.get("canonical_transcript_id"),
        "gnomad_gene_id": gene_data.get("gene_id"),
        "counts": {
            "gene_variants": int(len(gene_df)),
            "missense_variants": int(len(missense_df)),
            "converted_variants": int(converted_df["protein_change"].notna().sum()),
            "af_enriched_variants": int(len(af_df)),
        },
        "files": {key: str(path) for key, path in stage_files.items()},
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def load_manual_curated_variants(gene_name, use_manual=False):
    gene_lower = gene_name.lower()
    csv_path = stage_data_dir(gene_lower) / "manual" / "curated_variants.csv"
    xlsx_path = stage_data_dir(gene_lower) / "manual" / "curated_variants.xlsx"
    if not use_manual:
        return None, None, None
    if csv_path.exists():
        return pd.read_csv(csv_path), str(csv_path), "csv"
    if xlsx_path.exists():
        return pd.read_excel(xlsx_path, sheet_name="Sheet1"), str(xlsx_path), "xlsx"
    raise FileNotFoundError(
        f"--use-manual was requested but no manual curated file was found in "
        f"{display_path(stage_data_dir(gene_lower) / 'manual')}"
    )


def build_master_table(gene_name, use_manual=False):
    gene_lower = gene_name.lower()
    gene_upper = gene_name.upper()
    base_dir = stage_output_dir(gene_lower)
    clinvar_path = base_dir / "clinvar" / "05_final_deduplicated_rows.csv"
    kgp_path = base_dir / "1kgp" / "05_missense_gene_unique_prot_chan.csv"
    if not kgp_path.exists():
        legacy_kgp_path = base_dir / "1kgp" / "05_missense_cdkl5_unique_prot_chan.csv"
        if legacy_kgp_path.exists():
            kgp_path = legacy_kgp_path
    gnomad_path = base_dir / "gnomad" / "04_af_enriched.csv"
    output_dir = base_dir / "master"
    output_dir.mkdir(parents=True, exist_ok=True)

    gene_config = load_gene_config(gene_upper)
    reference_fasta = (TGVR_ROOT / gene_config["reference"]["fasta"]).resolve()
    reference_sequence = read_fasta_sequence(reference_fasta)
    if len(reference_sequence) != int(gene_config["reference"]["sequence_length"]):
        raise ValueError(
            f"Reference sequence length mismatch for {gene_upper}: "
            f"expected {gene_config['reference']['sequence_length']}, found {len(reference_sequence)}"
        )

    clinvar_df = pd.read_csv(clinvar_path)
    clinvar_standardized = pd.DataFrame(
        {
            "Name": clinvar_df["Name"],
            "Gene(s)": gene_upper,
            "Protein change": clinvar_df["Name"].apply(extract_protein_change),
            "Condition(s)": clinvar_df["PhenotypeList"],
            "Germline classification": clinvar_df["ClinicalSignificance"],
            "Source": "ClinVar (live)",
        }
    )
    clinvar_standardized = clinvar_standardized[clinvar_standardized["Protein change"].notna()].copy()
    clinvar_standardized = clinvar_standardized[
        clinvar_standardized["Protein change"].apply(is_simple_missense_change)
    ].copy()
    clinvar_standardized = clinvar_standardized.drop_duplicates(subset=["Protein change"]).copy()

    kgp_df = pd.read_csv(kgp_path)
    kgp_standardized = pd.DataFrame(
        {
            "Name": kgp_df["HGVSp"],
            "Gene(s)": gene_upper,
            "Protein change": kgp_df["Protein change"],
            "Condition(s)": "Healthy",
            "Germline classification": "Benign",
            "Source": "The 1000 Genomes Project",
        }
    )
    clinvar_set = set(clinvar_standardized["Protein change"].dropna().astype(str))
    kgp_unique = kgp_standardized[~kgp_standardized["Protein change"].astype(str).isin(clinvar_set)].copy()

    manual_standardized = pd.DataFrame(columns=CLINVAR_COLUMNS)
    manual_unique = pd.DataFrame(columns=CLINVAR_COLUMNS)
    manual_df, manual_source_path, manual_source_format = load_manual_curated_variants(gene_upper, use_manual=use_manual)
    manual_source_mode = "provided"

    if manual_df is not None:
        required_columns = {"Mutation", "Consequence", "Source"}
        missing_columns = sorted(required_columns - set(manual_df.columns))
        if missing_columns:
            raise ValueError(f"Manual curated file is missing required columns: {', '.join(missing_columns)}")

        manual_standardized = pd.DataFrame(
            {
                "Name": manual_df["Source"],
                "Gene(s)": gene_upper,
                "Protein change": manual_df["Mutation"],
                "Condition(s)": pd.NA,
                "Germline classification": manual_df["Consequence"],
                "Source": manual_df["Source"],
            }
        )
        manual_standardized = manual_standardized[
            manual_standardized["Protein change"].apply(is_simple_missense_change)
        ].copy()
        manual_standardized = manual_standardized.drop_duplicates(subset=["Protein change"]).copy()
        combined_pre_manual_set = set(pd.concat([clinvar_standardized["Protein change"], kgp_unique["Protein change"]]).dropna().astype(str))
        manual_unique = manual_standardized[~manual_standardized["Protein change"].astype(str).isin(combined_pre_manual_set)].copy()

    combined_df = pd.concat([clinvar_standardized, kgp_unique, manual_unique], ignore_index=True)
    combined_df = combined_df.drop_duplicates(subset=["Protein change"]).copy()

    gnomad_df = pd.read_csv(gnomad_path)
    gnomad_lookup = gnomad_df[["protein_change", "max_af"]].rename(
        columns={"protein_change": "Protein change", "max_af": "gnomAD Allele Frequency"}
    )
    combined_with_gaf_df = combined_df.merge(gnomad_lookup, on="Protein change", how="left")
    sequence_check_df = combined_with_gaf_df.copy()
    sequence_check_df["Mutation"] = sequence_check_df["Protein change"]
    sequence_check_df = add_mutation_components(sequence_check_df, mutation_col="Mutation")
    sequence_check_df["match_status"] = sequence_check_df.apply(
        lambda row: check_reference_match(row["wild"], row["position"], reference_sequence), axis=1
    )
    final_df = sequence_check_df[sequence_check_df["match_status"] == "Match"].copy()
    mismatches_df = sequence_check_df[sequence_check_df["match_status"] != "Match"].copy()

    feature_config = gene_config.get("features", {}).get("primary_domain")
    feature_label = feature_config["label"] if feature_config else None
    feature_start = feature_config["start"] if feature_config else None
    feature_end = feature_config["end"] if feature_config else None
    if feature_start is not None and feature_end is not None:
        feature_df = final_df[final_df["position"].between(feature_start, feature_end)].copy()
    else:
        feature_df = None

    labels = sorted(final_df["Germline classification"].dropna().astype(str).unique())
    full_length_summary_df = pd.DataFrame(
        {
            "Label": labels,
            "Count": [int((final_df["Germline classification"].astype(str) == label).sum()) for label in labels],
        }
    )
    totals_rows = [
        {
            "Scope": "Full length",
            "Residue range": f"1-{gene_config['reference']['sequence_length']}",
            "Variant count": int(len(final_df)),
        }
    ]
    domain_summary_df = pd.DataFrame(columns=["Label", "Count"])
    if feature_df is not None and feature_label is not None:
        domain_labels = sorted(feature_df["Germline classification"].dropna().astype(str).unique())
        domain_summary_df = pd.DataFrame(
            {
                "Label": domain_labels,
                "Count": [int((feature_df["Germline classification"].astype(str) == label).sum()) for label in domain_labels],
            }
        )
        totals_rows.append(
            {
                "Scope": feature_label,
                "Residue range": f"{feature_start}-{feature_end}",
                "Variant count": int(len(feature_df)),
            }
        )
    classification_totals_df = pd.DataFrame(totals_rows)

    stage_files = {
        "clinvar_standardized": output_dir / "01_clinvar_standardized.csv",
        "kgp_standardized": output_dir / "02_1kgp_standardized.csv",
        "kgp_unique": output_dir / "03_1kgp_unique_vs_clinvar.csv",
        "manual_standardized": output_dir / "04_manual_standardized.csv",
        "manual_unique": output_dir / "05_manual_unique_vs_clinvar_1kgp.csv",
        "combined": output_dir / "06_clinvar_1kgp_manual_combined.csv",
        "combined_with_gaf": output_dir / "07_clinvar_1kgp_manual_combined_with_gnomad_af.csv",
        "sequence_check": output_dir / "08_sequence_check.csv",
        "mismatches": output_dir / "09_sequence_mismatches.csv",
        "final": output_dir / "10_master_dataset_final.csv",
        "full_length_summary": output_dir / "11_full_length_classification_summary.csv",
        "domain_summary": output_dir / "12_domain_classification_summary.csv",
        "classification_totals": output_dir / "13_scope_totals.csv",
        "manual_note": output_dir / "MANUAL_STEP_NOTE.md",
    }
    clinvar_standardized.to_csv(stage_files["clinvar_standardized"], index=False)
    kgp_standardized.to_csv(stage_files["kgp_standardized"], index=False)
    kgp_unique.to_csv(stage_files["kgp_unique"], index=False)
    manual_standardized.to_csv(stage_files["manual_standardized"], index=False)
    manual_unique.to_csv(stage_files["manual_unique"], index=False)
    combined_df.to_csv(stage_files["combined"], index=False)
    combined_with_gaf_df.to_csv(stage_files["combined_with_gaf"], index=False)
    sequence_check_df.to_csv(stage_files["sequence_check"], index=False)
    mismatches_df.to_csv(stage_files["mismatches"], index=False)
    final_df.to_csv(stage_files["final"], index=False)
    full_length_summary_df.to_csv(stage_files["full_length_summary"], index=False)
    domain_summary_df.to_csv(stage_files["domain_summary"], index=False)
    classification_totals_df.to_csv(stage_files["classification_totals"], index=False)

    manual_note = "\n".join(
        [
            "# Manual Step Note",
            "",
            "This master-table workflow supports manually curated variants.",
            "",
            "Expected file location:",
            f"- workflow/tgvr/00_data/{gene_lower}/01_variant_curation/manual/curated_variants.csv",
            "",
            "What happens:",
            "- manual curated variants are added only if they are not already present in the ClinVar + 1KGP combined set",
            "- gnomAD AF is added after the combined table is assembled",
            "- separate summary files are written for full-length labels and, when configured, for a gene-specific domain",
        ]
    )
    stage_files["manual_note"].write_text(manual_note + "\n", encoding="utf-8")

    summary = {
        "gene": gene_upper,
        "manual_sources": {
            "manual_curated_variants": {
                "included": manual_df is not None,
                "mode": manual_source_mode if manual_df is not None else ("not_requested" if not use_manual else "missing"),
                "path": manual_source_path,
                "format": manual_source_format,
            }
        },
        "counts": {
            "clinvar_standardized": int(len(clinvar_standardized)),
            "kgp_standardized": int(len(kgp_standardized)),
            "kgp_unique_vs_clinvar": int(len(kgp_unique)),
            "manual_standardized": int(len(manual_standardized)),
            "manual_unique_vs_clinvar_1kgp": int(len(manual_unique)),
            "combined_clinvar_1kgp_manual": int(len(combined_df)),
            "combined_with_gnomad_af": int(len(combined_with_gaf_df)),
            "rows_with_gnomad_af": int(combined_with_gaf_df["gnomAD Allele Frequency"].notna().sum()),
            "sequence_check_rows": int(len(sequence_check_df)),
            "mismatch_rows": int(len(mismatches_df)),
            "final_post_sequence_check": int(len(final_df)),
        },
        "domain_ranges": {
            "full_length": f"1-{gene_config['reference']['sequence_length']}",
        },
        "files": {key: str(path) for key, path in stage_files.items()},
        "reference": {
            "fasta": str(reference_fasta),
            "source": gene_config["reference"]["source"],
            "source_url": gene_config["reference"]["source_url"],
            "sequence_length": int(gene_config["reference"]["sequence_length"]),
            "uniprot_accession": gene_config["reference"].get("uniprot_accession"),
        },
    }
    if feature_df is not None and feature_label is not None:
        summary["counts"][f"{feature_label.lower().replace(' ', '_')}_post_sequence_check"] = int(len(feature_df))
        summary["domain_ranges"][feature_label.lower().replace(" ", "_")] = f"{feature_start}-{feature_end}"
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def load_master_summary_scope(gene_name, scope):
    gene_upper = gene_name.upper()
    output_dir = stage_output_dir(gene_name.lower()) / "master"
    full_path = output_dir / "11_full_length_classification_summary.csv"
    domain_path = output_dir / "12_domain_classification_summary.csv"
    totals_path = output_dir / "13_scope_totals.csv"

    if scope == "full":
        if not full_path.exists():
            raise FileNotFoundError(f"Full-length summary not found: {display_path(full_path)}. Run the master stage first.")
        return pd.read_csv(full_path), full_path

    gene_config = load_gene_config(gene_upper)
    feature_config = gene_config.get("features", {}).get("primary_domain")
    if feature_config is None:
        raise ValueError(f"No special domain range is configured for {gene_upper}.")
    if not domain_path.exists():
        raise FileNotFoundError(f"Domain summary not found: {display_path(domain_path)}. Run the master stage first.")
    return pd.read_csv(domain_path), domain_path


def print_summary_table(gene_name, scope):
    summary_df, summary_path = load_master_summary_scope(gene_name, scope)
    totals_path = stage_output_dir(gene_name.lower()) / "master" / "13_scope_totals.csv"
    totals_df = pd.read_csv(totals_path) if totals_path.exists() else pd.DataFrame()
    gene_config = load_gene_config(gene_name.upper())

    if scope == "full":
        scope_label = "Full-length classification summary"
        totals_row = totals_df[totals_df["Scope"] == "Full length"]
    else:
        feature_label = gene_config["features"]["primary_domain"]["label"]
        scope_label = f"{feature_label} classification summary"
        totals_row = totals_df[totals_df["Scope"] == feature_label]

    RUN_LOGGER.log(scope_label + ":")
    if not totals_row.empty:
        row = totals_row.iloc[0]
        RUN_LOGGER.log(f"  - residue range: {row['Residue range']}")
        RUN_LOGGER.log(f"  - total variants in this scope: {row['Variant count']}")
    for _, row in summary_df.iterrows():
        RUN_LOGGER.log(f"  - {row['Label']}: {row['Count']}")
    RUN_LOGGER.log(f"Summary file: {display_path(summary_path)}")


def parse_residue_range(range_text):
    match = re.fullmatch(r"(\d+)-(\d+)", str(range_text).strip())
    if not match:
        raise ValueError("Range must be formatted as START-END, for example 13-297.")
    start, end = int(match.group(1)), int(match.group(2))
    if start < 1 or end < start:
        raise ValueError("Range must satisfy 1 <= START <= END.")
    return start, end


def print_custom_range_summary(gene_name, range_text, label=None):
    start, end = parse_residue_range(range_text)
    gene_upper = gene_name.upper()
    master_dir = stage_output_dir(gene_name.lower()) / "master"
    final_path = master_dir / "10_master_dataset_final.csv"
    if not final_path.exists():
        raise FileNotFoundError(f"Master dataset not found: {display_path(final_path)}. Run the master stage first.")

    gene_config = load_gene_config(gene_upper)
    final_df = pd.read_csv(final_path)
    if "position" not in final_df.columns:
        raise ValueError("Master dataset is missing the position column required for range summaries.")
    range_df = final_df[final_df["position"].between(start, end)].copy()
    labels = sorted(range_df["Germline classification"].dropna().astype(str).unique())
    summary_df = pd.DataFrame(
        {
            "Label": labels,
            "Count": [int((range_df["Germline classification"].astype(str) == item).sum()) for item in labels],
        }
    )

    summary_label = label or f"Residues {start}-{end}"
    safe_label = re.sub(r"[^a-z0-9]+", "_", summary_label.lower()).strip("_") or "range"
    range_summary_path = master_dir / f"range_summary_{safe_label}_{start}_{end}.csv"
    summary_df.to_csv(range_summary_path, index=False)

    protein_length = int(gene_config["reference"]["sequence_length"])
    if end > protein_length:
        RUN_LOGGER.log(
            f"Warning: requested range {start}-{end} extends beyond the configured protein length {protein_length}."
        )

    RUN_LOGGER.log(f"{summary_label} classification summary:")
    RUN_LOGGER.log(f"  - residue range: {start}-{end}")
    RUN_LOGGER.log(f"  - total variants in this range: {len(range_df)}")
    for _, row in summary_df.iterrows():
        RUN_LOGGER.log(f"  - {row['Label']}: {row['Count']}")
    RUN_LOGGER.log(f"Summary file: {display_path(range_summary_path)}")


def main():
    global RUN_LOGGER
    parser = argparse.ArgumentParser(
        description="Run the live variant-curation workflow for a gene and save outputs into workflow/tgvr/outputs/<gene>/01_variant_curation/."
    )
    parser.add_argument("gene", help="Gene symbol, for example CDKL5")
    parser.add_argument("uniprot_id", help="Canonical UniProt accession for the target protein, for example O76039")
    parser.add_argument(
        "--stage",
        choices=["clinvar", "1kgp", "gnomad", "master", "summary-full", "summary-range", "all"],
        default="all",
        help="Run a single stage or the whole workflow (default: all)",
    )
    parser.add_argument(
        "--clinvar-bulk",
        default=None,
        help=(
            "Path to downloaded ClinVar bulk file "
            "(default: workflow/tgvr/outputs/<gene>/01_variant_curation/clinvar/cache/clinvar_variant_summary.txt.gz)"
        ),
    )
    parser.add_argument(
        "--1kgp-mode",
        dest="kgp_mode",
        choices=["auto", "cached", "live", "vcf", "palmetto"],
        default="auto",
        help=(
            "How to build the 1000 Genomes stage. "
            "'cached' uses a stable workbook under workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/cache/, "
            "'live' queries Ensembl Phase 3 overlap directly, "
            "'vcf' uses the public phase3 chrX crossmap VCF, "
            "'palmetto' runs the legacy Palmetto-backed raw 1KGP flow directly from this script, "
            "and 'auto' prefers cached then falls back to live."
        ),
    )
    parser.add_argument(
        "--1kgp-workbook",
        dest="kgp_workbook",
        default=None,
        help=(
            "Optional path to a prepared 1KGP workbook with Sheet1 columns like the legacy "
            "1kgp_<gene>_grch38.xlsx export."
        ),
    )
    parser.add_argument(
        "--1kgp-vcf",
        dest="kgp_vcf",
        default=None,
        help=(
            "Optional path to a phase3 chrX crossmap VCF. "
            "If omitted in --1kgp-mode vcf, the public source is downloaded into "
            "workflow/tgvr/outputs/<gene>/01_variant_curation/1kgp/source/."
        ),
    )
    parser.add_argument(
        "--1kgp-palmetto-action",
        dest="kgp_palmetto_action",
        choices=["setup", "submit", "status", "log", "fetch"],
        default=None,
        help="Palmetto action to use when --1kgp-mode palmetto is selected.",
    )
    parser.add_argument("--palmetto-remote-root", default=None, help="Override the remote TGVR root used by the Palmetto-backed 1KGP flow.")
    parser.add_argument("--palmetto-target", default=None, help="Override the Palmetto SSH target used by the Palmetto-backed 1KGP flow.")
    parser.add_argument("--palmetto-socket", default=None, help="Override the SSH control socket path used by the Palmetto-backed 1KGP flow.")
    parser.add_argument("--palmetto-remote-run-dir", default=None, help="Explicit remote Palmetto run directory for fetch/log actions.")
    parser.add_argument("--palmetto-job-id", default=None, help="Optional Palmetto Slurm job id for status/log actions.")
    parser.add_argument("--gnomad-json", help="Path to prepared raw gnomAD gene JSON response")
    parser.add_argument(
        "--manual-file",
        help=(
            "Optional path to curated_variants.csv or curated_variants.xlsx. "
            "If provided, the file is copied into 00_data/<gene>/01_variant_curation/manual/ before the master stage."
        ),
    )
    parser.add_argument(
        "--use-manual",
        action="store_true",
        help=(
            "Use the standard TGVR manual curated file during the master stage: "
            "workflow/tgvr/00_data/<gene>/01_variant_curation/manual/curated_variants.csv. "
            "If this flag is not passed, the master stage runs without manual curated variants."
        ),
    )
    parser.add_argument(
        "--range",
        dest="residue_range",
        help="Residue range for --stage summary-range, formatted as START-END",
    )
    parser.add_argument(
        "--label",
        help="Optional label for --stage summary-range, for example 'Kinase domain'",
    )
    args = parser.parse_args()

    gene_upper = args.gene.upper()
    supplied_uniprot = args.uniprot_id.upper()
    if args.stage == "summary-range" and not args.residue_range:
        raise SystemExit("--range START-END is required for --stage summary-range")
    if args.manual_file:
        args.use_manual = True

    RUN_LOGGER = RunLogger(gene_upper, args.stage.replace("-", "_"))
    gene_dir = ensure_gene_workspace(gene_upper)
    ensure_gene_input_dirs(gene_upper)
    try:
        RUN_LOGGER.log(f"Gene workspace: {display_path(gene_dir)}")
        RUN_LOGGER.log(f"Run log: {display_path(RUN_LOGGER.run_log_path)}")

        config_path = gene_config_path(gene_upper)
        if not config_path.exists():
            RUN_LOGGER.log(
                f"No TGVR gene config found for {gene_upper}. Initializing a standalone TGVR reference "
                f"from UniProt accession {supplied_uniprot}."
            )
            config, created_config_path, created_fasta_path = initialize_gene_reference(gene_upper, supplied_uniprot)
            RUN_LOGGER.log(f"Created gene config: {display_path(created_config_path)}")
            RUN_LOGGER.log(f"Saved canonical reference FASTA: {display_path(created_fasta_path)}")
        else:
            config = load_gene_config(gene_upper)

        expected_uniprot = str(config["reference"].get("uniprot_accession", "")).upper()
        if expected_uniprot and supplied_uniprot != expected_uniprot:
            raise SystemExit(
                f"UniProt ID mismatch for {gene_upper}: expected {expected_uniprot}, received {supplied_uniprot}"
            )

        if args.manual_file:
            installed = maybe_install_manual_file(gene_upper, args.manual_file)
            RUN_LOGGER.log(f"Installed manual curated file: {display_path(installed)}")
            RUN_LOGGER.log("Manual curated variants: enabled for this run via --manual-file")
        elif args.use_manual and args.stage in {"master", "all"}:
            existing_manual = find_manual_file(gene_upper)
            if existing_manual is None:
                raise FileNotFoundError(
                    f"--use-manual was requested but no manual curated file was found in "
                    f"{display_path(stage_data_dir(gene_upper) / 'manual')}"
                )
            RUN_LOGGER.log(f"Using existing manual curated file: {display_path(existing_manual)}")
            RUN_LOGGER.log("Manual curated variants: enabled for this run via --use-manual")
        elif args.stage in {"master", "all"}:
            RUN_LOGGER.log("Manual curated variants: not requested for this run")

        if args.stage in {"clinvar", "all"}:
            RUN_LOGGER.log("Running ClinVar stage...")
            RUN_LOGGER.log("  - filtering latest ClinVar bulk data into gene, SNV, missense, and condition-specific tables")
            clinvar_bulk = ensure_default_clinvar_bulk(gene_upper, args.clinvar_bulk)
            clinvar_summary = build_clinvar_outputs(gene_upper, clinvar_bulk)
            print_stage_counts("ClinVar summary:", clinvar_summary["counts"])

        if args.stage in {"1kgp", "all"}:
            RUN_LOGGER.log("Running 1000 Genomes stage...")
            if args.kgp_mode == "palmetto":
                if args.stage != "1kgp":
                    raise SystemExit("--1kgp-mode palmetto is only supported with --stage 1kgp")
                if not args.kgp_palmetto_action:
                    raise SystemExit("--1kgp-palmetto-action is required when --1kgp-mode palmetto is selected")
                RUN_LOGGER.log(
                    f"  - running the legacy Palmetto-backed 1KGP path with action '{args.kgp_palmetto_action}'"
                )
                run_palmetto_1kgp_action(gene_upper, args)
                return 0

            if args.kgp_mode in {"auto", "cached"}:
                workbook_path = ensure_default_1kgp_workbook(gene_upper, args.kgp_workbook)
            else:
                workbook_path = None
            if args.kgp_mode == "vcf":
                vcf_path = ensure_default_1kgp_vcf(gene_upper, args.kgp_vcf)
            else:
                vcf_path = None

            if args.kgp_mode == "cached" and workbook_path is None:
                raise FileNotFoundError(
                    "No cached 1KGP workbook was found. Provide --1kgp-workbook or place "
                    f"{display_path(runtime_cache_dir(gene_upper, '1kgp') / f'1kgp_{gene_upper.lower()}_grch38.xlsx')} "
                    "before running with --1kgp-mode cached."
                )

            if workbook_path is not None:
                RUN_LOGGER.log("  - using a cached 1KGP workbook as the primary benign-source backend")
                RUN_LOGGER.log(f"  - cached workbook: {display_path(workbook_path)}")
                kgp_summary = build_1kgp_outputs_cached(gene_upper, workbook_path)
            elif vcf_path is not None:
                RUN_LOGGER.log("  - using the public phase3 chrX crossmap VCF as the benign-source backend")
                RUN_LOGGER.log(f"  - VCF source: {display_path(vcf_path)}")
                kgp_summary = build_1kgp_outputs_vcf(gene_upper, vcf_path)
            else:
                RUN_LOGGER.log("  - no cached workbook found; querying Ensembl Phase 3 live")
                RUN_LOGGER.log("  - keeping gene-specific missense changes from the live overlap/VEP path")
                kgp_summary = build_1kgp_outputs(gene_upper)
            print_stage_counts("1000 Genomes summary:", kgp_summary["counts"])

        if args.stage in {"gnomad", "all"}:
            RUN_LOGGER.log("Running gnomAD stage...")
            RUN_LOGGER.log("  - fetching gene variants, keeping missense rows, and preparing allele-frequency lookup tables")
            gnomad_json = ensure_default_gnomad_json(gene_upper, args.gnomad_json)
            gnomad_summary = build_gnomad_outputs(gene_upper, gnomad_json)
            print_stage_counts("gnomAD summary:", gnomad_summary["counts"])

        if args.stage in {"master", "all"}:
            RUN_LOGGER.log("Building master table...")
            if args.use_manual:
                RUN_LOGGER.log(
                    "  - merging ClinVar, unique 1KGP, requested manual curated variants, "
                    "and gnomAD AF before sequence checking"
                )
            else:
                RUN_LOGGER.log(
                    "  - merging ClinVar, unique 1KGP, and gnomAD AF before sequence checking "
                    "(manual curated variants not included)"
                )
            master_summary = build_master_table(gene_upper, use_manual=args.use_manual)
            print_stage_counts("Master table summary:", master_summary["counts"])
            final_path = gene_dir / "master" / "10_master_dataset_final.csv"
            RUN_LOGGER.log(f"Final master dataset: {display_path(final_path)}")

        if args.stage == "summary-full":
            RUN_LOGGER.log("Reading full-length classification summary...")
            print_summary_table(gene_upper, "full")

        if args.stage == "summary-range":
            RUN_LOGGER.log("Reading custom range classification summary...")
            print_custom_range_summary(gene_upper, args.residue_range, args.label)

        return 0
    except Exception as exc:
        RUN_LOGGER.log(f"ERROR: {exc}")
        RUN_LOGGER.log(traceback.format_exc().rstrip())
        return 1
    finally:
        RUN_LOGGER.close()


if __name__ == "__main__":
    raise SystemExit(main())
