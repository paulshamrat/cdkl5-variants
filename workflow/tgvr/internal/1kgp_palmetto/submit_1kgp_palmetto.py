import argparse
import os
import re
import subprocess
from datetime import datetime
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
TGVR_ROOT = SCRIPT_PATH.parents[2]
TEMPLATE_PATH = TGVR_ROOT / "palmetto_templates" / "1kgp_cdkl5_grch38_allvar_noid.sh.tmpl"

DEFAULT_REMOTE_ROOT = os.environ.get("TGVR_PALMETTO_REMOTE_ROOT", "/home/YOUR_USERNAME/tgvr")
DEFAULT_PALMETTO_TARGET = os.environ.get(
    "TGVR_PALMETTO_TARGET", "YOUR_USERNAME@slogin.palmetto.clemson.edu"
)
DEFAULT_LEGACY_1KGP_ROOT = "/project/ealexov/compbio/shamrat/250419_1000Genome"


def render_template(template_text, replacements):
    rendered = template_text
    for key, value in replacements.items():
        rendered = rendered.replace(key, value)
    return rendered


def run_ssh(sock_path, target, remote_command, stdin_text=None):
    return subprocess.run(
        ["ssh", "-S", sock_path, target, remote_command],
        input=stdin_text,
        text=True,
        capture_output=True,
        check=True,
    )


def extract_job_id(stdout_text):
    match = re.search(r"Submitted batch job (\d+)", stdout_text)
    return match.group(1) if match else None


def build_parser():
    parser = argparse.ArgumentParser(
        description="Submit the legacy-style CDKL5 GRCh38 1KGP Palmetto job into a safe writable home workspace through the Palmetto bridge."
    )
    parser.add_argument("gene", help="Currently only CDKL5 is supported")
    parser.add_argument(
        "--remote-root",
        default=DEFAULT_REMOTE_ROOT,
        help=f"Remote writable TGVR root on Palmetto (default: {DEFAULT_REMOTE_ROOT})",
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
    parser.add_argument(
        "--legacy-1kgp-root",
        default=DEFAULT_LEGACY_1KGP_ROOT,
        help=f"Read-only legacy 1KGP project root on Palmetto used only by the separate asset-setup step (default: {DEFAULT_LEGACY_1KGP_ROOT})",
    )
    parser.add_argument(
        "--asset-root",
        default=None,
        help="Remote TGVR-owned 1KGP asset root on Palmetto. Defaults to <remote-root>/resources/1kgp",
    )
    parser.add_argument(
        "--submit-only",
        action="store_true",
        help="Submit the Palmetto job and stop after printing the remote paths and job id",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    gene_upper = args.gene.upper()
    if gene_upper != "CDKL5":
        raise SystemExit("submit_1kgp_palmetto.py currently supports only CDKL5")

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    remote_run_dir = f"{args.remote_root}/outputs/cdkl5/01_variant_curation/1kgp/grch38_allvar_noid_{timestamp}"
    remote_script_path = f"{args.remote_root}/jobs/1kgp_cdkl5_grch38_allvar_noid_{timestamp}.sh"
    asset_root = args.asset_root or f"{args.remote_root}/resources/1kgp"

    template_text = TEMPLATE_PATH.read_text(encoding="utf-8")
    script_text = render_template(
        template_text,
        {
            "__REMOTE_ROOT__": args.remote_root,
            "__REMOTE_RUN_DIR__": remote_run_dir,
            "__ASSET_ROOT__": asset_root,
        },
    )

    setup_command = (
        f"mkdir -p {args.remote_root}/jobs "
        f"{args.remote_root}/outputs/cdkl5/01_variant_curation/1kgp "
        f"{remote_run_dir}/logs"
    )
    run_ssh(args.socket, args.palmetto_target, setup_command)
    run_ssh(
        args.socket,
        args.palmetto_target,
        (
            f"test -f {asset_root}/vep.sif "
            f"-a -d {asset_root}/vep_cache_GRCh38 "
            f"-a -f {asset_root}/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz "
            f"-a -f {asset_root}/GRCh38/phase3.chrX.GRCh38.GT.crossmap.vcf.gz.tbi"
        ),
    )
    run_ssh(args.socket, args.palmetto_target, f"cat > {remote_script_path}", stdin_text=script_text)
    run_ssh(args.socket, args.palmetto_target, f"chmod +x {remote_script_path} && bash -n {remote_script_path}")
    submit_result = run_ssh(args.socket, args.palmetto_target, f"sbatch {remote_script_path}")

    job_id = extract_job_id(submit_result.stdout)
    print(f"Remote script: {remote_script_path}")
    print(f"Remote run dir: {remote_run_dir}")
    print(f"Remote asset root: {asset_root}")
    if job_id:
        print(f"Submitted batch job {job_id}")
        print(
            "Check status with:\n"
            f"ssh -S {args.socket} {args.palmetto_target} \"squeue -j {job_id} -o '%i %T %M %R'\""
        )
        print(
            "Check log with:\n"
            f"ssh -S {args.socket} {args.palmetto_target} "
            f"\"sed -n '1,160p' {remote_run_dir}/logs/tgvr_1kgp_cdkl5_grch38_allvar_noid_{job_id}.log\""
        )
        print(
            "Fetch the completed run back into local TGVR with:\n"
            f"python workflow/tgvr/scripts/manage_1kgp_palmetto.py {gene_upper} fetch "
            f"--remote-run-dir {remote_run_dir}"
        )

    if submit_result.stderr.strip():
        print(submit_result.stderr.strip())

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
