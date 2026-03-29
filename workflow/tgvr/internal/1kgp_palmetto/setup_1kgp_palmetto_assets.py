import argparse
import os
import subprocess


DEFAULT_REMOTE_ROOT = os.environ.get("TGVR_PALMETTO_REMOTE_ROOT", "/home/YOUR_USERNAME/tgvr")
DEFAULT_PALMETTO_TARGET = os.environ.get(
    "TGVR_PALMETTO_TARGET", "YOUR_USERNAME@slogin.palmetto.clemson.edu"
)
DEFAULT_LEGACY_1KGP_ROOT = "/project/ealexov/compbio/shamrat/250419_1000Genome"


def run_ssh(sock_path, target, remote_command):
    return subprocess.run(
        ["ssh", "-S", sock_path, target, remote_command],
        text=True,
        capture_output=True,
        check=True,
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description="Copy the required legacy Palmetto 1KGP runtime assets into a TGVR-owned home workspace so the raw 1KGP backend no longer depends on the legacy project path at run time."
    )
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
        help=f"Read-only legacy 1KGP project root on Palmetto (default: {DEFAULT_LEGACY_1KGP_ROOT})",
    )
    parser.add_argument(
        "--asset-root",
        default=None,
        help="Remote TGVR-owned 1KGP asset root on Palmetto. Defaults to <remote-root>/resources/1kgp",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    asset_root = args.asset_root or f"{args.remote_root}/resources/1kgp"
    remote_cmd = f"""
set -euo pipefail
ASSET_ROOT="{asset_root}"
LEGACY_ROOT="{args.legacy_1kgp_root}"
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
    result = run_ssh(args.socket, args.palmetto_target, remote_cmd)
    print(f"Remote asset root: {asset_root}")
    print(result.stdout.strip())
    if result.stderr.strip():
        print(result.stderr.strip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
