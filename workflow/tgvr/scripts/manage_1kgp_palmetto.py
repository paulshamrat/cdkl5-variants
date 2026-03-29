import argparse
import os
import subprocess
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
SCRIPTS_DIR = SCRIPT_PATH.parent
INTERNAL_1KGP_DIR = SCRIPT_PATH.parents[1] / "internal" / "1kgp_palmetto"

DEFAULT_REMOTE_ROOT = os.environ.get("TGVR_PALMETTO_REMOTE_ROOT", "/home/YOUR_USERNAME/tgvr")
DEFAULT_PALMETTO_TARGET = os.environ.get(
    "TGVR_PALMETTO_TARGET", "YOUR_USERNAME@slogin.palmetto.clemson.edu"
)


def run_local_python(script_name, args):
    cmd = [sys.executable, str(INTERNAL_1KGP_DIR / script_name), *args]
    return subprocess.run(cmd, check=True)


def run_ssh(sock_path, target, remote_command):
    return subprocess.run(
        ["ssh", "-S", sock_path, target, remote_command],
        text=True,
        capture_output=True,
        check=True,
    )


def target_username(target):
    return target.split("@", 1)[0] if "@" in target else target


def build_parser():
    parser = argparse.ArgumentParser(
        description="Single user-facing manager for the TGVR Palmetto-backed raw 1KGP workflow."
    )
    parser.add_argument("gene", help="Currently only CDKL5 is supported")
    parser.add_argument(
        "action",
        choices=["setup", "submit", "status", "log", "fetch"],
        help="Which Palmetto-backed 1KGP action to run",
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
        "--remote-run-dir",
        default=None,
        help="Explicit remote run directory to use for fetch/log. If omitted, the latest matching run under the remote TGVR root is used.",
    )
    parser.add_argument(
        "--job-id",
        default=None,
        help="Optional specific Slurm job id to inspect with status/log. If omitted, status shows all your jobs and log uses the latest matching log in the selected run dir.",
    )
    return parser


def resolve_remote_run_dir(args):
    if args.remote_run_dir:
        return args.remote_run_dir
    latest_cmd = (
        f"ls -1dt {args.remote_root}/outputs/{args.gene.lower()}/01_variant_curation/1kgp/grch38_allvar_noid_* "
        "2>/dev/null | head -n 1"
    )
    result = run_ssh(args.socket, args.palmetto_target, latest_cmd)
    remote_run_dir = result.stdout.strip()
    if not remote_run_dir:
        raise SystemExit("No remote 1KGP Palmetto run directory was found.")
    return remote_run_dir


def do_status(args):
    username = target_username(args.palmetto_target)
    if args.job_id:
        cmd = f"squeue -j {args.job_id} -o '%i %T %M %R %j'"
    else:
        cmd = f"squeue -u {username} -o '%i %T %M %R %j'"
    result = run_ssh(args.socket, args.palmetto_target, cmd)
    print(result.stdout.strip())


def do_log(args):
    remote_run_dir = resolve_remote_run_dir(args)
    if args.job_id:
        log_path = f"{remote_run_dir}/logs/tgvr_1kgp_cdkl5_grch38_allvar_noid_{args.job_id}.log"
    else:
        latest_log_cmd = f"ls -1t {remote_run_dir}/logs/tgvr_1kgp_cdkl5_grch38_allvar_noid_*.log 2>/dev/null | head -n 1"
        result = run_ssh(args.socket, args.palmetto_target, latest_log_cmd)
        log_path = result.stdout.strip()
        if not log_path:
            raise SystemExit("No remote 1KGP Palmetto log file was found.")
    result = run_ssh(args.socket, args.palmetto_target, f"sed -n '1,200p' {log_path}")
    print(result.stdout.rstrip())


def main():
    parser = build_parser()
    args = parser.parse_args()

    action = args.action
    common_args = [
        "--remote-root",
        args.remote_root,
        "--palmetto-target",
        args.palmetto_target,
        "--socket",
        args.socket,
    ]

    if action == "setup":
        run_local_python("setup_1kgp_palmetto_assets.py", common_args)
        return 0

    if action == "submit":
        run_local_python("submit_1kgp_palmetto.py", [args.gene, *common_args])
        return 0

    if action == "fetch":
        fetch_args = [args.gene, *common_args]
        if args.remote_run_dir:
            fetch_args.extend(["--remote-run-dir", args.remote_run_dir])
        run_local_python("fetch_1kgp_palmetto.py", fetch_args)
        return 0

    if action == "status":
        do_status(args)
        return 0

    if action == "log":
        do_log(args)
        return 0

    raise SystemExit(f"Unsupported action: {action}")


if __name__ == "__main__":
    raise SystemExit(main())
