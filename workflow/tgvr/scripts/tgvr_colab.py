import argparse
import json
import shutil
import sys
from pathlib import Path


SCRIPT_PATH = Path(__file__).resolve()
REPO_ROOT = SCRIPT_PATH.parents[3]
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from workflow.tgvr.internal.colab_1kgp import build_colab_1kgp_workbook, build_drive_workspace


DEFAULT_DRIVE_ROOT = "/content/drive/MyDrive/tgvr_colab"


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

    copy_map = [
        (
            REPO_ROOT / "workflow" / "tgvr" / "scripts" / "tgvr_colab.py",
            drive_root / "scripts" / "tgvr_colab.py",
        ),
        (
            REPO_ROOT / "workflow" / "tgvr" / "internal" / "colab_1kgp.py",
            drive_root / "internal" / "colab_1kgp.py",
        ),
    ]

    for source, destination in copy_map:
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, destination)

    return {
        "root": str(drive_root),
        "scripts_dir": str((drive_root / "scripts").resolve()),
        "internal_dir": str((drive_root / "internal").resolve()),
        "runner": str((drive_root / "scripts" / "tgvr_colab.py").resolve()),
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
