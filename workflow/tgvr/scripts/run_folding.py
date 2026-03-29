import argparse
from pathlib import Path
import sys


SCRIPT_PATH = Path(__file__).resolve()
TGVR_ROOT = SCRIPT_PATH.parents[1]
if str(TGVR_ROOT) not in sys.path:
    sys.path.insert(0, str(TGVR_ROOT))

from src.folding import run_folding_stage
from src.folding.core import DDGUN_METHODS


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run the fresh-start TGVR 02_folding stage for a gene. Current focus: DDGun-first setup and input preparation."
    )
    parser.add_argument("gene", help="Gene symbol, for example CDKL5")
    parser.add_argument("uniprot_id", help="Canonical UniProt accession for the target protein, for example O76039")
    parser.add_argument(
        "--stage",
        choices=["status", "prepare-ddgun"],
        default="status",
        help="Current TGVR 02_folding stages after cleanup (default: status)",
    )
    parser.add_argument(
        "--range",
        dest="residue_range",
        help="Optional residue range for DDGun preparation, formatted as START-END",
    )
    parser.add_argument(
        "--label",
        help='Optional label for DDGun preparation, for example "Kinase domain"',
    )
    parser.add_argument(
        "--method",
        choices=sorted(DDGUN_METHODS.keys()),
        help="Current folding method choices after cleanup, for example ddgun_seq",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    if args.stage == "prepare-ddgun" and not args.method:
        raise SystemExit("--method is required for --stage prepare-ddgun")
    run_folding_stage(
        args.gene.upper(),
        args.uniprot_id.upper(),
        args.stage,
        range_text=args.residue_range,
        label=args.label,
        method=args.method,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
