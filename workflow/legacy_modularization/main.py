import argparse

from src.binding import run_binding_workflow
from src.data_cleaning import run_data_cleaning_workflow
from src.folding import run_folding_workflow


def build_parser():
    parser = argparse.ArgumentParser(description="Workflow entrypoint for modular CDKL5 analysis.")
    parser.add_argument(
        "command",
        nargs="?",
        default="data-cleaning",
        choices=["data-cleaning", "folding", "binding"],
        help="Workflow stage to run.",
    )
    parser.add_argument(
        "--step",
        choices=["all", "prepare", "analyze"],
        default="all",
        help="Sub-step for the folding or binding stage (default: all). Ignored for data-cleaning.",
    )
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "data-cleaning":
        summary = run_data_cleaning_workflow()
        print("Completed data-cleaning workflow.")
        print(f"Combined final rows: {summary['summary']['combined_rows_after_sequence_filter']}")
        print(f"Sequence mismatches removed: {summary['summary']['sequence_mismatches']}")
    elif args.command == "folding":
        result = run_folding_workflow(step=args.step)
        print(f"Completed folding workflow ({args.step}).")
        if "prepare" in result:
            print(f"Prepared mutation inputs: {result['prepare']['mutation_count']}")
        if "analyze" in result:
            print(f"DDG summary rows: {result['analyze']['rows']}")
    elif args.command == "binding":
        result = run_binding_workflow(step=args.step)
        print(f"Completed binding workflow ({args.step}).")
        if "prepare" in result:
            print(f"Prepared mutation inputs: {result['prepare']['mutation_count']}")
        if "analyze" in result:
            print(f"Binding DDG summary rows: {result['analyze']['rows']}")


if __name__ == "__main__":
    main()
