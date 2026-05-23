from __future__ import annotations

import argparse
from pathlib import Path

from crispr_phage_predictor.benchmark import (
    evaluate_benchmark_run,
    summarize_benchmark_evaluation,
)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Join a SABR run output folder to curated benchmark labels."
    )
    parser.add_argument("run_dir", type=Path, help="Saved SABR output folder")
    parser.add_argument(
        "--benchmark",
        type=Path,
        default=Path("data") / "curation" / "benchmark_labels.tsv",
        help="Benchmark label TSV",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional path for joined benchmark evaluation CSV",
    )
    args = parser.parse_args()

    evaluation = evaluate_benchmark_run(args.run_dir, args.benchmark)
    summary = summarize_benchmark_evaluation(evaluation)

    print("Benchmark evaluation summary")
    print(summary.to_string(index=False))
    print()
    print("Rows by expectation result")
    print(evaluation["score_expectation_result"].value_counts().to_string())

    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        evaluation.to_csv(args.output, index=False)
        print()
        print(f"Wrote joined evaluation to {args.output}")


if __name__ == "__main__":
    main()
