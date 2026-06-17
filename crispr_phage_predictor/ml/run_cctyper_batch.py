from __future__ import annotations

import argparse
import csv
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path


FASTA_EXTENSIONS = {".fa", ".faa", ".fasta", ".fna", ".ffn", ".frn"}
DEFAULT_GROUP_LABELS = {
    "resistant": "resistant",
    "susceptible": "susceptible",
}


@dataclass(frozen=True)
class CCTyperJob:
    fasta_path: Path
    output_dir: Path
    genome_id: str
    source_group: str
    phenotype_label: str


def discover_cctyper_jobs(
    input_root: str | Path,
    output_root: str | Path,
    group_labels: dict[str, str] | None = None,
) -> list[CCTyperJob]:
    input_root = Path(input_root)
    output_root = Path(output_root)
    group_labels = group_labels or DEFAULT_GROUP_LABELS

    if not input_root.exists():
        raise FileNotFoundError(f"Input folder does not exist: {input_root}")

    jobs: list[CCTyperJob] = []
    for fasta_path in sorted(input_root.rglob("*")):
        if not fasta_path.is_file() or fasta_path.suffix.lower() not in FASTA_EXTENSIONS:
            continue

        source_group = _source_group(input_root, fasta_path)
        phenotype_label = group_labels.get(source_group.lower(), "")
        genome_id = _genome_id(fasta_path)
        output_dir = output_root / source_group / f"{genome_id}_cctyper"
        jobs.append(
            CCTyperJob(
                fasta_path=fasta_path,
                output_dir=output_dir,
                genome_id=genome_id,
                source_group=source_group,
                phenotype_label=phenotype_label,
            )
        )

    return jobs


def run_cctyper_jobs(
    jobs: list[CCTyperJob],
    *,
    cctyper_command: str = "cctyper",
    threads: int = 4,
    db_path: str | Path | None = None,
    resume: bool = True,
    simplelog: bool = True,
    keep_tmp: bool = False,
) -> list[dict[str, str]]:
    results: list[dict[str, str]] = []

    for index, job in enumerate(jobs, start=1):
        job.output_dir.parent.mkdir(parents=True, exist_ok=True)
        expected_table = job.output_dir / "crisprs_near_cas.tab"

        if resume and expected_table.exists():
            results.append(_result_row(job, "skipped_existing", 0))
            print(f"[{index}/{len(jobs)}] skipped existing: {job.fasta_path}")
            continue

        command = [
            cctyper_command,
            "-t",
            str(threads),
            str(job.fasta_path),
            str(job.output_dir),
        ]
        if db_path:
            command.extend(["--db", str(db_path)])
        if simplelog:
            command.append("--simplelog")
        if keep_tmp:
            command.append("--keep_tmp")

        print(f"[{index}/{len(jobs)}] running: {job.fasta_path}")
        completed = subprocess.run(command, check=False)
        status = "completed" if completed.returncode == 0 else "failed"
        results.append(_result_row(job, status, completed.returncode))

    return results


def write_manifest(jobs: list[CCTyperJob], manifest_path: str | Path) -> None:
    rows = [
        {
            "cctyper_output_dir": str(job.output_dir),
            "genome_id": job.genome_id,
            "organism": "",
            "taxonomy": "",
            "assembly_level": "",
            "source_group": job.source_group,
            "phenotype_label": job.phenotype_label,
            "fasta_path": str(job.fasta_path),
        }
        for job in jobs
    ]
    _write_csv(manifest_path, rows)


def write_run_summary(rows: list[dict[str, str]], summary_path: str | Path) -> None:
    _write_csv(summary_path, rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Run CRISPRCasTyper/CCTyper over a folder of bacterial FASTA files "
            "and write a SABR-compatible cctyper_manifest.csv."
        )
    )
    parser.add_argument("input_root", type=Path, help="Folder containing bacterial FASTA files.")
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("data/curation/cctyper_runs"),
        help="Folder where CCTyper output directories will be written.",
    )
    parser.add_argument(
        "--manifest",
        type=Path,
        default=Path("data/curation/cctyper_manifest.csv"),
        help="CSV manifest for crispr_phage_predictor.ml.collect_cctyper_training.",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=Path("data/curation/cctyper_batch_summary.csv"),
        help="CSV run summary with completed/skipped/failed status.",
    )
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--cctyper-command", default="cctyper")
    parser.add_argument("--db", type=Path, default=None, help="Optional CCTyper database path.")
    parser.add_argument("--no-resume", action="store_true", help="Rerun jobs even if output exists.")
    parser.add_argument("--keep-tmp", action="store_true", help="Pass --keep_tmp to CCTyper.")
    parser.add_argument("--dry-run", action="store_true", help="Only write the manifest and summary.")
    args = parser.parse_args()

    jobs = discover_cctyper_jobs(args.input_root, args.output_root)
    if not jobs:
        raise SystemExit(f"No FASTA files found under {args.input_root}")

    write_manifest(jobs, args.manifest)
    print(f"Wrote manifest for {len(jobs)} FASTA files: {args.manifest}")

    if args.dry_run:
        rows = [_result_row(job, "dry_run", 0) for job in jobs]
    else:
        rows = run_cctyper_jobs(
            jobs,
            cctyper_command=args.cctyper_command,
            threads=args.threads,
            db_path=args.db,
            resume=not args.no_resume,
            keep_tmp=args.keep_tmp,
        )

    write_run_summary(rows, args.summary)
    print(f"Wrote run summary: {args.summary}")


def _source_group(input_root: Path, fasta_path: Path) -> str:
    relative = fasta_path.relative_to(input_root)
    if len(relative.parts) > 1:
        return _safe_name(relative.parts[0])
    return "un grouped".replace(" ", "_")


def _genome_id(fasta_path: Path) -> str:
    return _safe_name(fasta_path.stem)


def _safe_name(value: str) -> str:
    clean = re.sub(r"[^A-Za-z0-9_.-]+", "_", value.strip())
    clean = clean.strip("._-")
    return clean or "unnamed"


def _result_row(job: CCTyperJob, status: str, return_code: int) -> dict[str, str]:
    return {
        "fasta_path": str(job.fasta_path),
        "cctyper_output_dir": str(job.output_dir),
        "genome_id": job.genome_id,
        "source_group": job.source_group,
        "phenotype_label": job.phenotype_label,
        "status": status,
        "return_code": str(return_code),
    }


def _write_csv(path: str | Path, rows: list[dict[str, str]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
