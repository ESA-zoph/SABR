from __future__ import annotations

import argparse
import importlib.util
import os
import shutil
from dataclasses import dataclass
from pathlib import Path


REQUIRED_COMMANDS = ["prodigal", "hmmsearch", "minced", "blastn"]
EXPECTED_DB_FILES = ["Cas.profile.hmm", "Profiles", "type_dict.tab"]


@dataclass(frozen=True)
class EnvironmentCheck:
    name: str
    ok: bool
    detail: str


def check_cctyper_environment(db_path: str | Path | None = None) -> list[EnvironmentCheck]:
    cctyper_command = shutil.which("cctyper")
    cctyper_package = importlib.util.find_spec("cctyper")
    checks = [
        EnvironmentCheck(
            name="cctyper",
            ok=cctyper_command is not None or cctyper_package is not None,
            detail=cctyper_command or "python package importable",
        )
    ]
    for command in REQUIRED_COMMANDS:
        checks.append(
            EnvironmentCheck(
                name=command,
                ok=shutil.which(command) is not None,
                detail=shutil.which(command) or "not found on PATH",
            )
        )

    raw_db_path = str(db_path or os.environ.get("CCTYPER_DB", "")).strip()
    resolved_db_path = Path(raw_db_path) if raw_db_path else None
    db_ok = resolved_db_path is not None and resolved_db_path.exists()
    checks.append(
        EnvironmentCheck(
            name="CCTYPER_DB",
            ok=db_ok,
            detail=str(resolved_db_path) if db_ok else "not set or path does not exist",
        )
    )

    if db_ok and resolved_db_path is not None:
        for expected in EXPECTED_DB_FILES:
            candidate = resolved_db_path / expected
            checks.append(
                EnvironmentCheck(
                    name=f"db:{expected}",
                    ok=candidate.exists(),
                    detail=str(candidate) if candidate.exists() else "missing",
                )
            )
    return checks


def main() -> None:
    parser = argparse.ArgumentParser(description="Check local CCTyper runtime requirements.")
    parser.add_argument("--db", type=Path, default=None, help="Optional CCTyper database path")
    args = parser.parse_args()

    checks = check_cctyper_environment(args.db)
    for check in checks:
        status = "OK" if check.ok else "MISSING"
        print(f"{status}\t{check.name}\t{check.detail}")

    if not all(check.ok for check in checks):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
