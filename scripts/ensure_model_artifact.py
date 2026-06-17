from __future__ import annotations

import hashlib
import os
from pathlib import Path
from urllib.error import URLError
from urllib.request import urlopen


DEFAULT_MODEL_PATH = Path("models") / "cas_subtype_extratrees.joblib"
CHUNK_SIZE = 1024 * 1024


def main() -> int:
    model_path = Path(os.environ.get("SABR_MODEL_PATH", str(DEFAULT_MODEL_PATH)))
    model_url = os.environ.get("SABR_MODEL_URL", "").strip()
    expected_sha256 = os.environ.get("SABR_MODEL_SHA256", "").strip().lower()

    if model_path.exists():
        if expected_sha256 and file_sha256(model_path) != expected_sha256:
            print(f"Existing model artifact at {model_path} failed SHA-256 validation.")
            model_path.unlink()
        else:
            print(f"Using model artifact at {model_path}.")
            return 0

    if not model_url:
        print(
            "No runtime model artifact found and SABR_MODEL_URL is not set. "
            "SABR will start with its documented missing-artifact fallback."
        )
        return 0

    model_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = model_path.with_suffix(model_path.suffix + ".download")
    print(f"Downloading SABR model artifact to {model_path}...")

    try:
        download(model_url, tmp_path)
    except (OSError, URLError) as exc:
        print(f"Model artifact download failed: {exc}")
        return 1

    if expected_sha256:
        actual_sha256 = file_sha256(tmp_path)
        if actual_sha256 != expected_sha256:
            tmp_path.unlink(missing_ok=True)
            print(
                "Downloaded model artifact failed SHA-256 validation: "
                f"expected {expected_sha256}, got {actual_sha256}."
            )
            return 1

    tmp_path.replace(model_path)
    print(f"Model artifact ready at {model_path}.")
    return 0


def download(url: str, path: Path) -> None:
    with urlopen(url, timeout=60) as response, path.open("wb") as output:
        while True:
            chunk = response.read(CHUNK_SIZE)
            if not chunk:
                break
            output.write(chunk)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


if __name__ == "__main__":
    raise SystemExit(main())
