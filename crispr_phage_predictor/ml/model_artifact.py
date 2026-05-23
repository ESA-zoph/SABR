from __future__ import annotations

from dataclasses import dataclass
import hashlib
from pathlib import Path
from typing import Any

import joblib
import pandas as pd
from sklearn.ensemble import ExtraTreesClassifier

from crispr_phage_predictor.ml.dataset import load_repeat_cas_training_table
from crispr_phage_predictor.ml.features import build_repeat_feature_table, feature_columns
from crispr_phage_predictor.ml.train_classifier import _filter_min_class_count


DEFAULT_MODEL_PATH = Path("models") / "cas_subtype_extratrees.joblib"
DEFAULT_TRAINING_TABLE = (
    Path("data") / "training" / "repeats_cas_types_augmented_vink_genbank_targeted.csv"
)


@dataclass(frozen=True)
class CasSubtypeModelArtifact:
    model: Any
    feature_names: list[str]
    classes: list[str]
    metadata: dict[str, Any]


def train_extra_trees_artifact(
    training_table_path: Path = DEFAULT_TRAINING_TABLE,
    min_class_count: int = 20,
    random_state: int = 42,
    n_estimators: int = 400,
) -> CasSubtypeModelArtifact:
    table = load_repeat_cas_training_table(training_table_path)
    table = _filter_min_class_count(table, min_class_count=min_class_count)
    feature_table = build_repeat_feature_table(table)
    features = feature_columns(feature_table)
    model = ExtraTreesClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        class_weight="balanced",
        n_jobs=-1,
    )
    model.fit(feature_table[features], feature_table["cas_subtype"])
    metadata = {
        "method": "extra_trees",
        "training_table": str(training_table_path),
        "training_rows": int(len(table)),
        "min_class_count": int(min_class_count),
        "random_state": int(random_state),
        "n_estimators": int(n_estimators),
        "note": "Runtime artifact trained on all rows after min-class filtering; validation is documented separately.",
    }
    return CasSubtypeModelArtifact(
        model=model,
        feature_names=features,
        classes=[str(label) for label in model.classes_],
        metadata=metadata,
    )


def save_artifact(artifact: CasSubtypeModelArtifact, path: Path = DEFAULT_MODEL_PATH) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    joblib.dump(
        {
            "model": artifact.model,
            "feature_names": artifact.feature_names,
            "classes": artifact.classes,
            "metadata": artifact.metadata,
        },
        path,
    )
    return path


def load_artifact(path: Path = DEFAULT_MODEL_PATH) -> CasSubtypeModelArtifact:
    payload = joblib.load(path)
    return CasSubtypeModelArtifact(
        model=payload["model"],
        feature_names=list(payload["feature_names"]),
        classes=list(payload["classes"]),
        metadata=dict(payload["metadata"]),
    )


def model_artifact_metadata(path: Path = DEFAULT_MODEL_PATH) -> dict[str, Any]:
    """Return reproducibility metadata for a runtime model artifact."""
    metadata: dict[str, Any] = {
        "artifact_path": str(path),
        "artifact_exists": path.exists(),
        "artifact_sha256": "",
        "method": "",
        "training_table": "",
        "training_rows": None,
        "min_class_count": None,
        "random_state": None,
        "n_estimators": None,
        "classes": [],
        "load_error": "",
    }
    if not path.exists():
        return metadata

    metadata["artifact_sha256"] = _file_sha256(path)
    try:
        artifact = load_artifact(path)
    except Exception as exc:
        metadata["load_error"] = str(exc)
        return metadata

    artifact_metadata = artifact.metadata
    metadata.update(
        {
            "method": str(artifact_metadata.get("method", "")),
            "training_table": str(artifact_metadata.get("training_table", "")),
            "training_rows": artifact_metadata.get("training_rows"),
            "min_class_count": artifact_metadata.get("min_class_count"),
            "random_state": artifact_metadata.get("random_state"),
            "n_estimators": artifact_metadata.get("n_estimators"),
            "classes": artifact.classes,
        }
    )
    return metadata


def _file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prediction_table_for_arrays(arrays) -> pd.DataFrame:
    rows = []
    for array in arrays:
        rows.append(
            {
                "genome_id": array.genome_id,
                "contig_id": array.contig_id,
                "repeat_sequence": array.repeat_consensus,
                "repeat_length": array.repeat_length,
                "spacer_count": array.spacer_count,
                "mean_spacer_length": array.mean_spacer_length,
                "cas_type": "",
                "cas_subtype": "",
            }
        )
    return pd.DataFrame(rows)
