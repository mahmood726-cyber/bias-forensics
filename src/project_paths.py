from __future__ import annotations

import os
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTPUT_DIR = PROJECT_ROOT / "data" / "output"


def _nested_pairwise_dirs(parent: Path) -> list[Path]:
    if not parent.exists() or not parent.is_dir():
        return []
    candidates: list[Path] = []
    for child in parent.iterdir():
        if not child.is_dir():
            continue
        nested = child / "Pairwise70" / "data"
        candidates.append(nested)
    return candidates


def _candidate_pairwise_dirs() -> list[Path]:
    env_candidates = [
        os.environ.get("PAIRWISE70_DATA_DIR"),
        os.environ.get("PAIRWISE70_DIR"),
    ]
    candidates = [Path(candidate) for candidate in env_candidates if candidate]
    projects_root = PROJECT_ROOT.parent / "Projects"
    models_root = PROJECT_ROOT.parent / "Models"
    candidates.extend(
        [
            PROJECT_ROOT / "data" / "pairwise70",
            projects_root / "Pairwise70" / "data",
            models_root / "Pairwise70" / "data",
            PROJECT_ROOT.parent / "Pairwise70" / "data",
        ]
    )
    candidates.extend(_nested_pairwise_dirs(projects_root))
    candidates.extend(_nested_pairwise_dirs(models_root))
    return candidates


def resolve_pairwise_dir(pairwise_dir: str | Path | None = None, require_exists: bool = False) -> Path:
    if pairwise_dir is not None:
        resolved = Path(pairwise_dir)
    else:
        resolved = next((candidate for candidate in _candidate_pairwise_dirs() if candidate.exists()), _candidate_pairwise_dirs()[1])

    if require_exists and not resolved.exists():
        raise FileNotFoundError(
            f"Pairwise70 data directory not found: {resolved}. "
            "Set PAIRWISE70_DATA_DIR to the extracted Pairwise70 data directory."
        )
    return resolved
