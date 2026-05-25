import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.project_paths import resolve_pairwise_dir


def test_resolve_pairwise_dir_finds_nested_projects_pairwise_dir(tmp_path, monkeypatch):
    fake_project_root = tmp_path / "BiasForensics"
    nested_pairwise = tmp_path / "Projects" / "mahmood789" / "Pairwise70" / "data"
    nested_pairwise.mkdir(parents=True)

    monkeypatch.setattr("src.project_paths.PROJECT_ROOT", fake_project_root)

    resolved = resolve_pairwise_dir(require_exists=True)

    assert resolved == nested_pairwise


def test_explicit_pairwise_dir_still_wins(tmp_path):
    explicit = tmp_path / "custom" / "pairwise70"
    explicit.mkdir(parents=True)

    resolved = resolve_pairwise_dir(pairwise_dir=explicit, require_exists=True)

    assert resolved == explicit
