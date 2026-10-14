from pathlib import Path

from mkdocs.commands.build import build
from mkdocs.config import load_config


def test_docs_build_strict(tmp_path):
    config = load_config(
        config_file=str(Path(__file__).resolve().parents[1] / "mkdocs.yml"),
        site_dir=str(tmp_path / "site"),
        strict=True,
    )

    build(config)

    assert (tmp_path / "site" / "index.html").is_file()
