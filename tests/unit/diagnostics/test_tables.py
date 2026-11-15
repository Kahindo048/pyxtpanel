"""Embedded tables: manifest, checksum, and agreement with an independent implementation."""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

import pytest

from pyxtpanel._errors import DataError
from pyxtpanel.diagnostics._shared import tables
from pyxtpanel.diagnostics._shared.adf import AdfFit, adf_pvalue

GOLDEN = Path(__file__).resolve().parents[2] / "golden" / "mackinnon_statsmodels.json"
REQUIRED = ("title", "source", "usage", "transcribed", "verified_against", "sha256")


def test_every_table_file_is_in_the_manifest_with_its_provenance() -> None:
    files = {p.name for p in tables._directory().iterdir() if p.suffix == ".csv"}
    manifest = tables.manifest()
    assert files == set(manifest)
    for name, entry in manifest.items():
        assert all(entry.get(key) for key in REQUIRED), name
        assert (
            entry["sha256"]
            == hashlib.sha256(
                (tables._directory() / name).read_bytes().replace(b"\r\n", b"\n")
            ).hexdigest()
        )


def test_a_modified_table_is_refused(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    for item in tables._directory().iterdir():
        shutil.copy(item, tmp_path / item.name)
    monkeypatch.setattr(tables, "_directory", lambda: tmp_path)
    tables.manifest.cache_clear()
    tables._read.cache_clear()
    try:
        assert not tables.load_table("mackinnon_adf.csv").empty
        tables._read.cache_clear()
        csv = tmp_path / "mackinnon_adf.csv"
        csv.write_text(csv.read_text().replace("2.1659", "2.1660"), encoding="utf-8")
        with pytest.raises(DataError, match="does not match the checksum"):
            tables.load_table("mackinnon_adf.csv")
        with pytest.raises(DataError, match="not in the manifest"):
            tables.load_table("unknown.csv")
    finally:
        tables.manifest.cache_clear()
        tables._read.cache_clear()


def test_windows_line_endings_do_not_change_the_checksum(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    for item in tables._directory().iterdir():
        shutil.copy(item, tmp_path / item.name)
    csv = tmp_path / "mackinnon_adf.csv"
    csv.write_bytes(csv.read_bytes().replace(b"\n", b"\r\n"))
    monkeypatch.setattr(tables, "_directory", lambda: tmp_path)
    tables.manifest.cache_clear()
    tables._read.cache_clear()
    try:
        assert len(tables.load_table("mackinnon_adf.csv")) == 2
    finally:
        tables.manifest.cache_clear()
        tables._read.cache_clear()


def test_the_mackinnon_table_matches_statsmodels() -> None:
    reference = json.loads(GOLDEN.read_text(encoding="utf-8"))["cases"]
    table = tables.load_table("mackinnon_adf.csv").set_index("case")
    for case, expected in reference.items():
        row = table.loc[case]
        for key in ("tau_star", "tau_min", "tau_max"):
            assert row[key] == pytest.approx(expected[key], abs=1e-12)
        for i in range(4):
            assert row[f"small_g{i}"] == pytest.approx(expected["small"][i], abs=1e-12)
            assert row[f"large_g{i}"] == pytest.approx(expected["large"][i], abs=1e-12)


def test_the_mackinnon_p_values_match_statsmodels() -> None:
    reference = json.loads(GOLDEN.read_text(encoding="utf-8"))["cases"]
    for case, expected in reference.items():
        name = "constant" if case == "c" else "trend"
        for tau, p in expected["p"].items():
            fit = AdfFit(float(tau), 0.0, 20, 15, 1)
            assert adf_pvalue(fit, name) == pytest.approx(p, abs=1e-12)  # type: ignore[arg-type]


def test_returned_tables_are_copies() -> None:
    first = tables.load_table("mackinnon_adf.csv")
    first.loc[0, "tau_star"] = 99.0
    assert tables.load_table("mackinnon_adf.csv").loc[0, "tau_star"] == -1.61
