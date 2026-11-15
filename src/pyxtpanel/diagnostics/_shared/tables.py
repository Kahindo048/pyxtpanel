"""Numerical tables the tests read (L.1.3): published coefficients and critical values.

``scipy`` has no Dickey-Fuller distribution, no moments of the Im-Pesaran-Shin statistic and no
Levin-Lin-Chu adjustments.  They are small published tables, kept as CSV files in
``pyxtpanel/_data/tables/`` with a ``MANIFEST.json`` that records, for each file, its source (the
paper and the table), the date of the transcription, what it was verified against and the SHA-256
of its content.  :func:`load_table` refuses a file that is missing from the manifest or whose
content no longer matches it, so a table cannot change without the change being recorded.
"""

from __future__ import annotations

import hashlib
import io
import json
from functools import cache
from importlib import resources
from pathlib import Path

import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError

_i18n.register(
    "diag",
    {
        "diag.table.unknown": {
            "en": "table {name!r} is not in the manifest of the embedded tables.",
            "fr": "la table {name!r} ne figure pas dans le manifeste des tables embarquées.",
        },
        "diag.table.corrupt": {
            "en": "table {name!r} does not match the checksum recorded in the manifest "
            "(expected {expected}, found {found}); it was modified without updating the manifest.",
            "fr": "la table {name!r} ne correspond pas à la somme de contrôle du manifeste "
            "(attendue {expected}, trouvée {found}) ; elle a été modifiée sans mettre le "
            "manifeste à jour.",
        },
    },
)


def _directory() -> Path:
    """Return the folder of the embedded tables (the package is installed as plain files)."""
    return Path(str(resources.files("pyxtpanel"))) / "_data" / "tables"


def _digest(content: bytes) -> str:
    """Return the SHA-256 of ``content`` with Windows line endings read as Unix ones."""
    return hashlib.sha256(content.replace(b"\r\n", b"\n")).hexdigest()


@cache
def manifest() -> dict[str, dict[str, str]]:
    """Return the manifest: ``{file name: {title, source, usage, transcribed, ..., sha256}}``."""
    text = (_directory() / "MANIFEST.json").read_text(encoding="utf-8")
    data: dict[str, dict[str, str]] = json.loads(text)
    return data


def load_table(name: str) -> pd.DataFrame:
    """Return a copy of the table ``name`` (a file of ``_data/tables/``).

    Example: ``load_table("mackinnon_adf.csv")``.

    Raises
    ------
    DataError
        If the file is not in the manifest or its checksum does not match.
    """
    return _read(name).copy()


@cache
def _read(name: str) -> pd.DataFrame:
    """Read, check and cache the table ``name``."""
    entry = manifest().get(name)
    if entry is None:
        raise DataError(_i18n.message("diag.table.unknown", name=name), key="diag.table.unknown")
    content = (_directory() / name).read_bytes()
    found = _digest(content)
    if found != entry["sha256"]:
        raise DataError(
            _i18n.message(
                "diag.table.corrupt", name=name, expected=entry["sha256"][:12], found=found[:12]
            ),
            key="diag.table.corrupt",
        )
    return pd.read_csv(io.BytesIO(content))
