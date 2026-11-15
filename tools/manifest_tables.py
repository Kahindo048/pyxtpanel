"""Record the checksums of the embedded tables in ``_data/tables/MANIFEST.json`` (L.1.3).

    python -m tools.manifest_tables --write    # update the SHA-256 of the listed files
    python -m tools.manifest_tables --check    # fail if a file and the manifest disagree

The provenance fields (title, source, usage, transcribed, verified_against) are written by hand
in the manifest when a table is added; this tool only maintains the checksums and refuses a
table file that has no entry.
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

DIRECTORY = Path(__file__).resolve().parents[1] / "src" / "pyxtpanel" / "_data" / "tables"
REQUIRED = ("title", "source", "usage", "transcribed", "verified_against")


def digest(path: Path) -> str:
    """Return the SHA-256 of a file with Windows line endings read as Unix ones."""
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main(argv: list[str]) -> int:
    """Check or update the manifest; return the exit code."""
    manifest_path = DIRECTORY / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    problems = []
    files = sorted(p.name for p in DIRECTORY.glob("*.csv"))
    for name in files:
        entry = manifest.get(name)
        if entry is None:
            problems.append(f"{name}: no entry in the manifest (add its provenance by hand)")
            continue
        missing = [key for key in REQUIRED if not entry.get(key)]
        if missing:
            problems.append(f"{name}: missing fields {missing}")
        found = digest(DIRECTORY / name)
        if found != entry.get("sha256"):
            if "--write" in argv:
                entry["sha256"] = found
            else:
                problems.append(f"{name}: checksum differs from the manifest")
    for name in set(manifest) - set(files):
        problems.append(f"{name}: in the manifest but the file does not exist")
    if "--write" in argv and not problems:
        manifest_path.write_text(
            json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        print("manifest updated")
        return 0
    for problem in problems:
        print(problem)
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
