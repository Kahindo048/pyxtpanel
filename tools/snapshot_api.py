"""Photograph of the public API (A.5.7, G.1.4): names, signatures and the array contract.

    python -m tools.snapshot_api            # compare with tests/contract/api_snapshot.json
    python -m tools.snapshot_api --write    # accept the current API as the new snapshot

A change that removes or alters a name is a breaking change: it needs a major version (or, before
1.0, a note in CHANGELOG.md and a deprecation period).  Adding a name only needs the snapshot to
be rewritten.
"""

from __future__ import annotations

import dataclasses
import inspect
import json
import sys
from pathlib import Path

TARGET = Path(__file__).resolve().parents[1] / "tests" / "contract" / "api_snapshot.json"


def snapshot() -> dict[str, object]:
    """Describe the public API as plain data."""
    import pyxtpanel
    from pyxtpanel.contract import CONTRACT_VERSION, PanelArrays

    functions = {}
    for name in pyxtpanel.__all__:
        value = getattr(pyxtpanel, name)
        if inspect.isfunction(value):
            functions[name] = str(inspect.signature(value))
    methods = {
        name: str(inspect.signature(member))
        for name, member in inspect.getmembers(pyxtpanel.PanelData, inspect.isfunction)
        if not name.startswith("_")
    }
    return {
        "names": sorted(pyxtpanel.__all__),
        "functions": functions,
        "panel_data_methods": methods,
        "contract_version": CONTRACT_VERSION,
        "panel_arrays_fields": [f.name for f in dataclasses.fields(PanelArrays)],
    }


def main(argv: list[str]) -> int:
    current = json.dumps(snapshot(), indent=2, sort_keys=True) + "\n"
    if "--write" in argv:
        TARGET.write_text(current, encoding="utf-8", newline="\n")
        print(f"wrote {TARGET}")
        return 0
    if not TARGET.exists() or TARGET.read_text(encoding="utf-8") != current:
        print("the public API differs from tests/contract/api_snapshot.json")
        return 1
    print("the public API matches the snapshot")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
