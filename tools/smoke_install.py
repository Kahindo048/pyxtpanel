"""Smoke test of an INSTALLED package (J.3.3): run it from outside the repository, with ``-I``.

    python -I tools/smoke_install.py            # in an environment where only the wheel is installed

It checks that the package comes from ``site-packages`` (not from ``src/``), that the data files of
the wheel are there (the numerical tables of the unit-root tests), and that the main paths work with
numpy, pandas and scipy only.  Matplotlib may be present or not: the figures either work or refuse
with a message that names the extra.
"""

from __future__ import annotations

import importlib.util
import sys

import numpy as np
import pandas as pd

import pyxtpanel as xt


def main() -> None:
    location = str(xt.__file__)
    assert "site-packages" in location or "dist-packages" in location, f"not installed: {location}"
    print("pyxtpanel", xt.__version__, "from", location)

    rng = np.random.default_rng(0)
    frame = pd.DataFrame(
        {
            "id": np.repeat(np.arange(1, 21), 8),
            "t": np.tile(np.arange(2000, 2008), 20),
            "x": rng.normal(size=160),
            "y": rng.normal(size=160),
        }
    ).drop(index=[3, 4, 50, 51, 90])  # holes
    panel = xt.xtset(frame, "id", "t")
    assert panel.n_units == 20 and panel.has_gaps

    assert panel.xtdescribe().summary()
    assert panel.xtsum("x").to_frame().shape[0] == 3
    assert panel.lag("x").isna().sum() > 0
    assert panel.within("x").shape[0] == len(frame)
    assert panel.audit().summary().startswith("Panel audit")

    # the embedded numerical tables are shipped and verified against their checksums
    result = xt.run_test("fisher", panel, "x", lags=1)
    assert 0.0 <= float(result.pvalue) <= 1.0
    assert len(xt.list_tests()) >= 35

    if importlib.util.find_spec("matplotlib") is None:
        try:
            panel.plot("lines", variable="x")
        except ImportError as error:
            assert "pyxtpanel[plot]" in str(error), error
        else:
            raise AssertionError("plot should need matplotlib")
    print("python", sys.version.split()[0], "- smoke test OK")


if __name__ == "__main__":
    main()
