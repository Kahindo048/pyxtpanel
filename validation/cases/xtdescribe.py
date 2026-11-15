"""``xtdescribe`` against Stata (E.6.1): counts, distribution of T_i and the pattern table.

The patterns are not returned by Stata: the do-file writes the output of ``xtdescribe`` to a text
log and a Mata function parses the table (frequency, percent, cumulative percent, then one
numeric code per column: 0 for ``.``, the digit, 10 for ``X``; ``-1`` where the cell is empty).
"""

from __future__ import annotations

import warnings

import numpy as np

from validation import panels
from validation.registry import case, export, register_prelude

DOMAIN = "xtdescribe"
PANELS = (
    "balanced",
    "gaps",
    "entry_exit",
    "weak",
    "delta2",
    "quarterly",
    "monthly",
    "unsorted_ids",
    "one_obs_unit",
    "re_effects",
    "sparse_t",
    "wide",
)
OPTIONS = ((9, 100), (3, 100), (1000, 100), (0, 100), (9, 40))
STATS = [
    "N",
    "mean",
    "sd",
    "min",
    "p1",
    "p5",
    "p10",
    "p25",
    "p50",
    "p75",
    "p90",
    "p95",
    "p99",
    "max",
]

PARSER = """\
mata:
real matrix pyx_patterns(string scalar fname)
{
    real scalar fh, p, k, j, nT
    string scalar line, left, right, c
    string rowvector tok
    real matrix M
    real rowvector row
    fh = fopen(fname, "r")
    nT = 0
    while ((line = fget(fh)) != J(0,0,"")) {
        p = strpos(line, "|")
        if (p == 0) continue
        right = strtrim(substr(line, p + 1, .))
        if (substr(right, 1, 1) != "(" & right != "Pattern" & substr(right, 1, 7) != "Pattern") nT = max((nT, strlen(right)))
    }
    fclose(fh)
    fh = fopen(fname, "r")
    M = J(0, 3 + nT, .)
    while ((line = fget(fh)) != J(0,0,"")) {
        p = strpos(line, "|")
        if (p == 0) continue
        left = strtrim(substr(line, 1, p - 1))
        right = strtrim(substr(line, p + 1, .))
        tok = tokens(left)
        if (cols(tok) < 2 | cols(tok) > 3) continue
        if (missing(strtoreal(tok[1]))) continue
        row = J(1, 3 + nT, -1)
        for (k = 1; k <= cols(tok); k++) row[k] = strtoreal(tok[k])
        if (substr(right, 1, 1) != "(") {
            for (j = 1; j <= strlen(right); j++) {
                c = substr(right, j, 1)
                row[3 + j] = (c == "." ? 0 : (c == "X" ? 10 : strtoreal(c)))
            }
        }
        M = M \\ row
    }
    fclose(fh)
    return(M)
}
end
"""

register_prelude(DOMAIN, PARSER)


def _option(name: str) -> str:
    delta = panels.SPECS[name].delta
    return f", delta({delta})" if delta != 1 else ""


def _python(name: str, patterns: int, width: int) -> dict[str, np.ndarray]:
    from pyxtpanel import xtset
    from pyxtpanel.describe.patterns import participation_patterns

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # sparse_t: every period looks like a gap with delta=1
        panel = xtset(panels.load(name), "id", "t", delta=panels.SPECS[name].delta)
    result = panel.xtdescribe(patterns=patterns, width=width)
    ti = result.ti[STATS].to_numpy().reshape(-1, 1)
    out: dict[str, np.ndarray] = {"ti": ti}
    if patterns == 0:
        return out
    found = participation_patterns(panel.index, width)
    n_columns = found.counts.shape[1]
    shown = result.patterns
    rows = []
    for k, row in enumerate(shown.itertuples()):
        rows.append(
            [
                row.freq,
                round(row.percent, 2),
                round(row.cum_percent, 2),
                *found.counts[k].clip(max=10),
            ]
        )
    if result.other_freq:
        share = round(100 * result.other_freq / result.n_units, 2)
        rows.append([result.other_freq, share, 100.0, *([-1.0] * n_columns)])
    footer = [10 if on else 0 for on in found.footer]
    rows.append([result.n_units, 100.0, -1.0, *footer])
    out["patterns"] = np.array(rows, dtype=float)
    return out


def _register(name: str, patterns: int, width: int) -> None:
    tag = f"xtdescribe_{name}_p{patterns}_w{width}"
    log = f"../out/{tag}.log"
    scalars = [f"scalar s_{s} = r({s if s != 'N' else 'N'})" for s in STATS]
    lines = [
        f"xtset id t{_option(name)}",
        "set linesize 255",
        f'quietly log using "{log}", text replace name(pyxlog)',
        f"xtdescribe, patterns({patterns}) width({width})",
        *scalars,
        "quietly log close pyxlog",
        export(tag, "ti", "(" + " \\ ".join(f'st_numscalar("s_{s}")' for s in STATS) + ")"),
    ]
    if patterns != 0:
        lines.append(export(tag, "patterns", f'pyx_patterns("{log}")'))

    @case(tag, domain=DOMAIN, dataset=f"sim_{name}", stata=lines, atol=1e-9, rtol=1e-9)
    def _run(
        name: str = name, patterns: int = patterns, width: int = width
    ) -> dict[str, np.ndarray]:
        return _python(name, patterns, width)


for _name in PANELS:
    for _patterns, _width in OPTIONS:
        if _name == "wide" and _width == 100 and _patterns not in (9, 1000):
            continue
        _register(_name, _patterns, _width)
