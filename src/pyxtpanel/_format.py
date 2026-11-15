"""Plain-text, Markdown and LaTeX rendering of tables, shared by every result (layer 0)."""

from __future__ import annotations

import numpy as np
import pandas as pd

_LATEX_ESCAPES = {
    "\\": r"\textbackslash{}",
    "&": r"\&",
    "%": r"\%",
    "_": r"\_",
    "#": r"\#",
    "$": r"\$",
}


def number(value: float, digits: int = 4) -> str:
    """Format like Stata's ``%g``: ``.`` for missing, significant digits otherwise."""
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return "."
    return f"{value:.{digits}g}"


def _flat(frame: pd.DataFrame) -> pd.DataFrame:
    """Frame with the index turned into columns and one-line column names."""
    flat = frame.reset_index()
    flat.columns = [
        " ".join(str(part) for part in col if str(part)) if isinstance(col, tuple) else str(col)
        for col in flat.columns
    ]
    return flat


def _text(value: object, digits: int) -> str:
    if isinstance(value, (float, np.floating)):
        return number(float(value), digits)
    return str(value)


def markdown_table(frame: pd.DataFrame, digits: int = 4) -> str:
    """Render a DataFrame as a GitHub Markdown table (index columns first)."""
    flat = _flat(frame)
    rows = [[_text(v, digits) for v in row] for row in flat.itertuples(index=False)]
    header = "| " + " | ".join(flat.columns) + " |"
    rule = "|" + "|".join("---" for _ in flat.columns) + "|"
    return "\n".join([header, rule, *("| " + " | ".join(r) + " |" for r in rows)])


def latex_escape(text: str) -> str:
    """Escape the characters that are special in LaTeX."""
    return "".join(_LATEX_ESCAPES.get(c, c) for c in text)


def latex_table(frame: pd.DataFrame, digits: int = 4) -> str:
    """Render a DataFrame as a LaTeX ``tabular`` (index first, special characters escaped)."""
    flat = _flat(frame)

    def cell(value: object) -> str:
        return "".join(_LATEX_ESCAPES.get(c, c) for c in _text(value, digits))

    end = r" \\"
    head = " & ".join(cell(c) for c in flat.columns) + end
    rows = [" & ".join(cell(v) for v in row) + end for row in flat.itertuples(index=False)]
    begin = r"\begin{tabular}{" + "l" * flat.shape[1] + "}"
    return "\n".join([begin, r"\hline", head, r"\hline", *rows, r"\hline", r"\end{tabular}"])
