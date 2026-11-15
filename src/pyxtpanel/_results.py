"""Base of the result objects of the package (layer 0).

Results are frozen dataclasses.  They are read by the estimation packages, so fields are only
ever added (with a default), never renamed or removed within a minor version.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

from pyxtpanel._types import Lang

if TYPE_CHECKING:
    import pandas as pd


@dataclass(frozen=True, kw_only=True)
class Result:
    """Base class of result objects: immutable, with a table form and a text summary."""

    def to_frame(self) -> pd.DataFrame:
        """Return the result as a DataFrame."""
        raise NotImplementedError

    def to_markdown(self, digits: int = 4) -> str:
        """Return the table form as a Markdown table."""
        from pyxtpanel._format import markdown_table

        return markdown_table(self.to_frame(), digits)

    def to_latex(self, digits: int = 4) -> str:
        """Return the table form as a LaTeX ``tabular``."""
        from pyxtpanel._format import latex_table

        return latex_table(self.to_frame(), digits)

    def summary(self, lang: Lang | None = None) -> str:
        """Return a readable text report in ``lang`` (default: the global ``lang`` option)."""
        raise NotImplementedError

    def __str__(self) -> str:
        return self.summary()

    def _repr_html_(self) -> str:
        return self.to_frame().to_html()
