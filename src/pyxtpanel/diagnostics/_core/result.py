"""``TestResult``: the one result object of every statistical test (L.1.1).

A test returns a :class:`TestResult` whatever its family, so that a script, a report or an
estimation package reads ``statistic``, ``pvalue`` and ``df`` the same way for an Im-Pesaran-Shin
test, a Hansen test or a Hausman test.  The texts (title, null and alternative hypotheses) are
not stored in the object: they come from the message catalogue under ``diag.<name>.title``,
``diag.<name>.null`` and ``diag.<name>.alternative``, so that one result prints in French or
in English on demand.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._results import Result

if TYPE_CHECKING:
    from pyxtpanel._types import Lang

_i18n.register(
    "diag",
    {
        "diag.result.statistic": {"en": "statistic", "fr": "statistique"},
        "diag.result.pvalue": {"en": "p-value", "fr": "p-value"},
        "diag.result.df": {"en": "degrees of freedom", "fr": "degrés de liberté"},
        "diag.result.n_obs": {"en": "observations", "fr": "observations"},
        "diag.result.n_groups": {"en": "units", "fr": "unités"},
        "diag.result.null": {"en": "Ho", "fr": "H0"},
        "diag.result.alternative": {"en": "Ha", "fr": "H1"},
        "diag.result.options": {"en": "Options", "fr": "Options"},
        "diag.result.notes": {"en": "Notes", "fr": "Remarques"},
        "diag.result.details": {"en": "Details", "fr": "Détail"},
        "diag.result.sample": {
            "en": "{n_obs} observations, {n_groups} units",
            "fr": "{n_obs} observations, {n_groups} unités",
        },
        "diag.result.sample_obs": {"en": "{n_obs} observations", "fr": "{n_obs} observations"},
        "diag.result.rejected": {
            "en": "Ho is rejected at the {level:.0%} level.",
            "fr": "H0 est rejetée au seuil de {level:.0%}.",
        },
        "diag.result.not_rejected": {
            "en": "Ho is not rejected at the {level:.0%} level.",
            "fr": "H0 n'est pas rejetée au seuil de {level:.0%}.",
        },
    },
)


@dataclass(frozen=True)
class Note:
    """A remark a test makes about its own run, kept as a message key so that it translates.

    Attributes
    ----------
    key
        Message key in the catalogue (``diag.<test>.<note>``).
    fields
        Values the message template needs.
    """

    key: str
    fields: Mapping[str, object] = field(default_factory=dict)

    def text(self, lang: Lang | None = None) -> str:
        """Return the remark in ``lang``."""
        return _i18n.message(self.key, lang, **self.fields)


@dataclass(frozen=True, eq=False, kw_only=True)
class TestResult(Result):
    """Outcome of a statistical test.

    Attributes
    ----------
    name
        Registered name of the test (``"ips"``, ``"hansen"``); the key of its texts.
    statistic
        The test statistic (the one the reference software reports as the headline).
    pvalue
        The p-value of ``statistic`` under the null hypothesis.
    statistic_name
        How the statistic is called (``"t-bar"``, ``"chi2"``, ``"z"``), a symbol that is not
        translated.
    df
        Degrees of freedom: none, one (chi-square, t) or two (F).
    n_obs, n_groups
        Observations and units that entered the test (``None`` when it does not apply).
    variable
        Name of the tested variable for a test on data.
    options
        The options the test ran with (lags, trend, kernel...), read-only.
    details
        Per-unit or per-variant table (the ADF statistic of each unit, the seven statistics of
        Pedroni), or ``None``.
    notes
        Remarks the test makes about its own run (a fallback used, a unit dropped): :class:`Note`
        objects, which translate, or plain text.
    """

    # Not a pytest test class despite the name.
    __test__ = False

    name: str
    statistic: float
    pvalue: float
    statistic_name: str = "statistic"
    df: tuple[float, ...] = ()
    n_obs: int | None = None
    n_groups: int | None = None
    variable: str | None = None
    options: Mapping[str, object] = field(default_factory=dict)
    details: pd.DataFrame | None = None
    notes: tuple[Note | str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "options", MappingProxyType(dict(self.options)))
        object.__setattr__(self, "df", tuple(float(d) for d in self.df))
        object.__setattr__(self, "statistic", float(self.statistic))
        object.__setattr__(self, "pvalue", float(self.pvalue))

    # ------------------------------------------------------------------ texts

    def title(self, lang: Lang | None = None) -> str:
        """Return the name of the test in ``lang``."""
        return _i18n.message(f"diag.{self.name}.title", lang)

    def null(self, lang: Lang | None = None) -> str:
        """Return the null hypothesis in ``lang``."""
        return _i18n.message(f"diag.{self.name}.null", lang)

    def alternative(self, lang: Lang | None = None) -> str:
        """Return the alternative hypothesis in ``lang``."""
        return _i18n.message(f"diag.{self.name}.alternative", lang)

    def rejects(self, level: float = 0.05) -> bool:
        """Return ``True`` when the null hypothesis is rejected at ``level`` (``pvalue < level``).

        A ``NaN`` p-value never rejects.
        """
        return bool(not np.isnan(self.pvalue) and self.pvalue < level)

    # ------------------------------------------------------------------ forms

    def to_frame(self) -> pd.DataFrame:
        """Return one row: statistic, p-value, degrees of freedom, sample sizes."""
        row: dict[str, object] = {"statistic": self.statistic, "pvalue": self.pvalue}
        for i, value in enumerate(self.df, start=1):
            row["df" if len(self.df) == 1 else f"df{i}"] = value
        if self.n_obs is not None:
            row["n_obs"] = self.n_obs
        if self.n_groups is not None:
            row["n_groups"] = self.n_groups
        return pd.DataFrame(row, index=pd.Index([self.name], name="test"))

    def summary(self, lang: Lang | None = None) -> str:
        """Return the report: title, hypotheses, statistic, p-value, sample and details."""
        head = self.title(lang) + (f" ({self.variable})" if self.variable else "")
        label = _i18n.message
        stat = f"{self.statistic_name} = {self.statistic:.4f}"
        if self.df:
            shown = ", ".join(f"{d:g}" for d in self.df)
            stat += f"   {label('diag.result.df', lang)} = {shown}"
        stat += f"   {label('diag.result.pvalue', lang)} = {self.pvalue:.4f}"
        lines = [
            head,
            f"{label('diag.result.null', lang)}: {self.null(lang)}",
            f"{label('diag.result.alternative', lang)}: {self.alternative(lang)}",
            "",
            stat,
        ]
        if self.n_obs is not None and self.n_groups is not None:
            lines.append(
                label("diag.result.sample", lang, n_obs=self.n_obs, n_groups=self.n_groups)
            )
        elif self.n_obs is not None:
            lines.append(label("diag.result.sample_obs", lang, n_obs=self.n_obs))
        if self.options:
            shown_options = ", ".join(f"{k}={v}" for k, v in self.options.items())
            lines.append(f"{label('diag.result.options', lang)}: {shown_options}")
        lines.extend(
            f"{label('diag.result.notes', lang)}: {self.note_text(n, lang)}" for n in self.notes
        )
        if self.details is not None and len(self.details):
            lines += ["", self.details.to_string(float_format=lambda v: f"{v:.4f}")]
        return "\n".join(lines)

    @staticmethod
    def note_text(note: Note | str, lang: Lang | None = None) -> str:
        """Return a remark as text in ``lang``."""
        return note.text(lang) if isinstance(note, Note) else note

    def _repr_html_(self) -> str:
        table = self.to_frame().to_html(float_format=lambda v: f"{v:.4f}")
        if self.details is None or not len(self.details):
            return table
        return table + self.details.to_html(float_format=lambda v: f"{v:.4f}")
