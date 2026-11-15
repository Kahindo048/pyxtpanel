"""Variance inflation factors of the regressors, from bare arrays (category B, descriptive, L.8.9).

Not a test (there is no hypothesis and no p-value): it measures how much the variance of each
coefficient is inflated by the correlation of its regressor with the others.  Validated against
Stata 17: ``estat vif`` after ``regress`` (centered, with correlated regressors, with a dummy,
and ``uncentered``); the regressors are exported from Stata and given to the function.  The
cases are in ``validation/cases/diagnostics_specification_vif.py``.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import NumericalError, SpecificationError
from pyxtpanel._results import Result
from pyxtpanel._types import FloatArray, Lang
from pyxtpanel.diagnostics._core import as_array
from pyxtpanel.diagnostics._shared.ols import RANK_TOLERANCE, ols_fit

_i18n.register(
    "diag",
    {
        "diag.vif.summary": {
            "en": "Variance inflation factors ({kind}): {count} variables, mean {mean:.2f}, "
            "largest {largest:.2f} ({name}).",
            "fr": "Facteurs d'inflation de la variance ({kind}) : {count} variables, moyenne "
            "{mean:.2f}, plus grand {largest:.2f} ({name}).",
        },
        "diag.vif.centered": {"en": "centered", "fr": "centrés"},
        "diag.vif.uncentered": {"en": "uncentered", "fr": "non centrés"},
        "diag.vif.too_few": {
            "en": "the factors need at least 2 regressors, got {count}.",
            "fr": "les facteurs exigent au moins 2 régresseurs, reçu {count}.",
        },
        "diag.vif.names": {
            "en": "names must have one entry per column ({count}), got {got}.",
            "fr": "names doit avoir une entrée par colonne ({count}), reçu {got}.",
        },
        "diag.vif.constant": {
            "en": "regressor {name!r} is constant: its factor is not defined (the constant is "
            "added by the regression, leave it out).",
            "fr": "le régresseur {name!r} est constant : son facteur n'est pas défini (la "
            "constante est ajoutée par la régression, ne la fournissez pas).",
        },
        "diag.vif.collinear": {
            "en": "regressor {name!r} is an exact linear combination of the others: its factor "
            "is infinite (Stata omits such a variable from the model; drop it first).",
            "fr": "le régresseur {name!r} est une combinaison linéaire exacte des autres : son "
            "facteur est infini (Stata omet une telle variable du modèle ; retirez-la "
            "d'abord).",
        },
    },
)


@dataclass(frozen=True, kw_only=True)
class VarianceInflation(Result):
    """Result of :func:`vif`.

    Attributes
    ----------
    vif
        Variance inflation factor of each regressor, ``1 / (1 - R^2)`` of its regression on the
        others (sorted from the largest, unless ``sort=False``).
    uncentered
        Whether the regressions have no constant and the constant is itself a row (``"_cons"``).
    """

    vif: pd.Series
    uncentered: bool = False

    @property
    def tolerance(self) -> pd.Series:
        """``1 / VIF``: the share of the regressor not explained by the others."""
        result: pd.Series = 1.0 / self.vif
        return result.rename("1/VIF")

    @property
    def mean(self) -> float:
        """Mean of the factors."""
        return float(self.vif.mean())

    def to_frame(self) -> pd.DataFrame:
        """Return the factors and their inverses, one row per regressor."""
        return pd.DataFrame({"VIF": self.vif, "1/VIF": self.tolerance})

    def summary(self, lang: Lang | None = None) -> str:
        """Return the mean and the largest factor as text (the table is :meth:`to_frame`)."""
        position = int(np.argmax(self.vif.to_numpy()))
        kind = _i18n.message(
            "diag.vif.uncentered" if self.uncentered else "diag.vif.centered", lang
        )
        return _i18n.message(
            "diag.vif.summary",
            lang,
            kind=kind,
            count=len(self.vif),
            mean=self.mean,
            largest=float(self.vif.iloc[position]),
            name=self.vif.index[position],
        )


def vif(
    x: object, *, names: object = None, uncentered: bool = False, sort: bool = True
) -> VarianceInflation:
    """Return the variance inflation factor of each column of ``x``.

    Each regressor is regressed on the others (with a constant) and its factor is
    ``1 / (1 - R^2)`` (Stata's ``estat vif``).  With ``uncentered`` the constant is one of the
    regressors, the regressions have no constant and ``R^2`` is the uncentered one (Stata's
    ``estat vif, uncentered``), which gives a factor to the constant too.

    Parameters
    ----------
    x
        ``(n, k)`` regressors, ``k >= 2``, without the constant.  A DataFrame supplies the names.
    names
        Names of the columns (default: the DataFrame columns, else ``x1``, ``x2``, ...).
    uncentered
        Uncentered factors, with ``"_cons"`` as the last regressor.
    sort
        Sort from the largest factor (ties by name); ``False`` keeps the order of the columns
        (Stata's ``nosort``).

    Raises
    ------
    SpecificationError
        If there are fewer than 2 columns or the names do not match.
    NumericalError
        If a regressor is constant or an exact linear combination of the others.
    """
    # Stata: estat vif -> for each regressor, regress it on the others, vif = 1 / (1 - e(r2));
    # uncentered: the constant joins the regressors and the regressions have no constant
    labels = list(x.columns) if isinstance(x, pd.DataFrame) else None
    z = as_array("x", x, ndim=2)
    count = z.shape[1]
    if names is not None:
        labels = [str(v) for v in np.asarray(names).ravel()]
        if len(labels) != count:
            raise SpecificationError(
                _i18n.message("diag.vif.names", count=count, got=len(labels)),
                key="diag.vif.names",
            )
    labels = labels if labels is not None else [f"x{i + 1}" for i in range(count)]
    if count < 2:
        raise SpecificationError(
            _i18n.message("diag.vif.too_few", count=count), key="diag.vif.too_few"
        )
    if uncentered:
        z = np.column_stack([z, np.ones(z.shape[0])])
        labels = [*labels, "_cons"]
    factors = np.array([_factor(z, j, labels[j], uncentered) for j in range(z.shape[1])])
    series = pd.Series(factors, index=pd.Index(labels, name="variable"), name="VIF")
    if sort:
        order = sorted(range(len(labels)), key=lambda i: (-factors[i], labels[i]))
        series = series.iloc[order]
    return VarianceInflation(vif=series, uncentered=uncentered)


def _factor(z: FloatArray, j: int, name: str, uncentered: bool) -> float:
    """Return ``1 / (1 - R^2)`` of column ``j`` on the others."""
    y = z[:, j]
    others = np.delete(z, j, axis=1)
    if uncentered:
        total = float(y @ y)
        norm = np.sqrt(np.sum(others**2, axis=0))
        scaled = others / norm
        coefficients, *_ = np.linalg.lstsq(scaled, y, rcond=RANK_TOLERANCE)
        resid = y - scaled @ coefficients
        unexplained = float(resid @ resid) / total
    else:
        if np.ptp(y) == 0.0:
            raise NumericalError(
                _i18n.message("diag.vif.constant", name=name), key="diag.vif.constant"
            )
        fit = ols_fit(y, others)
        unexplained = fit.rss / fit.tss
    if unexplained <= 1e-12:
        raise NumericalError(
            _i18n.message("diag.vif.collinear", name=name), key="diag.vif.collinear"
        )
    return 1.0 / unexplained
