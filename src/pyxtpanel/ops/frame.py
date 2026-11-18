"""Operators applied to the variables of a panel, results aligned on the user's rows (layer 3)."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from typing import Literal

import numpy as np
import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError
from pyxtpanel._types import FloatArray
from pyxtpanel.ops._spec import OPERATORS
from pyxtpanel.ops.diff import difference, seasonal
from pyxtpanel.ops.groups import group_mean
from pyxtpanel.ops.lags import lag
from pyxtpanel.ops.parse import Term, parse_terms
from pyxtpanel.ops.within import twoway_demean, within
from pyxtpanel.panel.arrays import sorted_columns
from pyxtpanel.panel.base import PanelBase
from pyxtpanel.panel.index import PanelIndex

_i18n.register(
    "ops",
    {
        "ops.negative_amount": {
            "en": "{kind} needs a non-negative number of periods, got {k}.",
            "fr": "{kind} exige un nombre de périodes positif ou nul, reçu {k}.",
        },
    },
)

OperatorKind = Literal["lag", "lead", "diff", "seasonal"]


def evaluate(index: PanelIndex, values: FloatArray, term: Term) -> FloatArray:
    """Values of the operators of ``term`` applied to ``values`` (sorted order).

    The operators do not commute when periods are missing.  Stata applies them from the variable
    outwards (seasonal differences, then differences, then the shift), which is the canonical
    spelling read from the right.
    """
    out = values
    for period in term.seasonal:
        out = seasonal(index, out, period)
    if term.diff:
        out = difference(index, out, term.diff)
    if term.shift:
        out = lag(index, out, term.shift)
    return out


def operate(panel: PanelBase, terms: Iterable[str | Term]) -> pd.DataFrame:
    """Evaluate Stata-style terms (``"L2.x"``, ``"D.y"``, ``"L(1/3).z"``) on a panel.

    Parameters
    ----------
    panel
        The declared panel.
    terms
        Terms such as ``"L2.x"``, ``"D.y"`` or ``"L(1/3).z"``, as strings or parsed terms.

    Returns
    -------
    pandas.DataFrame
        One column per distinct term, named in canonical form (``"L2.x"``), on the index and in
        the row order of the declared frame.  Rows outside the panel are missing.

    Raises
    ------
    SpecificationError
        If a term is malformed.
    DataError
        If a variable is unknown or not numeric.
    """
    parsed: list[Term] = []
    for term in parse_terms(terms):
        if term.name not in {t.name for t in parsed}:
            parsed.append(term)
    bases = list(dict.fromkeys(t.base for t in parsed))
    values = sorted_columns(panel, bases)
    column = {name: values[:, j] for j, name in enumerate(bases)}
    out = np.empty((panel.n_obs, len(parsed)), dtype=np.float64)
    for j, term in enumerate(parsed):
        out[:, j] = evaluate(panel.index, column[term.base], term)
    frame = panel.align(out)
    assert isinstance(frame, pd.DataFrame)
    frame.columns = [t.name for t in parsed]
    return frame


def _term(kind: OperatorKind, name: str, k: int) -> Term:
    if kind in ("diff", "seasonal") and k < 0:
        raise SpecificationError(
            _i18n.message("ops.negative_amount", kind=kind, k=k), key="ops.negative_amount"
        )
    if kind == "lag":
        return Term(name, shift=k)
    if kind == "lead":
        return Term(name, shift=-k)
    if kind == "diff":
        return Term(name, diff=k)
    return Term(name, seasonal=(k,) if k else ())


def apply_operator(
    panel: PanelBase,
    kind: OperatorKind,
    variables: str | Sequence[str],
    amount: int | Sequence[int] = 1,
) -> pd.Series | pd.DataFrame:
    """``lag``, ``lead``, ``diff`` or ``seasonal`` of one or several variables.

    ``amount`` is the number of periods (``lag``, ``lead``, ``seasonal``) or the order
    (``diff``); a list gives one column per value.  The result is a Series when ``variables``
    and ``amount`` are both single, otherwise a DataFrame whose columns run over the variables
    and, within each, over the amounts.
    """
    names = [variables] if isinstance(variables, str) else list(variables)
    amounts = [amount] if isinstance(amount, (int, np.integer)) else [int(k) for k in amount]
    terms = [_term(kind, name, int(k)) for name in names for k in amounts]
    frame = operate(panel, terms)
    if isinstance(variables, str) and isinstance(amount, (int, np.integer)):
        series: pd.Series = frame.iloc[:, 0]
        return series
    return frame


def _selected(
    panel: PanelBase, variables: str | Sequence[str], listwise: bool
) -> tuple[list[str], FloatArray]:
    names = [variables] if isinstance(variables, str) else list(variables)
    values = sorted_columns(panel, names)
    if listwise:  # Stata's xtdata: a row with any missing variable is dropped for all
        values[np.isnan(values).any(axis=1)] = np.nan
    return names, values


def _present(
    panel: PanelBase, variables: str | Sequence[str], names: list[str], out: FloatArray
) -> pd.Series | pd.DataFrame:
    frame = panel.align(out.reshape(panel.n_obs, -1))
    assert isinstance(frame, pd.DataFrame)
    frame.columns = names
    if isinstance(variables, str):
        series: pd.Series = frame.iloc[:, 0]
        return series
    return frame


def transform(
    panel: PanelBase,
    name: str,
    variables: str | Sequence[str],
    *,
    listwise: bool = False,
    **params: object,
) -> pd.Series | pd.DataFrame:
    """Apply the registered transformation ``name`` to variables; results keep the row alignment.

    Parameters
    ----------
    panel
        The declared panel.
    name
        A name of :data:`pyxtpanel.ops._spec.OPERATORS` (``"within"``, ``"between"``, ``"fod"``,
        ``"rolling"``...).
    variables
        One name (a Series is returned) or several (a DataFrame with the same column names).
    listwise
        Drop, for all the variables, the rows where any of them is missing, before computing
        (what Stata's ``xtdata`` does).  Default: each variable uses its own observed rows.
    **params
        Parameters of the transformation.

    Returns
    -------
    pandas.Series or pandas.DataFrame
        A Series for one variable, a DataFrame for several, on the index of the declared frame.
    """
    func = OPERATORS.get(name)
    names, values = _selected(panel, variables, listwise)
    out = func(panel.index, values, **params)
    assert isinstance(out, np.ndarray)
    return _present(panel, variables, names, out)


def eliminate_effects(
    panel: PanelBase,
    variables: str | Sequence[str],
    effects: Literal["unit", "time", "twoway"] = "unit",
    *,
    add_mean: bool = False,
    listwise: bool = False,
    tol: float = 1e-12,
    maxiter: int = 10_000,
) -> pd.Series | pd.DataFrame:
    """Within transformation removing unit effects, period effects or both.

    ``"twoway"`` puts the convergence report in ``result.attrs["convergence"]``.
    """
    names, values = _selected(panel, variables, listwise)
    report = None
    if effects == "unit":
        out = within(panel.index, values, add_mean=add_mean)
    elif effects == "time":
        data = np.asarray(values, dtype=np.float64)
        out = data - group_mean(panel.index, data, "time")
        if add_mean:
            out = out + np.nanmean(data, axis=0, keepdims=True)
    else:
        out, report = twoway_demean(panel.index, values, tol=tol, maxiter=maxiter)
        if add_mean:
            out = out + np.nanmean(values, axis=0, keepdims=True)
    result = _present(panel, variables, names, out)
    if report is not None:
        result.attrs["convergence"] = report.as_dict()
    return result


def term_columns(
    panel: PanelBase, items: str | Sequence[str | Term]
) -> tuple[list[str], FloatArray]:
    """Columns for variable names or Stata-style terms, as a sorted ``(n_obs, k)`` array.

    A name that is a column of the frame is used as it is (so names such as ``"gdp.pc"`` work);
    anything else is read as a term (``"L2.x"``, ``"D.x"``).  Returns the names (canonical for
    terms) and the values in the sorted order of the panel.
    """
    listed = [items] if isinstance(items, (str, Term)) else list(items)
    names: list[str] = []
    columns: list[FloatArray] = []
    for item in listed:
        if isinstance(item, str) and item in panel.frame.columns:
            names.append(item)
            columns.append(sorted_columns(panel, [item])[:, 0])
            continue
        for term in parse_terms([item]):
            names.append(term.name)
            columns.append(evaluate(panel.index, sorted_columns(panel, [term.base])[:, 0], term))
    return names, np.column_stack(columns) if columns else np.empty((panel.n_obs, 0))
