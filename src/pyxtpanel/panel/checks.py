"""Checks run when a panel is declared (what Stata's ``xtset`` verifies).

Stata refuses repeated (unit, time) pairs and classifies the panel as strongly balanced, weakly
balanced or unbalanced in ``r(balanced)``.
"""

from __future__ import annotations

import difflib
from collections.abc import Iterable
from dataclasses import dataclass
from typing import Literal

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, DeclarationError
from pyxtpanel._types import IntArray

Balance = Literal["strongly balanced", "weakly balanced", "unbalanced"]

#: Columns listed, and close names suggested, when a variable is not found.
MAX_LISTED, MAX_SUGGESTED = 10, 3
#: Rows and periods named in a message.
MAX_ROWS_LISTED = 3
MAX_PERIODS_LISTED = 5

_i18n.register(
    "panel",
    {
        "panel.variable_not_found": {
            "en": "Variable '{name}' not found.",
            "fr": "Variable « {name} » introuvable.",
        },
        "panel.variable_suggest": {
            "en": "Did you mean {suggestions}?",
            "fr": "Vouliez-vous dire {suggestions} ?",
        },
        "panel.variable_available": {
            "en": "Available columns: {columns}.",
            "fr": "Colonnes disponibles : {columns}.",
        },
        "panel.keys_not_found": {
            "en": "The panel identifier and the time variable were not given and the index does "
            "not have two levels; pass id= and time=.",
            "fr": "L'identifiant du panel et la variable de temps ne sont pas fournis et l'index "
            "n'a pas deux niveaux ; passez id= et time=.",
        },
        "panel.one_key_given": {
            "en": "Pass both id= and time=, or neither (to use a two-level index).",
            "fr": "Passez id= et time= ensemble, ou aucun des deux (pour utiliser un index à "
            "deux niveaux).",
        },
        "panel.duplicate_column": {
            "en": "Column '{name}' appears more than once; the panel keys must be unique columns.",
            "fr": "La colonne « {name} » apparaît plusieurs fois ; les clés du panel doivent être "
            "des colonnes uniques.",
        },
        "panel.key_missing": {
            "en": "The {what} variable '{name}' has {n} missing value(s) (first rows: {rows}). "
            "Remove them, or pass missing_keys='drop' to ignore these observations as Stata "
            "does.",
            "fr": "La variable {what} « {name} » a {n} valeur(s) manquante(s) (premières lignes : "
            "{rows}). Retirez-les, ou passez missing_keys='drop' pour ignorer ces observations "
            "comme Stata.",
        },
        "panel.missing_keys_option": {
            "en": "missing_keys must be 'error' or 'drop', got {value!r}.",
            "fr": "missing_keys doit valoir 'error' ou 'drop', reçu {value!r}.",
        },
        "panel.duplicates": {
            "en": "Repeated observations: {n} repeated pair(s) of ({id}, {time}), for example "
            "{examples}. A unit can have only one observation per period.",
            "fr": "Observations répétées : {n} couple(s) ({id}, {time}) répété(s), par exemple "
            "{examples}. Une unité ne peut avoir qu'une observation par période.",
        },
        "panel.id_unsortable": {
            "en": "The identifier '{name}' mixes types that cannot be sorted ({error}).",
            "fr": "L'identifiant « {name} » mélange des types qui ne se trient pas ({error}).",
        },
        "panel.no_observations": {
            "en": "No usable observation: the data are empty or every row has a missing key.",
            "fr": "Aucune observation utilisable : les données sont vides ou chaque ligne a une "
            "clé manquante.",
        },
        "panel.balance_example": {
            "en": "First incomplete unit: {id}={unit}, which lacks {time} {periods}.",
            "fr": "Première unité incomplète : {id}={unit}, qui n'a pas {time} {periods}.",
        },
    },
)


def unknown_variable(name: str, available: Iterable[object]) -> DataError:
    """Return the error for a variable that is not in the data, with what could have been meant.

    The message names the variable, suggests the columns whose name is close to it (ignoring the
    case) and lists the available columns, so that a typo is found without opening the data.
    """
    columns = [str(column) for column in available]
    lowered = {column.lower(): column for column in columns}
    close = difflib.get_close_matches(str(name).lower(), list(lowered), n=MAX_SUGGESTED, cutoff=0.6)
    text = _i18n.message("panel.variable_not_found", name=name)
    if close:
        suggestions = ", ".join(f"'{lowered[c]}'" for c in close)
        text += " " + _i18n.message("panel.variable_suggest", suggestions=suggestions)
    if columns:
        listed = ", ".join(columns[:MAX_LISTED]) + (" ..." if len(columns) > MAX_LISTED else "")
        text += " " + _i18n.message("panel.variable_available", columns=listed)
    return DataError(text, key="panel.variable_not_found")


def check_unique_pairs(
    unit: IntArray,
    rank: IntArray,
    id_name: str,
    time_name: str,
    show: tuple[np.ndarray, np.ndarray],
) -> IntArray:
    """Raise if a (unit, period) pair is repeated (Stata ``xtset`` error 451).

    The check sorts the rows by (unit, period); that order is returned so that the caller sorts
    only once.

    Parameters
    ----------
    unit, rank
        Integer unit code and calendar rank of each row.
    id_name, time_name
        Names used in the message.
    show
        Identifier and time value of each row as they are in the data, for the examples.

    Returns
    -------
    numpy.ndarray
        Positions that put the rows in (unit, period) order.

    Raises
    ------
    DeclarationError
        If a pair appears twice; the message counts them and gives the first three.
    """
    span = int(np.max(rank)) + 1 if rank.size else 1
    key = unit.astype(np.int64) * span + rank
    order = np.argsort(key, kind="stable")
    sorted_key = key[order]
    repeated = np.flatnonzero(sorted_key[1:] == sorted_key[:-1])
    if repeated.size == 0:
        return order
    rows = order[repeated[:MAX_ROWS_LISTED] + 1]
    examples = ", ".join(f"({id_name}={show[0][i]}, {time_name}={show[1][i]})" for i in rows)
    raise DeclarationError(
        _i18n.message(
            "panel.duplicates", n=int(repeated.size), id=id_name, time=time_name, examples=examples
        ),
        stata_code=451,
        key="panel.duplicates",
    )


@dataclass(frozen=True)
class BalanceSummary:
    """Balance of a panel as ``xtset`` reports it, with the counts behind the verdict.

    Attributes
    ----------
    kind
        ``r(balanced)`` of Stata: strongly balanced (every unit observed in every period),
        weakly balanced (same number of periods in each unit, but not the same periods) or
        unbalanced.
    n_units
        Number of units.
    n_periods
        Number of distinct periods observed in the whole panel.
    n_incomplete
        Number of units not observed in every period.
    first_incomplete
        Code of the first incomplete unit (``-1`` when there is none).
    missing_ranks
        Calendar ranks, among the observed periods, that the first incomplete unit lacks.
    """

    kind: Balance
    n_units: int
    n_periods: int
    n_incomplete: int
    first_incomplete: int = -1
    missing_ranks: tuple[int, ...] = ()


def balance_summary(offsets: IntArray, unit: IntArray, rank: IntArray) -> BalanceSummary:
    """Classify a panel like ``xtset`` does in ``r(balanced)``.

    ``offsets``, ``unit`` and ``rank`` describe rows sorted by unit then period, with no
    repeated pair.
    """
    n_units = int(offsets.size - 1)
    if rank.size == 0:
        return BalanceSummary("strongly balanced", 0, 0, 0)
    periods = np.unique(rank)
    counts = np.diff(offsets)
    incomplete = np.flatnonzero(counts != periods.size)
    if incomplete.size == 0:
        kind: Balance = "strongly balanced"
    elif np.all(counts == counts[0]):
        kind = "weakly balanced"
    else:
        kind = "unbalanced"
    if incomplete.size == 0:
        return BalanceSummary(kind, n_units, int(periods.size), 0)
    first = int(incomplete[0])
    present = rank[offsets[first] : offsets[first + 1]]
    missing = np.setdiff1d(periods, present)
    return BalanceSummary(
        kind,
        n_units,
        int(periods.size),
        int(incomplete.size),
        first,
        tuple(int(m) for m in missing),
    )
