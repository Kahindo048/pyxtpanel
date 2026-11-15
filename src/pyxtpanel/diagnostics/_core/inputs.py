"""Normalised inputs of the tests (L.1.5, L.1.7).

Category A tests read a declared panel and variable names: :func:`panel_variables` returns the
sorted ``float64`` arrays of the contract, with the rows where a variable is missing removed, and
refuses, with a message that names the first offending unit, what the test cannot handle (an
unbalanced panel, a gap inside a unit, a unit too short).  Category B tests read bare arrays:
:func:`as_array` checks their shape and that they are finite.

Stata refuses the same cases for the same tests (``xtunitroot llc`` wants a strongly balanced
panel, ``xtserial`` wants no gap); each test states what it needs and the refusals are the
same everywhere.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import replace
from typing import TYPE_CHECKING, Literal

import numpy as np

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.panel.arrays import PanelArrays, build_arrays

if TYPE_CHECKING:
    from pyxtpanel._types import FloatArray
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "diag",
    {
        "diag.input.no_data": {
            "en": "{test}: no observation is left once missing values are removed.",
            "fr": "{test} : aucune observation ne reste une fois les valeurs manquantes retirées.",
        },
        "diag.input.not_balanced": {
            "en": "{test} needs a strongly balanced panel: unit {unit!r} is observed {have} "
            "times out of {need} periods. Balance the panel first (PanelData.balance(), or "
            "select complete units).",
            "fr": "{test} exige un panel cylindré : l'unité {unit!r} est observée {have} fois "
            "sur {need} périodes. Cylindrez d'abord le panel (PanelData.balance() ou sélection "
            "des unités complètes).",
        },
        "diag.input.gap": {
            "en": "{test} needs consecutive observations: unit {unit!r} has a gap after "
            "period {period}. Use PanelData.tsfill() to see the missing rows and gap_types() "
            "to classify the gaps.",
            "fr": "{test} exige des observations consécutives : l'unité {unit!r} a un trou après "
            "la période {period}. Utilisez PanelData.tsfill() pour voir les lignes manquantes et "
            "gap_types() pour classer les trous.",
        },
        "diag.input.too_short": {
            "en": "{test} needs at least {need} observations per unit: unit {unit!r} has {have}.",
            "fr": "{test} exige au moins {need} observations par unité : l'unité {unit!r} "
            "en a {have}.",
        },
        "diag.input.too_few_units": {
            "en": "{test} needs at least {need} units, got {have}.",
            "fr": "{test} exige au moins {need} unités, reçu {have}.",
        },
        "diag.input.array_shape": {
            "en": "{name} must be {expected}, got an array of shape {got}.",
            "fr": "{name} doit être {expected}, reçu un tableau de forme {got}.",
        },
        "diag.input.not_finite": {
            "en": "{name} contains {count} missing or infinite values; the test needs finite "
            "numbers.",
            "fr": "{name} contient {count} valeurs manquantes ou infinies ; le test exige des "
            "nombres finis.",
        },
    },
)


def panel_variables(
    panel: PanelBase,
    variables: Sequence[str],
    *,
    test: str,
    balanced: Literal["required", "both"] = "both",
    consecutive: bool = False,
    min_obs: int = 1,
    min_units: int = 1,
) -> PanelArrays:
    """Return the arrays of ``variables`` for a category A test, after the checks it needs.

    Rows where any of ``variables`` is missing are removed, and units left without a row
    disappear.  Then:

    * ``balanced="required"``: every unit must be observed in every period that remains;
    * ``consecutive=True``: the periods of a unit must follow each other without a gap;
    * every unit must keep at least ``min_obs`` rows, and at least ``min_units`` units remain.

    Parameters
    ----------
    panel
        The declared panel.
    variables
        Names of the numeric variables the test reads (missing in any one drops the row).
    test
        Name shown in the messages.
    balanced
        What the test needs (``"required"`` refuses an unbalanced panel).
    consecutive
        Refuse a unit with a gap between two of its observations.
    min_obs, min_units
        Smallest numbers of observations per unit and of units.

    Raises
    ------
    DataError
        When a requirement is not met; the message names the first offending unit.
    """
    arrays = _drop_missing(build_arrays(panel, list(variables)))
    if arrays.n_obs == 0:
        raise DataError(_i18n.message("diag.input.no_data", test=test), key="diag.input.no_data")
    counts = arrays.counts
    labels = arrays.unit_labels
    if balanced == "required":
        need = int(np.unique(arrays.time_index).size)
        short = np.flatnonzero(counts != need)
        if short.size:
            first = int(short[0])
            raise DataError(
                _i18n.message(
                    "diag.input.not_balanced",
                    test=test,
                    unit=plain_label(labels[first]),
                    have=int(counts[first]),
                    need=need,
                ),
                stata_code=459,
                key="diag.input.not_balanced",
            )
    if consecutive:
        step = np.diff(arrays.time_index)
        inside = np.ones(step.size, dtype=bool)
        inside[arrays.offsets[1:-1] - 1] = False  # the step between two units is not a step
        bad = np.flatnonzero(inside & (step != 1))
        if bad.size:
            row = int(bad[0])
            raise DataError(
                _i18n.message(
                    "diag.input.gap",
                    test=test,
                    unit=plain_label(labels[int(arrays.unit_of_row[row])]),
                    period=panel.axis.labels(np.array([arrays.time_index[row]]))[0],
                ),
                stata_code=451,
                key="diag.input.gap",
            )
    small = np.flatnonzero(counts < min_obs)
    if small.size:
        first = int(small[0])
        raise DataError(
            _i18n.message(
                "diag.input.too_short",
                test=test,
                unit=plain_label(labels[first]),
                have=int(counts[first]),
                need=min_obs,
            ),
            key="diag.input.too_short",
        )
    if arrays.n_units < min_units:
        raise DataError(
            _i18n.message(
                "diag.input.too_few_units", test=test, have=arrays.n_units, need=min_units
            ),
            key="diag.input.too_few_units",
        )
    return arrays


def plain_label(label: object) -> object:
    """Return a numpy scalar as a Python value: messages show ``2``, not ``np.int64(2)``."""
    return label.item() if isinstance(label, np.generic) else label


def _drop_missing(arrays: PanelArrays) -> PanelArrays:
    """Keep the rows without a missing value; units left without a row disappear."""
    if arrays.valid.all():
        return arrays
    keep = arrays.valid
    unit = arrays.unit_of_row[keep]
    counts = np.bincount(unit, minlength=arrays.n_units)
    present = counts > 0
    code = np.cumsum(present) - 1
    return replace(
        arrays,
        values=np.ascontiguousarray(arrays.values[keep]),
        offsets=np.concatenate([[0], np.cumsum(counts[present])]).astype(np.int64),
        unit_labels=arrays.unit_labels[present],
        unit_of_row=code[unit].astype(np.int64),
        time_index=arrays.time_index[keep],
        source_row=arrays.source_row[keep],
        valid=np.ones(int(keep.sum()), dtype=bool),
    )


def as_array(
    name: str,
    value: object,
    *,
    ndim: int,
    shape: tuple[int | None, ...] | None = None,
    finite: bool = True,
) -> FloatArray:
    """Return ``value`` as a checked ``float64`` array of ``ndim`` dimensions (category B input).

    Parameters
    ----------
    name
        Name of the argument, shown in the messages.
    value
        Anything ``numpy.asarray`` reads as numbers (list, array, Series, DataFrame).
    ndim
        Required number of dimensions.
    shape
        Required shape, with ``None`` for a free dimension (``(None, 3)`` = any number of rows,
        3 columns).
    finite
        Refuse ``NaN`` and infinities.

    Raises
    ------
    SpecificationError
        If the dimensions or the shape are wrong.
    DataError
        If ``finite`` and a value is missing or infinite.
    """
    try:
        array = np.asarray(value, dtype=np.float64)
    except (TypeError, ValueError):
        raise SpecificationError(
            _i18n.message(
                "diag.input.array_shape", name=name, expected="numeric", got=type(value).__name__
            ),
            key="diag.input.array_shape",
        ) from None
    wanted = shape if shape is not None else (None,) * ndim
    if array.ndim != ndim or any(
        w is not None and w != s for w, s in zip(wanted, array.shape, strict=False)
    ):
        expected = " x ".join("n" if w is None else str(w) for w in wanted)
        raise SpecificationError(
            _i18n.message(
                "diag.input.array_shape", name=name, expected=expected, got=str(array.shape)
            ),
            key="diag.input.array_shape",
        )
    if finite and not np.isfinite(array).all():
        raise DataError(
            _i18n.message(
                "diag.input.not_finite", name=name, count=int((~np.isfinite(array)).sum())
            ),
            key="diag.input.not_finite",
        )
    return np.ascontiguousarray(array)
