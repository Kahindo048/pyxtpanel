"""Exceptions and warnings of the package (layer 0).

The errors that Stata reports with a return code keep it in :attr:`PyxtpanelError.stata_code`
(451 repeated time values, 459 missing key, 198 invalid option, 2000 no observations...), so
that users coming from Stata recognise them.  :attr:`PyxtpanelError.key` is the key of the
message in the catalogue and identifies the rule that failed.

The estimation packages derive their own root error from :class:`PyxtpanelError` and reuse
:class:`DataError`, :class:`SpecificationError` and :class:`NumericalError`, so that a single
``except DataError`` catches an error of the foundation as well as one of an estimator.
"""

from __future__ import annotations


class PyxtpanelError(Exception):
    """Base class of every error raised by pyxtpanel."""

    #: Return code of the equivalent Stata error, when there is one.
    stata_code: int | None = None
    #: Key of the message in the catalogue (identifies the rule that failed).
    key: str | None = None

    def __init__(
        self, message: str, *, stata_code: int | None = None, key: str | None = None
    ) -> None:
        super().__init__(message)
        if stata_code is not None:
            self.stata_code = stata_code
        self.key = key


class DataError(PyxtpanelError, ValueError):
    """The data do not satisfy a requirement (balance, duplicates, missing values...)."""


class DeclarationError(DataError):
    """The panel cannot be declared (``xtset``): repeated pairs, missing keys, bad time."""


class TimeAxisError(DataError):
    """The time variable does not define a regular axis (frequency, step, irregular dates)."""


class SpecificationError(PyxtpanelError, ValueError):
    """An option or an argument is invalid or inconsistent."""

    stata_code = 198


class NumericalError(PyxtpanelError, ArithmeticError):
    """A numerical step failed (singular matrix, no convergence, empty sample)."""


class PyxtpanelWarning(UserWarning):
    """Non-fatal condition worth reporting (infinite value read as missing, fallback used)."""


#: Stata return codes the package keeps, with the rule they stand for.  A test checks that no
#: error carries a code that is not listed here.
STATA_CODES: dict[int, str] = {
    110: "a variable that must be new already exists",
    134: "too many distinct values to tabulate",
    198: "invalid syntax, option or argument",
    451: "repeated time values within a panel",
    459: "something that should not be missing is missing (or too few distinct values)",
    2000: "no observations",
}
