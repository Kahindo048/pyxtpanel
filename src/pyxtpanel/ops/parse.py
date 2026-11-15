"""Stata time-series notation: ``"L2.x"``, ``"LD.x"``, ``"S12.x"``, ``"L(1/3).x"`` (layer 3).

Pure text, no data.  A term names a variable and a chain of operators: ``L`` lag, ``F`` lead,
``D`` difference, ``S`` seasonal difference, each with an optional number (``L2`` second lag,
``D2`` second difference, ``S12`` twelve-period difference) or a list in parentheses
(``"L(1/2).x"`` is ``"L.x"`` and ``"L2.x"``).  A chain reduces to a net shift, a difference order
and the seasonal periods; the canonical spelling is the one Stata writes (``"ll.x"`` and
``"l2.x"`` are ``"L2.x"``, ``"dl.x"`` is ``"LD.x"``, ``"lf.x"`` is ``"x"``), which makes ``"x"``
and ``"L0.x"`` the same term.

Origin: ``pyxtdpthresh.data.operators`` (52 spellings validated against Stata).
"""

from __future__ import annotations

import itertools
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field

from pyxtpanel import _i18n
from pyxtpanel._errors import SpecificationError

_i18n.register(
    "ops",
    {
        "ops.term_invalid": {
            "en": "invalid term {term!r}: expected a variable name, optionally preceded by "
            "operators such as L., L2., F., D., S12., LD. or L(1/3).",
            "fr": "terme invalide {term!r} : un nom de variable, éventuellement précédé "
            "d'opérateurs tels que L., L2., F., D., S12., LD. ou L(1/3).",
        },
    },
)

# A letter (Unicode: Stata 14+ accepts "année") or an underscore, then word characters
_NAME = re.compile(r"^[^\W\d]\w*$")
# a number, or a list in parentheses that may hold numlist steps such as 1(2)5
_TOKEN = re.compile(r"([lfdsLFDS])(\((?:[^()]|\([^()]*\))*\)|\d*)")


@dataclass(frozen=True)
class Term:
    """A variable with time-series operators, in canonical form.

    Attributes
    ----------
    base
        Name of the variable.
    shift
        Net lag (``L`` minus ``F``); negative for a lead.
    diff
        Order of differencing (``D``).
    seasonal
        Periods of the seasonal differences (``S``).
    raw
        The text the term was written as (kept for messages; not compared).
    """

    base: str
    shift: int = 0
    diff: int = 0
    seasonal: tuple[int, ...] = ()
    raw: str = field(default="", compare=False)

    @property
    def is_operated(self) -> bool:
        """Whether any operator remains after reduction."""
        return bool(self.shift or self.diff or self.seasonal)

    @property
    def netlag(self) -> int:
        """Net lag of the term (negative for a lead)."""
        return self.shift

    @property
    def has_dsf(self) -> bool:
        """Whether the term has a ``D``, ``S`` or ``F`` operator (a current or future value)."""
        return self.diff > 0 or bool(self.seasonal) or self.shift < 0

    @property
    def name(self) -> str:
        """Canonical spelling, e.g. ``"L2D.x"``."""
        if not self.is_operated:
            return self.base
        ops = ""
        if self.shift:
            ops += ("L" if self.shift > 0 else "F") + (
                str(abs(self.shift)) if abs(self.shift) > 1 else ""
            )
        if self.diff:
            ops += "D" + (str(self.diff) if self.diff > 1 else "")
        ops += "".join("S" + (str(p) if p > 1 else "") for p in self.seasonal)
        return f"{ops}.{self.base}"

    def __str__(self) -> str:
        return self.name


def _invalid(term: str) -> SpecificationError:
    return SpecificationError(
        _i18n.message("ops.term_invalid", term=term), stata_code=198, key="ops.term_invalid"
    )


def _numbers(text: str, term: str) -> list[int]:
    """Numbers of an operator suffix: ``""``, ``"2"``, ``"(1/3)"``, ``"(1 3)"``, ``"(1(2)5)"``."""
    if not text:
        return [1]
    if not text.startswith("("):
        return [int(text)]
    values: list[int] = []
    for part in re.split(r"[\s,]+", text[1:-1].strip()):
        if m := re.fullmatch(r"(\d+)/(\d+)", part):
            values += range(int(m.group(1)), int(m.group(2)) + 1)
        elif m := re.fullmatch(r"(\d+)\((\d+)\)(\d+)", part):
            values += range(int(m.group(1)), int(m.group(3)) + 1, max(int(m.group(2)), 1))
        elif part.isdigit():
            values.append(int(part))
        else:
            raise _invalid(term)
    if not values:
        raise _invalid(term)
    return values


def _reduce(base: str, tokens: list[tuple[str, int]], raw: str) -> Term:
    shift = diff_order = 0
    periods: list[int] = []
    for letter, n in tokens:
        if letter == "l":
            shift += n
        elif letter == "f":
            shift -= n
        elif letter == "d":
            diff_order += n
        elif n > 0:
            periods.append(n)
    return Term(base, shift, diff_order, tuple(periods), raw)


def parse_term(item: str | Term) -> list[Term]:
    """Parse a term; a range such as ``"L(1/2).x"`` gives several terms.

    Raises
    ------
    SpecificationError
        If the text is not a variable name with valid operators (Stata code 198).
    """
    # Stata: syntax varlist(ts), canonical spelling (pyxtdpthresh: _xdpt_tsterm)
    if isinstance(item, Term):
        return [item]
    text = str(item).strip()
    *prefixes, base = text.split(".")
    if not _NAME.match(base):
        raise _invalid(text)
    if not prefixes:
        return [Term(base, raw=text)]
    chain = "".join(prefixes)
    position, groups = 0, []
    while position < len(chain):
        match = _TOKEN.match(chain, position)
        if match is None or match.end() == position:
            raise _invalid(text)
        groups.append((match.group(1).lower(), _numbers(match.group(2), text)))
        position = match.end()
    if not groups:
        raise _invalid(text)
    return [
        _reduce(base, list(zip((g[0] for g in groups), combo, strict=True)), text)
        for combo in itertools.product(*(g[1] for g in groups))
    ]


def parse_terms(items: Iterable[str | Term]) -> list[Term]:
    """Parse several terms, expanding ranges, in order."""
    return [term for item in items for term in parse_term(item)]


def names(terms: Sequence[Term]) -> list[str]:
    """Canonical spelling of each term."""
    return [t.name for t in terms]
