"""Reading and writing panel files, including the ``xtset`` Stata records in a ``.dta``.

When Stata saves a dataset after ``xtset`` it stores the panel variable, the time variable and
the time step as *characteristics* of the dataset (``_dta[_TSpanel]``, ``_dta[_TStvar]``,
``_dta[_TSdelta]``).  :func:`stata_xtset` reads them without Stata, and :func:`write_stata`
writes them back, which pandas cannot do, so that Stata finds the panel declared.
"""

from __future__ import annotations

import contextlib
import re
import struct
from collections.abc import Hashable
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any

import pandas as pd

from pyxtpanel import _i18n
from pyxtpanel._errors import DataError, SpecificationError
from pyxtpanel.time.formats import freq_from_format

if TYPE_CHECKING:
    from pyxtpanel.panel.base import PanelBase

_i18n.register(
    "io",
    {
        "io.unsupported_file": {
            "en": "Unsupported file type '{suffix}' (use .dta, .csv or .parquet).",
            "fr": "Type de fichier non géré « {suffix} » (utilisez .dta, .csv ou .parquet).",
        },
        "io.no_xtset": {
            "en": "The file '{path}' carries no xtset declaration: pass id= and time=.",
            "fr": "Le fichier « {path} » ne porte pas de déclaration xtset : passez id= et time=.",
        },
        "io.delta_not_integer": {
            "en": "The time step recorded in the file ({delta}) is not a whole number of periods.",
            "fr": "Le pas de temps enregistré dans le fichier ({delta}) n'est pas un entier de "
            "périodes.",
        },
        "io.id_not_numeric": {
            "en": "Stata needs a numeric panel variable, but '{name}' is {dtype}. Encode the units "
            "as integers before writing.",
            "fr": "Stata exige une variable de panel numérique, mais « {name} » est {dtype}. "
            "Codez les unités en entiers avant d'écrire.",
        },
        "io.week_unsupported": {
            "en": "Weekly data cannot be written with a Stata format (Stata weeks are not "
            "calendar weeks).",
            "fr": "Les données hebdomadaires ne peuvent pas être écrites avec un format Stata "
            "(les semaines de Stata ne sont pas des semaines calendaires).",
        },
    },
)

#: Width of the names stored in the characteristics section, by ``.dta`` release.
_NAME_WIDTH = {117: 33, 118: 129, 119: 129}
#: Index, in the ``<map>`` of a release 117-119 file, of the characteristics and the data sections.
_MAP_CHARACTERISTICS, _MAP_DATA = 8, 9
#: Entries of the map that follow the characteristics section and move when it grows.
_MAP_AFTER = range(_MAP_DATA, 14)


@dataclass(frozen=True)
class XtsetInfo:
    """What Stata recorded in a ``.dta`` file when ``xtset`` or ``tsset`` was run.

    Attributes
    ----------
    panel, time
        Panel and time variables (``None`` when the file was saved without ``xtset``, or without
        a panel variable after ``tsset``).
    delta
        The time step in units of the time format (``1.0`` for consecutive periods).
    """

    panel: str | None = None
    time: str | None = None
    delta: float | None = None


def _stata_delta(text: str) -> float | None:
    """Read ``_dta[_TSdelta]``, which Stata stores in its hexadecimal format (``%21x``)."""
    match = re.fullmatch(r"([+-]?)(\d)\.([0-9a-fA-F]*)X([+-])([0-9a-fA-F]+)", text.strip())
    if match is None:
        return None
    sign, lead, fraction, exp_sign, exponent = match.groups()
    mantissa = float.fromhex(f"0x{lead}.{fraction or '0'}")
    value = mantissa * 2.0 ** (int(exponent, 16) * (-1 if exp_sign == "-" else 1))
    return -value if sign == "-" else value


def _hex_delta(delta: int) -> str:
    """Write a time step in Stata's ``%21x`` format (``+1.0000000000000X+000`` for 1)."""
    mantissa, exponent = float(delta).hex()[2:].split("p")
    lead, _, fraction = mantissa.partition(".")
    sign = "-" if exponent.startswith("-") else "+"
    return f"+{lead}.{fraction.ljust(13, '0')}X{sign}{abs(int(exponent)):03x}"


def _characteristics(raw: bytes, release: int, order: str) -> dict[tuple[str, str], str]:
    """Return the characteristics ``(variable, name) -> text`` of the section ``raw``."""
    width = _NAME_WIDTH[release]
    found: dict[tuple[str, str], str] = {}
    position = raw.find(b"<ch>")
    while position >= 0:
        (size,) = struct.unpack(f"{order}I", raw[position + 4 : position + 8])
        body = raw[position + 8 : position + 8 + size]
        names = (body[:width], body[width : 2 * width], body[2 * width :])
        variable, name, text = (part.split(b"\0")[0].decode("latin-1") for part in names)
        found[(variable, name)] = text
        position = raw.find(b"<ch>", position + 8 + size)
    return found


def _read_header(head: bytes) -> tuple[int, str, int] | None:
    release_tag = re.search(rb"<release>(\d+)</release>", head)
    order_tag = re.search(rb"<byteorder>(MSF|LSF)</byteorder>", head)
    map_at = head.find(b"<map>")
    if release_tag is None or order_tag is None or map_at < 0:
        return None
    release = int(release_tag.group(1))
    if release not in _NAME_WIDTH:
        return None
    return release, "<" if order_tag.group(1) == b"LSF" else ">", map_at + 5


def stata_xtset(path: str | Path) -> XtsetInfo:
    """Read the ``xtset`` panel and time variables recorded in a ``.dta`` file.

    Releases 117 to 119 (Stata 13 and later) are read; any other file, or a file that cannot be
    parsed, gives an empty result.
    """
    try:
        with Path(path).open("rb") as handle:
            head = handle.read(4096)
            header = _read_header(head)
            if header is None:
                return XtsetInfo()
            release, order, map_start = header
            offsets = struct.unpack(f"{order}14Q", head[map_start : map_start + 14 * 8])
            start, stop = offsets[_MAP_CHARACTERISTICS], offsets[_MAP_DATA]
            handle.seek(start)
            section = handle.read(stop - start)
        facts = _characteristics(section, release, order)
    except (OSError, struct.error, ValueError, UnicodeDecodeError):
        return XtsetInfo()
    delta = facts.get(("_dta", "_TSdelta"))
    return XtsetInfo(
        panel=facts.get(("_dta", "_TSpanel")) or None,
        time=facts.get(("_dta", "_TStvar")) or None,
        delta=_stata_delta(delta) if delta else None,
    )


def read_table(path: str | Path) -> pd.DataFrame:
    """Read a ``.dta``, ``.csv`` or ``.parquet`` file into a DataFrame.

    Stata value labels are not applied, so labelled variables stay numeric as they are in Stata.
    """
    file = Path(path)
    suffix = file.suffix.lower()
    if suffix == ".dta":
        return read_stata(file)
    if suffix == ".csv":
        return pd.read_csv(file)
    if suffix == ".parquet":
        return pd.read_parquet(file)
    raise DataError(
        _i18n.message("io.unsupported_file", suffix=suffix or "(none)"), key="io.unsupported_file"
    )


def read_stata(path: str | Path) -> pd.DataFrame:
    """Read a Stata ``.dta`` file and keep what it says about its panel declaration.

    Date variables are converted to ``datetime64``.  ``frame.attrs`` holds ``id``, ``time``,
    ``delta`` (an integer, or ``None``) and ``time_format`` (the ``%t`` format of the time
    variable, or ``None``); all are ``None`` when the file was saved without ``xtset``.  Names
    are **never guessed** from the columns.
    """
    file = Path(path)
    with pd.read_stata(file, convert_categoricals=False, iterator=True) as reader:
        frame = reader.read()
        formats = list(getattr(reader, "_fmtlist", []))
    info = stata_xtset(file)
    usable = info.panel in frame.columns and info.time in frame.columns
    delta = None
    if usable and info.delta is not None and float(info.delta).is_integer():
        delta = int(info.delta)
    time_format = None
    if usable and len(formats) == len(frame.columns):
        time_format = dict(zip(frame.columns, formats, strict=True)).get(info.time)
    frame.attrs.update(
        id=info.panel if usable else None,
        time=info.time if usable else None,
        delta=delta,
        delta_recorded=info.delta if usable else None,
        time_format=time_format,
    )
    return frame


def panel_from_stata(cls: type[PanelBase], path: str | Path, **options: Any) -> PanelBase:
    """Declare the panel recorded in a ``.dta`` file (see :meth:`PanelBase.from_stata`)."""
    frame = read_stata(path)
    attrs = frame.attrs
    id_name = options.pop("id", None) or attrs.get("id")
    time_name = options.pop("time", None) or attrs.get("time")
    if id_name is None or time_name is None:
        raise SpecificationError(
            _i18n.message("io.no_xtset", path=str(path)), stata_code=459, key="io.no_xtset"
        )
    from_file = time_name == attrs.get("time")
    recorded = attrs.get("delta_recorded")
    if (
        from_file
        and options.get("delta") is None
        and recorded is not None
        and attrs["delta"] is None
    ):
        raise DataError(
            _i18n.message("io.delta_not_integer", delta=recorded), key="io.delta_not_integer"
        )
    options.setdefault("delta", attrs.get("delta") if from_file else None)
    fmt = attrs.get("time_format")
    if fmt and from_file and pd.api.types.is_datetime64_any_dtype(frame[time_name]):
        with contextlib.suppress(SpecificationError):  # %tw, %tc: infer the frequency instead
            options.setdefault("freq", freq_from_format(fmt))
    return cls.from_frame(frame, id_name, time_name, **options)


# ---------------------------------------------------------------- writing


def _characteristic_block(variable: str, name: str, text: str, release: int, order: str) -> bytes:
    width = _NAME_WIDTH[release]
    payload = (
        variable.encode("latin-1").ljust(width, b"\0")
        + name.encode("latin-1").ljust(width, b"\0")
        + text.encode("latin-1")
        + b"\0"
    )
    return b"<ch>" + struct.pack(f"{order}I", len(payload)) + payload + b"</ch>"


def inject_xtset(path: str | Path, panel: str, time: str, delta: int) -> None:
    """Add the ``xtset`` characteristics to a ``.dta`` file written by pandas.

    The file is rewritten in place: the characteristics section is replaced and the map offsets of
    the sections that follow are shifted by the growth.

    Raises
    ------
    DataError
        If the file is not a release 117-119 ``.dta``.
    """
    file = Path(path)
    raw = bytearray(file.read_bytes())
    header = _read_header(bytes(raw[:4096]))
    if header is None:
        raise DataError(f"'{file}' is not a Stata 13-15 .dta file", key="io.not_dta")
    release, order, map_start = header
    offsets = list(struct.unpack(f"{order}14Q", bytes(raw[map_start : map_start + 14 * 8])))
    start, stop = offsets[_MAP_CHARACTERISTICS], offsets[_MAP_DATA]
    open_tag, close_tag = b"<characteristics>", b"</characteristics>"
    assert raw[start : start + len(open_tag)] == open_tag
    blocks = b"".join(
        _characteristic_block("_dta", name, text, release, order)
        for name, text in (
            ("_TStvar", time),
            ("_TSpanel", panel),
            ("_TSdelta", _hex_delta(delta)),
            ("_TSitrvl", "1"),
            ("tis", time),
            ("iis", panel),
        )
    )
    section = open_tag + blocks + close_tag
    growth = len(section) - (stop - start)
    raw[start:stop] = section
    for entry in _MAP_AFTER:
        offsets[entry] += growth
    raw[map_start : map_start + 14 * 8] = struct.pack(f"{order}14Q", *offsets)
    file.write_bytes(bytes(raw))


_STATA_DATE_CODE = {"D": "td", "M": "tm", "Q": "tq", "H": "th", "Y": "ty"}


def write_stata(panel: PanelBase, path: str | Path) -> None:
    """Write the panel as a Stata 14 ``.dta`` that keeps its ``xtset`` declaration.

    The data are written sorted by unit then period.  Dates are written with the matching ``%t``
    format (a date that is not a period start is written as the start of its period).  Weekly
    data are refused.

    Raises
    ------
    SpecificationError
        If the identifier is not numeric (Stata's ``xtset`` needs a numeric panel variable) or
        the data are weekly.
    """
    axis = panel.axis
    data = panel.to_frame(order="sorted")
    keys = [n for n in (panel.id_name, panel.time_name) if n in data.index.names]
    if keys:
        data = data.reset_index(level=keys)
    if not pd.api.types.is_numeric_dtype(data[panel.id_name]):
        raise SpecificationError(
            _i18n.message("io.id_not_numeric", name=panel.id_name, dtype=data[panel.id_name].dtype),
            key="io.id_not_numeric",
        )
    dates: dict[Hashable, str] = {}
    if axis.kind in ("datetime", "period", "stata"):
        if axis.freq == "W":
            raise SpecificationError(
                _i18n.message("io.week_unsupported"), key="io.week_unsupported"
            )
        assert axis.freq is not None
        data[panel.time_name] = axis.datetimes(panel.index.rank)
        dates[panel.time_name] = _STATA_DATE_CODE[axis.freq]
    file = Path(path)
    data.to_stata(file, write_index=False, convert_dates=dates or None, version=118)  # type: ignore[arg-type]
    inject_xtset(file, panel.id_name, panel.time_name, axis.delta)


def _stata_dates(axis: Any, panel: PanelBase) -> pd.DatetimeIndex:
    """Dates (period starts) of a panel whose time was Stata-coded integers."""
    from pyxtpanel.time.axis import _start_dates  # same layer; private helper

    return pd.DatetimeIndex(_start_dates(axis.ordinal(panel.index.rank), axis.freq))
