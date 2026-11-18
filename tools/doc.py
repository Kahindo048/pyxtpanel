"""Generated annexes and coverage of ``docs/documentation.md`` (J.2).

    python -m tools.doc            # check: annexes up to date, every name covered by its part
    python -m tools.doc --write    # rewrite the annexes between their markers
    python -m tools.doc --latex    # docs/documentation.tex from the Markdown (needs pandoc)
    python -m tools.doc --pdf      # and docs/documentation.pdf (needs xelatex)

The documentation is one hand-written document (parts 0 to 9) followed by four annexes that are
generated from the registries, so they cannot drift from the code:

* annex A: every function with its Stata command, its R or Python equivalent and its summary;
* annex B: one card per statistical test (signature, family, category, reference, validation);
* annex C: the Stata return codes kept and the error classes;
* annex D: the behaviours of Stata reproduced on purpose, from the ``# STATA-QUIRK:`` comments.

The check also fails when a registered function, operator, plot, exploration tool or test is not
named (in backticks) in the part of the document that is supposed to teach it, so a feature added
without documentation is caught.  The Python blocks of parts 0 to 8 are run by
``tests/contract/test_documentation.py``.
"""

from __future__ import annotations

import inspect
import json
import re
import shutil
import subprocess
import sys
import tempfile
import textwrap
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "docs" / "documentation.md"
TEX = ROOT / "docs" / "documentation.tex"
PDF = ROOT / "docs" / "documentation.pdf"
HEADER = ROOT / "docs" / "latex" / "header.tex"

MARK = "<!-- annex:{key}:{edge} -->"
ANNEXES = ("a", "b", "c", "d", "e")


def cell(text: object) -> str:
    """A table cell: one line, no pipe."""
    return (
        str(text if text not in (None, "") else "—").replace("|", "\\|").replace("\n", " ").strip()
    )


# ------------------------------------------------------------------ annex A


#: Nearest function elsewhere for the registries that do not carry one (from ETAT_DES_LIEUX.md).
EQUIVALENTS = {
    "xtset(frame, id, time, delta=)": "plm::pdata.frame; linearmodels PanelData",
    "PanelData.tsfill()": "plm::make.pbalanced, make.pconsecutive; panelr balance_panel",
    "PanelData.tsfill(full=True)": "plm::make.pbalanced",
    "duplicates(frame, id, time)": "plm (refuses duplicates)",
    "PanelData.xtdescribe()": "plm::pdim (size only)",
    "PanelData.xtsum()": "collapse::qsu(by = pid) (equivalence to xtsum not confirmed)",
    "PanelData.xtcorr()": "plm::cortab; collapse::pwcor, pwnobs",
    "PanelData.missing_report()": "collapse::fnobs",
    "transform:lag": "plm::lag; collapse::flag",
    "transform:lead": "plm::lead; collapse::flag",
    "transform:difference": "plm::diff; collapse::fdiff",
    "transform:seasonal": "plm::diff (lag argument)",
    "transform:growth": "collapse::fgrowth",
    "transform:within": "plm::Within; collapse::fwithin; linearmodels demean",
    "transform:between": "plm::Between; collapse::fbetween; linearmodels mean",
    "transform:group_count": "collapse::fnobs; linearmodels count",
    "transform:group_mean": "collapse::fbetween",
    "transform:twoway": "linearmodels general_demean",
}


def annex_a() -> str:
    """Every function, one table per area."""
    from pyxtpanel.describe._spec import DESCRIPTIONS
    from pyxtpanel.diagnostics._core import TESTS
    from pyxtpanel.explore._spec import TOOLS
    from pyxtpanel.ops._spec import OPERATORS
    from pyxtpanel.plot._spec import KINDS
    from tools.stata_map import STRUCTURAL

    head = ["| Call | Stata | R or Python equivalent | What it does |", "|---|---|---|---|"]
    out: list[str] = []

    def table(
        title: str, rows: list[tuple[str, object, object, object]], *, raw: bool = False
    ) -> None:
        out.extend([f"### {title}", "", *head])
        for call, stata, other, what in rows:
            first = call if raw else f"`{call}`"
            other = (
                other
                or EQUIVALENTS.get(call)
                or EQUIVALENTS.get(call.replace('PanelData.transform("', "transform:").rstrip('")'))
            )
            out.append(f"| {first} | {cell(stata)} | {cell(other)} | {cell(what)} |")
        out.append("")

    table(
        "Declaring, selecting, reading",
        [(python, stata, EQUIVALENTS.get(python.strip("`")), None) for stata, python in STRUCTURAL],
        raw=True,
    )
    table(
        "Describing",
        [
            (f"PanelData.{n}()", s.stata, None, s.summary)
            for n, s in sorted(DESCRIPTIONS.specs(), key=lambda p: p[0])
        ],
    )
    table(
        "Operators and transformations",
        [
            (f'PanelData.transform("{n}")', s.stata, None, s.summary)
            for n, s in sorted(OPERATORS.specs(), key=lambda p: p[0])
        ],
    )
    table(
        "Exploration tools",
        [
            (f"PanelData.{n}()" if s.on_panel else f"xt.{n}()", s.stata, s.similar, s.summary)
            for n, s in sorted(TOOLS.specs(), key=lambda p: p[0])
        ],
    )
    table(
        "Figures",
        [
            (f'PanelData.plot("{n}")', s.stata, None, s.summary)
            for n, s in sorted(KINDS.specs(), key=lambda p: p[0])
        ],
    )
    table(
        "Statistical tests",
        [
            (f'xt.run_test("{n}")', s.stata, s.r, s.summary)
            for n, s in sorted(TESTS.specs(), key=lambda p: p[0])
        ],
    )
    return "\n".join(out).rstrip()


# ------------------------------------------------------------------ annex B


def signature_text(function: object) -> str:
    """``(a: int, *, b: str = 'x')`` without the quotes the postponed annotations add."""
    parameters = list(inspect.signature(function).parameters.values())  # type: ignore[arg-type]
    pieces: list[str] = []
    star = False
    for parameter in parameters:
        if parameter.kind is parameter.KEYWORD_ONLY and not star:
            pieces.append("*")
            star = True
        text = parameter.name
        if parameter.annotation is not parameter.empty:
            text += f": {parameter.annotation}"
        if parameter.default is not parameter.empty:
            text += f" = {parameter.default!r}"
        pieces.append(text)
    return "(" + ", ".join(pieces) + ")"


def annex_b() -> str:
    """One card per test, by family."""
    from pyxtpanel.diagnostics._core import TESTS

    cards: list[str] = []
    by_family: dict[str, list[str]] = {}
    for name, spec in sorted(TESTS.specs(), key=lambda p: (p[1].family, p[0])):
        signature = signature_text(TESTS.get(name))
        kind = "panel (category A)" if spec.category == "A" else "arrays (category B)"
        balanced = {
            "required": "balanced panel required",
            "both": "balanced or unbalanced",
            "na": "no panel structure",
        }[spec.balanced]
        lines = [
            f"#### `{name}`",
            "",
            spec.summary,
            "",
            "```python",
            textwrap.fill(
                f"{name}{signature}",
                width=88,
                subsequent_indent="    ",
                break_long_words=False,
                break_on_hyphens=False,
            ),
            "```",
            "",
            f"- **Input**: {kind}; {balanced}",
            f"- **Family**: {spec.family}",
            f"- **Stata**: {cell(spec.stata)}",
            f"- **R**: {cell(spec.r)}",
            f"- **Reference**: {cell(spec.reference)}",
            f"- **Validation**: {len(spec.validation)} cases run in Stata 17 and frozen",
        ]
        if spec.aliases:
            lines.append(f"- **Aliases**: {', '.join(f'`{a}`' for a in spec.aliases)}")
        by_family.setdefault(spec.family, []).append("\n".join(lines))
    for family, family_cards in by_family.items():
        cards.append(f"### {family.capitalize()}\n\n" + "\n\n".join(family_cards))
    return "\n\n".join(cards)


# ------------------------------------------------------------------ annex C


def annex_c() -> str:
    """Return codes kept and error classes."""
    import pyxtpanel._errors as errors

    rows = ["| Code | Meaning |", "|---:|---|"]
    rows += [f"| {code} | {cell(text)} |" for code, text in sorted(errors.STATA_CODES.items())]
    classes = ["| Class | Raised when |", "|---|---|"]
    for name, obj in sorted(vars(errors).items()):
        if inspect.isclass(obj) and obj.__module__ == errors.__name__ and name != "STATA_CODES":
            first = (inspect.getdoc(obj) or "").splitlines()[0] if inspect.getdoc(obj) else ""
            classes.append(f"| `{name}` | {cell(first)} |")
    return (
        "### Stata return codes kept\n\n"
        + "\n".join(rows)
        + "\n\n### Error classes\n\n"
        + "\n".join(classes)
    )


# ------------------------------------------------------------------ annex D


def quirks() -> list[tuple[str, str]]:
    """``(module, text)`` of each ``# STATA-QUIRK:`` comment, with its continuation lines joined."""
    found: list[tuple[str, str]] = []
    for path in sorted((ROOT / "src" / "pyxtpanel").rglob("*.py")):
        lines = path.read_text(encoding="utf-8").splitlines()
        module = path.relative_to(ROOT / "src").with_suffix("").as_posix().replace("/", ".")
        i = 0
        while i < len(lines):
            match = re.match(r"(\s*)# STATA-QUIRK:\s*(.+)", lines[i])
            if match:
                indent, text = match.group(1), match.group(2).strip()
                while i + 1 < len(lines) and re.match(
                    rf"{indent}#\s+(?!STATA-QUIRK|Stata:)\S", lines[i + 1]
                ):
                    i += 1
                    text += " " + lines[i].strip().lstrip("#").strip()
                found.append((module, text))
            i += 1
    return found


def annex_d() -> str:
    """Behaviours of Stata that the package reproduces on purpose, by module."""
    by_module: dict[str, list[str]] = {}
    for module, text in quirks():
        by_module.setdefault(module, []).append(text)
    out = [
        "Each item is a place where the package follows what Stata does rather than the textbook, "
        "and where a test compares the result with Stata 17.",
        "",
    ]
    for module, texts in by_module.items():
        out.append(f"**`{module}`**")
        out.append("")
        out.extend(f"- {text}" for text in texts)
        out.append("")
    return "\n".join(out).rstrip()


# ------------------------------------------------------------------ annex E


def annex_e() -> str:
    """Validation journal: each frozen reference file, the Stata that produced it, its cases."""
    rows = ["| Domain | Stata | Frozen on | Cases |", "|---|---|---|---:|"]
    total = 0
    for path in sorted((ROOT / "tests" / "golden").glob("*_golden.json")):
        master = json.loads(path.read_text(encoding="utf-8"))
        meta = master.get("meta", {})
        cases = master.get("cases", {})
        total += len(cases)
        rows.append(
            f"| `{meta.get('domain', path.stem)}` | {cell(meta.get('stata'))} | "
            f"{cell(meta.get('frozen_on'))} | {len(cases)} |"
        )
    rows.append(f"| **Total** | | | **{total}** |")
    return "\n".join(rows)


RENDER = {"a": annex_a, "b": annex_b, "c": annex_c, "d": annex_d, "e": annex_e}


# ------------------------------------------------------------------ the document


def replace_annex(text: str, key: str, body: str) -> str:
    """``text`` with the region of annex ``key`` replaced by ``body``."""
    start, end = MARK.format(key=key, edge="start"), MARK.format(key=key, edge="end")
    a, b = text.index(start) + len(start), text.index(end)
    return text[:a] + "\n\n" + body + "\n\n" + text[b:]


def rendered(text: str) -> str:
    for key in ANNEXES:
        text = replace_annex(text, key, RENDER[key]())
    return text


def parts(text: str) -> dict[int, str]:
    """Text of each numbered part (``## 7. Statistical tests``), by number."""
    found: dict[int, str] = {}
    pieces = re.split(r"^## (\d+)\. ", text, flags=re.M)
    for number, body in zip(pieces[1::2], pieces[2::2], strict=True):
        found[int(number)] = body
    return found


def python_blocks(text: str) -> list[tuple[int, str]]:
    """``(part, code)`` of the Python blocks of parts 0 to 8 that are meant to run.

    A block whose info string is ``python skip`` is shown but not run.
    """
    blocks: list[tuple[int, str]] = []
    for number, body in parts(text).items():
        if number > 8:
            continue
        for match in re.finditer(r"^```python\n(.*?)^```", body, flags=re.M | re.S):
            blocks.append((number, match.group(1)))
    return blocks


def coverage(text: str) -> list[str]:
    """Names of registered features that the part meant to teach them does not mention."""
    from pyxtpanel.describe._spec import DESCRIPTIONS
    from pyxtpanel.diagnostics._core import TESTS
    from pyxtpanel.explore._spec import TOOLS
    from pyxtpanel.ops._spec import OPERATORS
    from pyxtpanel.plot._spec import KINDS

    sections = parts(text)
    wanted = {
        3: [n for n, _ in OPERATORS.specs()],
        4: [n for n, _ in DESCRIPTIONS.specs()],
        5: [n for n, _ in TOOLS.specs()],
        6: [n for n, _ in KINDS.specs()],
        7: [n for n, _ in TESTS.specs()],
    }
    missing = []
    for number, names in wanted.items():
        body = sections.get(number, "")
        for name in names:
            mentioned = re.search(rf"`{name}[`(]|\"{name}\"|{name}\(", body)
            if not mentioned:
                missing.append(f"part {number}: `{name}` is not mentioned")
    return missing


# ------------------------------------------------------------------ LaTeX


def _weights(rows: list[list[str]]) -> list[int]:
    """Relative width of each column of a table, from the length of its cells."""
    widths = []
    for column in zip(*rows, strict=True):
        lengths = sorted(len(cell_text.strip()) for cell_text in column)
        typical = lengths[int(0.9 * (len(lengths) - 1))]  # a long outlier does not decide
        # a code span cannot be broken: its column must be wide enough for the longest one
        code = max(
            (len(span) for cell_text in column for span in re.findall(r"`([^`]+)`", cell_text)),
            default=0,
        )
        widths.append(min(max(typical, int(code * 1.7) + 2, 6), 90))
    return widths


def weighted_tables(text: str) -> str:
    """Pipe tables whose separator row carries the relative widths pandoc uses for long cells."""
    lines = text.split("\n")
    out: list[str] = []
    i = 0
    while i < len(lines):
        if (
            lines[i].startswith("|")
            and i + 1 < len(lines)
            and re.fullmatch(r"\|[-:| ]+\|", lines[i + 1])
        ):
            j = i
            while j < len(lines) and lines[j].startswith("|"):
                j += 1
            rows = [
                [c for c in re.split(r"(?<!\\)\|", line.strip())[1:-1]]
                for line in lines[i:j]
                if not re.fullmatch(r"\|[-:| ]+\|", line)
            ]
            weights = _weights(rows)
            aligns = [c.strip() for c in lines[i + 1].strip().strip("|").split("|")]
            seps = []
            for weight, align in zip(weights, aligns, strict=True):
                dashes = "-" * weight
                left = ":" if align.startswith(":") else ""
                right = ":" if align.endswith(":") else ""
                seps.append(left + dashes + right)
            out.append(lines[i])
            out.append("|" + "|".join(seps) + "|")
            out.extend(lines[i + 2 : j])
            i = j
        else:
            out.append(lines[i])
            i += 1
    return "\n".join(out)


def latex_markdown(text: str) -> str:
    """The Markdown prepared for pandoc: titles, anchors and table widths."""
    text = re.sub(r"\A# .*\n", "", text)  # the title goes in the metadata
    text = text.replace("\n# Annexes\n", "\n")
    text = text.replace("```python skip", "```python")  # pandoc takes one word as the language
    text = re.sub(r"\]\(#\d+-", "](#", text)  # pandoc drops the number from a heading id
    return weighted_tables(text)


def build_latex() -> Path:
    """Write ``docs/documentation.tex`` from the Markdown with pandoc."""
    if shutil.which("pandoc") is None:
        raise SystemExit("pandoc is needed for --latex")
    from pyxtpanel import __version__

    with tempfile.TemporaryDirectory() as folder:
        source = Path(folder) / "documentation.md"
        source.write_text(latex_markdown(TARGET.read_text(encoding="utf-8")), encoding="utf-8")
        subprocess.run(
            [
                "pandoc", str(source), "-f", "markdown+tex_math_dollars", "-t", "latex", "-s",
                "-o", str(TEX), "--toc", "--toc-depth=2", "--top-level-division=section",
                "-H", str(HEADER),
                "-M", "title=pyxtpanel documentation", "-M", f"subtitle=version {__version__}",
                "-V", "documentclass=article", "-V", "papersize=a4", "-V", "fontsize=10pt",
                "-V", "geometry:margin=2.3cm", "-V", "colorlinks=true", "-V", "linkcolor=blue",
                "-V", "urlcolor=blue", "-V", "mainfont=TeX Gyre Pagella",
                "-V", "monofont=TeX Gyre Cursor", "-V", "monofontoptions=Scale=0.88",
                "-V", "lang=en",
            ],
            check=True,
        )  # fmt: skip
    return TEX


def build_pdf() -> Path:
    """Compile ``docs/documentation.tex`` with xelatex (twice, for the table of contents)."""
    engine = shutil.which("xelatex")
    if engine is None:
        raise SystemExit("xelatex is needed for --pdf")
    build_latex()
    with tempfile.TemporaryDirectory() as folder:
        for _ in range(2):
            run = subprocess.run(
                [
                    engine,
                    "-interaction=nonstopmode",
                    "-halt-on-error",
                    f"-output-directory={folder}",
                    str(TEX),
                ],
                capture_output=True,
                text=True,
                cwd=TEX.parent,
            )
            if run.returncode:
                print(run.stdout[-3000:])
                raise SystemExit("xelatex failed")
        shutil.copy(Path(folder) / "documentation.pdf", PDF)
        log = (Path(folder) / "documentation.log").read_text(encoding="utf-8", errors="replace")
    missing = sorted(set(re.findall(r"Missing character: There is no (.+?) in font", log)))
    overfull = len(re.findall(r"Overfull \\hbox", log))
    print(f"pdf: {PDF}; missing characters: {missing or 'none'}; overfull boxes: {overfull}")
    return PDF


def main() -> int:
    if "--pdf" in sys.argv:
        build_pdf()
        return 0
    if "--latex" in sys.argv:
        print(f"wrote {build_latex()}")
        return 0
    text = TARGET.read_text(encoding="utf-8")
    if "--write" in sys.argv:
        TARGET.write_text(rendered(text), encoding="utf-8")
        print(f"wrote {TARGET}")
        return 0
    problems = []
    if rendered(text) != text:
        problems.append("the annexes are out of date: python -m tools.doc --write")
    problems += coverage(text)
    for problem in problems:
        print(problem)
    print("documentation: up to date" if not problems else f"{len(problems)} problem(s)")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
