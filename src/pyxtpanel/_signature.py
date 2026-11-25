"""Check a call against the signature of the function it reaches (M.2.1).

The dispatchers (``run_panel_test`` for the tests, ``plot`` for the figures) call a registered
function with whatever the user typed.  Python's own ``TypeError`` ("missing a required argument")
says neither which function was meant nor what it accepts, so the dispatchers bind the call to the
signature first and report the problem themselves.  The two helpers here are the part they share.
"""

from __future__ import annotations

import inspect
import re

__all__: list[str] = []


def reason(signature: inspect.Signature, error: TypeError) -> str:
    """Return Python's reason for the refused call, worded the same on every Python version.

    Python 3.12 says "missing a required keyword-only argument" for an option that can only be
    given by name; 3.10 and 3.11 say "missing a required argument" for it.  The message is shown
    to users and compared in tests, so the 3.12 wording is used everywhere.

    Parameters
    ----------
    signature
        Signature of the function the call was bound to.
    error
        The ``TypeError`` raised by ``signature.bind``.

    Returns
    -------
    str
        The reason, to put in an error message.
    """
    text = str(error)
    match = re.fullmatch(r"missing a required argument: '(\w+)'", text)
    parameter = signature.parameters.get(match.group(1)) if match else None
    if parameter is not None and parameter.kind is inspect.Parameter.KEYWORD_ONLY:
        return f"missing a required keyword-only argument: '{parameter.name}'"
    return text


def accepted(signature: inspect.Signature) -> str:
    """Return the arguments after the first (the panel), with the default of each option.

    Parameters
    ----------
    signature
        Signature of the function.

    Returns
    -------
    str
        Comma-separated list such as ``"variable, n_units=12, mean=True"``.
    """
    parts = []
    for index, (name, parameter) in enumerate(signature.parameters.items()):
        if index == 0:
            continue  # the panel
        if parameter.default is inspect.Parameter.empty:
            parts.append(name)
        else:
            parts.append(f"{name}={parameter.default!r}")
    return ", ".join(parts)
