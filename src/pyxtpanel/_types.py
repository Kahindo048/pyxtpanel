"""Type aliases shared by the whole package (layer 0)."""

from __future__ import annotations

from typing import Literal

import numpy as np
import numpy.typing as npt

FloatArray = npt.NDArray[np.float64]
IntArray = npt.NDArray[np.int64]
BoolArray = npt.NDArray[np.bool_]

#: Language of messages and reports.
Lang = Literal["en", "fr"]
#: Languages the package is translated into.
LANGUAGES: tuple[Lang, ...] = ("en", "fr")
