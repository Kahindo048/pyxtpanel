"""Settings shared by the whole test suite."""

from __future__ import annotations

# hypothesis is a test dependency; if it is missing, the property tests are not collected anyway.
try:
    from hypothesis import settings
except ImportError:
    pass
else:
    # No time limit per example: a shared CI runner (a slow Windows machine, a busy Linux one) can
    # be several times slower than a laptop, and a property that is correct would fail with
    # DeadlineExceeded.  What the properties check is the result, not the speed (the speed guards
    # are in tests/benchmarks).
    settings.register_profile("pyxtpanel", deadline=None)
    settings.load_profile("pyxtpanel")
