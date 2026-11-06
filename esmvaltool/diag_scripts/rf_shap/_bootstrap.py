"""Make the rf_shap scientific kernel importable for the diagnostic.

Preference order:
1. Vendored copy at ``_vendor/`` (upstream ESMValTool PR layout)
2. Development repo ``src/`` next to an editable checkout
"""

from __future__ import annotations

import sys
from pathlib import Path

_DIAG_DIR = Path(__file__).resolve().parent
_VENDOR = _DIAG_DIR / "_vendor"
_DEV_SRC = _DIAG_DIR.parents[3] / "src"


def ensure_rf_shap_on_path() -> Path | None:
    for candidate in (_VENDOR, _DEV_SRC):
        if (candidate / "rf_shap").is_dir():
            path = str(candidate)
            if path not in sys.path:
                sys.path.insert(0, path)
            return candidate
    return None
