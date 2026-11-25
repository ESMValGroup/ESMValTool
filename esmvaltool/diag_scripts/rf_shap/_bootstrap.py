"""Make the rf_shap scientific kernel importable for the diagnostic.

Preference order (no vendored copy):
1. Already-installed ``rf_shap`` package (pip/conda dependency)
2. Development repo ``src/`` next to an editable overlay checkout
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

_DIAG_DIR = Path(__file__).resolve().parent
_DEV_SRC = _DIAG_DIR.parents[3] / "src"


def ensure_rf_shap_on_path() -> Path | None:
    if importlib.util.find_spec("rf_shap") is not None:
        return None
    if (_DEV_SRC / "rf_shap").is_dir():
        path = str(_DEV_SRC)
        if path not in sys.path:
            sys.path.insert(0, path)
        return _DEV_SRC
    raise ModuleNotFoundError(
        "rf_shap is not installed. Install the standalone package into the "
        "ESMValTool environment, e.g. `pip install -e /path/to/RF_SHAP_ESMValTool[explain]` "
        "(see issue ESMValGroup/ESMValTool#4609)."
    )
