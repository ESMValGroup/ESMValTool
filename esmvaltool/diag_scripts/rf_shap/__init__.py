"""rf_shap as an ESMValTool diagnostic package.

Import the scientific API through ESMValTool::

    import esmvaltool.diag_scripts.rf_shap as rf
    from esmvaltool.diag_scripts.rf_shap import DiagnosisPipeline
    from esmvaltool.diag_scripts.rf_shap.models import build_model
    from esmvaltool.diag_scripts.rf_shap.cube import align_xy

The ESMValTool recipe entry point remains ``diagnostic.py``.
"""

import sys

from rf_shap.namespace import install_namespace

install_namespace(sys.modules[__name__])
