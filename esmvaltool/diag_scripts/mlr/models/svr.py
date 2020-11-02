"""Support Vector Regression model.

Use ``mlr_model_type: svr`` to use this MLR model in the recipe.

"""

import logging
import os

from sklearn.svm import SVR

from esmvaltool.diag_scripts.mlr.models import MLRModel

logger = logging.getLogger(os.path.basename(__file__))


@MLRModel.register_mlr_model('svr')
class SVRModel(MLRModel):
    """Support Vector Regression model."""

    _CLF_TYPE = SVR
