"""Huber Regression model.

Use ``mlr_model_type: huber`` to use this MLR model in the recipe.

"""

import logging
import os

from sklearn.linear_model import HuberRegressor

from esmvaltool.diag_scripts.mlr.models import MLRModel
from esmvaltool.diag_scripts.mlr.models.linear_base import LinearModel

logger = logging.getLogger(os.path.basename(__file__))


@MLRModel.register_mlr_model('huber')
class HuberRegressionModel(LinearModel):
    """Huber Regression model."""

    _CLF_TYPE = HuberRegressor
