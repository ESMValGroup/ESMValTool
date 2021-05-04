"""Ridge Regression model.

Use ``mlr_model_type: ridge`` to use this MLR model in the recipe.

"""

import logging
import os

from sklearn.linear_model import Ridge

from esmvaltool.diag_scripts.mlr.models import MLRModel
from esmvaltool.diag_scripts.mlr.models.linear_base import LinearModel

logger = logging.getLogger(os.path.basename(__file__))


@MLRModel.register_mlr_model('ridge')
class RidgeModel(LinearModel):
    """Ridge Regression model."""

    _CLF_TYPE = Ridge
