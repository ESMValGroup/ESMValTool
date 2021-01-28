"""Random Forest Regression model.

Use ``mlr_model_type: rfr`` to use this MLR model in the recipe.

"""

import logging
import os

from sklearn.ensemble import RandomForestRegressor

from esmvaltool.diag_scripts.mlr.models import MLRModel

logger = logging.getLogger(os.path.basename(__file__))


@MLRModel.register_mlr_model('rfr')
class RFRModel(MLRModel):
    """Random Forest Regression model."""

    _CLF_TYPE = RandomForestRegressor
