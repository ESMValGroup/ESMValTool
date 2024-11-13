"""Ridge Regression model with built-in CV.

Use ``mlr_model_type: ridge_cv`` to use this MLR model in the recipe.

"""

import logging
import os

from sklearn.linear_model import RidgeCV

from esmvaltool.diag_scripts.mlr.models import MLRModel
from esmvaltool.diag_scripts.mlr.models.linear_base import LinearModel

logger = logging.getLogger(os.path.basename(__file__))


@MLRModel.register_mlr_model('ridge_cv')
class RidgeCVModel(LinearModel):
    """Ridge Regression model with built-in CV."""

    _CLF_TYPE = RidgeCV

    def fit(self):
        """Print final ``alpha`` after successful fitting."""
        super().fit()
        logger.info("Optimal alpha of Ridge model: α = %.5f",
                    self._clf.steps[-1][1].regressor_.alpha_)
