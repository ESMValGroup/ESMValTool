"""Lasso Regression model with built-in CV using LARS algorithm.

Use ``mlr_model_type: lasso_lars_cv`` to use this MLR model in the recipe.

"""

import logging
import os

from sklearn.linear_model import LassoLarsCV

from esmvaltool.diag_scripts.mlr.models import MLRModel
from esmvaltool.diag_scripts.mlr.models.linear_base import LinearModel

logger = logging.getLogger(os.path.basename(__file__))


@MLRModel.register_mlr_model('lasso_lars_cv')
class LassoLarsCVModel(LinearModel):
    """Lasso Regression model with built-in CV using LARS algorithm."""

    _CLF_TYPE = LassoLarsCV

    def fit(self):
        """Print final ``alpha`` after successful fitting."""
        super().fit()
        logger.info("Optimal alpha of Lasso model: α = %.5f",
                    self._clf.steps[-1][1].regressor_.alpha_)
