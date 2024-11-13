"""Gaussian Process Regression model (using :mod:`sklearn`).

Use ``mlr_model_type: gpr_sklearn`` to use this MLR model in the recipe.

"""

# pylint: disable=arguments-differ

import logging
import os

from sklearn.gaussian_process import GaussianProcessRegressor

from esmvaltool.diag_scripts.mlr.models import MLRModel

logger = logging.getLogger(os.path.basename(__file__))


class AdvancedGaussianProcessRegressor(GaussianProcessRegressor):
    """Expand :class:`sklearn.gaussian_process.GaussianProcessRegressor`."""

    def predict(self, x_data, return_var=False, return_cov=False):
        """Expand :meth:`predict` to accept ``return_var``."""
        pred = super().predict(x_data, return_std=return_var,
                               return_cov=return_cov)
        if return_var:
            return (pred[0], pred[1]**2)
        return pred


@MLRModel.register_mlr_model('gpr_sklearn')
class SklearnGPRModel(MLRModel):
    """Gaussian Process Regression model (:mod:`sklearn` implementation)."""

    _CLF_TYPE = AdvancedGaussianProcessRegressor

    def print_kernel_info(self):
        """Print information of the fitted kernel of the GPR model."""
        self._check_fit_status('Printing kernel')
        kernel = self._clf.steps[-1][1].regressor_.kernel_
        logger.info("Fitted kernel: %s", kernel)
        logger.info("All fitted log-hyperparameters:")
        for (idx, hyper_param) in enumerate(kernel.hyperparameters):
            logger.info("%s: %s", hyper_param, kernel.theta[idx])
