"""Gradient Boosting Regression model (using :mod:`xgboost`).

Use ``mlr_model_type: gbr_xgboost`` to use this MLR model in the recipe.

"""

import logging
import os

from xgboost import XGBRegressor

from esmvaltool.diag_scripts.mlr.models import MLRModel
from esmvaltool.diag_scripts.mlr.models.gbr_base import GBRModel

logger = logging.getLogger(os.path.basename(__file__))


@MLRModel.register_mlr_model('gbr_xgboost')
class XGBoostGBRModel(GBRModel):
    """Gradient Boosting Regression model (:mod:`xgboost` implementation)."""

    _CLF_TYPE = XGBRegressor

    def plot_training_progress(self, filename=None):
        """Plot training progress for training and (if possible) test data.

        Parameters
        ----------
        filename : str, optional (default: 'training_progress')
            Name of the plot file.

        """
        clf = self._clf.steps[-1][1].regressor_
        if not hasattr(clf, 'evals_result_'):
            raise AttributeError(
                "Plotting training progress for XGBRegressor model is not "
                "possible, necessary attribute 'evals_result_' is missing. "
                "This is usually cause by calling MLRModel.rfecv()")
        evals_result = clf.evals_result()
        train_score = evals_result['validation_0']['rmse']
        test_score = None
        if 'test' in self.data:
            test_score = evals_result['validation_1']['rmse']
        self._plot_training_progress(train_score, test_score=test_score,
                                     filename=filename)

    def _update_fit_kwargs(self, fit_kwargs):
        """Add transformed training and test data as fit kwargs."""
        fit_kwargs = super()._update_fit_kwargs(fit_kwargs)

        # Fit all transformers
        x_train = self.get_x_array('train')
        y_train = self.get_y_array('train')
        self._clf.fit_transformers_only(x_train, y_train, **fit_kwargs)
        self._clf.fit_target_transformer_only(y_train, **fit_kwargs)

        # Transform input data
        x_train = self._clf.transform_only(x_train)
        y_train = self._clf.transform_target_only(y_train)
        eval_set = [(x_train, y_train)]
        sample_weights = [self._get_sample_weights('train')]
        if 'test' in self.data:
            x_test = self._clf.transform_only(self.get_x_array('test'))
            y_test = self._clf.transform_target_only(self.get_y_array('test'))
            eval_set.append((x_test, y_test))
            sample_weights.append(self._get_sample_weights('test'))
        if self._get_sample_weights('all') is None:
            sample_weights = None

        # Update kwargs
        fit_kwargs.update({
            f'{self._clf.steps[-1][0]}__regressor__eval_metric':
            'rmse',
            f'{self._clf.steps[-1][0]}__regressor__eval_set':
            eval_set,
            f'{self._clf.steps[-1][0]}__regressor__sample_weight_eval_set':
            sample_weights,
        })
        logger.debug(
            "Updated keyword arguments of final regressor's fit() function "
            "with training and (if possible) test datasets for evaluation of "
            "prediction errors")
        return fit_kwargs
