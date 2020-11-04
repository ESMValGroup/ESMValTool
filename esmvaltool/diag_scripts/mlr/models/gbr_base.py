"""Base class for Gradient Boosting Regression model."""

import logging
import os

import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.mlr.models import MLRModel

logger = logging.getLogger(os.path.basename(__file__))


class GBRModel(MLRModel):
    """Base class for Gradient Boosting Regression models."""

    _CLF_TYPE = None

    def plot_feature_importance(self, filename=None, color_coded=True):
        """Plot feature importance.

        This function uses properties of the GBR model based on the number of
        appearances of that feature in the regression trees and the
        improvements made by the individual splits (see Friedman, 2001).

        Note
        ----
        The features plotted here are not necessarily the real input features,
        but the ones after preprocessing.

        Parameters
        ----------
        filename : str, optional (default: 'feature_importance')
           Name of the plot file.
        color_coded : bool, optional (default: True)
            If ``True``, mark positive (linear) correlations with red bars and
            negative (linear) correlations with blue bars. If ``False``, all
            bars are blue.

        """
        if not self._is_ready_for_plotting():
            return

        # Get plot path
        if filename is None:
            filename = 'feature_importance'
        new_filename = filename + '.' + self._cfg['output_file_type']
        plot_path = os.path.join(self._cfg['mlr_plot_dir'], new_filename)

        # Get feature importance dictionary and colors for bars
        feature_importance_dict = dict(zip(self.features_after_preprocessing,
                                           self._clf.feature_importances_))
        colors = self._get_colors_for_features(color_coded=color_coded)

        # Plot
        self._plot_feature_importance(feature_importance_dict, colors,
                                      plot_path)

    def _plot_training_progress(self,
                                train_score,
                                test_score=None,
                                filename=None):
        """Plot training progress during fitting."""
        if not self._is_ready_for_plotting():
            return
        logger.info("Plotting training progress for GBR model")
        if filename is None:
            filename = 'training_progress'
        (_, axes) = plt.subplots()
        x_values = np.arange(len(train_score), dtype=np.float64) + 1.0
        x_values_all = []
        scores_all = []
        data_types = []

        # Plot train score
        axes.plot(x_values,
                  train_score,
                  color='b',
                  linestyle='-',
                  label='train data')
        x_values_all.append(x_values)
        scores_all.append(train_score)
        data_types.append(np.full(x_values.shape, 'train'))

        # Plot test score if possible
        if test_score is not None:
            axes.plot(x_values,
                      test_score,
                      color='g',
                      linestyle='-',
                      label='test data')
            x_values_all.append(x_values)
            scores_all.append(test_score)
            data_types.append(np.full(x_values.shape, 'test'))

        # Appearance
        ylim = axes.get_ylim()
        axes.set_ylim(0.0, ylim[1])
        title = f"Training progress ({self._cfg['mlr_model_name']})"
        axes.set_title(title)
        axes.set_xlabel('Boosting iterations')
        axes.set_ylabel('Normalized RMSE')
        axes.legend(loc='upper right')
        new_filename = filename + '.' + self._cfg['output_file_type']
        plot_path = os.path.join(self._cfg['mlr_plot_dir'], new_filename)
        plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save provenance
        cube = mlr.get_1d_cube(
            np.concatenate(x_values_all),
            np.concatenate(scores_all),
            x_kwargs={'var_name': 'iteration',
                      'long_name': 'Boosting Iteration',
                      'units': 'no unit'},
            y_kwargs={'var_name': 'rmse',
                      'long_name': 'Normalized RMSE',
                      'units': '1',
                      'attributes': {'project': '', 'dataset': ''}},
        )
        cube.add_aux_coord(
            self._get_data_type_coord(np.concatenate(data_types)), 0)
        self._write_plot_provenance(
            cube, plot_path, ancestors=self.get_ancestors(prediction_names=[]),
            caption=title + '.', plot_types=['line'])
