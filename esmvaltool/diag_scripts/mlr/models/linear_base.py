"""Base class for linear Machine Learning Regression models."""

import logging
import os

import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.mlr.models import MLRModel

logger = logging.getLogger(os.path.basename(__file__))


class LinearModel(MLRModel):
    """Base class for linear Machine Learning models."""

    _CLF_TYPE = None

    def plot_coefs(self, filename=None):
        """Plot linear coefficients of models.

        Note
        ----
        The features plotted here are not necessarily the real input features,
        but the ones after preprocessing.

        Parameters
        ----------
        filename : str, optional (default: 'coefs')
           Name of the plot file.

        """
        if not self._is_ready_for_plotting():
            return
        logger.info("Plotting linear coefficients")
        if filename is None:
            filename = 'coefs'
        (_, axes) = plt.subplots()

        # Plot
        coefs = self._clf.coef_
        sorted_idx = np.argsort(coefs)
        pos = np.arange(sorted_idx.shape[0]) + 0.5
        axes.barh(pos, coefs[sorted_idx], align='center')

        # Plot appearance
        axes.tick_params(axis='y', which='minor', left=False, right=False)
        axes.tick_params(axis='y', which='major', left=True, right=False)
        y_tick_labels = self.features_after_preprocessing[sorted_idx]
        title = f"Linear coefficients ({self._cfg['mlr_model_name']})"
        axes.set_title(title)
        axes.set_yticks(pos)
        axes.set_yticklabels(y_tick_labels)
        axes.set_xlim(-np.max(np.abs(axes.get_xlim())),
                      np.max(np.abs(axes.get_xlim())))
        axes.axvline(0.0, color='k')

        # Save plot
        new_filename = filename + '.' + self._cfg['output_file_type']
        plot_path = os.path.join(self._cfg['mlr_plot_dir'], new_filename)
        plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save provenance
        cube = mlr.get_1d_cube(
            y_tick_labels,
            coefs[sorted_idx],
            x_kwargs={'var_name': 'feature',
                      'long_name': 'Feature name',
                      'units': 'no unit'},
            y_kwargs={'var_name': 'coef',
                      'long_name': '(Normalized) Linear Coefficients',
                      'units': '1',
                      'attributes': {'project': '', 'dataset': ''}},
        )
        self._write_plot_provenance(
            cube, plot_path, ancestors=self.get_ancestors(prediction_names=[]),
            caption=title + '.', plot_types=['bar'])

    def plot_feature_importance(self, filename=None, color_coded=True):
        """Plot feature importance given by linear coefficients.

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
        coefs = self._clf.coef_
        feature_importances = np.abs(coefs) / np.sum(np.abs(coefs))
        feature_importance_dict = dict(zip(self.features_after_preprocessing,
                                           feature_importances))
        colors = self._get_colors_for_features(color_coded=color_coded)

        # Plot
        self._plot_feature_importance(feature_importance_dict, colors,
                                      plot_path)
