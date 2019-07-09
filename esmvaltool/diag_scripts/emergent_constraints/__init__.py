"""Convenience functions for emergent constraints diagnostics."""
import logging

import numpy as np
from scipy import integrate, stats

logger = logging.getLogger(__name__)


def _check_input_arrays(*arrays):
    """Check the shapes of multiple arrays."""
    shape = None
    for array in arrays:
        if shape is None:
            shape = array.shape
        else:
            if array.shape != shape:
                raise ValueError("Expected input arrays with identical shapes")


def standard_prediction_error(x_data, y_data):
    """Return function to calculate standard prediction error.

    The standard prediction error of a linear regression is the error when
    predicting a new value which is not in the original data.

    Parameters
    ----------
    x_data : numpy.array
        x coordinates of the points.
    y_data : numpy.array
        y coordinates of the points.

    Returns
    -------
    callable
        Standard prediction error function for new x values.

    """
    _check_input_arrays(x_data, y_data)
    reg = stats.linregress(x_data, y_data)
    y_estim = reg.slope * x_data + reg.intercept
    n_data = x_data.shape[0]
    see = np.sqrt(np.sum(np.square(y_data - y_estim)) / (n_data - 2))
    x_mean = np.mean(x_data)
    ssx = np.sum(np.square(x_data - x_mean))

    def spe(x_new):
        """Return standard prediction error."""
        return see * np.sqrt(1.0 + 1.0 / n_data + (x_new - x_mean)**2 / ssx)

    return np.vectorize(spe)


def regression_line(x_data, y_data, n_points=100):
    """Return x and y coordinates of the regression line (mean and error).

    Parameters
    ----------
    x_data : numpy.array
        x coordinates of the points.
    y_data : numpy.array
        y coordinates of the points.
    n_points : int, optional (default: 100)
        Number of points for the regression lines.

    Returns
    -------
    dict
        `numpy.array`s for the keys `x`, `y_best_estim`, `y_minus_err`,
        `y_plus_err', 'rvalue', 'slope' and 'intercept'.

    """
    _check_input_arrays(x_data, y_data)
    spe = standard_prediction_error(x_data, y_data)
    out = {}
    reg = stats.linregress(x_data, y_data)
    x_range = max(x_data) - min(x_data)
    x_lin = np.linspace(min(x_data) - x_range, max(x_data) + x_range, n_points)
    out['y_best_estim'] = reg.slope * x_lin + reg.intercept
    out['y_minus_err'] = out['y_best_estim'] - spe(x_lin)
    out['y_plus_err'] = out['y_best_estim'] + spe(x_lin)
    out['x'] = x_lin
    out['rvalue'] = reg.rvalue
    out['slope'] = reg.slope
    out['intercept'] = reg.intercept
    return out


def gaussian_pdf(x_data, y_data, obs_mean, obs_std, n_points=100):
    """Calculate Gaussian probability densitiy function for target variable.

    Parameters
    ----------
    x_data : numpy.array
        x coordinates of the points.
    y_data : numpy.array
        y coordinates of the points.
    obs_mean : float
        Mean of observational data.
    obs_std : float
        Standard deviation of observational data.
    n_points : int, optional (default: 100)
        Number of points for the regression lines.

    Returns
    -------
    tuple of numpy.array
        x and y values for the PDF.

    """
    _check_input_arrays(x_data, y_data)
    spe = standard_prediction_error(x_data, y_data)
    reg = stats.linregress(x_data, y_data)

    def obs_pdf(x_new):
        """Return PDF of observations P(x)."""
        norm = np.sqrt(2.0 * np.pi * obs_std**2)
        return np.exp(-(x_new - obs_mean)**2 / 2.0 / obs_std**2) / norm

    def cond_pdf(x_new, y_new):
        """Return conditional PDF P(y|x)."""
        y_estim = reg.slope * x_new + reg.intercept
        norm = np.sqrt(2.0 * np.pi * spe(x_new)**2)
        return np.exp(-(y_new - y_estim)**2 / 2.0 / spe(x_new)**2) / norm

    def comb_pdf(x_new, y_new):
        """Return combined PDF P(y,x)."""
        return obs_pdf(x_new) * cond_pdf(x_new, y_new)

    # PDF of target variable P(y)
    y_range = max(y_data) - min(y_data)
    y_lin = np.linspace(min(y_data) - y_range, max(y_data) + y_range, n_points)
    y_pdf = [
        integrate.quad(comb_pdf, -np.inf, +np.inf, args=(y, ))[0]
        for y in y_lin
    ]
    return (y_lin, np.array(y_pdf))


def cdf(data, pdf):
    """Calculate cumulative distribution function for a PDF.

    Parameters
    ----------
    data : numpy.array
        Data points (x axis).
    pdf : numpy.array
        Corresponding probability density function (PDF).

    Returns
    -------
    numpy.array
        Corresponding cumulative distribution function (CDF).

    """
    idx_range = range(1, len(data) + 1)
    cum_dens = [integrate.simps(pdf[:idx], data[:idx]) for idx in idx_range]
    return np.array(cum_dens)
