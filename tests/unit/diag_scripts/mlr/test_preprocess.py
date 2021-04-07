"""Unit tests for the module :mod:`esmvaltool.diag_scripts.mlr.preprocess`."""

import numpy as np
import pytest

import esmvaltool.diag_scripts.mlr.preprocess as preprocess

X_ARR = np.arange(5)
TEST_GET_SLOPE = [
    (X_ARR, 3.14 * X_ARR, 3.14),
    (np.arange(1.0), np.arange(1.0), np.nan),
    (X_ARR, np.ma.masked_invalid([np.nan, 1.0, np.nan, 1.5, np.nan]), 0.25),
]


@pytest.mark.parametrize('x_arr,y_arr,output', TEST_GET_SLOPE)
def test_get_slope(x_arr, y_arr, output):
    """Test calculation of slope."""
    out = preprocess._get_slope(x_arr, y_arr)
    assert ((out == output) | (np.isnan(output) & np.isnan(output))).all()


Y_ARR_1 = np.ma.masked_invalid([np.nan, 1.0, 0.0, np.nan, -0.5])
Y_ARR_2 = np.ma.masked_invalid([np.nan, np.nan, np.nan, np.nan, -0.5])
Y_ARR_2x2 = np.ma.masked_invalid(
    [[2.1 * X_ARR, -3.14 * X_ARR, 0.8 * X_ARR],
     [Y_ARR_1.filled(np.nan),
      Y_ARR_1.filled(2.0),
      Y_ARR_2.filled(np.nan)]])
TEST_GET_SLOPE_VECTORIZED = [
    (X_ARR, Y_ARR_2x2,
     np.array([[2.1, -3.14, 0.8], [-0.46428571428571436, -0.4, np.nan]])),
    (X_ARR, np.ma.array([X_ARR, Y_ARR_2]), np.array([1.0, np.nan])),
]


@pytest.mark.parametrize('x_arr,y_arr,output', TEST_GET_SLOPE_VECTORIZED)
def test_get_slope_vectorized(x_arr, y_arr, output):
    """Test vectorized calculation of slope."""
    get_slope = np.vectorize(preprocess._get_slope, excluded=['x_arr'],
                             signature='(n),(n)->()')
    out = get_slope(x_arr, y_arr)
    assert (np.isclose(out, output) | (np.isnan(out) & np.isnan(output))).all()
