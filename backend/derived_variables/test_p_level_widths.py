from __future__ import print_function
from __future__ import division
import numpy as np
import pytest

from calculate_variables import _p_level_widths


def test_col_is_not_monotonic():
    sp = 1000
    top_limit = 5
    col = np.array([1, 2, 3, 2, 1])
    col = np.insert(col, 0, sp)
    col = np.append(col, top_limit)
    with pytest.raises(ValueError):
        _p_level_widths(col)


def test__p_level_widths_keeps_columns_length():
    sp = 1000
    top_limit = 5
    col = np.array([1000, 900, 800])
    col = np.insert(col, 0, sp)
    col = np.append(col, top_limit)
    assert len(_p_level_widths(col)) == len(col) - 2


def test_lowest_level_is_surface_pressure():
    sp = 1000
    top_limit = 5
    col = np.array([1000, 900, 800])
    col = np.insert(col, 0, sp)
    col = np.append(col, top_limit)
    result = np.array([50, 100, 845])
    assert all(_p_level_widths(col) == result)


def test_lowest_level_is_above_surface_pressure():
    sp = 1020
    top_limit = 5
    col = np.array([1000, 900, 800])
    col = np.insert(col, 0, sp)
    col = np.append(col, top_limit)
    result = np.array([70, 100, 845])
    assert all(_p_level_widths(col) == result)


def test_lowest_level_is_below_surface_pressure():
    sp = 970
    top_limit = 5
    col = np.array([np.NaN, 900, 800])
    col = np.insert(col, 0, sp)
    col = np.append(col, top_limit)
    result = np.array([0, 120, 845])
    assert all(_p_level_widths(col) == result)

    col = np.array([np.NaN, np.NaN, 900, 800])
    col = np.insert(col, 0, sp)
    col = np.append(col, top_limit)
    result = np.array([0, 0, 120, 845])
    assert all(_p_level_widths(col) == result)


def test_highest_level_is_top_limit():
    sp = 1020
    top_limit = 5
    col = np.array([1000, 900, 5])
    col = np.insert(col, 0, sp)
    col = np.append(col, top_limit)
    result = np.array([70, 50 + 895/2, 895/2])
    assert all(_p_level_widths(col) == result)


def test_highest_level_above_top_limit():
    sp = 1020
    top_limit = 5
    col = np.array([1000, 900, 3])
    col = np.insert(col, 0, sp)
    col = np.append(col, top_limit)
    with pytest.raises(ValueError):
        _p_level_widths(col)


