"""Tests for toz variable derivation functions."""

import numpy as np
import pytest

from esmvalcore.preprocessor._derive.toz import _p_level_widths


def test_col_is_not_monotonic():
    """Test for non-monotonic column."""
    plev = 1000
    top_limit = 5
    col = np.array([1, 2, 3, 2, 1])
    col = np.insert(col, 0, plev)
    col = np.append(col, top_limit)
    with pytest.raises(ValueError):
        _p_level_widths(col)


def test_keeping_column_length():
    """Test for level widths keeping column lenght."""
    plev = 1000
    top_limit = 5
    col = np.array([1000, 900, 800])
    col = np.insert(col, 0, plev)
    col = np.append(col, top_limit)
    assert len(_p_level_widths(col)) == len(col) - 2


def test_low_lev_surf_press():
    """Test for lowest level equal to surface pressure."""
    plev = 1000
    top_limit = 5
    col = np.array([1000, 900, 800])
    col = np.insert(col, 0, plev)
    col = np.append(col, top_limit)
    result = np.array([50, 100, 845])
    assert all(_p_level_widths(col) == result)


def test_low_lev_above_surf_press():
    """Test for lowest level above surface pressure."""
    plev = 1020
    top_limit = 5
    col = np.array([1000, 900, 800])
    col = np.insert(col, 0, plev)
    col = np.append(col, top_limit)
    result = np.array([70, 100, 845])
    assert all(_p_level_widths(col) == result)


def test_low_lev_below_surf_press():
    """Test for lowest level below surface pressure."""
    plev = 970
    top_limit = 5
    col = np.array([np.NaN, 900, 800])
    col = np.insert(col, 0, plev)
    col = np.append(col, top_limit)
    result = np.array([0, 120, 845])
    assert all(_p_level_widths(col) == result)

    col = np.array([np.NaN, np.NaN, 900, 800])
    col = np.insert(col, 0, plev)
    col = np.append(col, top_limit)
    result = np.array([0, 0, 120, 845])
    assert all(_p_level_widths(col) == result)


def test_high_level_top_limit():
    """Test for highest level equal to top limit."""
    plev = 1020
    top_limit = 5
    col = np.array([1000, 900, 5])
    col = np.insert(col, 0, plev)
    col = np.append(col, top_limit)
    result = np.array([70, 50 + 895 / 2, 895 / 2])
    assert all(_p_level_widths(col) == result)


def test_high_level_above_top_limit():
    """Test for highest level above top limit."""
    plev = 1020
    top_limit = 5
    col = np.array([1000, 900, 3])
    col = np.insert(col, 0, plev)
    col = np.append(col, top_limit)
    with pytest.raises(ValueError):
        _p_level_widths(col)
