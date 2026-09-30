"""Tests for Russell et al. (2018) diagnostic calculations."""

import numpy as np

from esmvaltool.diag_scripts.russell18jgr.russell_common import (
    southern_ocean_flux_sum,
)


def test_southern_ocean_flux_sum_preserves_input_masks():
    """Masked land cells must not contribute to integrated ocean flux."""
    flux = np.ma.array([[1.0, 2.0], [3.0, 4.0]], mask=[[0, 1], [0, 0]])
    area = np.ma.array([[10.0, 10.0], [10.0, 10.0]], mask=[[0, 0], [1, 0]])

    result = southern_ocean_flux_sum(flux, area, [-60.0, -30.0], 1.0)

    assert result == 50.0
