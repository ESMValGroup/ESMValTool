"""Tests for the zero emissions commitment diagnostic."""

from types import SimpleNamespace

import numpy as np

from esmvaltool.diag_scripts.climate_metrics.zec import calc_zec_x


def test_calc_zec_x_uses_twenty_annual_samples_for_each_model():
    """Each requested year covers 20 samples, with the upper end included."""
    first = np.arange(100, dtype=float)
    second = first**2
    zec = {
        "first": SimpleNamespace(data=first),
        "second": SimpleNamespace(data=second),
    }

    assert calc_zec_x(zec, 20) == {
        "first": np.mean(first[10:30]),
        "second": np.mean(second[10:30]),
    }
    assert calc_zec_x(zec, 50) == {
        "first": np.mean(first[40:60]),
        "second": np.mean(second[40:60]),
    }
