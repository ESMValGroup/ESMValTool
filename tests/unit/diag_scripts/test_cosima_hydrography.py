"""Analytic checks for the hydrographic and density diagnostics."""

import numpy as np

from esmvaltool.diag_scripts.cosima_cmip import density_compensation as dc
from esmvaltool.diag_scripts.cosima_cmip import hydrographic_benchmark as hb


def test_reference_metrics_use_paired_cells():
    """Bias, RMSE, and coverage share a paired model/atlas mask."""
    model = np.ma.array([[[[12.0, 13.0], [12.0, 12.0]]]])
    model[0, 0, 0, 1] = np.ma.masked
    reference = np.ma.array([[[[10.0, 10.0], [10.0, 10.0]]]])
    weights = np.ones((1, 2, 2))

    bias, rmse, coverage = hb._reference_metrics(  # noqa: SLF001
        model, reference, weights
    )

    np.testing.assert_allclose(bias, [[2.0]])
    np.testing.assert_allclose(rmse, [[2.0]])
    np.testing.assert_allclose(coverage, [[0.75]])


def test_control_drift_is_per_century():
    """A 0.01-degree yearly warming gives 1 degree per century."""
    years = np.arange(950.0, 990.0)
    values = (2.0 + 0.01 * (years - years[0]))[:, None, None, None]
    control = np.ma.array(np.broadcast_to(values, (40, 1, 2, 2)))
    weights = np.ones((1, 2, 2))

    means, anomaly, slope, coverage = hb._control_metrics(  # noqa: SLF001
        control, years, weights
    )

    np.testing.assert_allclose(means[:, 0, 0], values[:, 0, 0, 0])
    np.testing.assert_allclose(anomaly[0], 0.0)
    np.testing.assert_allclose(slope, [[1.0]])
    np.testing.assert_allclose(coverage, [[1.0]])


def test_density_components_close_with_opposing_effects():
    """Symmetric thermal and haline contributions sum to total density."""
    shape = (2, 2, 2)
    ref_t = np.ma.array(np.full(shape, 10.0))
    ref_s = np.ma.array(np.full(shape, 35.0))
    model_t = ref_t + 2.0
    model_s = ref_s + 0.1
    latitude = np.array([[-45.0, -45.0], [45.0, 45.0]])
    longitude = np.array([[10.0, 11.0], [10.0, 11.0]])

    thermal, haline, total = dc.density_components(
        model_t,
        model_s,
        ref_t,
        ref_s,
        np.array([100.0, 1000.0]),
        latitude,
        longitude,
    )

    np.testing.assert_allclose(total, thermal + haline, atol=1e-10)
    assert np.all(thermal < 0.0)
    assert np.all(haline > 0.0)
