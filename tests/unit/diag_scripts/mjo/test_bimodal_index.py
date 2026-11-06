# Copyright (C) 2026 ESMValTool development team
"""Scientific invariants of the MJO/BSISO EEOF calculation."""

from datetime import UTC, datetime, timedelta

import numpy as np
import pytest

from esmvaltool.diag_scripts.mjo.bimodal_index import (
    _bandpass,
    _fit_season,
    _lagged_matrix,
)


def test_filtered_lags_remain_on_continuous_daily_axis():
    """Filtering and adding lags must not join separate seasons."""
    days = np.arange("2001-01-01", "2003-01-01", dtype="datetime64[D]")
    signal = np.sin(2 * np.pi * np.arange(len(days)) / 45)
    values = signal[:, np.newaxis, np.newaxis]
    filtered, filtered_days = _bandpass(values, days, 25, 90, 139)
    matrix, dates = _lagged_matrix(
        filtered,
        filtered_days,
        [-10, -5, 0],
        [0],
    )

    assert len(filtered_days) == len(days) - 138
    assert np.all(np.diff(dates) == np.timedelta64(1, "D"))
    assert dates[0] == filtered_days[10]
    np.testing.assert_allclose(matrix[:, 0], filtered[:-10, 0, 0])
    np.testing.assert_allclose(matrix[:, 2], filtered[10:, 0, 0])


def test_pc_sign_flip_preserves_reconstruction_and_unit_variance():
    """Paired PC/EEOF sign changes retain the physical anomaly fields."""
    days = np.array(
        [
            datetime(2001, 1, 1, tzinfo=UTC) + timedelta(days=day)
            for day in range(730)
        ],
    )
    time = np.arange(len(days))
    matrix = np.column_stack(
        (
            np.sin(2 * np.pi * time / 45),
            np.cos(2 * np.pi * time / 55),
        ),
    )
    ordinary = _fit_season(matrix, days, [6, 7, 8], (1, 1, 2))
    flipped = _fit_season(
        matrix,
        days,
        [6, 7, 8],
        (1, 1, 2),
        flip_pc2=True,
    )
    selected = ordinary["season"]

    np.testing.assert_allclose(flipped["training"].std(axis=0), [1, 1])
    np.testing.assert_allclose(
        flipped["projected"][selected],
        flipped["training"],
    )
    np.testing.assert_allclose(
        ordinary["projected"][:, 0],
        flipped["projected"][:, 0],
    )
    np.testing.assert_allclose(
        ordinary["projected"][:, 1],
        -flipped["projected"][:, 1],
    )
    for result in (ordinary, flipped):
        patterns = result["eeofs"].reshape(2, 2)
        seasonal_mean = matrix[selected].mean(axis=0)
        reconstruction = result["projected"] @ patterns + seasonal_mean
        np.testing.assert_allclose(reconstruction, matrix, atol=1e-10)


def test_filter_rejects_even_window():
    days = np.arange("2001-01-01", "2002-01-01", dtype="datetime64[D]")
    with pytest.raises(ValueError, match="odd integer"):
        _bandpass(np.zeros((len(days), 1, 1)), days, 25, 90, 140)
