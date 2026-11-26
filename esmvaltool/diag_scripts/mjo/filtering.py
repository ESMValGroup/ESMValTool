"""Lanczos filtering utilities for MJO diagnostics."""

import numpy as np


def lanczos_weights(low_period, high_period, weights):
    """Return Lanczos band-pass weights for daily data."""
    if weights < 3 or weights % 2 == 0:
        msg = "weights must be an odd integer greater than 1"
        raise ValueError(msg)
    if low_period <= 0 or high_period <= low_period:
        msg = "Expected 0 < low_period < high_period, e.g. 20 and 100 days"
        raise ValueError(msg)

    half_width = weights // 2
    n = np.arange(-half_width, half_width + 1, dtype=float)

    # Cutoff frequencies in cycles per day.
    low_frequency = 1.0 / high_period
    high_frequency = 1.0 / low_period

    coefficients = np.empty(weights, dtype=float)
    coefficients[half_width] = 2.0 * (high_frequency - low_frequency)

    nonzero = n != 0
    coefficients[nonzero] = (
        np.sin(2.0 * np.pi * high_frequency * n[nonzero])
        - np.sin(2.0 * np.pi * low_frequency * n[nonzero])
    ) / (np.pi * n[nonzero])

    # Lanczos sigma window.
    coefficients *= np.sinc(n / (half_width + 1.0))

    # Do not normalize by coefficients.sum(): a band-pass filter has
    # approximately zero response at zero frequency.
    return coefficients


def _filter_series(series, coefficients):
    """Filter one time series and mask invalid endpoints."""
    series = np.ma.asarray(series, dtype=float)
    values = series.filled(np.nan)

    valid = np.isfinite(values)
    numerator = np.convolve(
        np.where(valid, values, 0.0),
        coefficients,
        mode="same",
    )

    # Require every value in the filter window to be valid.
    valid_count = np.convolve(
        valid.astype(int),
        np.ones(coefficients.size, dtype=int),
        mode="same",
    )
    output = np.ma.masked_where(valid_count < coefficients.size, numerator)

    half_width = coefficients.size // 2
    output[:half_width] = np.ma.masked
    output[-half_width:] = np.ma.masked
    return output


def lanczos_bandpass(
    cube,
    low_period=20,
    high_period=100,
    weights=91,
):
    """Apply a Lanczos band-pass filter along the cube time dimension."""
    if not cube.coords("time"):
        msg = "Input cube has no time coordinate"
        raise ValueError(msg)

    time_axis = cube.coord_dims("time")[0]
    coefficients = lanczos_weights(low_period, high_period, weights)

    data = np.moveaxis(np.ma.asarray(cube.data), time_axis, 0)
    flattened = data.reshape(data.shape[0], -1)

    filtered = np.ma.empty(flattened.shape, dtype=float)
    for column in range(flattened.shape[1]):
        filtered[:, column] = _filter_series(
            flattened[:, column],
            coefficients,
        )

    filtered = filtered.reshape(data.shape)
    filtered = np.moveaxis(filtered, 0, time_axis)

    result = cube.copy(data=filtered)
    result.long_name = (
        f"{cube.name()} {low_period}-{high_period} day filtered anomalies"
    )
    result.attributes["temporal_filter"] = (
        f"Lanczos band-pass; periods={low_period}-{high_period} days; "
        f"weights={weights}"
    )
    return result
