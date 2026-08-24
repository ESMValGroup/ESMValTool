"""Lagged regression utilities for MJO diagnostics."""

import iris
import numpy as np


def _paired_regression(field, index):
    """Regress each field column against an index using valid pairs."""
    field = np.ma.asarray(field, dtype=float)
    index = np.ma.asarray(index, dtype=float)

    output = np.ma.masked_all(field.shape[1], dtype=float)

    for longitude_index in range(field.shape[1]):
        y = np.ma.asarray(field[:, longitude_index])
        valid = (
            ~np.ma.getmaskarray(y)
            & ~np.ma.getmaskarray(index)
            & np.isfinite(y.filled(np.nan))
            & np.isfinite(index.filled(np.nan))
        )

        if valid.sum() < 3:
            continue

        x_valid = np.asarray(index[valid], dtype=float)
        y_valid = np.asarray(y[valid], dtype=float)

        # Remove small residual means caused by filtering and lag truncation.
        x_valid -= x_valid.mean()
        y_valid -= y_valid.mean()

        denominator = np.dot(x_valid, x_valid)
        if denominator > 0:
            output[longitude_index] = np.dot(x_valid, y_valid) / denominator

    return output


def lag_regression(cube, index, max_lag=60, *, standardize_index=True):
    """Regress a time-longitude field against an index at multiple lags.

    Negative lag means that the field leads the index. Positive lag means
    that the field follows the index.
    """
    data = np.ma.asarray(cube.data)
    time_axis = cube.coord_dims("time")[0]
    data = np.moveaxis(data, time_axis, 0)

    if data.ndim != 2:
        msg = (
            "Expected a two-dimensional time-longitude cube, "
            f"got shape {data.shape}"
        )
        raise ValueError(msg)

    index = np.ma.asarray(index, dtype=float).squeeze()
    if index.ndim != 1 or index.size != data.shape[0]:
        msg = "The index must be one-dimensional and match the time axis"
        raise ValueError(msg)

    if standardize_index:
        index_mean = np.ma.mean(index)
        index_std = np.ma.std(index)
        if not np.isfinite(index_std) or index_std == 0:
            msg = "The MJO index has zero or invalid variance"
            raise ValueError(msg)
        index = (index - index_mean) / index_std

    lags = np.arange(-max_lag, max_lag + 1)
    regression = np.ma.masked_all((lags.size, data.shape[1]), dtype=float)

    for lag_index, lag in enumerate(lags):
        if lag < 0:
            # Field precedes the reference-index event.
            field_lagged = data[:lag]
            index_lagged = index[-lag:]
        elif lag > 0:
            # Field follows the reference-index event.
            field_lagged = data[lag:]
            index_lagged = index[:-lag]
        else:
            field_lagged = data
            index_lagged = index

        regression[lag_index] = _paired_regression(
            field_lagged,
            index_lagged,
        )

    return regression, lags


def regression_cube(
    regression,
    lags,
    longitude_coord,
    units,
    latitude_coord=None,
):
    """Create an Iris cube containing the lag-regression result."""
    lag_coord = iris.coords.DimCoord(
        np.asarray(lags),
        long_name="lag",
        var_name="lag",
        units="days",
    )

    result = iris.cube.Cube(
        regression,
        long_name="MJO lagged regression",
        var_name="mjo_lag_regression",
        units=units,
        dim_coords_and_dims=[
            (lag_coord, 0),
            (longitude_coord.copy(), 1),
        ],
    )
    result.attributes["lag_convention"] = (
        "Negative lag: field leads index; positive lag: field follows index"
    )
    result.attributes["index_normalization"] = (
        "Reference index standardized to unit standard deviation"
    )
    if latitude_coord is not None:
        result.add_aux_coord(latitude_coord.copy())
    return result
