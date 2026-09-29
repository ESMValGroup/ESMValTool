"""Shared helpers for the Python port of the russell18jgr diagnostics.

Russell, J.L., et al., 2018, J. Geophysical Research - Oceans, 123,
3120-3143. https://doi.org/10.1002/2017JC013461

These helpers replicate the behaviour of the original NCL scripts in
esmvaltool/diag_scripts/russell18jgr/ using numpy/iris/matplotlib only.
"""

import logging

import iris
import numpy as np

logger = logging.getLogger(__name__)

# Constants used throughout the original NCL scripts
EARTH_RADIUS = 6.37e06  # m (value used in the NCL scripts)
DEG2RAD = 0.0174533  # value used in the NCL scripts


def load_cube(filename, short_name=None):
    """Load a single cube from an ESMValTool preprocessed file."""
    if short_name is not None:
        try:
            return iris.load_cube(
                filename, iris.NameConstraint(var_name=short_name)
            )
        except Exception:  # noqa: BLE001 - fall through to plain load
            pass
    return iris.load_cube(filename)


def coord_1d_or_2d(cube, axis):
    """Return the latitude/longitude points array (1D or 2D).

    axis is "latitude" or "longitude".  ESMValTool preprocessed ocean
    files carry either 1D dimension coordinates or 2D auxiliary
    coordinates; this hides the difference (the NCL scripts did the same
    by reading lat/lon variables from the input file).
    """
    coord = cube.coord(axis)
    return np.asarray(coord.points)


def has_regular_grid(cube):
    """True if the cube has 1D lat/lon coordinates."""
    return (
        cube.coord("latitude").ndim == 1 and cube.coord("longitude").ndim == 1
    )


def lat_1d(cube):
    """A representative 1D latitude array.

    For regular grids this is the latitude coordinate itself; for
    curvilinear grids the first column of the 2D array is used (the NCL
    scripts used column 0 or 10 of the 2D lat array interchangeably).
    """
    lat = coord_1d_or_2d(cube, "latitude")
    if lat.ndim == 2:
        return lat[:, 0]
    return lat


def lon_1d(cube):
    """A representative 1D longitude array (row 0 for 2D grids)."""
    lon = coord_1d_or_2d(cube, "longitude")
    if lon.ndim == 2:
        return lon[0, :]
    return lon


def time_mean(cube):
    """Average over the time dimension (if present)."""
    if any(c.name() == "time" for c in cube.coords()):
        if cube.coord("time").shape[0] > 1:
            return cube.collapsed("time", iris.analysis.MEAN)
        return cube[0] if cube.coord_dims("time") else cube
    return cube


def monthly_climatology(cube):
    """12-month climatology of a monthly time series (like clmMonTLL)."""
    import iris.coord_categorisation

    if not any(c.name() == "month_number" for c in cube.coords()):
        iris.coord_categorisation.add_month_number(cube, "time")
    clim = cube.aggregated_by("month_number", iris.analysis.MEAN)
    # Sort January..December
    order = np.argsort(clim.coord("month_number").points)
    return clim[order]


def surface_field(cube):
    """The surface (shallowest) level, if the cube has a depth axis.

    ``ph`` is a full-depth field in CMIP5 (time, lev, lat, lon) while
    russell18jgr figure 8 is the *surface* pH, so the shallowest level
    is selected.  2D fields such as tauu and fgco2 are returned
    unchanged.
    """
    zcoords = [
        c
        for c in cube.coords(dim_coords=True)
        if c.attributes.get("positive") in ("down", "up")
        or c.name() in ("depth", "lev", "olevel")
    ]
    if not zcoords:
        return cube
    zcoord = zcoords[0]
    index = [slice(None)] * cube.ndim
    index[cube.coord_dims(zcoord)[0]] = int(np.argmin(np.abs(zcoord.points)))
    return cube[tuple(index)]


def sea_ice_percent(cube):
    """Sea ice concentration in percent, with fill values masked.

    Some models publish a fraction (0-1) rather than a percentage, and
    some store missing data as 1e20 rather than as a masked/NaN value.
    The NCL scripts decided between the two with ``max(sic) < 5``,
    which a single 1e20 fill value defeats - the fraction is then never
    scaled, nothing exceeds the 15% threshold and the plot comes out
    empty.  Fill values are therefore masked before the units are
    decided, and the cube units are trusted when they are meaningful.
    """
    data = np.ma.masked_invalid(np.ma.asarray(cube.data, dtype=float))
    data = np.ma.masked_greater(data, 1.0e3)
    unit = str(cube.units).strip().lower()
    if "%" in unit or "percent" in unit:
        is_fraction = False
    elif unit in ("1", "1.0", "", "unknown", "none", "fraction"):
        is_fraction = True
    else:
        is_fraction = data.count() > 0 and data.max() < 5.0
    if is_fraction:
        data = data * 100.0
    return np.ma.masked_outside(data, 0.0, 100.0)


def zonal_mean(masked_data):
    """Mean over the last (longitude/i) axis honouring the mask."""
    data = np.ma.masked_invalid(masked_data)
    return np.ma.mean(data, axis=-1)


def style_for(dataset, styleset="CMIP5"):
    """Return the plotting style dict for a dataset.

    Keys: color, dash, thick, mark (matplotlib-compatible).  ``styleset``
    is the recipe's ``styleset`` option and selects the style file
    (CMIP5, CMIP6, ...).  The ESMValTool style file is used when
    available (i.e. when running as a diagnostic); notebooks that only
    have ESMValCore installed fall back to the matplotlib colour cycle.
    """
    try:
        from esmvaltool.diag_scripts.shared.plot import get_dataset_style

        return get_dataset_style(dataset, str(styleset).lower())
    except Exception:  # noqa: BLE001 - be permissive about style lookup
        logger.debug("No ESMValTool style for %s, using defaults", dataset)
        return {
            "color": None,
            "dash": "-",
            "thick": 1.5,
            "mark": "o",
            "avgstd": 0,
            "facecolor": "none",
        }


def closest_index(value, array):
    """Index of the array element closest to value (NCL closest_val)."""
    array = np.asarray(array, dtype=float)
    return int(np.nanargmin(np.abs(array - value)))


def ncl_cell_area(lat, lon, shape):
    """Cell areas for a regular lat-lon grid, as computed in NCL.

    Replicates the manual areacello calculation of the NCL scripts
    (russell18jgr-fig5g / fig7i / fig9b / fig9c): constant dlat/dlon
    taken from grid points 19 and 20, area = dx * dy with
    dx = R * cos(lat) * dlon and dy = R * dlat.

    Returns a 2D array broadcast to `shape` (lat, lon).
    """
    lat = np.asarray(lat, dtype=float)
    lon = np.asarray(lon, dtype=float)
    dlat = abs(lat[20] - lat[19])
    dlon = abs(lon[20] - lon[19])
    clat = np.cos(lat * DEG2RAD)
    dx = dlon * EARTH_RADIUS * DEG2RAD * clat
    dy = dlat * EARTH_RADIUS * DEG2RAD
    dxdy = dx * dy  # (nlat,)
    return np.broadcast_to(dxdy[:, np.newaxis], shape)


def depth_to_dbar(depth):
    """Approximate pressure (decibars) from depth (m), as in NCL/POP."""
    depth = np.asarray(depth, dtype=float)
    bars = (
        0.059808 * (np.exp(-0.025 * depth) - 1.0)
        + 0.100766 * depth
        + 2.28405e-7 * depth**2
    )
    return 10.0 * bars


def rho_mwjf(theta, salt, depth):
    """Potential density from McDougall, Wright, Jackett & Feistel (2003).

    Port of NCL's ``rho_mwjf`` (as used by the POP ocean model).

    Parameters
    ----------
    theta : array
        Potential temperature in degrees Celsius.
    salt : array
        Salinity (psu).
    depth : float
        Reference depth in metres (0, 1977 and 3948 m give reference
        pressures of approximately 0, 2000 and 4000 dbar).

    Returns
    -------
    array
        Density in g/cm^3 (i.e. ~1.0265); the NCL scripts convert this
        to sigma with ``1000 * (rho - 1)``.
    """
    t = np.ma.masked_invalid(np.ma.asarray(theta, dtype=float))
    s = np.ma.masked_invalid(np.ma.asarray(salt, dtype=float))
    s = np.ma.where(s < 0.0, 0.0, s)
    p = float(depth_to_dbar(depth))

    p001 = 0.001  # scales the numerator so that rho is in g/cm^3

    # Numerator coefficients
    mwjfnp0s0t0 = 9.99843699e2 * p001
    mwjfnp0s0t1 = 7.35212840e0 * p001
    mwjfnp0s0t2 = -5.45928211e-2 * p001
    mwjfnp0s0t3 = 3.98476704e-4 * p001
    mwjfnp0s1t0 = 2.96938239e0 * p001
    mwjfnp0s1t1 = -7.23268813e-3 * p001
    mwjfnp0s2t0 = 2.12382341e-3 * p001
    mwjfnp1s0t0 = 1.04004591e-2 * p001
    mwjfnp1s0t2 = 1.03970529e-7 * p001
    mwjfnp1s1t0 = 5.18761880e-6 * p001
    mwjfnp2s0t0 = -3.24041825e-8 * p001
    mwjfnp2s0t2 = -1.23869360e-11 * p001

    # Denominator coefficients
    mwjfdp0s0t0 = 1.0e0
    mwjfdp0s0t1 = 7.28606739e-3
    mwjfdp0s0t2 = -4.60835542e-5
    mwjfdp0s0t3 = 3.68390573e-7
    mwjfdp0s0t4 = 1.80809186e-10
    mwjfdp0s1t0 = 2.14691708e-3
    mwjfdp0s1t1 = -9.27062484e-6
    mwjfdp0s1t3 = -1.78343643e-10
    mwjfdp0sqt0 = 4.76534122e-6
    mwjfdp0sqt1 = 1.63410736e-9
    mwjfdp1s0t0 = 5.30848875e-6
    mwjfdp2s0t3 = -3.03175128e-16
    mwjfdp3s0t1 = -1.27934137e-17

    # Pressure-dependent combined coefficients
    nums0t0 = mwjfnp0s0t0 + p * (mwjfnp1s0t0 + p * mwjfnp2s0t0)
    nums0t1 = mwjfnp0s0t1
    nums0t2 = mwjfnp0s0t2 + p * (mwjfnp1s0t2 + p * mwjfnp2s0t2)
    nums0t3 = mwjfnp0s0t3
    nums1t0 = mwjfnp0s1t0 + p * mwjfnp1s1t0
    nums1t1 = mwjfnp0s1t1
    nums2t0 = mwjfnp0s2t0

    dens0t0 = mwjfdp0s0t0 + p * mwjfdp1s0t0
    dens0t1 = mwjfdp0s0t1 + (p**3) * mwjfdp3s0t1
    dens0t2 = mwjfdp0s0t2
    dens0t3 = mwjfdp0s0t3 + (p**2) * mwjfdp2s0t3
    dens0t4 = mwjfdp0s0t4
    dens1t0 = mwjfdp0s1t0
    dens1t1 = mwjfdp0s1t1
    dens1t3 = mwjfdp0s1t3
    densqt0 = mwjfdp0sqt0
    densqt1 = mwjfdp0sqt1

    numerator = (
        nums0t0
        + t * (nums0t1 + t * (nums0t2 + nums0t3 * t))
        + s * (nums1t0 + nums1t1 * t + nums2t0 * s)
    )
    denominator = (
        dens0t0
        + t * (dens0t1 + t * (dens0t2 + t * (dens0t3 + dens0t4 * t)))
        + s * (dens1t0 + dens1t1 * t + dens1t3 * t**3)
        + s * np.ma.sqrt(s) * (densqt0 + densqt1 * t**2)
    )
    return numerator / denominator


def westerly_band_width(tauu_zonal, lat):
    """Latitudinal width of the Southern Hemisphere westerly wind band.

    Replicates the NCL fig9a/9b logic: starting from the latitude
    closest to 50S, find the first zero crossing of the zonal-mean
    zonal wind stress equatorward (+ to -) and the last zero crossing
    poleward of 50S (- to +), both linearly interpolated.
    """
    tauu = np.ma.masked_invalid(np.asarray(tauu_zonal, dtype=float))
    lat = np.asarray(lat, dtype=float)
    a1 = closest_index(-50.0, lat)
    a2 = closest_index(-75.0, lat)

    final_lat = None
    for i in range(a1, len(lat) // 2):
        if tauu[i] >= 0 and tauu[i + 1] < 0:
            lat1, lat2 = lat[i], lat[i + 1]
            final_lat = lat1 - (
                (lat1 - lat2) * tauu[i] / (tauu[i] - tauu[i + 1])
            )
            break
    lower_lat = None
    for i in range(a2, a1 + 1):
        if tauu[i] < 0 and tauu[i + 1] >= 0:
            lat1, lat2 = lat[i], lat[i + 1]
            lower_lat = lat1 - (
                (lat1 - lat2) * tauu[i] / (tauu[i] - tauu[i + 1])
            )
    if final_lat is None or lower_lat is None:
        raise ValueError(
            "Could not find the zero crossings of the zonal wind stress "
            "around 50S needed to compute the westerly band width."
        )
    return float(final_lat - lower_lat)


def southern_ocean_flux_sum(flux2d, areacello, lat, factor):
    """Integrated flux south of 30S: sum(flux * area) * factor.

    Replicates the NCL fig9 scripts: flux per cell is summed along
    longitude and then summed from the pole up to the latitude row
    closest to 30S (assumes south-to-north ordering, which ESMValTool
    preprocessed CMIP data has).
    """
    flux = np.ma.masked_invalid(np.asarray(flux2d, dtype=float))
    area = np.ma.masked_invalid(np.asarray(areacello, dtype=float))
    per_cell = flux * area * factor
    per_lat = per_cell.sum(axis=-1)
    a = closest_index(-30.0, lat)
    return float(per_lat[: a + 1].sum())


def extract_isoline(lon2d, lat2d, field, level):
    """Extract the (lon, lat) points of an isoline of a 2D field.

    Replicates the NCL get_isolines-based front extraction of
    russell18jgr-fig3b(-2).ncl: contour the field against its lat/lon
    coordinates, take the longest contour segment as the main front and
    append any further segments whose start is within 20 degrees
    longitude of the current end (stray points are discarded).

    Returns (x, y): longitudes and latitudes of the isoline, ordered
    with increasing longitude.
    """
    import matplotlib.pyplot as plt

    fig = plt.figure()
    try:
        contours = plt.contour(
            lon2d, lat2d, np.ma.masked_invalid(field), [level]
        )
        segments = [seg for seg in contours.allsegs[0] if len(seg) > 1]
    finally:
        plt.close(fig)
    if not segments:
        raise ValueError(f"No isoline found at level {level}")
    segments.sort(key=len, reverse=True)
    main = segments[0]
    x = list(main[:, 0])
    y = list(main[:, 1])
    for seg in segments[1:]:
        if abs(seg[0, 0] - x[-1]) < 20.0:
            x.extend(seg[:, 0])
            y.extend(seg[:, 1])
    x = np.asarray(x)
    y = np.asarray(y)
    if x[min(11, len(x) - 1)] < x[min(1, len(x) - 1)]:
        x = x[::-1]
        y = y[::-1]
    return x, y


def replicate_isoline(x, y):
    """Repeat the isoline three times to span -360..360 degrees."""
    x_all = np.concatenate([x - 360.0, x, x + 360.0])
    y_all = np.concatenate([y, y, y])
    return x_all, y_all


def to_kelvin(data):
    """Make sure temperatures are in Kelvin (as done in the NCL code).

    Values of exactly 0 are masked when the field is already in Kelvin
    (some models use 0 as a fill value).
    """
    data = np.ma.masked_invalid(data)
    if data.max() > 273.0:
        return np.ma.masked_values(data, 0.0) if data.min() == 0 else data
    return data + 273.0


def provenance_record(
    caption,
    ancestors,
    plot_types=("geo",),
    statistics=("mean",),
    domains=("sh",),
):
    """Standard provenance record for the russell18jgr diagnostics."""
    return {
        "caption": caption,
        "statistics": list(statistics),
        "domains": list(domains),
        "plot_types": list(plot_types),
        "authors": ["russell_joellen", "pandde_amarjiit"],
        "references": ["russell18jgr"],
        "ancestors": list(ancestors),
    }


def select_by_short_name(metadata, short_name):
    """List of metadata dicts with the given short_name."""
    return [m for m in metadata if m["short_name"] == short_name]


def match_by_dataset(metadata, dataset):
    """First metadata dict for the given dataset name, or None."""
    for meta in metadata:
        if meta["dataset"] == dataset:
            return meta
    return None


def year_range_str(metas):
    """'YYYY-YYYY' over all datasets (used in output file names)."""
    start = min(int(m["start_year"]) for m in metas)
    end = max(int(m["end_year"]) for m in metas)
    return f"{start:04d}-{end:04d}"
