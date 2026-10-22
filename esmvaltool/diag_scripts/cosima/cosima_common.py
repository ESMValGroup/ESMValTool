"""Shared helpers for the COSIMA-to-CMIP diagnostics.

The COSIMA notebooks lean on three libraries that ESMValTool does not
depend on: ``xgcm`` (staggered-grid differencing), ``xhistogram``
(weighted binning) and ``metpy`` (cross sections).  Intensive fields can
be regridded to a regular latitude/longitude grid, while face-integrated
transports remain on their native grid.  The replacements live here, next to the
TEOS-10 wrappers, the contour extraction that replaces
``Create_contour.ipynb`` and the binning that replaces ``xhistogram``.

Nothing in this module touches the ESMValTool diagnostic interface, so
it can be imported from a plain script or a notebook.
"""

import contextlib
import logging
import warnings

import numpy as np

import iris
import iris.analysis
import iris.analysis.cartography
import iris.coord_systems
import iris.coords
import iris.cube
import iris.util
from esmvalcore.preprocessor import climate_statistics

try:  # gsw is required for every density-space diagnostic
    import gsw
except ImportError:  # pragma: no cover - environments without TEOS-10
    gsw = None

LOGGER = logging.getLogger(__name__)

EARTH_RADIUS = 6371000.0  # m, matches iris.analysis.cartography
RHO0 = 1035.0  # kg m-3, Boussinesq reference density used by MOM5/MOM6
CP_OCEAN = 3992.1  # J kg-1 degC-1, TEOS-10 heat capacity used by MOM6

# ---------------------------------------------------------------------------
# cube basics
# ---------------------------------------------------------------------------


#: Values at or above this are treated as missing.  Covers CMOR's
#: ``missing_value`` of 1e20 and netCDF's default fill of 9.9692e36,
#: while staying well above any physical value in CMOR units (the
#: largest is hfbasin, a few times 1e15 W).
FILL_VALUE_THRESHOLD = 1.0e19


def masked_data(cube, dtype=float):
    """Realised, fully masked data of a cube.

    Always go through this rather than ``cube.core_data()``.  The latter
    returns a *lazy* dask array, on which ``numpy.ma`` operations are
    silently no-ops: the mask is not applied, and fill values flow
    straight into the arithmetic.  This realises the data, applies any
    mask the file carries, and additionally masks non-finite values and
    unconverted fill values, which real CMOR files do contain.
    """
    data = np.ma.masked_invalid(np.ma.asarray(cube.data, dtype=dtype))
    return np.ma.masked_where(
        np.abs(np.ma.filled(data, 0.0)) >= FILL_VALUE_THRESHOLD, data)


def load_cube(filename, short_name):
    """Load a single cube, tolerating missing ``var_name``."""
    cubes = iris.load(filename)
    for cube in cubes:
        if cube.var_name == short_name:
            return cube
    return cubes[0]


def time_mean(cube):
    """Average over time if the cube still carries a time coordinate."""
    if cube.coords("time") and cube.coord("time").shape[0] > 1:
        return climate_statistics(cube, operator="mean", period="full")
    if cube.coords("time"):
        return iris.util.squeeze(cube)
    return cube


def year_range_str(metadata):
    """``start_year-end_year`` spanning a list of metadata dicts."""
    starts = [m["start_year"] for m in metadata if "start_year" in m]
    ends = [m["end_year"] for m in metadata if "end_year" in m]
    if not starts:
        return ""
    return f"{min(starts)}-{max(ends)}"


def _coord_points(cube, name):
    """1D points for ``name``, taking column/row 0 of a 2D aux coord."""
    coord = cube.coord(name)
    if coord.ndim == 1:
        return coord.points
    points = coord.points
    # On a curvilinear grid latitude varies down a column and longitude
    # along a row; take the slice that actually varies.
    if name == "latitude":
        return points[:, 0]
    return points[0, :]


def lat_1d(cube):
    """Representative 1D latitude, also for 2D (curvilinear) coordinates."""
    return _coord_points(cube, "latitude")


def lon_1d(cube):
    """Representative 1D longitude, also for 2D coordinates."""
    return _coord_points(cube, "longitude")


def lat_2d(cube):
    """2D latitude field, broadcasting a 1D coordinate if needed."""
    coord = cube.coord("latitude")
    if coord.ndim == 2:
        return coord.points
    lon = cube.coord("longitude").points
    return np.broadcast_to(coord.points[:, None], (coord.shape[0], lon.size))


def lon_2d(cube):
    """2D longitude field, broadcasting a 1D coordinate if needed."""
    coord = cube.coord("longitude")
    if coord.ndim == 2:
        return coord.points
    lat = cube.coord("latitude").points
    return np.broadcast_to(coord.points[None, :], (lat.size, coord.shape[0]))


def guess_bounds(cube, coords=("latitude", "longitude")):
    """Add bounds in place where they are missing."""
    for name in coords:
        if cube.coords(name) and not cube.coord(name).has_bounds():
            try:
                cube.coord(name).guess_bounds()
            except ValueError:  # scalar or single-point coordinate
                pass
    return cube


def depth_coord_name(cube):
    """Name of the vertical coordinate, or None for a 2D field."""
    for name in ("depth", "ocean_sigma_z", "olevel", "lev"):
        if cube.coords(name):
            return name
    for coord in cube.coords(dim_coords=True):
        if coord.units.is_convertible("m") and coord.name() != "latitude":
            return coord.name()
    return None


def surface_field(cube):
    """Shallowest level of a 3D field; a no-op for 2D fields."""
    name = depth_coord_name(cube)
    if name is None:
        return cube
    depth = cube.coord(name)
    index = int(np.argmin(np.abs(depth.points)))
    constraint = iris.Constraint(**{name: depth.points[index]})
    return cube.extract(constraint)


# ---------------------------------------------------------------------------
# grid metrics (the xgcm replacement)
# ---------------------------------------------------------------------------


def cell_lengths(cube):
    """Zonal and meridional cell lengths in metres.

    Returns ``(dx, dy)`` broadcast to the horizontal shape of ``cube``.
    ``dx`` shrinks with the cosine of latitude, ``dy`` does not.
    """
    lat = lat_1d(cube)
    lon = lon_1d(cube)
    dlon = np.gradient(lon) * np.pi / 180.0
    dlat = np.gradient(lat) * np.pi / 180.0
    coslat = np.cos(np.deg2rad(lat))
    dx = EARTH_RADIUS * coslat[:, None] * dlon[None, :]
    dy = EARTH_RADIUS * np.broadcast_to(dlat[:, None], dx.shape)
    return dx, np.array(dy)


def horizontal_gradient(field, cube, periodic=True):
    """Central-difference gradient of ``field`` in m-1 units.

    ``field`` is a 2D masked array on the horizontal grid of ``cube``.
    Returns ``(d/dx, d/dy)``.  The zonal direction wraps when
    ``periodic`` is set, which is correct for a global CMOR grid.
    """
    dxm, dym = cell_lengths(cube)
    data = np.ma.masked_invalid(np.ma.asarray(field, dtype=float))
    filled = data.filled(np.nan)

    if periodic:
        padded = np.concatenate(
            [filled[:, -1:], filled, filled[:, :1]], axis=1)
        dfdx = (padded[:, 2:] - padded[:, :-2]) / (2.0 * dxm)
    else:
        dfdx = np.full_like(filled, np.nan)
        dfdx[:, 1:-1] = (filled[:, 2:] - filled[:, :-2]) / (2.0 * dxm[:, 1:-1])

    dfdy = np.full_like(filled, np.nan)
    dfdy[1:-1, :] = (filled[2:, :] - filled[:-2, :]) / (2.0 * dym[1:-1, :])

    return np.ma.masked_invalid(dfdx), np.ma.masked_invalid(dfdy)


#: First baroclinic gravity wave speed in m/s.  A canonical
#: mid-latitude value; Chelton et al. (1998) find roughly 2.8 m/s in
#: the tropics falling to ~1 m/s at high latitude, so this is good to
#: well within a factor of two, which is enough for a resolution test.
GRAVITY_WAVE_SPEED = 2.0

#: L_d / grid spacing above which a grid is eddy-permitting, and above
#: which it is eddy-resolving.  Two points per deformation radius lets
#: eddies exist; four begins to represent them.
EDDY_PERMITTING_RATIO = 2.0
EDDY_RESOLVING_RATIO = 4.0


def rossby_radius(lat, wave_speed=GRAVITY_WAVE_SPEED):
    """First baroclinic Rossby deformation radius in metres.

    ``L_d = c / |f|``, the length scale of mesoscale eddies.  It
    diverges at the equator, where the equatorial deformation radius
    ``sqrt(c / 2 beta)`` takes over, so the value is capped there.
    """
    lat = np.asarray(lat, dtype=float)
    fcor = np.abs(coriolis(lat))
    beta = 2.0 * 7.292115e-5 / EARTH_RADIUS
    equatorial = np.sqrt(wave_speed / (2.0 * beta))
    with np.errstate(divide="ignore", invalid="ignore"):
        radius = wave_speed / fcor
    return np.minimum(np.where(np.isfinite(radius), radius, equatorial),
                      equatorial)


def eddy_resolution_ratio(cube, wave_speed=GRAVITY_WAVE_SPEED):
    """Deformation radius divided by grid spacing, per latitude.

    This is the question a vorticity or EKE diagnostic has to ask
    before it means anything: *can this grid support an eddy here?*
    Below :data:`EDDY_PERMITTING_RATIO` it cannot, and the vorticity
    field is the shear of the resolved large-scale currents rather
    than anything mesoscale.

    A 1-degree ocean -- what most CMIP6 models run -- does not reach 2
    anywhere outside the deep tropics.
    """
    lat = lat_1d(cube)
    dxm, dym = cell_lengths(cube)
    # The COARSER direction binds: an eddy has to be resolved in both.
    # Taking the finer one flatters a grid that is anisotropic at high
    # latitude, which a uniform lat/lon grid always is -- dx shrinks
    # with cos(lat) while dy does not.
    spacing = np.maximum(np.abs(dxm).max(axis=1), np.abs(dym).max(axis=1))
    return lat, rossby_radius(lat, wave_speed) / spacing


def relative_vorticity(u_field, v_field, cube):
    """Vertical component of relative vorticity, ``dv/dx - du/dy``.

    This is the regular-grid equivalent of the ``xgcm`` curl in
    ``Relative_Vorticity.ipynb``.  Both velocity components must already
    sit on the same (tracer) points, which is what ``regrid`` delivers.
    """
    dvdx, _ = horizontal_gradient(v_field, cube)
    _, dudy = horizontal_gradient(u_field, cube)
    return dvdx - dudy


def coriolis(lat):
    """Coriolis parameter in s-1."""
    if gsw is not None:
        return np.asarray(gsw.f(np.asarray(lat, dtype=float)))
    omega = 7.292115e-5
    return 2.0 * omega * np.sin(np.deg2rad(np.asarray(lat, dtype=float)))


def gravity(lat, depth=0.0):
    """Gravitational acceleration in m s-2."""
    if gsw is not None:
        return np.asarray(gsw.grav(np.asarray(lat, dtype=float), depth))
    return np.full(np.shape(lat), 9.7963)


def geostrophic_velocities(zos_field, cube, min_lat=5.0):
    """Surface geostrophic velocities from sea level.

    ``u = -(g/f) deta/dy`` and ``v = (g/f) deta/dx``.  Geostrophic
    balance is not valid close to the equator, where ``f`` vanishes.  The
    equatorial band through and including ``min_lat`` degrees either side
    of the equator is therefore masked (5 degrees by default).

    This finite-difference implementation requires the rectilinear grid
    produced by the recipe's ``regrid`` preprocessor.  Applying the
    one-dimensional spherical metrics to a native curvilinear MOM grid
    would silently use one arbitrary column of latitude and give a wrong
    derivative and equatorial mask.
    """
    if (cube.coord("latitude").ndim != 1
            or cube.coord("longitude").ndim != 1):
        raise ValueError(
            "geostrophic_velocities requires a rectilinear latitude/"
            "longitude grid; regrid zos before calculating derivatives")
    detadx, detady = horizontal_gradient(zos_field, cube)
    lat = lat_1d(cube)
    equator = np.abs(lat) <= min_lat
    # Mask f *before* dividing.  Masking the result afterwards still
    # evaluates 1/f on the equator, which can leak infinities or enormous
    # near-equatorial values into a reduction if a backend drops the mask.
    equator = np.broadcast_to(equator[:, None], detadx.shape)
    fcor = np.ma.masked_where(
        equator, np.broadcast_to(coriolis(lat)[:, None], detadx.shape))
    grav = np.broadcast_to(gravity(lat)[:, None], detadx.shape)
    u_g = -(grav / fcor) * detady
    v_g = (grav / fcor) * detadx
    return np.ma.masked_invalid(u_g), np.ma.masked_invalid(v_g)


def time_interval_days(cube, use_bounds=True):
    """Length of each time step in days, in the cube's own calendar.

    Do not convert a time coordinate to a hard-coded reference unit to
    do this.  ``cf_units`` refuses to convert between calendars, so
    ``coord.units.convert(points, "days since 1850-01-01")`` raises on
    any cube whose calendar is not ``standard`` -- and CMIP6 is full of
    them: ACCESS-CM2 is ``proleptic_gregorian``, UKESM1 and HadGEM3 are
    ``360_day``, several others are ``noleap``.

    Converting is also the wrong idea even where it works.  Month
    lengths are calendar-dependent, and they are exactly what a
    flux-weighted average needs: every month is 30 days in a
    ``360_day`` model and February never has 29 in a ``noleap`` one.
    Going through ``num2date`` lets the model's own calendar answer.

    Returns per-step lengths from the coordinate bounds, or the spacing
    between points where there are no bounds.
    """
    coord = cube.coord("time")
    if use_bounds and coord.has_bounds():
        dates = coord.units.num2date(coord.bounds)
        return np.array([(end - start).total_seconds() / 86400.0
                         for start, end in dates])

    dates = np.atleast_1d(coord.units.num2date(coord.points))
    if dates.size < 2:
        return np.array([30.0])
    spacing = np.array([(dates[i + 1] - dates[i]).total_seconds() / 86400.0
                        for i in range(dates.size - 1)])
    return np.concatenate([spacing, spacing[-1:]])


def time_step_days(cube):
    """Median spacing of the time coordinate, in days.

    Returns ``None`` when the cube has no usable time axis.
    """
    if not cube.coords("time") or cube.coord("time").shape[0] < 2:
        return None
    return float(np.median(time_interval_days(cube, use_bounds=False)))


def eddy_kinetic_energy(u_series, v_series):
    """Surface EKE from velocity time series, ``0.5*(u'^2 + v'^2)``."""
    u_anom = u_series - np.ma.mean(u_series, axis=0)
    v_anom = v_series - np.ma.mean(v_series, axis=0)
    return 0.5 * (np.ma.mean(u_anom**2, axis=0) + np.ma.mean(v_anom**2, axis=0))


def horizontal_slice(cube):
    """First slice of a cube over its two horizontal dimensions.

    ``cube.slices(["latitude", "longitude"])`` cannot be used here: on a
    curvilinear grid the two coordinates share both dimensions and iris
    rejects them as non-orthogonal.  Slicing by dimension index works
    for regular and curvilinear grids alike.
    """
    horizontal_dims = sorted(
        set(cube.coord_dims("latitude")) | set(cube.coord_dims("longitude")))
    if cube.ndim == len(horizontal_dims):
        return cube
    index = tuple(slice(None) if dim in horizontal_dims else 0
                  for dim in range(cube.ndim))
    return cube[index]


def spherical_polygon_area(lon_vertices, lat_vertices, radius=EARTH_RADIUS):
    """Area in m2 of spherical polygons given their corner coordinates.

    ``lon_vertices`` and ``lat_vertices`` are in degrees with the corners
    along the last axis, which is the ``(..., 4)`` layout CMOR uses for
    ``vertices_longitude`` / ``vertices_latitude``.  Any number of
    corners works as long as they are ordered around the cell.

    The area is the spherical excess, accumulated over a triangle fan
    from the first corner, using the Van Oosterom & Strackee (1983)
    form:

        tan(E/2) = |v1 . (v2 x v3)|
                   / (1 + v1.v2 + v2.v3 + v3.v1)

    with ``v`` the unit vectors of the corners.  This is *exact* for
    great-circle edges and therefore invariant under rotation of the
    sphere, which matters because a tripolar grid is exactly a rotated
    and stretched one: a formula that is only correct for cells aligned
    with latitude circles would degrade precisely where the grid stops
    being aligned.

    The obvious alternative, the spherical shoelace formula
    ``R^2/2 sum (dlon)(sin lat_i + sin lat_{i+1})``, treats each edge as
    linear in ``(lon, sin lat)`` rather than as a great circle.  It is
    exact for a latitude-longitude box and accurate to second order for
    small cells, but it is not rotation invariant and it fails badly for
    cells spanning large longitude ranges, such as those beside a
    displaced pole.

    The two disagree by about 2e-4 relative on a 2-degree cell, and
    neither is "wrong": they assume different cell edges.  A cell whose
    north and south edges follow *parallels* has the shoelace area,
    which is what ``iris.analysis.cartography.area_weights`` returns; a
    cell with great-circle edges has this one.  Both tile the sphere
    exactly.  The discrepancy is far smaller than the difference
    between either and the model's own ``areacello``, which is why
    ``areacello`` is still preferred when it is published.
    """
    lon = np.deg2rad(np.asarray(lon_vertices, dtype=float))
    lat = np.deg2rad(np.asarray(lat_vertices, dtype=float))
    unit = np.stack([np.cos(lat) * np.cos(lon),
                     np.cos(lat) * np.sin(lon),
                     np.sin(lat)], axis=-1)

    first = unit[..., 0, :]
    excess = np.zeros(unit.shape[:-2])
    for corner in range(1, unit.shape[-2] - 1):
        second = unit[..., corner, :]
        third = unit[..., corner + 1, :]
        triple = np.abs(np.einsum("...i,...i->...",
                                  first, np.cross(second, third)))
        denominator = (1.0
                       + np.einsum("...i,...i->...", first, second)
                       + np.einsum("...i,...i->...", second, third)
                       + np.einsum("...i,...i->...", third, first))
        excess += 2.0 * np.arctan2(triple, denominator)
    return excess * radius**2


def area_from_bounds(cube, radius=EARTH_RADIUS):
    """Cell area in m2 from 2D coordinate bounds, or ``None``.

    This is what makes a tripolar grid workable without ``areacello``.
    CMOR requires bounds on ``latitude`` and ``longitude``, and on a
    curvilinear grid those bounds are the four corner vertices of every
    cell, so the area follows exactly from
    :func:`spherical_polygon_area`.

    ESMValCore does not do this: ``try_adding_calculated_cell_area``
    handles regular and rotated-pole grids and raises
    ``CoordinateMultiDimError`` for everything else.  The information is
    in the file regardless.

    Returns ``None`` when the cube has no usable 2D bounds, so callers
    can fall back.
    """
    horizontal = horizontal_slice(cube)
    lat_coord = horizontal.coord("latitude")
    lon_coord = horizontal.coord("longitude")
    if lat_coord.ndim != 2 or lon_coord.ndim != 2:
        return None
    if lat_coord.bounds is None or lon_coord.bounds is None:
        return None
    if lat_coord.bounds.shape[-1] < 3:
        return None
    return np.ma.masked_invalid(
        spherical_polygon_area(lon_coord.bounds, lat_coord.bounds,
                               radius=radius))


def cell_area(cube):
    """Cell area in m2, preferring the model's own ``areacello``.

    Order of preference:

    1. an ``areacello`` cell measure attached to the cube, which is what
       ESMValCore builds from a ``supplementary_variables`` entry.  This
       is the model's own area and is always the best answer: it
       accounts for the actual grid generation, including any partial
       cells at the coast.
    2. a spherical polygon area from the 2D corner vertices, via
       :func:`area_from_bounds`.  CMOR mandates bounds on latitude and
       longitude, so this works on a tripolar grid with no
       ``areacello`` at all.  It assumes great-circle cell edges on a
       sphere, so it differs from the model's own areas by well under a
       percent for an ocean grid.
    3. spherical areas from 1D latitude/longitude bounds, for a regular
       grid.

    Only if all three fail does this raise, and then with instructions.
    """
    horizontal = horizontal_slice(cube)

    for measure in horizontal.cell_measures():
        if measure.var_name in ("areacello", "areacella", "areacell"):
            return np.ma.masked_invalid(
                np.ma.asarray(
                    horizontal.cell_measure(measure.name()).core_data(),
                    dtype=float).reshape(horizontal.shape))

    from_bounds = area_from_bounds(horizontal)
    if from_bounds is not None:
        return from_bounds

    if horizontal.coord("latitude").ndim > 1:
        raise ValueError(
            "this cube has 2D (curvilinear) latitude/longitude and no cell "
            "corner bounds, so no cell area can be obtained. Either attach "
            "areacello -- in a recipe with\n"
            "    supplementary_variables:\n"
            "      - short_name: areacello\n"
            "        mip: Ofx\n"
            "(add 'exp: piControl' if the experiment you want does not "
            "publish it), in a notebook with "
            "dataset.add_supplementary(short_name='areacello', mip='Ofx') "
            "-- or keep the coordinate bounds, which CMOR requires and "
            "which are enough to compute the area exactly.")

    # iris defaults to an Earth radius of 6367470 m where the cube has
    # no coordinate system, while cell_lengths uses the mean radius of
    # 6371000 m. That is a 0.11% difference in area, small but free to
    # avoid: pin the same sphere so lengths and areas agree.
    horizontal = horizontal.copy()
    guess_bounds(horizontal)
    if horizontal.coord("latitude").coord_system is None:
        sphere = iris.coord_systems.GeogCS(EARTH_RADIUS)
        horizontal.coord("latitude").coord_system = sphere
        horizontal.coord("longitude").coord_system = sphere
    return iris.analysis.cartography.area_weights(horizontal)


# ---------------------------------------------------------------------------
# TEOS-10 wrappers
# ---------------------------------------------------------------------------


@contextlib.contextmanager
def quiet_gsw():
    """Silence gsw's masked-array warning for the enclosed call.

    ``gsw`` emits ``'where' used without 'out', expect uninitialized
    memory in output`` every time it is handed a masked array, which is
    every time here.  The uninitialised values sit only underneath the
    mask and never reach the result -- asserted in
    ``validation/test_cosima_common.py`` -- but the warning fires once
    per call and buries everything else a notebook prints.

    Scoped to this one message so genuine gsw warnings still surface.
    """
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore", message=".*'where' used without 'out'.*",
            category=UserWarning)
        yield


def _require_gsw():
    if gsw is None:
        raise ImportError(
            "gsw (TEOS-10) is required for the density diagnostics; "
            "install it with `pixi add gsw`")


def pressure_from_depth(depth, lat):
    """Sea pressure in dbar from depth in metres (positive down)."""
    _require_gsw()
    depth = np.asarray(depth, dtype=float)
    lat = np.asarray(lat, dtype=float)
    return gsw.p_from_z(-depth, lat)


def absolute_salinity(practical_salinity, pressure, lon, lat):
    """Practical salinity (PSU) to absolute salinity (g/kg)."""
    _require_gsw()
    with quiet_gsw():
        return gsw.SA_from_SP(practical_salinity, pressure, lon, lat)


def conservative_temperature(abs_salinity, potential_temperature):
    """Potential temperature (degC) to conservative temperature (degC)."""
    _require_gsw()
    with quiet_gsw():
        return gsw.CT_from_pt(abs_salinity, potential_temperature)


def _broadcast_to_field(depth, lat, lon, shape):
    """Broadcast depth/lat/lon onto a ``(depth, lat, lon)`` field shape."""
    depth_b = np.broadcast_to(np.asarray(depth)[:, None, None], shape)
    lat_b = np.broadcast_to(np.asarray(lat)[None, :, None], shape)
    lon_b = np.broadcast_to(np.asarray(lon)[None, None, :], shape)
    return depth_b, lat_b, lon_b


def teos10_fields(thetao_cube, so_cube, temperature_is_conservative=False):
    """Absolute salinity, conservative temperature and pressure.

    Accepts fields with optional leading dimensions, such as time:
    ``(lat, lon)``, ``(depth, lat, lon)`` or
    ``(time, depth, lat, lon)``.  It returns masked arrays broadcast to
    the cube shape.  CMIP ``thetao`` is potential temperature, so it is
    converted unless the caller says otherwise (MOM5 ``temp`` is already
    conservative).
    """
    _require_gsw()
    salinity = masked_data(so_cube)
    temperature = masked_data(thetao_cube)

    if salinity.shape != temperature.shape:
        raise ValueError("thetao and so must have the same shape")

    def broadcast_1d(points, axis):
        """Broadcast a coordinate along its corresponding data axis."""
        reshape = [1] * salinity.ndim
        reshape[axis] = len(points)
        return np.broadcast_to(np.asarray(points).reshape(reshape),
                               salinity.shape)

    lat = lat_1d(so_cube)
    lon = lon_1d(so_cube)
    lat_axis = so_cube.coord_dims("latitude")[0]
    lon_axis = so_cube.coord_dims("longitude")[-1]
    depth_name = depth_coord_name(so_cube)
    if depth_name is None:
        depth_b = np.zeros(salinity.shape)
    else:
        depth = so_cube.coord(depth_name).points
        depth_axis = so_cube.coord_dims(depth_name)[0]
        depth_b = broadcast_1d(depth, depth_axis)
    lat_b = broadcast_1d(lat, lat_axis)
    lon_b = broadcast_1d(lon, lon_axis)
    pressure = pressure_from_depth(depth_b, lat_b)
    abs_sal = absolute_salinity(salinity, pressure, lon_b, lat_b)
    if temperature_is_conservative:
        cons_temp = temperature
    else:
        cons_temp = conservative_temperature(abs_sal, temperature)

    return abs_sal, cons_temp, pressure


def potential_density_anomaly(thetao_cube, so_cube, p_ref=2000,
                              temperature_is_conservative=False):
    """Potential density anomaly sigma_n in kg m-3.

    ``p_ref`` is the reference pressure in dbar; 0, 1000, 2000, 3000 and
    4000 map onto the ``gsw.sigma0``...``gsw.sigma4`` family.
    """
    _require_gsw()
    sigma_functions = {
        0: gsw.sigma0,
        1000: gsw.sigma1,
        2000: gsw.sigma2,
        3000: gsw.sigma3,
        4000: gsw.sigma4,
    }
    if p_ref not in sigma_functions:
        raise ValueError(
            f"p_ref must be one of {sorted(sigma_functions)}, got {p_ref}")
    abs_sal, cons_temp, _ = teos10_fields(
        thetao_cube, so_cube,
        temperature_is_conservative=temperature_is_conservative)
    with quiet_gsw():
        return np.ma.masked_invalid(
            sigma_functions[p_ref](abs_sal, cons_temp))


#: Reference pressures and the depth ranges they are patched over, after
#: Reid (1994) and Lynn & Reid (1968).  Used by
#: :func:`patched_potential_density`.
PATCH_LEVELS = (
    (0, 0.0, 500.0),
    (1000, 500.0, 1500.0),
    (2000, 1500.0, 2500.0),
    (3000, 2500.0, 3500.0),
    (4000, 3500.0, np.inf),
)


def patched_potential_density(thetao_cube, so_cube,
                              temperature_is_conservative=False):
    """Continuous patched potential density, an approximation to gamma_n.

    ``Neutral_density.ipynb`` uses ``pygamma_n`` (Jackett & McDougall
    1997), which is not installable from conda-forge.  This is the
    classical alternative: potential density referenced to 0, 1000,
    2000, 3000 and 4000 dbar, each patch offset so the composite is
    continuous across the joins, following Reid (1994).

    It is *not* neutral density.  It tracks neutral surfaces well within
    each patch and is discontinuous in slope at the joins; use it for
    water-mass identification, not for quantitative neutral-surface
    work.  When ``pygamma_n`` is importable the diagnostics use it
    instead.
    """
    _require_gsw()
    abs_sal, cons_temp, _ = teos10_fields(
        thetao_cube, so_cube,
        temperature_is_conservative=temperature_is_conservative)
    depth_name = depth_coord_name(so_cube)
    if depth_name is None:
        return potential_density_anomaly(
            thetao_cube, so_cube, p_ref=0,
            temperature_is_conservative=temperature_is_conservative)

    depth = np.asarray(so_cube.coord(depth_name).points, dtype=float)
    sigma_functions = {
        0: gsw.sigma0, 1000: gsw.sigma1, 2000: gsw.sigma2,
        3000: gsw.sigma3, 4000: gsw.sigma4,
    }

    result = np.ma.masked_all(abs_sal.shape)
    offset = 0.0
    previous = None
    with quiet_gsw():
        for p_ref, top, bottom in PATCH_LEVELS:
            levels = np.where((depth >= top) & (depth < bottom))[0]
            if levels.size == 0:
                continue
            sigma = sigma_functions[p_ref](abs_sal[levels],
                                           cons_temp[levels])
            if previous is not None:
                # Match the new patch to the old one at the shallowest
                # level of the new patch, so there is no jump.
                top_level = levels[0]
                overlap_previous = sigma_functions[previous](
                    abs_sal[top_level], cons_temp[top_level])
                overlap_new = sigma_functions[p_ref](
                    abs_sal[top_level], cons_temp[top_level])
                offset += float(
                    np.ma.median(overlap_previous - overlap_new + offset)
                    - offset)
            result[levels] = sigma + offset
            previous = p_ref
    return np.ma.masked_invalid(result)


def neutral_density(thetao_cube, so_cube, temperature_is_conservative=False):
    """Neutral density gamma_n, via ``pygamma_n`` when it is available.

    Falls back to :func:`patched_potential_density`.  The return value
    is ``(density, method)`` so callers can label their output honestly.
    """
    try:
        import pygamma_n
    except ImportError:
        return (
            patched_potential_density(
                thetao_cube, so_cube,
                temperature_is_conservative=temperature_is_conservative),
            "patched potential density (gamma_n approximation)",
        )

    _require_gsw()
    abs_sal, cons_temp, pressure = teos10_fields(
        thetao_cube, so_cube,
        temperature_is_conservative=temperature_is_conservative)
    with quiet_gsw():
        in_situ = gsw.t_from_CT(abs_sal, cons_temp, pressure)
    lat = lat_1d(so_cube)
    lon = lon_1d(so_cube)

    gamma = np.ma.masked_all(abs_sal.shape)
    for j in range(abs_sal.shape[1]):
        for i in range(abs_sal.shape[2]):
            column_s = abs_sal[:, j, i]
            column_t = in_situ[:, j, i]
            column_p = pressure[:, j, i]
            # pygamma_n's Fortran solver does not treat NaN as a missing
            # observation.  Passing a filled masked array therefore emits a
            # solver warning for partially wet columns and can contaminate
            # their result.  Give it only the physically observed levels;
            # the original vertical pressure values are retained.
            valid = (~np.ma.getmaskarray(column_s)
                     & ~np.ma.getmaskarray(column_t)
                     & ~np.ma.getmaskarray(column_p)
                     & np.isfinite(column_s)
                     & np.isfinite(column_t)
                     & np.isfinite(column_p))
            if np.count_nonzero(valid) < 2:
                continue
            values, _, _ = pygamma_n.gamma_n(
                np.asarray(column_s[valid]), np.asarray(column_t[valid]),
                np.asarray(column_p[valid]), lon[i], lat[j])
            gamma[valid, j, i] = values
    return np.ma.masked_invalid(gamma), "gamma_n (Jackett & McDougall 1997)"


def thermal_expansion(abs_sal, cons_temp, pressure):
    """Thermal expansion coefficient alpha in K-1."""
    _require_gsw()
    with quiet_gsw():
        return gsw.alpha(abs_sal, cons_temp, pressure)


def haline_contraction(abs_sal, cons_temp, pressure):
    """Haline contraction coefficient beta in kg/g."""
    _require_gsw()
    with quiet_gsw():
        return gsw.beta(abs_sal, cons_temp, pressure)


# ---------------------------------------------------------------------------
# binning (the xhistogram replacement)
# ---------------------------------------------------------------------------


def bin_sum(values, bin_variable, bins, axis=None):
    """Sum ``values`` into ``bins`` of ``bin_variable``.

    ``values`` and ``bin_variable`` must have the same shape.  Returns
    an array whose first axis runs over the bin centres and whose
    remaining axes are the axes of ``values`` *not* listed in ``axis``.
    With ``axis=None`` everything is collapsed into a 1D profile.

    This is ``xhistogram.xarray.histogram(bin_variable, bins=bins,
    weights=values, dim=axis)`` without the xarray dependency.
    """
    values = np.ma.masked_invalid(np.ma.asarray(values, dtype=float))
    bin_variable = np.ma.masked_invalid(
        np.ma.asarray(bin_variable, dtype=float))
    bins = np.asarray(bins, dtype=float)
    n_bins = bins.size - 1

    if axis is None:
        axis = tuple(range(values.ndim))
    elif isinstance(axis, int):
        axis = (axis,)
    axis = tuple(a % values.ndim for a in axis)
    keep = tuple(a for a in range(values.ndim) if a not in axis)

    moved = np.ma.transpose(values, keep + axis)
    moved_bin = np.ma.transpose(bin_variable, keep + axis)
    kept_shape = moved.shape[:len(keep)]
    flat = moved.reshape(int(np.prod(kept_shape)) if keep else 1, -1)
    flat_bin = moved_bin.reshape(flat.shape)

    out = np.zeros((n_bins,) + kept_shape)
    flat_out = out.reshape(n_bins, -1)
    # np.digitize returns 0 below the first edge and n_bins+1 above the
    # last, so a bin index of 1..n_bins is inside the range.
    for column in range(flat.shape[0]):
        data = flat[column]
        coordinate = flat_bin[column]
        valid = ~np.ma.getmaskarray(data) & ~np.ma.getmaskarray(coordinate)
        if not valid.any():
            continue
        index = np.digitize(
            np.ma.filled(coordinate, np.nan)[valid], bins)
        inside = (index >= 1) & (index <= n_bins)
        flat_out[:, column] = np.bincount(
            index[inside] - 1,
            weights=np.ma.filled(data, 0.0)[valid][inside],
            minlength=n_bins)
    return out


def bin_centres(bins):
    """Midpoints of a monotonic array of bin edges."""
    bins = np.asarray(bins, dtype=float)
    return 0.5 * (bins[1:] + bins[:-1])


def rebin_to_density(transport, density, bins, depth_axis=0):
    """Rebin a z-coordinate transport onto density bins.

    ``Transformation_from_Depth_to_Potential_Density.ipynb`` does this
    online in MOM; CMIP6 publishes only z-coordinate ``vmo``, so the
    rebinning is offline and uses monthly-mean density.  That
    underestimates the eddy contribution to the transport, which is why
    the diagnostic reports it as a separate method rather than as an
    equivalent of the online diagnostic.
    """
    return bin_sum(transport, density, bins, axis=depth_axis)


#: Bin widths, matching the COSIMA original's
#: ``np.arange(floor(min) - 0.5, ceil(max) + 0.6, 0.5)`` for temperature
#: and the equivalent 0.1 for salinity.  The *ranges* are derived from
#: the data, as the original derives them, rather than fixed.
#:
#: Width matters more here than it looks.  A time-averaged section has
#: only as many distinct T and S values as it has grid cells, and the
#: vertical levels set that: a 50-level column gives ~50 distinct
#: salinities. Bins finer than that spacing leave empty bins *between*
#: populated ones, which the eye reads as structure but is aliasing.
DEFAULT_TEMPERATURE_BIN_WIDTH = 0.5
DEFAULT_SALINITY_BIN_WIDTH = 0.1


def ts_bin_edges(values, width):
    """Bin edges spanning the data, at the given width.

    Reproduces the original's construction: floor the minimum, ceil the
    maximum, pad by one bin either side.
    """
    finite = np.ma.compressed(np.ma.masked_invalid(values))
    if finite.size == 0:
        return np.array([0.0, width])
    low = np.floor(finite.min()) - width
    high = np.ceil(finite.max()) + width
    return np.arange(low, high + width, width)


def report_aliasing(histogram):
    """Warn when the bins are finer than the data can populate.

    Counts empty bins that sit *between* populated ones on the same
    temperature row.  Those gaps are not water that is absent, they are
    bins narrower than the spacing of the available values, and on the
    figure they appear as white streaks through the middle of the
    distribution.
    """
    filled = histogram > 0
    gaps = 0
    for row in filled:
        populated = np.where(row)[0]
        if populated.size > 1:
            gaps += (populated[-1] - populated[0] + 1) - populated.size
    if filled.sum() and gaps > 0.25 * filled.sum():
        LOGGER.warning(
            "%d of %d populated bins have an empty bin beside them within "
            "the occupied range; the bins are finer than the data can fill "
            "and the figure will show white streaks. Widen "
            "temperature_bin_width / salinity_bin_width, or keep the time "
            "axis so more values contribute.", gaps, int(filled.sum()))
    return gaps


def volume_weighted_ts_histogram(cons_temp, abs_sal, volume,
                                 temperature_bins, salinity_bins):
    """Volume-weighted joint T-S histogram.

    Returns an array of shape ``(len(temperature_bins) - 1,
    len(salinity_bins) - 1)`` holding the total volume in each cell.
    """
    cons_temp = np.ma.masked_invalid(np.ma.asarray(cons_temp, dtype=float))
    abs_sal = np.ma.masked_invalid(np.ma.asarray(abs_sal, dtype=float))
    volume = np.ma.masked_invalid(np.ma.asarray(volume, dtype=float))
    volume = np.ma.masked_array(
        np.broadcast_to(volume, cons_temp.shape),
        mask=np.broadcast_to(np.ma.getmaskarray(volume), cons_temp.shape))

    valid = (~np.ma.getmaskarray(cons_temp) & ~np.ma.getmaskarray(abs_sal)
             & ~np.ma.getmaskarray(volume))
    if not valid.any():
        return np.zeros((len(temperature_bins) - 1, len(salinity_bins) - 1))
    hist, _, _ = np.histogram2d(
        np.ma.filled(cons_temp, np.nan)[valid],
        np.ma.filled(abs_sal, np.nan)[valid],
        bins=[np.asarray(temperature_bins), np.asarray(salinity_bins)],
        weights=np.ma.filled(volume, 0.0)[valid])
    return hist


#: Coordinate names that carry a latitude, in order of preference.
#: ``msftyz`` is the reason the list is not just ``("latitude",)``: its
#: CMOR dimension is ``gridlatitude`` (iris ``grid_latitude``) because
#: the streamfunction is defined on the model's own y index rather than
#: on true latitude.  ``msftmz`` and ``hfbasin`` use ``latitude``.
LATITUDE_COORD_NAMES = ("latitude", "grid_latitude")


def latitude_coord(cube):
    """The cube's latitude coordinate, whatever CMOR called it.

    Raises with the names that *are* present, since a missing latitude
    is nearly always a variable whose CMOR dimension differs rather
    than genuinely unlabelled data.
    """
    for name in LATITUDE_COORD_NAMES:
        if cube.coords(name):
            return cube.coord(name)
    for coord in cube.coords():
        if str(coord.units) in ("degrees_north", "degrees_N", "degree_north"):
            return coord
    raise ValueError(
        f"no latitude-like coordinate on this cube; it has "
        f"{sorted(c.name() for c in cube.coords())}")


def row_latitude(cube):
    """Representative latitude of each grid row.

    The honest label for a quantity summed along the ``i`` index of a
    curvilinear grid, such as an overturning streamfunction.
    :func:`lat_1d` takes one arbitrary column instead, which south of a
    tripolar fold is close to the same thing and north of it is not.

    Reduces to the latitude coordinate itself on a regular grid.
    """
    coord = latitude_coord(cube)
    if coord.ndim == 1:
        return coord.points
    return np.ma.mean(np.ma.masked_invalid(coord.points), axis=1)


def latitude_bin_edges(cube, n_bins=None):
    """Bin edges spanning the cube's *whole* latitude field.

    Do not build these from :func:`lat_1d`.  On a CMORised tripolar grid
    there is no 1D nominal latitude to fall back on -- ``latitude`` is a
    2D auxiliary coordinate -- and ``lat_1d`` returns one arbitrary
    column of it.  North of the tripolar fold the grid rows stop
    following latitude circles, so that column's maximum depends on
    which meridian it happens to sit on: for a displaced-pole Arctic it
    can stop dead at the fold (~65N) while the field itself reaches
    89.5N.  Edges built from it silently discard every cell beyond,
    which is the entire Arctic.

    The upstream COSIMA notebook does not hit this because MOM5 output
    carries ``yt_ocean``, a genuine 1D nominal latitude spanning the
    full grid, alongside the 2D ``geolat_t``.  CMOR keeps only the 2D
    field, so the edges have to be constructed instead of borrowed.

    Returns ``n_bins + 1`` edges from -90 to 90, defaulting to one bin
    per grid row.
    """
    lat_field = lat_2d(cube)
    if n_bins is None:
        n_bins = lat_field.shape[0]
    return np.linspace(-90.0, 90.0, int(n_bins) + 1)


def latitude_bin_sum(values, lat_field, lat_bins):
    """Sum ``values`` into latitude bins of a possibly curvilinear grid.

    This is the ``groupby_bins('geolat_t', ...)`` of
    ``Meridional_heat_transport.ipynb``.  On a regular grid it reduces
    to a zonal sum.
    """
    lat_field = np.broadcast_to(np.asarray(lat_field), np.shape(values))
    return bin_sum(values, lat_field, lat_bins)


# ---------------------------------------------------------------------------
# contours (the Create_contour.ipynb replacement)
# ---------------------------------------------------------------------------


def extract_contour(field, level, longest_only=True):
    """Contour of ``field`` at ``level``, in index space.

    ``field`` is a 2D array, typically ``deptho``.  Returns ``(x, y)``
    fractional index coordinates.  As in ``Create_contour.ipynb`` only
    the longest path is kept by default, which discards contours around
    isolated seamounts and depressions.
    """
    import matplotlib

    matplotlib.use("Agg", force=False)
    import matplotlib.pyplot as plt

    data = np.ma.masked_invalid(np.ma.asarray(field, dtype=float))
    figure = plt.figure()
    try:
        axes = figure.add_subplot(111)
        contour_set = axes.contour(np.ma.filled(data, 0.0), levels=[level])
        paths = contour_set.get_paths()
        segments = []
        for path in paths:
            vertices = path.vertices
            if vertices.shape[0] < 2:
                continue
            # A single Path may hold several disjoint polylines; split
            # them where the step between consecutive vertices jumps.
            steps = np.abs(np.diff(vertices[:, 0]).astype(int))
            breaks = np.where(steps > 1)[0]
            start = 0
            for end in list(breaks) + [vertices.shape[0] - 1]:
                if end + 1 - start >= 2:
                    segments.append(vertices[start:end + 1])
                start = end + 1
    finally:
        plt.close(figure)

    if not segments:
        return np.array([]), np.array([])
    if longest_only:
        chosen = max(segments, key=lambda seg: seg.shape[0])
        return chosen[:, 0], chosen[:, 1]
    stacked = np.concatenate(segments, axis=0)
    return stacked[:, 0], stacked[:, 1]


def contour_cells(x_contour, y_contour, shape):
    """Snap a contour to unique, adjacently-connected grid cells.

    Returns ``(rows, cols)`` integer index arrays ordered along the
    contour.  Diagonal steps are filled in with an intermediate cell,
    which ``Cross-contour_transport.ipynb`` requires for the transport
    masks to be exact.
    """
    rows = np.clip(np.rint(y_contour).astype(int), 0, shape[0] - 1)
    cols = np.clip(np.rint(x_contour).astype(int), 0, shape[1] - 1)

    out_rows, out_cols = [], []
    for row, col in zip(rows, cols):
        if out_rows and out_rows[-1] == row and out_cols[-1] == col:
            continue
        if out_rows:
            drow = row - out_rows[-1]
            dcol = col - out_cols[-1]
            if abs(drow) == 1 and abs(dcol) == 1:
                # Insert the cell that makes the step adjacent rather
                # than diagonal; step in x first, as the notebook does.
                out_rows.append(out_rows[-1])
                out_cols.append(col)
        out_rows.append(row)
        out_cols.append(col)
    return np.array(out_rows), np.array(out_cols)


def cross_contour_masks(rows, cols, shape):
    """Masks selecting the faces crossed by a contour, with sign.

    Returns ``(mask_u, mask_v)``, each of shape ``shape``, holding +1 or
    -1 on the cells whose eastward (``mask_u``) or northward
    (``mask_v``) face is crossed by the contour, and 0 elsewhere.  The
    sign is set so that the summed transport is positive *outward*, away
    from the shallow side of the contour.

    Multiplying the mask by ``umo``/``vmo`` and summing gives the
    cross-contour transport of ``Cross-contour_transport.ipynb``.
    """
    rows = np.asarray(rows, dtype=int)
    cols = np.asarray(cols, dtype=int)
    if (rows.size != cols.size or not rows.size
            or np.any(rows < 0) or np.any(rows >= shape[0])
            or np.any(cols < 0) or np.any(cols >= shape[1])):
        raise ValueError(
            "contour indices do not fit the transport grid "
            f"{tuple(shape)}. deptho and umo/vmo must cover the same "
            "horizontal region; do not clip or interpolate mass transport.")

    mask_u = np.zeros(shape)
    mask_v = np.zeros(shape)
    for index in range(len(rows) - 1):
        row, col = rows[index], cols[index]
        drow = rows[index + 1] - row
        dcol = cols[index + 1] - col
        # Walking the contour with the shallow side on the left, a step
        # in +x crosses the northward face and a step in +y crosses the
        # eastward face.
        if dcol == 1:
            mask_v[row, col] += 1
        elif dcol == -1:
            mask_v[row, col] -= 1
        elif drow == 1:
            mask_u[row, col] -= 1
        elif drow == -1:
            mask_u[row, col] += 1
    return mask_u, mask_v


def common_face_shape(*shapes):
    """Overlap on which native staggered transport faces can be indexed.

    MOM's ``umo`` and ``vmo`` live on different faces.  After a geographic
    ``extract_region`` their latitude coordinates can select one fewer edge
    row from one component: ACCESS-CM2, for example, yields ``(45, 360)``
    and ``(44, 360)``.  That is a normal C-grid edge, not a reason to
    interpolate either extensive transport.

    The common index region is the shorter latitude dimension.  It is only
    valid for this one-edge-row stagger; different zonal dimensions or a
    larger mismatch are evidence that the fields do not cover the same
    region, and must fail rather than silently discard transports.
    """
    horizontal = [tuple(shape[-2:]) for shape in shapes]
    if not horizontal or any(len(shape) != 2 or min(shape) <= 0
                             for shape in horizontal):
        raise ValueError("transport fields need non-empty horizontal grids")
    rows = [shape[0] for shape in horizontal]
    columns = [shape[1] for shape in horizontal]
    if len(set(columns)) != 1:
        raise ValueError(
            f"transport grids have different zonal dimensions: {horizontal}")
    if max(rows) - min(rows) > 1:
        raise ValueError(
            "transport grids differ by more than one latitude edge row: "
            f"{horizontal}; check that extract_region selected the same "
            "native domain")
    return min(rows), columns[0]


def cross_contour_faces(rows, cols, shape):
    """Ordered signed faces crossed by a contour.

    Returns ``(rows, cols, components, signs)`` in contour order, where
    ``components`` is ``"u"`` or ``"v"``.  Unlike the 2-D masks this keeps
    the path ordering needed for the original recipe's cumulative
    along-contour transport plot.
    """
    # Reuse the mask constructor's strict grid validation.  A silent clip
    # would turn a mass-transport integral into a different contour.
    cross_contour_masks(rows, cols, shape)
    face_rows, face_cols, components, signs = [], [], [], []
    for index in range(len(rows) - 1):
        row, col = rows[index], cols[index]
        drow = rows[index + 1] - row
        dcol = cols[index + 1] - col
        if dcol == 1:
            component, sign = "v", 1.0
        elif dcol == -1:
            component, sign = "v", -1.0
        elif drow == 1:
            component, sign = "u", -1.0
        elif drow == -1:
            component, sign = "u", 1.0
        else:
            continue
        face_rows.append(row)
        face_cols.append(col)
        components.append(component)
        signs.append(sign)
    return (np.asarray(face_rows, dtype=int), np.asarray(face_cols, dtype=int),
            np.asarray(components), np.asarray(signs, dtype=float))


def describe_contour(rows, cols, mask_u, mask_v, shape):
    """Classify a contour and apply the right closure check.

    There are two kinds of contour here and they close differently.

    A *loop* -- around a seamount or a basin -- returns to its start,
    so every face it crosses one way it crosses back and both mask
    sums are zero.

    A *circumpolar* contour, which is what the 1000 m isobath around
    Antarctica is, does not close in index space at all: it wraps
    through the periodic zonal boundary. Walking it eastward crosses
    one northward face per column, so ``mask_v`` sums to roughly the
    number of columns rather than to zero. That is correct, and a
    check demanding zero would report the Antarctic case as broken --
    which is exactly what the first Gadi run did, printing
    "v-face: +359 (0 means closed)".

    Returns ``(kind, message)``.
    """
    n_columns = shape[-1]
    spans_globe = (cols.max() - cols.min()) >= 0.9 * (n_columns - 1)
    v_sum, u_sum = float(mask_v.sum()), float(mask_u.sum())

    if spans_globe and abs(abs(v_sum) - n_columns) <= 0.1 * n_columns:
        return "circumpolar", (
            f"circumpolar contour: {rows.size} cells spanning the full "
            f"zonal range, {v_sum:+.0f} net northward-face crossings for "
            f"{n_columns} columns -- one per column, as it should be")
    if not spans_globe and abs(v_sum) < 2 and abs(u_sum) < 2:
        return "loop", (
            f"closed loop: {rows.size} cells, face crossings balance "
            f"({u_sum:+.0f} zonal, {v_sum:+.0f} meridional)")
    return "open", (
        f"contour does not close: {rows.size} cells, {u_sum:+.0f} zonal "
        f"and {v_sum:+.0f} meridional net crossings over {n_columns} "
        f"columns. The transport across it is not a well-defined "
        f"exchange; check contour_depth and the region.")


def contour_indices_to_coordinates(rows, cols, lat, lon):
    """Latitude/longitude of integer contour cells."""
    lat = np.asarray(lat)
    lon = np.asarray(lon)
    if lat.ndim == 2:
        return lat[rows, cols], lon[rows, cols]
    return lat[rows], lon[cols]


def wrap_longitude(values, reference_lon):
    """Express longitudes in the same convention as ``reference_lon``.

    CMIP6 grids come both as 0-360 and as -180-180, and a section
    longitude written for one is silently wrong on the other: a naive
    ``interp`` clamps everything out of range onto the first column.
    Mapping into ``[lon[0], lon[0] + 360)`` makes the request
    convention-independent.
    """
    reference_lon = np.asarray(reference_lon, dtype=float)
    origin = float(reference_lon[0])
    return origin + (np.asarray(values, dtype=float) - origin) % 360.0


def sample_section(cube, start, end, n_points=60):
    """Bilinearly sample a 3D cube along a straight lon/lat path.

    ``Cross-slope_section.ipynb`` uses
    ``metpy.interpolate.cross_section``; metpy is not an ESMValTool
    dependency, so this does the same job with the grid's own indices.
    The zonal seam wraps, and the longitudes are normalised with
    :func:`wrap_longitude` so the section lands where it was asked for
    whichever convention the grid uses.

    Returns ``(section, lons, lats)`` with ``section`` shaped
    ``(depth, n_points)``.

    The data is realised once, before the per-point loop.  Indexing a
    lazy cube inside the loop makes every point its own dask task and
    ships the whole field to the scheduler again each time, which on
    Gadi produced repeated "Sending large graph of size 33 MiB"
    warnings and dominated the runtime.
    """
    lon = lon_1d(cube)
    lat = lat_1d(cube)
    lons = np.linspace(start[0], end[0], n_points)
    lats = np.linspace(start[1], end[1], n_points)

    data = masked_data(cube)
    if data.ndim == 4:
        data = np.ma.mean(data, axis=0)
    # Realised, contiguous, and never touched lazily again below.
    data = np.ma.masked_array(np.ascontiguousarray(np.ma.filled(data, np.nan)),
                              mask=np.ma.getmaskarray(data))
    data = np.ma.masked_invalid(data)

    col = np.interp(wrap_longitude(lons, lon), lon, np.arange(lon.size))
    row = np.interp(lats, lat, np.arange(lat.size))

    section = np.ma.masked_all((data.shape[0], n_points))
    for point in range(n_points):
        j0 = int(np.floor(row[point]))
        i0 = int(np.floor(col[point]))
        j1 = min(j0 + 1, lat.size - 1)
        i1 = (i0 + 1) % lon.size  # the zonal direction wraps
        weight_j = row[point] - j0
        weight_i = col[point] - i0
        corners = [
            ((1 - weight_j) * (1 - weight_i), data[:, j0, i0]),
            ((1 - weight_j) * weight_i, data[:, j0, i1]),
            (weight_j * (1 - weight_i), data[:, j1, i0]),
            (weight_j * weight_i, data[:, j1, i1]),
        ]
        accumulated = np.ma.zeros(data.shape[0])
        weights = np.zeros(data.shape[0])
        for weight, values in corners:
            valid = ~np.ma.getmaskarray(values)
            accumulated[valid] += weight * values[valid]
            weights[valid] += weight
        section[:, point] = np.ma.masked_where(
            weights == 0, accumulated / np.where(weights == 0, 1.0, weights))
    return section, lons, lats


def bathymetry_on_grid(deptho_cube, reference_cube):
    """Bathymetry on the same horizontal grid as ``reference_cube``.

    ``deptho`` lives on the tracer grid while ``umo``/``vmo`` live on
    cell faces.  On a native (unregridded) grid those are genuinely
    different arrays: ``extract_region`` selects by latitude *value*,
    and because the staggered latitudes sit half a cell away it can
    return one row fewer for the transports than for ``deptho``.  A
    contour traced on the tracer grid then cannot index the transports
    at all -- shapes ``(50, 45, 360)`` against ``(50, 44, 360)``, which
    is exactly how this failed on ACCESS-CM2.

    Where the shapes already agree, ``deptho`` is used as it is.  Where
    they do not, the bathymetry comes from the reference field's own
    deepest wet cell, which is on the right grid by construction.  That
    is the same substitution ``Along_Isobath_Average.ipynb`` makes for
    its own reasons, and it costs the continuous depth: the contour is
    then quantised to the model's vertical levels.

    Returns ``(bathymetry, source)``.
    """
    reference_shape = reference_cube.shape[-2:]
    if deptho_cube is not None:
        depth = masked_data(deptho_cube)
        if depth.shape == reference_shape:
            return depth, "deptho"
        LOGGER.warning(
            "deptho is %s but the field it must index is %s -- these are "
            "staggered grids that extract_region has trimmed differently. "
            "Falling back to the bottom depth implied by the field itself, "
            "which is quantised to the model levels. To use deptho, put "
            "both on a common grid -- but never regrid a transport.",
            depth.shape, reference_shape)

    depth_name = depth_coord_name(reference_cube)
    if depth_name is None:
        raise ValueError(
            "no deptho on the right grid and no depth axis to derive one")
    levels = reference_cube.coord(depth_name).points
    data = masked_data(reference_cube)
    while data.ndim > 3:
        data = data[0]
    marker = np.ma.masked_array(
        np.broadcast_to(levels[:, None, None], data.shape).astype(float),
        mask=np.ma.getmaskarray(data))
    valid = ~np.ma.getmaskarray(marker)
    any_valid = valid.any(axis=0)
    deepest = valid.shape[0] - 1 - np.argmax(valid[::-1], axis=0)
    bathymetry = np.take_along_axis(
        np.ma.filled(marker, np.nan), deepest[None, ...], axis=0)[0]
    return (np.ma.masked_array(bathymetry, mask=~any_valid),
            "deepest wet level of the transport field")


def report_section_coverage(section, label="section"):
    """Log how much of a sampled section actually has water.

    A section that misses the ocean -- wrong endpoints, wrong longitude
    convention, a path over land -- produces an all-masked array, and a
    blank figure looks much like a failed one.  Saying so up front is
    the difference between a bug and a result.
    """
    total = section.size
    valid = int(np.ma.count(section))
    if valid == 0:
        LOGGER.warning(
            "the %s is entirely masked: no ocean along that path. Check "
            "the endpoints and that the section longitude is inside the "
            "extracted region.", label)
    elif valid < 0.1 * total:
        LOGGER.warning(
            "only %d of %d %s points have data (%.0f%%); the path is "
            "mostly over land or outside the region",
            valid, total, label, 100 * valid / total)
    else:
        LOGGER.info("%s has data at %d of %d points (%.0f%%)",
                    label, valid, total, 100 * valid / total)
    return valid


def isobath_mask(deptho, minimum, maximum):
    """Boolean mask of cells whose depth lies in ``[minimum, maximum)``."""
    depth = np.ma.masked_invalid(np.ma.asarray(deptho, dtype=float))
    return (depth >= minimum) & (depth < maximum)


def topographic_unit_vectors(deptho, cube):
    """Unit vectors tangent and normal to the bathymetry contours.

    Returns ``(t_x, t_y, n_x, n_y)``.  The normal points towards deeper
    water; the tangent is the normal rotated so that shallow water lies
    to the left, matching the sketch in
    ``Along-slope-velocities.ipynb``:

        u_along = u * h_y/|grad h| - v * h_x/|grad h|
        v_cross = u * h_x/|grad h| + v * h_y/|grad h|
    """
    dhdx, dhdy = horizontal_gradient(deptho, cube)
    magnitude = np.ma.sqrt(dhdx**2 + dhdy**2)
    magnitude = np.ma.masked_where(magnitude == 0, magnitude)
    n_x = dhdx / magnitude
    n_y = dhdy / magnitude
    return n_y, -n_x, n_x, n_y


def project_along_slope(u_field, v_field, deptho, cube):
    """Along- and cross-slope velocity components."""
    t_x, t_y, n_x, n_y = topographic_unit_vectors(deptho, cube)
    along = u_field * t_x + v_field * t_y
    cross = u_field * n_x + v_field * n_y
    return along, cross


def nearest_distance(mask, cube):
    """Great-circle distance in metres from each cell to a True cell.

    The ``Nearest_Neighbour_Distance.ipynb`` utility; used for distance
    to the sea-ice edge.  Cells that are themselves True get 0.
    """
    mask = np.asarray(mask, dtype=bool)
    lat = np.deg2rad(lat_2d(cube))
    lon = np.deg2rad(lon_2d(cube))
    targets = np.where(mask)
    if targets[0].size == 0:
        return np.ma.masked_all(mask.shape)

    target_lat = lat[targets]
    target_lon = lon[targets]
    distance = np.empty(mask.shape)
    for j in range(mask.shape[0]):
        # Haversine from every cell of row j to every target at once.
        dlat = target_lat[None, :] - lat[j][:, None]
        dlon = target_lon[None, :] - lon[j][:, None]
        haversine = (np.sin(dlat / 2.0)**2
                     + np.cos(lat[j])[:, None] * np.cos(target_lat)[None, :]
                     * np.sin(dlon / 2.0)**2)
        distance[j] = EARTH_RADIUS * 2.0 * np.arcsin(
            np.sqrt(np.clip(haversine, 0.0, 1.0))).min(axis=1)
    return distance


# ---------------------------------------------------------------------------
# vertical structure
# ---------------------------------------------------------------------------


def bottom_cell(cube):
    """Value in the deepest unmasked cell of every water column.

    ``Extract_Variables_at_Ocean_Bottom.ipynb`` uses MOM's ``kmt``;
    CMIP6 has no equivalent, so the bottom cell is found from the mask
    of the field itself.
    """
    name = depth_coord_name(cube)
    if name is None:
        return masked_data(cube)
    axis = cube.coord_dims(name)[0]
    data = masked_data(cube)
    data = np.ma.masked_array(np.moveaxis(data, axis, 0),
                              mask=np.moveaxis(
                                  np.ma.getmaskarray(data), axis, 0))
    valid = ~np.ma.getmaskarray(data)
    any_valid = valid.any(axis=0)
    # The deepest valid level is the last True along the depth axis.
    deepest = valid.shape[0] - 1 - np.argmax(valid[::-1], axis=0)
    result = np.take_along_axis(
        np.ma.filled(data, np.nan), deepest[None, ...], axis=0)[0]
    return np.ma.masked_array(result, mask=~any_valid)


def cumulative_from_bottom(values, axis=0):
    """Cumulative sum from the bottom up along ``axis``.

    ``numpy`` sums from index 0 down, so the bottom-up cumulative sum is
    the total minus the top-down sum, exactly as
    ``Overturning_Circulation.ipynb`` writes it.
    """
    values = np.ma.masked_invalid(np.ma.asarray(values, dtype=float))
    filled = np.ma.filled(values, 0.0)
    return np.cumsum(filled, axis=axis) - np.sum(
        filled, axis=axis, keepdims=True)


def check_transport_grid(cube, name="transport"):
    """Warn if an extensive transport looks like it has been regridded.

    ``umo`` and ``vmo`` are *extensive*: kg s-1 through one cell face,
    not a flux per unit area.  Interpolating them and then summing
    scales the total with the number of cells, so a 1x1 to 0.5x0.5
    regrid inflates a transport by roughly four.  Nothing in the
    interpolation complains -- the numbers stay plausible -- so this
    checks the one thing that is visible after the fact: whether the
    cube still sits on the model's own (curvilinear) grid.

    A regular lat/lon grid is not proof of regridding, since a few
    models publish ocean data on one, hence a warning rather than an
    error.
    """
    if cube.coord("latitude").ndim == 1:
        LOGGER.warning(
            "%s is on a regular lat/lon grid. If it was regridded, its "
            "total is wrong: %s is kg s-1 per cell face, and interpolation "
            "does not rescale it for the new cell size. Compute transports "
            "on the native grid, or convert to a flux per unit area before "
            "regridding and multiply by the new face area afterwards.",
            name, name)
        return False
    return True


def overturning_streamfunction(vmo, axis_depth=0, axis_lon=-1, rho0=RHO0):
    """Meridional overturning streamfunction in Sv from ``vmo``.

    ``vmo`` is northward *mass* transport in kg s-1.  The streamfunction
    is the zonal sum, converted to volume transport and cumulatively
    integrated from the sea floor upwards.
    """
    transport = np.ma.masked_invalid(np.ma.asarray(vmo, dtype=float))
    zonal = np.ma.sum(transport, axis=axis_lon)
    zonal = (zonal / rho0) / 1.0e6
    depth_axis = axis_depth if axis_depth < zonal.ndim else 0
    return cumulative_from_bottom(zonal, axis=depth_axis)


# ---------------------------------------------------------------------------
# sea ice seasonality (Massom et al. 2013)
# ---------------------------------------------------------------------------


def seaice_seasonality(concentration, threshold=0.15, consecutive_days=5):
    """Day of advance, day of retreat and season duration.

    ``concentration`` is ``(time, lat, lon)`` daily sea-ice
    concentration as a *fraction*, covering one sea-ice year (15
    February to 14 February).  Following Massom et al. (2013):

    - advance is the first day of the first run of ``consecutive_days``
      days at or above ``threshold``.  Note that the COSIMA notebook
      instead reports the day the run *reaches* five, which is four
      days later; Massom's definition is the first day, and that is
      what this returns;
    - retreat is the last day at or above ``threshold``;
    - duration is the difference.

    Cells that never reach the threshold are masked; cells that never
    drop below it get day 1 and the full year.
    """
    data = np.ma.masked_invalid(np.ma.asarray(concentration, dtype=float))
    above = np.ma.filled(data >= threshold, False)
    n_days = above.shape[0]
    days_above = above.sum(axis=0)

    no_ice = days_above == 0
    no_advance = days_above < consecutive_days
    always_ice = days_above == n_days

    # Length of the current run of "above threshold" days, reset to 0
    # whenever the concentration drops below the threshold.
    run = np.zeros(above.shape, dtype=int)
    counter = np.zeros(above.shape[1:], dtype=int)
    for day in range(n_days):
        counter = np.where(above[day], counter + 1, 0)
        run[day] = counter

    reached = run == consecutive_days
    advance = np.argmax(reached, axis=0).astype(float)
    # argmax returns 0 when nothing is True; those cells are masked by
    # no_advance below.  The advance day is the *start* of the run.
    advance = advance - (consecutive_days - 1) + 1  # 1-based day of year

    reversed_above = above[::-1]
    last_above = n_days - np.argmax(reversed_above, axis=0)
    retreat = last_above.astype(float)

    advance = np.ma.masked_array(advance, mask=no_advance)
    retreat = np.ma.masked_array(retreat, mask=no_ice)
    advance[always_ice] = 1.0
    retreat[always_ice] = float(n_days)
    duration = retreat - advance
    return advance, retreat, duration


#: Concentration below which a cell is not considered ice covered.
#: The 15% cutoff is the NSIDC convention that both the area and the
#: extent definitions are built on.
ICE_THRESHOLD = 0.15


def seaice_area(concentration, area, threshold=ICE_THRESHOLD,
                weighted=True):
    """Total sea-ice area or extent in m2, NSIDC conventions.

    Both quantities count only cells at or above ``threshold``:

    - *area* (``weighted=True``) is ``sum(concentration * cell_area)``;
    - *extent* (``weighted=False``) is ``sum(cell_area)``.

    Applying the threshold to the area as well as to the extent is what
    makes extent >= area, which is the relation these two numbers are
    always quoted with.  Weighting every cell regardless of
    concentration, as a naive reading of "area" suggests, breaks that
    relation as soon as the marginal ice zone is wide.

    Pass ``threshold=None`` to weight every cell, for the rare case
    where that is genuinely wanted.
    """
    conc = np.ma.masked_invalid(np.ma.asarray(concentration, dtype=float))
    area = np.ma.masked_invalid(np.ma.asarray(area, dtype=float))
    if threshold is not None:
        conc = np.ma.masked_where(conc < threshold, conc)
    if not weighted:
        return float(np.ma.sum(np.ma.masked_where(
            np.ma.getmaskarray(conc), area)))
    return float(np.ma.sum(conc * area))


def fraction_from_percent(cube):
    """Sea-ice concentration as a fraction, whatever the input units.

    CMOR ``siconc`` is a percentage but several models publish a
    fraction.  The units are trusted when they are meaningful and the
    data range is used only as a fallback -- masking the fill values
    first, which is the bug that made the Russell sea-ice figures come
    out empty.
    """
    data = masked_data(cube)
    data = np.ma.masked_greater(data, 1.0e15)
    units = str(cube.units).strip().lower()
    if units in ("%", "percent"):
        return data / 100.0
    if units in ("1", "", "unknown", "none"):
        return data / 100.0 if np.ma.max(data) > 5.0 else data
    return data / 100.0 if np.ma.max(data) > 5.0 else data


# ---------------------------------------------------------------------------
# cube construction and plotting
# ---------------------------------------------------------------------------


def make_cube(data, coords, var_name, units, long_name):
    """Build a cube from ``(points, standard_name, units)`` tuples."""
    cube = iris.cube.Cube(data, var_name=var_name, units=units,
                          long_name=long_name)
    for dim, (points, name, coord_units) in enumerate(coords):
        coord = iris.coords.DimCoord(
            np.asarray(points, dtype=float), units=coord_units)
        try:
            coord.standard_name = name
        except ValueError:
            coord.long_name = name
        cube.add_dim_coord(coord, dim)
    return cube


def south_polar_axes(figure, position, max_lat=-30.0):
    """Circular South Polar Stereographic axes."""
    import cartopy.crs as ccrs
    import matplotlib.path as mpath

    axes = figure.add_subplot(*position, projection=ccrs.SouthPolarStereo())
    axes.set_extent([-180, 180, -90, max_lat], crs=ccrs.PlateCarree())
    theta = np.linspace(0, 2 * np.pi, 100)
    verts = np.vstack([np.sin(theta), np.cos(theta)]).T
    axes.set_boundary(mpath.Path(verts * 0.5 + [0.5, 0.5]),
                      transform=axes.transAxes)
    return axes


def style_for(dataset, index=0):
    """Deterministic colour/linestyle for a dataset name."""
    colours = ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e",
               "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf"]
    dashes = ["-", "--", "-.", ":"]
    key = index if index else abs(hash(dataset)) % len(colours)
    return {"color": colours[key % len(colours)],
            "dash": dashes[(key // len(colours)) % len(dashes)],
            "thick": 1.8}
