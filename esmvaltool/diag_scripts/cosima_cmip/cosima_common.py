"""Shared array, geometry, and cube helpers for ocean hydrography.

This module contains only the helpers used by the hydrographic benchmark and
its density-compensation companion. Its functions are copied from the COSIMA
CMIP diagnostic collection so both diagnostics use identical area weights.
"""

import iris
import iris.analysis.cartography
import iris.coord_systems
import iris.coords
import iris.cube
import numpy as np

EARTH_RADIUS = 6371000.0
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
        np.abs(np.ma.filled(data, 0.0)) >= FILL_VALUE_THRESHOLD, data
    )


def load_cube(filename, short_name):
    """Load a single cube, tolerating missing ``var_name``."""
    cubes = iris.load(filename)
    for cube in cubes:
        if cube.var_name == short_name:
            return cube
    return cubes[0]


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


def horizontal_slice(cube):
    """First slice of a cube over its two horizontal dimensions.

    ``cube.slices(["latitude", "longitude"])`` cannot be used here: on a
    curvilinear grid the two coordinates share both dimensions and iris
    rejects them as non-orthogonal.  Slicing by dimension index works
    for regular and curvilinear grids alike.
    """
    horizontal_dims = sorted(
        set(cube.coord_dims("latitude")) | set(cube.coord_dims("longitude"))
    )
    if cube.ndim == len(horizontal_dims):
        return cube
    index = tuple(
        slice(None) if dim in horizontal_dims else 0
        for dim in range(cube.ndim)
    )
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
    unit = np.stack(
        [np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)],
        axis=-1,
    )

    first = unit[..., 0, :]
    excess = np.zeros(unit.shape[:-2])
    for corner in range(1, unit.shape[-2] - 1):
        second = unit[..., corner, :]
        third = unit[..., corner + 1, :]
        triple = np.abs(
            np.einsum("...i,...i->...", first, np.cross(second, third))
        )
        denominator = (
            1.0
            + np.einsum("...i,...i->...", first, second)
            + np.einsum("...i,...i->...", second, third)
            + np.einsum("...i,...i->...", third, first)
        )
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
        spherical_polygon_area(
            lon_coord.bounds, lat_coord.bounds, radius=radius
        )
    )


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
                    dtype=float,
                ).reshape(horizontal.shape)
            )

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
            "which are enough to compute the area exactly."
        )

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


def make_cube(data, coords, var_name, units, long_name):
    """Build a cube from ``(points, standard_name, units)`` tuples."""
    cube = iris.cube.Cube(
        data, var_name=var_name, units=units, long_name=long_name
    )
    for dim, (points, name, coord_units) in enumerate(coords):
        coord = iris.coords.DimCoord(
            np.asarray(points, dtype=float), units=coord_units
        )
        try:
            coord.standard_name = name
        except ValueError:
            coord.long_name = name
        cube.add_dim_coord(coord, dim)
    return cube
