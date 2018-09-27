# -*- coding: utf-8 -*-
"""Provides regridding for irregular grids"""

import ESMF
import iris
import numpy as np
from iris.exceptions import CoordinateNotFoundError

from ._mapping import get_empty_data, map_slices, ref_to_dims_index

ESMF_MANAGER = ESMF.Manager(debug=False)

ESMF_LON, ESMF_LAT = 0, 1

ESMF_REGRID_METHODS = {
    'linear': ESMF.RegridMethod.BILINEAR,
    'area_weighted': ESMF.RegridMethod.CONSERVE,
    'nearest': ESMF.RegridMethod.NEAREST_STOD,
}

# ESMF_REGRID_METHODS = {
#     'bilinear': ESMF.RegridMethod.BILINEAR,
#     'patch': ESMF.RegridMethod.PATCH,
#     'conserve': ESMF.RegridMethod.CONSERVE,
#     'nearest_stod': ESMF.RegridMethod.NEAREST_STOD,
#     'nearest_dtos': ESMF.RegridMethod.NEAREST_DTOS,
# }


def coords_iris_to_esmpy(lat, lon, circular):
    """Build ESMF compatible coordinate information from iris coords"""
    dim = len(lat.shape)
    if dim == 1:
        for coord in [lat, lon]:
            if not coord.has_bounds():
                coord.guess_bounds()
        esmpy_lat, esmpy_lon = np.meshgrid(lat.points, lon.points)
        lat_corners = np.concatenate([lat.bounds[:, 0], lat.bounds[-1:, 1]])
        if circular:
            lon_corners = lon.bounds[:, 0]
        else:
            lon_corners = np.concatenate([lon.bounds[:, 0],
                                          lon.bounds[-1:, 1]])
        esmpy_lat_corners, esmpy_lon_corners = np.meshgrid(lat_corners,
                                                           lon_corners)
    elif dim == 2:
        esmpy_lat, esmpy_lon = lat.points.T.copy(), lon.points.T.copy()
        lat_corners = np.concatenate([lat.bounds[:, :, 0],
                                      lat.bounds[-1:, :, 1]])
        lon_corners = np.concatenate([lon.bounds[:, :, 0],
                                      lon.bounds[-1:, :, 1]])
        esmpy_lat_corners = lat_corners.T
        esmpy_lon_corners = lon_corners.T
    else:
        raise RuntimeError('Coord dimension is {}. Expected 1 or 2.'
                           ''.format(dim))
    return esmpy_lat, esmpy_lon, esmpy_lat_corners, esmpy_lon_corners


def get_grid(esmpy_lat, esmpy_lon,
             esmpy_lat_corners, esmpy_lon_corners, circular):
    """Build EMSF grid from given coordinate information"""
    if circular:
        num_peri_dims = 1
    else:
        num_peri_dims = 0
    grid = ESMF.Grid(np.array(esmpy_lat.shape),
                     num_peri_dims=num_peri_dims,
                     staggerloc=[ESMF.StaggerLoc.CENTER])
    grid.get_coords(ESMF_LON)[...] = esmpy_lon
    grid.get_coords(ESMF_LAT)[...] = esmpy_lat
    grid.add_coords([ESMF.StaggerLoc.CORNER])
    grid_lon_corners = grid.get_coords(ESMF_LON,
                                       staggerloc=ESMF.StaggerLoc.CORNER)
    grid_lat_corners = grid.get_coords(ESMF_LAT,
                                       staggerloc=ESMF.StaggerLoc.CORNER)
    grid_lon_corners[...] = esmpy_lon_corners
    grid_lat_corners[...] = esmpy_lat_corners
    if grid.mask[0] is None:
        grid.add_item(ESMF.GridItem.MASK, ESMF.StaggerLoc.CENTER)
    return grid


def get_empty_field(cube, grid, remove_mask=False):
    """Build empty ESMF field from cube with given grid"""
    field = ESMF.Field(grid,
                       name=cube.long_name,
                       staggerloc=ESMF.StaggerLoc.CENTER)
    center_mask = grid.get_item(ESMF.GridItem.MASK, ESMF.StaggerLoc.CENTER)
    if np.ma.isMaskedArray(cube.data):
        field.data[...] = cube.data.data.T
        if remove_mask:
            center_mask[...] = 0
        else:
            center_mask[...] = cube.data.mask.T
    else:
        field.data[...] = cube.data.T
        center_mask[...] = 0
    return field


def is_lon_circular(lat, lon):
    """Determine if longitudes are circular"""
    if isinstance(lon, iris.coords.DimCoord):
        circular = lon.circular
    elif isinstance(lon, iris.coords.AuxCoord):
        if len(lon.shape) == 1:
            seam = lon.bounds[-1, 1] - lon.bounds[0, 0]
        elif len(lon.shape) == 2:
            n_to_s = lat.points[0, 0] > lat.points[-1, 0]
            if n_to_s:
                seam = (lon.bounds[1:-1, -1, (0, 1)]
                        - lon.bounds[1:-1, 0, (3, 2)])
            else:
                seam = (lon.bounds[1:-1, -1, (1, 2)]
                        - lon.bounds[1:-1, 0, (0, 3)])
        else:
            raise RuntimeError('AuxCoord longitude is higher dimensional'
                               'than 2d. Giving up.')
        circular = np.alltrue(abs(seam) % 360. < 1.e-3)
    else:
        raise RuntimeError('longitude is neither DimCoord nor AuxCoord.'
                           'Giving up.')
    return circular


def cube_to_empty_field(cube, circular_lon=None, remove_mask=False):
    """Build an empty ESMF field from a cube"""
    lat = cube.coord('latitude')
    lon = cube.coord('longitude')
    if circular_lon is None:
        circular = is_lon_circular(lat, lon)
    else:
        circular = circular_lon
    esmpy_coords = coords_iris_to_esmpy(lat, lon, circular)
    grid = get_grid(*esmpy_coords, circular=circular)
    field = get_empty_field(cube, grid, remove_mask)
    return field


def get_representant(cube, ref_to_slice):
    """Get a representative slice from a cube"""
    slice_dims = ref_to_dims_index(cube, ref_to_slice)
    rep_ind = [0] * cube.ndim
    for dim in slice_dims:
        rep_ind[dim] = slice(None, None)
    rep_ind = tuple(rep_ind)
    return cube[rep_ind]


def build_regridder_2d(src_rep, dst_rep, regrid_method, mask_threshold=.0):
    """Build regridder for 2d regridding"""
    dst_field = cube_to_empty_field(dst_rep)
    src_field = cube_to_empty_field(src_rep)
    regridding_arguments = {
        'srcfield': src_field,
        'dstfield': dst_field,
        'regrid_method': regrid_method,
        'unmapped_action': ESMF.UnmappedAction.IGNORE,
        'ignore_degenerate': True,
    }
    if np.ma.is_masked(src_rep.data):
        mask_regridder = ESMF.Regrid(src_mask_values=np.array([]),
                                     **regridding_arguments)
        src_field.data[...] = src_rep.data.mask.T
        regr_field = mask_regridder(src_field, dst_field)
        dst_mask = regr_field.data[...].T > mask_threshold
        center_mask = dst_field.grid.get_item(ESMF.GridItem.MASK,
                                              ESMF.StaggerLoc.CENTER)
        center_mask[...] = dst_mask.T
    else:
        dst_mask = False
    field_regridder = ESMF.Regrid(src_mask_values=np.array([1]),
                                  dst_mask_values=np.array([1]),
                                  **regridding_arguments)

    def regridder(src):
        """Regrid 2d for irregular grids"""
        res = get_empty_data(dst_rep.shape)
        data = src.data
        if np.ma.is_masked(data):
            data = data.data
        src_field.data[...] = data.T
        regr_field = field_regridder(src_field, dst_field)
        res.data[...] = regr_field.data[...].T
        res.mask[...] = dst_mask
        return res

    return regridder


def build_regridder_3d(src_rep, dst_rep, regrid_method, mask_threshold=.0):
    # pylint: disable=too-many-locals
    # The necessary refactoring will be done for the full 3d regridding.
    """Build regridder for 2.5d regridding"""
    dst_field = cube_to_empty_field(dst_rep[0], remove_mask=True)
    src_fields = []
    esmf_regridders = []
    dst_masks = []
    no_levels = src_rep.shape[0]
    regridding_arguments = {
        'dstfield': dst_field,
        'regrid_method': regrid_method,
        'unmapped_action': ESMF.UnmappedAction.IGNORE,
        'ignore_degenerate': True,
    }
    for level in range(no_levels):
        src_field = cube_to_empty_field(src_rep[level])
        src_fields.append(src_field)
        center_mask = dst_field.grid.get_item(ESMF.GridItem.MASK,
                                              ESMF.StaggerLoc.CENTER)
        if np.ma.is_masked(src_rep.data):
            mask_regridder = ESMF.Regrid(srcfield=src_field,
                                         src_mask_values=np.array([]),
                                         **regridding_arguments)
            src_field.data[...] = src_rep[level].data.mask.T
            regr_field = mask_regridder(src_field, dst_field)
            dst_mask = regr_field.data[...].T > mask_threshold
            dst_masks.append(dst_mask)
            center_mask[...] = dst_mask.T
        else:
            dst_masks.append(False)
            center_mask[...] = 0
        esmf_regridders.append(
            ESMF.Regrid(srcfield=src_field,
                        src_mask_values=np.array([1]),
                        dst_mask_values=np.array([1]),
                        **regridding_arguments)
        )

    def regridder(src):
        """Regrid 2.5d for irregular grids"""
        res = get_empty_data(dst_rep.shape)
        for i, (src_field, esmf_regridder, dst_mask) \
                in enumerate(zip(src_fields, esmf_regridders, dst_masks)):
            data = src[i].data
            if np.ma.is_masked(data):
                data = data.data
            src_field.data[...] = data.T
            regr_field = esmf_regridder(src_field, dst_field)
            res.data[i, ...] = regr_field.data[...].T
            res.mask[i, ...] = dst_mask
        return res

    return regridder


def build_regridder(src_rep, dst_rep, method, mask_threshold=.0):
    """Build regridders from representants"""
    regrid_method = ESMF_REGRID_METHODS[method]
    if src_rep.ndim == 2:
        regridder = build_regridder_2d(src_rep, dst_rep,
                                       regrid_method, mask_threshold)
    elif src_rep.ndim == 3:
        regridder = build_regridder_3d(src_rep, dst_rep,
                                       regrid_method, mask_threshold)
    return regridder


def get_grid_representant(cube, horizontal_only=False):
    """Extract the spatial grid from a cube"""
    horizontal_slice = ['latitude', 'longitude']
    if horizontal_only:
        ref_to_slice = horizontal_slice
    else:
        try:
            cube_z_coord = cube.coord(axis='Z')
            ref_to_slice = [cube_z_coord] + horizontal_slice
        except CoordinateNotFoundError:
            ref_to_slice = horizontal_slice
    return get_representant(cube, ref_to_slice)


def get_grid_representants(src, dst):
    """
    Construct cubes representing the source and destination grid

    This method constructs two new cubes that representant the grids,
    i.e. the spatial dimensions of the given cubes.

    Parameters
    ----------
    src: :class:`iris.cube.Cube`
        Cube to be regridded. Typically a time series of 2d or 3d slices.
    dst: :class:`iris.cube.Cube`
        Cube defining the destination grid. Usually just a 2d or 3d cube.

    Returns
    -------
    tuple of :class:`iris.cube.Cube`:
        A tuple containing two cubes, representing the source grid and the
        destination grid, respectively.
    """
    src_rep = get_grid_representant(src)
    dst_horiz_rep = get_grid_representant(dst, horizontal_only=True)
    if src_rep.ndim == 3:
        dst_shape = (src_rep.shape[0],)
        dim_coords = [src_rep.coord(dimensions=[0], dim_coords=True)]
    else:
        dst_shape = tuple()
        dim_coords = []
    dst_shape += dst_horiz_rep.shape
    dim_coords += dst_horiz_rep.coords(dim_coords=True)
    dim_coords_and_dims = [(c, i) for i, c in enumerate(dim_coords)]
    dst_rep = iris.cube.Cube(
        data=get_empty_data(dst_shape),
        standard_name=src.standard_name,
        long_name=src.long_name,
        var_name=src.var_name,
        units=src.units,
        attributes=src.attributes,
        cell_methods=src.cell_methods,
        dim_coords_and_dims=dim_coords_and_dims,
    )
    return src_rep, dst_rep


def regrid(src, dst, method='linear'):
    """
    Regrid src_cube to the grid defined by dst_cube.

    Regrid the data in src_cube onto the grid defined by dst_cube.

    Parameters
    ----------
    src_cube: :class:`iris.cube.Cube`
        Source data. Must have latitude and longitude coords.
        These can be 1d or 2d and should have bounds.
    dst_cube: :class:`iris.cube.Cube`
        Defines the target grid.
    regrid_method:
        Selects the regridding method.
        Can be 'linear', 'area_weighted',
        or 'nearest'. See ESMPy_.

    Returns
    -------
    :class:`iris.cube.Cube`:
        The regridded cube.


    .. _ESMPy: http://www.earthsystemmodeling.org/
       esmf_releases/non_public/ESMF_7_0_0/esmpy_doc/html/
       RegridMethod.html#ESMF.api.constants.RegridMethod
    """
    src_rep, dst_rep = get_grid_representants(src, dst)
    regridder = build_regridder(src_rep, dst_rep, method)
    res = map_slices(src, regridder, src_rep, dst_rep)
    return res
