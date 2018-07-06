#!/usr/bin/env python
# -*- coding: utf-8 -*-

import ESMF
import iris
from ._mapping import ref_to_dims_index, get_slice_spec, map_slices
import numpy as np

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
    dim = len(lat.shape)
    assert(dim == len(lon.shape))
    if dim == 1:
        for c in [lat, lon]:
            if not c.has_bounds():
                c.guess_bounds()
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
        if circular:
            lat_corners = np.concatenate([lat.bounds[:, :, 0],
                                          lat.bounds[-1:, :, 1]])
            lon_corners = np.concatenate([lon.bounds[:, :, 0],
                                          lon.bounds[-1:, :, 1]])
        else:
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
    if circular:
        num_peri_dims = 1
    else:
        num_peri_dims = 0
    grid = ESMF.Grid(np.array(esmpy_lat.shape), num_peri_dims=num_peri_dims,
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
    field = ESMF.Field(grid, name=cube.long_name,
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
    lat = cube.coord('latitude')
    lon = cube.coord('longitude')
    if circular_lon is None:
        circular = is_lon_circular(lat, lon)
    else:
        circular = circular_lon
    esmpy_coords = coords_iris_to_esmpy(lat, lon, circular)
    grid = get_grid(*esmpy_coords, circular)
    field = get_empty_field(cube, grid, remove_mask)
    return field


def get_representant(cube, ref_to_slice):
    slice_dims = ref_to_dims_index(cube, ref_to_slice)
    rep_ind = [0]*cube.ndim
    for d in slice_dims:
        rep_ind[d] = slice(None, None)
    rep_ind = tuple(rep_ind)
    return cube[rep_ind]


def build_regridder(src_rep, dst_rep, method):
    regrid_method = ESMF_REGRID_METHODS[method]
    src_field = cube_to_empty_field(src_rep)
    dst_field = cube_to_empty_field(dst_rep, remove_mask=True)
    esmf_regridder = ESMF.Regrid(src_field, dst_field,
                                 regrid_method=regrid_method,
                                 src_mask_values=np.array([1]),
                                 unmapped_action=ESMF.UnmappedAction.IGNORE,
                                 ignore_degenerate=True)

    def regridder(src):
        src_field.data[...] = src.data.data.T
        regr_field = esmf_regridder(src_field, dst_field)
        data = regr_field.data[...].T
        res = np.ma.masked_array(data, data == 0.)
        return res
    return regridder


def correct_metadata(cube):
    pass


def regrid(src, dst, method='linear'):
    """
    Regrids src_cube to the grid defined by dst_cube.

    Regrids the data in src_cube onto the grid defined by dst_cube.

    Args:

    * src_cube (:class:`iris.cube.Cube`)
        Source data. Must have latitude and longitude coords.
        These can be 1d or 2d and should have bounds.

    * dst_cube (:class:`iris.cube.Cube`)
        Defines the target grid.

    * regrid_method
        Selects the regridding method.
        Can be 'linear', 'area_weighted',
        or 'nearest'. See ESMPy_

    .. _ESMPy:
    http://www.earthsystemmodeling.org/esmf_releases/non_public/ESMF_7_0_0/
    esmpy_doc/html/RegridMethod.html#ESMF.api.constants.RegridMethod

    """
    ref_to_slice = ['latitude', 'longitude']
    src_rep = get_representant(src, ref_to_slice)
    dst_rep = get_representant(dst, ref_to_slice)
    regridder = build_regridder(src_rep, dst_rep, method)
    dst_slice_spec = get_slice_spec(dst, ref_to_slice)
    res = map_slices(src, regridder, ref_to_slice, *dst_slice_spec)
    correct_metadata(res)
    return res
