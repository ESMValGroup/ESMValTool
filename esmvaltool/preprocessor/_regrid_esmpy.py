# -*- coding: utf-8 -*-

import ESMF
ESMF_MANAGER = ESMF.Manager(debug=False)
import iris  # noqa
import numpy as np  # noqa
from scipy.interpolate import NearestNDInterpolator  # noqa
from scipy.ndimage.morphology import binary_closing  # noqa


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


def get_grid(esmpy_lat, esmpy_lon, esmpy_lat_corners, esmpy_lon_corners, circular):
    if circular:
        num_peri_dims = 1
    else:
        num_peri_dims = 0
    grid = ESMF.Grid(np.array(esmpy_lat.shape), num_peri_dims=num_peri_dims,
                     staggerloc=[ESMF.StaggerLoc.CENTER])
    grid.get_coords(ESMF_LON)[...] = esmpy_lon
    grid.get_coords(ESMF_LAT)[...] = esmpy_lat
    grid.add_coords([ESMF.StaggerLoc.CORNER])
    gridLonCorner = grid.get_coords(ESMF_LON, staggerloc=ESMF.StaggerLoc.CORNER)
    gridLatCorner = grid.get_coords(ESMF_LAT, staggerloc=ESMF.StaggerLoc.CORNER)
    gridLonCorner[...] = esmpy_lon_corners
    gridLatCorner[...] = esmpy_lat_corners
    return grid


def get_field(cube, grid):
    field = ESMF.Field(grid, name=cube.long_name,
                       staggerloc=ESMF.StaggerLoc.CENTER)
    if grid.mask[0] is None:
        center_mask = grid.add_item(ESMF.GridItem.MASK, ESMF.StaggerLoc.CENTER)
    else:
        center_mask = grid.get_item(ESMF.GridItem.MASK, ESMF.StaggerLoc.CENTER)
    if np.ma.isMaskedArray(cube.data):
        field.data[...] = cube.data.data.T
        mask = cube.data.mask.T
        center_mask[...] = mask
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


def cube_to_field(cube, circular_lon=None):
    lat = cube.coord('latitude')
    lon = cube.coord('longitude')
    if circular_lon is None:
        circular = is_lon_circular(lat, lon)
    else:
        circular = circular_lon
    try:
        esmpy_coords = coords_iris_to_esmpy(lat, lon, circular)
    except:
        print(cube)
        raise
    grid = get_grid(*esmpy_coords, circular)
    field = get_field(cube, grid)
    return field


def build_regridded_cube(src_cube, dst_cube, regridded_field):
    regridded_cube = dst_cube.copy()
    regridded_cube.metadata = src_cube.metadata
    regridded_cube.data.data[...] = regridded_field.data[...].T
    regridded_cube.data.mask[...] = regridded_field.grid.mask[0][...].T
    for ac in src_cube.aux_coords:
        if len(src_cube.coord_dims(ac)) == 0:
            regridded_cube.add_aux_coord(ac)
    return regridded_cube


def prepare_fields(src_field, dst_field, regrid_method):
    src_mask_field = ESMF.Field(src_field.grid)
    src_mask = src_field.grid.get_item(ESMF.GridItem.MASK)
    src_mask_field.data[...] = src_mask[...]
    dst_field.data[...] = 1
    regridder = ESMF.Regrid(src_mask_field, dst_field,
                            regrid_method=regrid_method,
                            src_mask_values=np.array([-1]),
                            dst_mask_values=np.array([-1]),
                            unmapped_action=ESMF.UnmappedAction.IGNORE,
                            ignore_degenerate=True)
    regridded_field = regridder(src_mask_field, dst_field)
    mask = regridded_field.data[...].copy()
    closed_mask = binary_closing(mask)
    closed_mask[0, :] = mask[0, :]
    closed_mask[:, -1] = mask[:, -1]
    closed_mask[-1, :] = mask[-1, :]
    closed_mask[:, 0] = mask[:, 0]
    mask = closed_mask
    mask[mask > 0.] = 1
    mask = mask.astype(int)
    dst_field.grid.get_item(ESMF.GridItem.MASK)[...] = mask
    src_mask_values = np.array([1])
    return src_mask_values


def regrid_esmpy(src_cube, dst_cube,
                 regrid_method='linear'):
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

    .. _ESMPy: http://www.earthsystemmodeling.org/esmf_releases/non_public/ESMF_7_0_0/esmpy_doc/html/RegridMethod.html#ESMF.api.constants.RegridMethod

    """
    regrid_method = ESMF_REGRID_METHODS[regrid_method]
    src_field = cube_to_field(src_cube)
    dst_field = cube_to_field(dst_cube)
    src_mask_values = prepare_fields(src_field, dst_field, regrid_method)
    regridder = ESMF.Regrid(src_field, dst_field, regrid_method=regrid_method,
                            src_mask_values=src_mask_values,
                            dst_mask_values=np.array([1]),
                            unmapped_action=ESMF.UnmappedAction.IGNORE,
                            ignore_degenerate=True)
    regridded_field = regridder(src_field, dst_field)
    regridded_cube = build_regridded_cube(src_cube, dst_cube, regridded_field)
    return regridded_cube


def regrid(src_cube, dst_cube, method):
    regrid_method = ESMF_REGRID_METHODS[method]
    src_representant = src_cube[0]
    dst_representant = dst_cube[0]
    src_field = cube_to_field(src_representant)
    dst_field = cube_to_field(dst_representant)
    src_mask_values = prepare_fields(src_field, dst_field, regrid_method)
    regridder = ESMF.Regrid(src_field, dst_field, regrid_method=regrid_method,
                            src_mask_values=src_mask_values,
                            dst_mask_values=np.array([1]),
                            unmapped_action=ESMF.UnmappedAction.IGNORE,
                            ignore_degenerate=True)

    def do_regrid(src):
        src_field = cube_to_field(src)
        regridded_field = regridder(src_field, dst_field)
        regridded_cube = build_regridded_cube(src,
                                              dst_representant,
                                              regridded_field)
        return regridded_cube
    lat_dims = src_cube.coord_dims('latitude')
    lon_dims = src_cube.coord_dims('longitude')
    slice_dims = set(lat_dims + lon_dims)
    cl = iris.cube.CubeList([do_regrid(src)
                             for src in src_cube.slices(slice_dims)])
    print(cl)
    return cl.merge_cube()
