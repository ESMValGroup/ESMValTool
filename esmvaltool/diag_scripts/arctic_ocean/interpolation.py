# -*- coding: utf-8 -*-
import logging
import os
import ESMF
import numpy as np
import pyresample
from cartopy.util import add_cyclic_point
from netCDF4 import Dataset

from esmvaltool.diag_scripts.arctic_ocean.getdata import load_meta

logger = logging.getLogger(os.path.basename(__file__))


def closest_depth(depths, depth):
    """Find closest depth.

    From vector of depths finds target depth,
    that is closest to desired depth, Also returns
    an index for the desired depth.
    """
    target_level = abs(abs(depths) - abs(depth)).argmin()
    target_depth = depths[target_level]
    logger.debug('target_depth: %s', target_depth)
    return target_depth, target_level


def interpolate_vert(depth_model, target_depth, data_model):
    """Vertical linear interpolation.

    Very simple linear interpolation of the model data to the
    desired depth. Can't extrapolate, so the limitation is that model
    data should have at least one level => and one level <= than the
    target depth.
    """
    # Simple vertical interpolation
    dep_up = [z for z in abs(depth_model) if z <= target_depth][-1]
    dep_lo = [z for z in abs(depth_model) if z > target_depth][0]
    i_up = 1 - abs(target_depth - dep_up) / (dep_lo - dep_up)
    i_lo = 1 - abs(target_depth - dep_lo) / (dep_lo - dep_up)

    iz_up = abs(abs(depth_model) - abs(dep_up)).argmin()
    iz_lo = abs(abs(depth_model) - abs(dep_lo)).argmin()

    data_up = data_model[iz_up, :, :]
    data_lo = data_model[iz_lo, :, :]
    if not isinstance(data_up, np.ma.MaskedArray):
        data_up = np.ma.masked_equal(data_up, 0)
        data_lo = np.ma.masked_equal(data_lo, 0)
    data = i_up * data_up
    data = data + i_lo * data_lo
    return data


def weighting(distance):
    """Weighting function for pyresample."""
    weight = 1 / distance**2
    return weight


def interpolate_pyresample(obs_file, mod_file, depth, cmor_var):
    """The 2d interpolation with pyresample.

    Simple realisation of horizontal 2d interpolation with
    pyresample. Before spatial interpolation data are linearly
    interpolated to the `obs_file` depth closest
    to the desired `depth`.

    Parameters
    ----------
    obs_file : str
        path to the file with target grid
    mod_file : str
        path to the file with source grid.
    depth : float
        desired depth.
    cmor_var : str
        cmor name of the variable to interpolate.
        In case of ocean it's isially 'thetao' or 'so'

    Returns
    -------
    lonc, latc, target_depth, data_onlev_obs_cyc, interpolated
    lonc : 2d np array
        numpy array with longitudes of target grid
    latc : 2d np array
        numpy array with longitudes of target grid
    target_depth : float
        the `obs_file` depth closest to the desired `depth`
        that was actually used for interpolation.
    interpolated : 2d np array
        Field with result of interpolation
    """
    # load observations data
    obs = Dataset(obs_file)
    data_obs = obs.variables[cmor_var][:]
    # salt_phc  = phc.variables['salt'][:]
    lon_obs = obs.variables['lon'][:]
    lat_obs = obs.variables['lat'][:]
    depth_obs = obs.variables['lev'][:]

    # Select depth in climatology that is closest to the desired depth
    target_depth, level_depth = closest_depth(depth_obs, depth)

    # climatology data on the level
    data_onlev_obs = data_obs[level_depth, :, :]

    # add cyclic point to data and coordinates
    data_onlev_obs_cyc, lon_obs_cyc = add_cyclic_point(data_onlev_obs,
                                                       lon_obs[0, :])
    lonc, latc = np.meshgrid(lon_obs_cyc, lat_obs[:, 0])

    # define target grid
    targ_def = pyresample.geometry.SwathDefinition(lons=lonc - 180, lats=latc)

    # Now working with the model
    model = Dataset(mod_file)
    data_model = model.variables['thetao'][:]
    lon_model = model.variables['lon'][:]
    lat_model = model.variables['lat'][:]

    lat_model[lat_model > 90] = 90
    depth_model = model.variables['lev'][:]

    # some ocean models still have 1d coordinates
    if lon_model.ndim == 2:
        lon2d, lat2d = lon_model, lat_model
    elif lon_model.ndim == 1:
        lon2d, lat2d = np.meshgrid(lon_model, lat_model)

    # Simple vertical interpolation
    data = interpolate_vert(depth_model, target_depth, data_model[0, :, :, :])

    # interpolation (weighted nearest neighbor)
    # here one can implement other methods available in pyresample
    # two main control parameters are radius_of_influence and neighbours,
    # that can be made available for the user too

    # define original (model) grid
    orig_def = pyresample.geometry.SwathDefinition(lons=lon2d - 180,
                                                   lats=lat2d)

    interpolated = pyresample.kd_tree.resample_custom(
        orig_def,
        data,
        targ_def,
        radius_of_influence=150000,
        neighbours=10,
        weight_funcs=weighting,
        fill_value=None)

    return lonc, latc, target_depth, data_onlev_obs_cyc, interpolated


def interpolate_esmf(obs_file, mod_file, depth, cmor_var):
    """The 2d interpolation with ESMF.

    Parameters
    ----------
    obs_file: str
        path to file with observations/climatology,
        will be used to extract the grid to interpolate on to.
    mod_file: str
        path to the file with model data.
    depth: int
        depth to interpolate to. First the closest depth from the
        observations will be selected and then.
    """
    metadata = load_meta(obs_file, fxpath=None)
    # obs, lon_obs, lat_obs, depth_obs, time, areacello = metadata
    obs = metadata['datafile']
    lon_obs = metadata['lon2d']
    lat_obs = metadata['lat2d']
    depth_obs = metadata['lev']

    data_obs = obs.variables[cmor_var][:]

    # Select depth in climatology that is closest to the desired depth
    target_depth, level_depth = closest_depth(depth_obs, depth)

    # climatology data on the level
    data_onlev_obs = data_obs[0, level_depth, :, :]

    # Setting ESMF
    grid_obs = ESMF.Grid(filename=obs_file, filetype=ESMF.FileFormat.GRIDSPEC)
    mask_obs = grid_obs.add_item(ESMF.GridItem.MASK)
    mask_obs[:] = data_onlev_obs.mask.astype('int').T
    distfield = ESMF.Field(
        grid_obs,
        staggerloc=ESMF.StaggerLoc.CENTER,
        name='OBS',
    )
    distfield.data[:] = 0.0

    # Now working with the model
    metadata = load_meta(mod_file, fxpath=None)
    model = metadata['datafile']
    depth_model = metadata['lev']

    data_model = model.variables[cmor_var][:]

    # Simple vertical interpolation
    data = interpolate_vert(depth_model, target_depth, data_model[0, :, :, :])

    # set the model grid
    grid_model = ESMF.Grid(filename=mod_file,
                           filetype=ESMF.FileFormat.GRIDSPEC)

    model_mask = grid_model.add_item(ESMF.GridItem.MASK)
    model_mask[:] = data.mask.astype('int').T

    # define asource field
    sourcefield = ESMF.Field(
        grid_model,
        staggerloc=ESMF.StaggerLoc.CENTER,
        name='Model',
    )
    # define a destination field
    distfield = ESMF.Field(
        grid_obs,
        staggerloc=ESMF.StaggerLoc.CENTER,
        name='Model_interp',
    )
    sourcefield.data[...] = data.T

    # define the regrider
    regrid = ESMF.Regrid(
        sourcefield,
        distfield,
        regrid_method=ESMF.RegridMethod.NEAREST_STOD,
        # regrid_method=ESMF.RegridMethod.BILINEAR,
        unmapped_action=ESMF.UnmappedAction.IGNORE,
        dst_mask_values=np.array([1]),
        src_mask_values=np.array([1]))
    # actual regriding
    distfield = regrid(sourcefield, distfield)
    # reshape the data and convert to masked array
    data_interpolated = distfield.data[:].T
    interpolated = np.ma.masked_equal(data_interpolated, 0)
    # add cyclic points
    data_onlev_obs_cyc, lon_obs_cyc = add_cyclic_point(data_onlev_obs,
                                                       coord=lon_obs[0, :])
    lonc, latc = np.meshgrid(lon_obs_cyc, lat_obs[:, 0])

    interpolated_cyc, lon_obs_cyc = add_cyclic_point(interpolated,
                                                     coord=lon_obs[0, :])

    return lonc, latc, target_depth, data_onlev_obs_cyc, interpolated_cyc
