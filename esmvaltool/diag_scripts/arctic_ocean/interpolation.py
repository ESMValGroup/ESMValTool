# -*- coding: utf-8 -*-
"""Part of the ESMValTool Arctic Ocean diagnostics.

This module contains functions for data interpolation.
"""
import logging
import os
import ESMF
import numpy as np
# import pyresample
from cartopy.util import add_cyclic_point
# from netCDF4 import Dataset

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


def define_esmf_field(ifile, data_onlevel, name):
    """Define ESMF field from netCDF file."""

    grid_obs = ESMF.Grid(filename=ifile, filetype=ESMF.FileFormat.GRIDSPEC)
    mask_obs = grid_obs.add_item(ESMF.GridItem.MASK)
    mask_obs[:] = data_onlevel.mask.astype('int').T
    esmf_field = ESMF.Field(
        grid_obs,
        staggerloc=ESMF.StaggerLoc.CENTER,
        name=name,
    )
    return esmf_field


def add_esmf_cyclic(metadata_obs, data_onlevel, interpolated):
    """Add cyclic points to interpolated data."""

    data_onlevel_cyc, lon_obs_cyc = add_cyclic_point(
        data_onlevel, coord=metadata_obs['lon2d'][0, :])

    lonc, latc = np.meshgrid(lon_obs_cyc, metadata_obs['lat2d'][:, 0])

    interpolated_cyc, lon_obs_cyc = add_cyclic_point(
        interpolated, coord=metadata_obs['lon2d'][0, :])
    return lonc, latc, data_onlevel_cyc, interpolated_cyc


def esmf_regriding(sourcefield, distfield, metadata_obs, data_onlev_obs):
    """Use ESMF fields to do the regriding."""
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
    data_interpolated = np.ma.masked_equal(data_interpolated, 0)
    lonc, latc, data_onlevel_cyc, interpolated_cyc = add_esmf_cyclic(
        metadata_obs, data_onlev_obs, data_interpolated)
    return lonc, latc, data_onlevel_cyc, interpolated_cyc


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
    metadata_obs = load_meta(obs_file, fxpath=None)
    metadata_mod = load_meta(mod_file, fxpath=None)

    data_obs = metadata_obs['datafile'].variables[cmor_var][:]
    data_model = metadata_mod['datafile'].variables[cmor_var][:]

    # Select depth in climatology that is closest to the desired depth
    target_depth, level_depth = closest_depth(metadata_obs['lev'], depth)

    # climatology and model data on the level
    data_onlev_obs = data_obs[0, level_depth, :, :]
    data_onlev_mod = interpolate_vert(metadata_mod['lev'], target_depth,
                                      data_model[0, :, :, :])

    # prepear interpolation fields
    distfield = define_esmf_field(obs_file, data_onlev_obs, 'OBS')
    distfield.data[:] = 0.0

    sourcefield = define_esmf_field(mod_file, data_onlev_mod, 'Model')
    sourcefield.data[...] = data_onlev_mod.T

    lonc, latc, data_onlev_obs_cyc, data_interpolated_cyc = esmf_regriding(
        sourcefield, distfield, metadata_obs, data_onlev_obs)

    return lonc, latc, target_depth, data_onlev_obs_cyc, data_interpolated_cyc
