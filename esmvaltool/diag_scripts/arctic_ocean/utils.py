# -*- coding: utf-8 -*-
import logging
import os
import shutil

import cmocean.cm as cmo
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
import palettable
import pyproj
import seawater as sw
from cdo import Cdo

mpl.use('agg')
logger = logging.getLogger(os.path.basename(__file__))


class DiagnosticError(Exception):
    """Error in diagnostic."""


def genfilename(basedir,
                variable=None,
                mmodel=None,
                region=None,
                data_type=None,
                extension=None,
                basis='arctic_ocean'):
    """Generate file name for the output data.

    Parameters
    ----------
    basedir: str
        base directory
    variable: str
        name of the variable
    mmodel: str
        name of the model
    region: str
        name of the region
    data_type: str
        type of the data, for example `timmean`
    extention: str
        fiel extention, for example `nc`
    basis: str
        basis name that can be used for series of
        diagnostics

    Returns
    -------
    ifilename: str
        path to the file
    """
    nname = [basis, region, mmodel, variable, data_type]
    nname_nonans = []
    for i in nname:
        if i:
            nname_nonans.append(i)
    basename = "_".join(nname_nonans)
    if extension:
        basename = basename + extension
    ifilename = os.path.join(basedir, basename)
    return ifilename


def timmean(model_filenames, mmodel, cmor_var, diagworkdir,
            observations='PHC'):
    """Create time mean of input data with cdo.

    Parameters
    ----------
    model_filenames: OrderedDict
        OrderedDict with model names as keys and input files as values.
    mmodel: str
        model name that will be processed
    cmor_var: str
        name of the CMOR variable
    diagworkdir: str
        path to the working directory
    observations: str
        name of observational/climatology data set.

    Returns
    -------
    None
    """
    logger.info("Calculate timmean %s for %s", cmor_var, mmodel)
    cdo = Cdo()
    ofilename = genfilename(diagworkdir,
                            cmor_var,
                            mmodel,
                            data_type='timmean',
                            extension='.nc')
    if mmodel != observations:
        cdo.timmean(input=model_filenames[mmodel], output=ofilename)
    else:
        shutil.copy2(model_filenames[mmodel], ofilename)


def get_clim_model_filenames(config, variable):
    """Extract model filenames from the configuration."""
    model_filenames = {}
    for key, value in config['input_data'].items():
        if value['short_name'] == variable:
            model_filenames[value['dataset']] = key
    return model_filenames


def get_fx_filenames(config, variable, fx_var):
    """Extract fx file names."""
    areacello_fxdataset = {}
    for _, value in config['input_data'].items():
        if value['short_name'] == variable:
            areacello_fxdataset[value['dataset']] = value['fx_files'][fx_var]
    return areacello_fxdataset


def find_observations_name(config):
    """Find "model name" of the observations data set.

    Assumes that there is only one observational data set.
    """
    obsname = []
    for _, value in config['input_data'].items():
        if value['project'] == "OBS":
            obsname = value['dataset']
            print(obsname)
    if not obsname:
        logger.info('Can\'t find observational (climatology) data')

    return obsname


def shiftedcolormap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    """Offset the "center" of a colormap.

    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.
    Source: https://stackoverflow.com/questions/7404116/
    defining-the-midpoint-of-a-colormap-in-matplotlib

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax/(vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          `midpoint` and 1.0.
    """
    cdict = {'red': [], 'green': [], 'blue': [], 'alpha': []}

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False),
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for regi, shii in zip(reg_index, shift_index):
        red, gren, blue, alpha = cmap(regi)

        cdict['red'].append((shii, red, red))
        cdict['green'].append((shii, gren, gren))
        cdict['blue'].append((shii, blue, blue))
        cdict['alpha'].append((shii, alpha, alpha))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    cm.register_cmap(cmap=newcmap)

    return newcmap


def dens_back(smin, smax, tmin, tmax):
    """Calculate density for TS diagram."""
    xdim = round((smax - smin) / 0.1 + 1, 0)
    ydim = round((tmax - tmin) + 1, 0)

    dens = np.zeros((int(ydim), int(xdim)))

    ti = np.linspace(tmin, tmax, ydim * 10)
    si = np.linspace(smin, smax, xdim * 10)

    si2, ti2 = np.meshgrid(si, ti)
    dens = sw.dens0(si2, ti2) - 1000
    return si2, ti2, dens


def get_cmap(cmap_name):
    """Return matplotlib colormap object.

    From matplotlib.cm or cmocean.
    Additional custom colormap for salinity is provided:
    - "custom_salinity1"
    """
    if cmap_name in cmo.cmapnames:
        colormap = cmo.cmap_d[cmap_name]
    elif cmap_name in cm.datad:
        colormap = cm.get_cmap(cmap_name)
    elif cmap_name == "custom_salinity1":
        colormap = shiftedcolormap(
            palettable.cubehelix.cubehelix3_16.mpl_colormap,
            start=0,
            midpoint=0.89,
            stop=0.9,
            name='shiftedcmap')
    else:
        raise ValueError('Get unrecognised name for the colormap `{}`.\
                            Colormaps should be from standard matplotlib \
                            set or from cmocean package.'.format(cmap_name))
    return colormap


def point_distance(lon_s4new, lat_s4new):
    """Calculate distance between points of the section.

    Parameters
    ----------
    lon_s4new: numpy array
        1d array of longitudes
    lat_s4new: numpy array
        1d array of latitudes

    Returns
    -------
    dist: numpy array
        1d array of distances between points in km.
    """
    g = pyproj.Geod(ellps='WGS84')
    (_, _, dist) = g.inv(lon_s4new[0:-1], lat_s4new[0:-1], lon_s4new[1:],
                         lat_s4new[1:])
    dist = dist.cumsum() / 1000
    dist = np.insert(dist, 0, 0)
    return dist
