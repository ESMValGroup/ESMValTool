"""
*********************************************************************
APPLICATE/TRR Ocean Diagnostics
*********************************************************************
"""
import logging
import os
import joblib
from collections import OrderedDict
import iris
import shutil
from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

import inspect
from netCDF4 import Dataset
import numpy as np
import os
import matplotlib as mpl
mpl.use('agg') #noqa
import matplotlib.pylab as plt
import math
from matplotlib import cm
from netCDF4 import num2date
#import seawater as sw
from collections import OrderedDict
from cdo import Cdo
import cmocean.cm as cmo
import matplotlib.cm as cm
from mpl_toolkits.basemap import Basemap
from mpl_toolkits.basemap import addcyclic
import pandas as pd
import pyresample
from scipy.interpolate import interp1d
import ESMF
import pyproj
#import seaborn as sns
import palettable
#plt.style.context('seaborn-talk')
#sns.set_context("paper")
import seawater as sw

class DiagnosticError(Exception):
    """Error in diagnostic"""
    
def timmean(model_filenames, mmodel,
            cmor_var, diagworkdir, observations = 'PHC'):
    logger.info("Calculate timmean {} for {}".format(cmor_var, mmodel))
    cdo = Cdo()
    ofilename = genfilename(diagworkdir, cmor_var,
                             mmodel,  data_type='timmean', extension='.nc')
    # ofilename = os.path.join(diagworkdir,
    #                         'arctic_ocean_{}_{}_timmean.nc'.format(cmor_var,
    #                                                               mmodel))
    if mmodel != observations:
        cdo.timmean(input=model_filenames[mmodel],
                output=ofilename)
    else:
        shutil.copy2(model_filenames[mmodel], ofilename)


def genfilename(basedir, variable=None,
                mmodel=None, region=None, 
                data_type=None, extension=None, basis='arctic_ocean'):
    nname = [basis,  region, mmodel, variable, data_type]
    nname_nonans = []
    for i in nname:
        if i:
            nname_nonans.append(i)
    #print(nname_nonans)
    basename = "_".join(nname_nonans)
    #print(basename)
    if extension:
        basename = basename+extension
    #print(basename)
    ifilename = os.path.join(basedir, basename)
    #print(ifilename)
    return ifilename

def get_clim_model_filenames(config, variable):
    model_filenames = {}
    for key, value in (config['input_data'].items()):
        if value['short_name'] == variable:
            model_filenames[value['dataset']] = key
    return(model_filenames)

def get_fx_filenames(config, variable, fx_var):
    areacello_fxdataset = {}
    for key, value in (config['input_data'].items()):
        if value['short_name'] == variable:
             areacello_fxdataset[value['dataset']] = value['fx_files'][fx_var]
    return(areacello_fxdataset)

def find_observations_name(config, variable):
    ''' Find "model name" of the observations data set (climatology in our case)
    Assumes that there is only one observational data set.
    '''
    obsname = []
    for key, value in (config['input_data'].items()):
        if value['project'] == "OBS":
            obsname = value['dataset']
            print(obsname)
    
    if not obsname:
        logger.info('Can\'t find observational (climatology) data')
    
    return obsname

def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.
    Source: https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib

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
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = mpl.colors.LinearSegmentedColormap(name, cdict)
    cm.register_cmap(cmap=newcmap)

    return newcmap


def dens_back(smin, smax, tmin, tmax):

    xdim = round((smax-smin)/0.1+1,0)
    ydim = round((tmax-tmin)+1,0)

    dens = np.zeros((int(ydim),int(xdim)))

    ti = np.linspace(tmin,tmax,ydim*10)
    si = np.linspace(smin,smax,xdim*10)

    si2, ti2 = np.meshgrid(si, ti)
    dens = sw.dens0(si2, ti2)-1000
    return si2, ti2, dens

def get_cmap(cmap_name):
    '''Return matplotlib colormap object 
    from matplotlib.cm or cmocean.
    Additional custom colormap for salinity is provided:
    - "custom_salinity1"
    ''' 
    #hack to support different versions of cmocean
    try:
        cmo_names = cmo.cm.cmapnames
    except:
        cmo_names = cmo.cmapnames

    if cmap_name in cmo.cmapnames:
        colormap = cmo.cmap_d[cmap_name]
    elif cmap_name in cm.datad:
        colormap = cm.get_cmap(cmap_name)
    elif cmap_name=="custom_salinity1":
        colormap = shiftedColorMap(palettable.cubehelix.cubehelix3_16.mpl_colormap, 
                                   start=0,
                                   midpoint=0.89, 
                                   stop=0.9,
                                   name='shiftedcmap')
    else:
        raise ValueError('Get unrecognised name for the colormap `{}`.\
                            Colormaps should be from standard matplotlib \
                            set or from cmocean package.'.format(cmap_name))
    return(colormap)
    