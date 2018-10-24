import logging
import os
import joblib
from collections import OrderedDict
import iris

from esmvaltool.diag_scripts.shared import run_diagnostic
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))

import inspect
from netCDF4 import Dataset
import numpy as np
import os
import matplotlib as mpl
mpl.use('agg')
import matplotlib.pylab as plt
import math
from matplotlib import cm
from netCDF4 import num2date
#import seawater as sw
from collections import OrderedDict
from cdo import Cdo
import cmocean.cm as cmo
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



def hofm_regions(region, lon2d, lat2d):
    if region == 'EB':
        # Eurasian Basin of the Arctic Ocean
        indi, indj = np.where((lon2d > 300) & (lat2d > 80))
        indi2, indj2 = np.where((lon2d < 100) & (lat2d > 80))
        indi3, indj3 = np.where((lon2d > 100) & (lon2d < 140) & (lat2d > 66))

        indexesi = np.hstack((indi, indi2, indi3))
        indexesj = np.hstack((indj, indj2, indj3))
    elif region == 'AB':
        # Amerasian basin of the Arctic Ocean
        indi, indj = np.where((lon2d >= 260) & (lon2d <= 300) & (lat2d >= 80))
        indi2, indj2 = np.where((lon2d >= 140) & (lon2d < 260) & (lat2d > 66))

        indexesi = np.hstack((indi, indi2))
        indexesj = np.hstack((indj, indj2))
    elif region == 'Barents_sea':

        indi, indj = np.where((lon2d >= 20) & (lon2d <= 55) & (lat2d >= 70) & (lat2d <= 80))

        indexesi = indi
        indexesj = indj
    elif region == 'North_sea':
        # Amerasian basin of the Arctic Ocean
        indi, indj = np.where((lon2d >= 355) & (lon2d <= 360) & (lat2d >= 50)& (lat2d <= 60))
        indi2, indj2 = np.where((lon2d >= 0) & (lon2d <= 10) & (lat2d >= 50)& (lat2d <= 60))

        indexesi = np.hstack((indi, indi2))
        indexesj = np.hstack((indj, indj2))
    else:
        print('Region {} is not recognized'.format(region))
    return indexesi, indexesj

def transect_points(transect, mult = 2):
    if transect == 'AWpath':
        lon_s4=np.array([17.6, 16.5, 16.05, 15.6, 15.1, 14.1, 13.0, 12.0, 10.0, 8.0, 4.0
             , 4.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0
             , 110.0, 120.0, 130.0, 140.0])
        lat_s4=np.array([69.0, 70.6, 71.3, 72.02, 72.8, 73.8, 75.0, 76.0, 77.0, 78.0, 79.0, 80.0
    , 81.0, 81.8, 81.8, 82.6, 83.0, 83.2, 83.1, 82.8, 82.5, 81.8, 79.7, 78.2
    , 78.7, 79.7])
    elif transect == 'Fram':
        lon_s4=np.array([-22.3732, -21.4796, -19.6479, -18.1074, -16.7828, -15.504 ,
            -14.2042, -12.9771, -11.6642, -10.1892,  -8.7414,  -7.719 ,
                -6.3646,  -4.4803,  -3.4232,  -2.4435,  -1.615 ,  -0.6752,
                0.343 ,   1.6947,   2.7157,   3.7374,   4.6099,   5.5097,
                6.3754,   7.2394,   7.9238,   8.7029,   9.7338,  10.4462,
                11.0559,  12.0102,  13.3313])
        lat_s4=np.array([ 78.9373,  78.9276,  78.9183,  78.9356,  78.9346,  78.9334,
                78.9425,  78.9434,  78.9274,  78.9392,  78.9287,  78.9262,
                78.9392,  78.95  ,  78.9405,  78.9347,  78.9334,  78.922 ,
                78.9287,  78.9131,  78.919 ,  78.9215,  78.9242,  78.909 ,
                78.8995,  78.8874,  78.8865,  78.9026,  78.8992,  78.8841,
                78.8793,  78.8715,  78.9012])

    else:
        print('Transect {} is not recognized'.format(transect))
        
    x = np.linspace(1, lon_s4.shape[0], num=lon_s4.shape[0], endpoint=True)
    f = interp1d(x,lon_s4)
    g = interp1d(x,lat_s4)
    xnew = np.linspace(1, lon_s4.shape[0], num=mult*lon_s4.shape[0],
                    endpoint=True)

    lon_s4new = f(xnew)
    lat_s4new = g(xnew)
    return lon_s4new, lat_s4new
