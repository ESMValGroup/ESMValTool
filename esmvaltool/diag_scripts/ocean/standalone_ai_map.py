"""
This is an independent one off script to make a bathymetry map
for the Ascension Islands plots.

Ascensionb island region is:
    central_longitude = -14.25  +/-3 # Northh -11.25
    central_latitude = -7.56 +/3 # East

"""
import matplotlib
# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')

import iris

import numpy as np

from shelve import open as shopen
import os
from glob import glob
from netCDF4 import Dataset, num2date
from matplotlib import pyplot

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvalcore.preprocessor._time import extract_time, annual_statistics, regrid_time
from esmvalcore.preprocessor._regrid import extract_levels, regrid
from esmvalcore.preprocessor._area import extract_region

# Want the figure to be a Orthographic
def main():
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.Orthographic(-10, -20))

    ax.add_feature(cfeature.OCEAN, zorder=0)
    ax.add_feature(cfeature.LAND, zorder=0, edgecolor='black')

    # lat_bnd = 20.
    # lon_bnd = 30.
    # central_longitude = -14.25 #W #-160.+3.5
    # central_latitude = -7.56
    # ax0.set_extent([central_longitude-lon_bnd,
    #                central_longitude+lon_bnd,
    #                central_latitude-lat_bnd,
    #                central_latitude+lat_bnd, ])

    ax.set_global()
    ax.gridlines()


    plt.show()

if __name__ == '__main__':
    main()
