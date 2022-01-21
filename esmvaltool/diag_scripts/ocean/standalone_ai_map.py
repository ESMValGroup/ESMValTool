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
import os, shutil
from glob import glob
from netCDF4 import Dataset, num2date
from matplotlib import pyplot
import matplotlib.patches as mpatches

import cartopy.crs as ccrs
import cartopy.feature as cfeature


from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvalcore.preprocessor._time import extract_time, annual_statistics, regrid_time
from esmvalcore.preprocessor._regrid import extract_levels, regrid
from esmvalcore.preprocessor._area import extract_region

# Want the figure to be a Orthographic

def compute_radius(ortho, radius_degrees, proj= ccrs.PlateCarree(), lat=0, lon=0):
    """
    catlculate the correct radius:
    from:
    https://stackoverflow.com/questions/52105543/drawing-circles-with-cartopy-in-orthographic-projection
    """
    phi1 = lat + radius_degrees if lat <= 0 else lat - radius_degrees
    _, y1 = ortho.transform_point(lon, phi1, proj)
    return abs(y1)

def add_shadow(fig,ax, proj,):

    num_shadows = 10
    for i,w in enumerate(np.linspace(0., 2., num_shadows+1)): 
        h = w/8.
        alpha = 1./num_shadows
        ellipse = mpatches.Ellipse(xy=(0, -0.1), width=w, height=h, fc=(0,0,0,alpha), transform=ax.transAxes,zorder = 40)
        ax.add_artist(ellipse)

    return fig, ax


def main():
    path = diagtools.folder('images/')+'AI_map_4.png'
    #if os.path.exists(path):
        
    fig = pyplot.figure()
    #proj = ccrs.Orthographic(-10, -20)
    proj = ccrs.PlateCarree()
    ax = fig.add_subplot(1, 1, 1, projection=proj)
#   ax.stock_img()

    #ax.add_feature(cfeature.OCEAN, zorder=2)
    #ax.add_feature(cfeature.LAND, zorder=2, edgecolor='black')


    # lat_bnd = 20.
    # lon_bnd = 30.
    central_longitude = -14.25 #W #-160.+3.5
    central_latitude = -7.56
    # ax.set_extent([central_longitude-lon_bnd,
    #                central_longitude+lon_bnd,
    #                central_latitude-lat_bnd,
    #                central_latitude+lat_bnd, ])

    # Add transect transect at longitude 23°W, latitude -3°:3° and depth from 0 to 400m

    r_ortho = compute_radius(proj, 3., proj=proj, lat = central_latitude, lon=central_longitude,)
    ax.add_patch(mpatches.Circle(xy=[central_longitude, central_latitude], radius=r_ortho, color='black', alpha=0.3, transform=proj, zorder=30))

    ax.plot([-23, -23], [-3, 3], c='green', lw=2., transform=proj, ) #transform=ccrs.PlateCarree(), zorder=20,)

    #ax.set_global()
    ax.gridlines()
    # ax.background_img(name='50-natural-earth-1-downsampled.png')
    ax.stock_img()


    fig, ax = add_shadow(fig, ax, proj)
    #pyplot.show()
    print('saving figure to:', path)
    pyplot.savefig(path)
    pyplot.close()
    
    copy_path = diagtools.folder('/home/users/ldemora/Master_CRACAB_plots/figure1')+os.path.basename(path)
    print('copying to ', copy_path)
    shutil.copy(path, copy_path)
    
    print('Now run on pmpc:')
    print('rsync -avP ldemora@xfer1.jasmin.ac.uk:/home/users/ldemora/Master_CRACAB_plots /users/modellers/ledm/ImagesFromJasmin/.')


if __name__ == '__main__':
    main()
