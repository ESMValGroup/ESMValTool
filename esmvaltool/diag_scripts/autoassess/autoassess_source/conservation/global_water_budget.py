'''
Module with routines to calculate global water budgets
'''

import datetime
import os

import iris
import numpy as np

from .loaddata import load_run_ss
from .area_utils import area_average


def fluxes_submodel(run, stash_f, f_mult):
    '''
    function to calculate time-average global water fluxes from a UM submodel
    (atmosphere, ocean, sea ice, soil moisture, land snow, TRIP).

    It assumes that units are mks:  Kg/m2s.  If this is not the case, use
    f_mult multiplicative factors to achieve that.
    Use also f_mult to multiply fields by fractional masks, etc.
    '''
    # Initialize array of global budgets with zeros:
    fval = np.zeros(len(stash_f)+1, dtype=np.float)

    # TODO use NaN
    # Set value of data-missing indicator:
    mdi = -10000.0

    # Read data and calculate global budgets

    # read data:
    try:
        for (i, stash) in enumerate(stash_f):
            # VPREDOI
            # need the correct files
            # ppfy = load_run_ss(run, 'annual', stash)
            ppfy = load_run_ss(run, 'monthly', 'eastward_wind', lbproc=192)

            ppfy *= f_mult[i]
            # Calculate global mean, time series of global values
            ppfyg = area_average(ppfy, weighted=True,
                                 aggregator=iris.analysis.SUM)
            # multi-year mean of global values
            ppfy_total = ppfyg.collapsed('time', iris.analysis.MEAN)

            # VPREDOI
            # reducing the dim; this should not be done
            # in the real diag
            # fval[i] = ppfy_total.data/1.0e9  # TODO magic number
            fval[i] = ppfy_total.data[i]/1.0e9  # TODO magic number

        fval[-1] = fval.sum()

    except iris.exceptions.ConstraintMismatchError:
        print("ERROR:  Missing data!!!")
        print("An MDI will be assigned to all water fluxes of this sub-model")
        fval[:] = mdi

    return fval


def fluxes_ocean_submodel(expid, mesh, areas, opath, wfpath, yi, yf):
    '''
    Function to calculate time-average global water fluxes that go into
    NEMO ocean model
    '''
    # load land-fraction-mask an# array with global values of Rd area cubes:
    mask = iris.load_cube(mesh)   # NEMO ocean-fraction mask
    area = iris.load_cube(areas)  # what NEMO  uses for grid cell areas

    # Initialize array of global budgets with zeros:
    r = np.zeros(yf-yi+1, dtype=np.float)    # global values of R
    net = np.zeros(yf-yi+1, dtype=np.float)  # global net flux into the ocean
    pme = np.zeros(yf-yi+1, dtype=np.float)  # global P-E into the ocean

    # initialize data missing error flag:
    errf = 1

    # TODO use NaN
    # Set value of data-missing indicator:
    mdi = -10000.0

    # Calculate water fluxes

    # Constant iceberg-calving flux:
    # Read pp-file and calculate global budget:

    wfix = iris.load_cube(wfpath)
    iceberg = area_average(wfix, weighted=True, aggregator=iris.analysis.SUM)
    iceberg /= 1.0e9  # TODO magic number

    # calculate  R, P-E and net flux:

    # TODO local path
    # read NEMO netcdf files:
    filetmp = opath + '/' + expid + 'o_1y_{0}1201_{1}1130_grid_T.nc'
    nemofiles = [filetemp.format(x, x+1) for x in range(yi, yf+1, 1)]
    for (i, model_file) in enumerate(nemofiles):
        if os.path.exists(model_file):
            sorunoff_name = 'water_flux_into_sea_water_from_rivers'
            sorunoff = iris.load_cube(model_file, sorunoff_name)
            sowaflup_name = 'water_flux_out_of_sea_ice_and_sea_water'
            sowaflup = iris.load_cube(model_file, sowaflup_name)

            # change mdi values to zeros:
            sorunoff.data = np.where(sorunoff.data < 1.0, sorunoff.data, 0.0)
            sowaflup.data = np.where(sowaflup.data < 1.0, sowaflup.data, 0.0)

            # get rid of duplicate points (tripolar grid) and set to
            #  zero over land:
            rflux = sorunoff*mask*area
            nflux = sowaflup*mask*area

            # TODO magic number
            # calculate global budgets:
            r[i] = rflux.data.sum() / 1.0e9 - iceberg.data
            net[i] = nflux.data.sum() / (-1.0e9)
            # P-E is obtained as a residual from net, R and iceberg-calving
            pme[i] = net[i]-r[i]-iceberg.data

        else:
            errf = 0

    if errf:
        fval = np.array([np.mean(pme), np.mean(r), iceberg.data, np.mean(net)])
    else:
        print("ERROR:  Missing data!!!")
        print("An MDI value will be assigned to all ocean water fluxes")
        fval = np.array([mdi, mdi, mdi, mdi])

    return fval
