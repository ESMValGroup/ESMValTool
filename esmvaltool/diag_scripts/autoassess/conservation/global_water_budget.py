"""Module with routines to calculate global water budgets"""

import os

import iris
import numpy as np

from esmvaltool.diag_scripts.autoassess.loaddata import load_run_ss
from .area_utils import area_average


def fluxes_submodel(run, stash_f, f_mult):
    """
    Calculate time avg fluxes

    function to calculate time-average global water fluxes from a UM submodel
    (atmosphere, ocean, sea ice, soil moisture, land snow, TRIP).

    It assumes that units are mks:  Kg/m2s.  If this is not the case, use
    f_mult multiplicative factors to achieve that.
    Use also f_mult to multiply fields by fractional masks, etc.
    """
    ########################
    # Stash codes used here:
    ########################
    # m01s08i234: surface_runoff_flux; CMOR-Lmon: mrros
    # m01s08i235: subsurface_runoff_flux; CMOR-Lmon: mrro ?
    # m01s08i245: no std name; ~ Inland basin runoff; nope CMOR
    # m01s26i004: water_flux_into_sea_water_from_rivers; nope CMOR
    # m01s04i204: stratiform_snowfall_flux -> prsn: snowfall_flux
    # m01s05i206: convective_snowfall_flux -> prsn: snowfall_flux
    # m01s03i298: water_sublimation_flux
    # -> sbl: surface_snow_and_ice_sublimation_flux -- hard to find
    # m01s03i353: water_sublimation_flux
    # (?, from sea-ice over sea, masked over land)
    # -> sbl: surface_snow_and_ice_sublimation_flux -- hard to find
    # m01s08i231: surface_snow_melt_flux (over land, masked over sea) -> ???
    # m01s04i203: stratiform_rainfall_flux -> ???
    # m01s05i205: convective_rainfall_flux -> ???
    # m01s03i223: surface_upward_water_flux-> evspsbl: water_evaporation_flux
    # m01s03i232: surface_upward_water_flux
    # (from sea, masked over land) -> evspsbl: water_evaporation_flux
    # m01s05i216: precipitation_flux -> pr: precipitation_flux

    # Data loading by standard_name rather than by stash
    stash_dict = {
        'm01s08i234': 'surface_runoff_flux',
        'm01s08i235': 'surface_runoff_flux',
        'm01s08i245': 'surface_runoff_flux',
        'm01s26i004': 'surface_runoff_flux',
        'm01s04i204': 'snowfall_flux',
        'm01s05i206': 'snowfall_flux',
        'm01s03i298': 'water_evaporation_flux',
        'm01s03i353': 'water_evaporation_flux',
        'm01s08i231': 'water_evaporation_flux',
        'm01s04i203': 'surface_runoff_flux',
        'm01s05i205': 'surface_runoff_flux',
        'm01s03i223': 'water_evaporation_flux',
        'm01s03i232': 'water_evaporation_flux',
        'm01s05i216': 'precipitation_flux'
    }

    # Initialize array of global budgets with zeros:
    fval = np.zeros(len(stash_f) + 1, dtype=np.float)

    # TODO use NaN
    # Set value of data-missing indicator:
    mdi = -10000.0

    # Read data and calculate global budgets

    # read data:
    try:
        for (i, stash) in enumerate(stash_f):

            ppfy = load_run_ss(run, 'annual', stash_dict[stash])

            # VPREDOI TODO
            # this fails unless using the correct variables
            # ppfy *= f_mult[i]

            # Calculate global mean, time series of global values
            ppfyg = area_average(
                ppfy, weighted=True, aggregator=iris.analysis.SUM)
            # multi-year mean of global values
            ppfy_total = ppfyg.collapsed('time', iris.analysis.MEAN)
            fval[i] = ppfy_total.data / 1.0e9  # TODO magic number

        fval[-1] = fval.sum()

    except iris.exceptions.ConstraintMismatchError:
        print("ERROR:  Missing data!!!")
        print("An MDI will be assigned to all water fluxes of this sub-model")
        fval[:] = mdi

    return fval


def fluxes_ocean_submodel(expid, mesh, areas, opath, wfpath, y_i, y_f):
    """
    Compute ocean fluxes

    Function to calculate time-average global water fluxes that go into
    NEMO ocean model
    """
    # load land-fraction-mask an# array with global values of Rd area cubes:
    mask = iris.load_cube(mesh)  # NEMO ocean-fraction mask
    area = iris.load_cube(areas)  # what NEMO  uses for grid cell areas

    # Initialize array of global budgets with zeros:
    r_arr = np.zeros(y_f - y_i + 1, dtype=np.float)  # global values of R
    net_arr = np.zeros(
        y_f - y_i + 1, dtype=np.float)  # global net flux into the ocean
    pme = np.zeros(y_f - y_i + 1, dtype=np.float)  # global P-E into the ocean

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
    filetemp = opath + '/' + expid + 'o_1y_{0}1201_{1}1130_grid_T.nc'
    nemofiles = [filetemp.format(x, x + 1) for x in range(y_i, y_f + 1, 1)]
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
            rflux = sorunoff * mask * area
            nflux = sowaflup * mask * area

            # TODO magic number
            # calculate global budgets:
            r_arr[i] = rflux.data.sum() / 1.0e9 - iceberg.data
            net_arr[i] = nflux.data.sum() / (-1.0e9)
            # P-E is obtained as a residual from net, R and iceberg-calving
            pme[i] = net_arr[i] - r_arr[i] - iceberg.data

        else:
            errf = 0

    if errf:
        fval = np.array([np.mean(pme), np.mean(r_arr),
                         iceberg.data, np.mean(net_arr)])
    else:
        print("ERROR:  Missing data!!!")
        print("An MDI value will be assigned to all ocean water fluxes")
        fval = np.array([mdi, mdi, mdi, mdi])

    return fval
