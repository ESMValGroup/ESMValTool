"""
Autoassess Conservation

Module with routines to estimate conservation of water in all sub-models
in a GC configuration. Presently it calculates conservation as long-term
water fluxes across various sub-models, using annual mean data. It is
expected that in the future it will include exact calculations.
"""

import subprocess
import sys
import matplotlib.pyplot as plt

import iris

from . import global_water_budget as gwb
from .loaddata import load_run_ss
from .matplotlib_table import render_mpl_table


def resolution(cube):
    """
    Get resolution

    Get resolution for cube and whether it is on ENDGame grid
    Assume on full p grid
    Different algorithm required if longitude or latitude missing
    """
    # TODO: Extend to C grid variables?
    resol = "n" + str(cube.coord('longitude').points.size / 2)
    endgame = (cube.coord('latitude').points.size % 2 == 0)
    return resol, endgame


def global_freshwater_fluxes_over_various_GC_cubmodels(run):
    """
    Calculate fluxes

    Function to calculate long-term water fluxes in varios sub-models
    The conservation will be measured using as units 1e9 Kg/m2s~ Sv
    It is assumed that data comes from a model configuration of GA5
    and GO5 or higher.
    """
    metrics = dict()

    # Preliminaries to produce table with fluxes:
    expid = run['runid']
    fluxes_table = expid + '_global_freshwater_fluxes_table'
    # filename for storing fluxes table

    table = []
    val_format = '{:.3f}'
    table.append([expid, ''])
    table.append(['GLOBAL FRESHWATER FLUXES', '(1e9 Kg/s ~ Sv)'])

    # Location of various masks, ocean-data location, ocean-area values, etc.
    # ATMOSPHERE:
    # TODO: Get fields from experiment data

    # determine atmospheric horizontal model resolution:
    # VPREDOI
    # need the right data file
    # pptest = load_run_ss(run, 'seasonal', 'precipitation_flux')
    # m01s05i216: total precipitation rate; arbitrary cube
    pptest = load_run_ss(run, 'monthly', 'eastward_wind', lbproc=192)
    resol, endgame = resolution(pptest)

    # TODO local paths
    # land fraction mask:
    if endgame:
        lfpath = run['ancil_root'] + '/masks/qrparm.landfrac_' + resol + 'e.pp'
    else:
        lfpath = run['ancil_root'] + '/masks/qrparm.landfrac_' + resol + '.pp'

    # glacial mask:
    # TODO :  put mask in central directory
    if endgame:
        gmpath = run['ancil_root'] + '/conservation/glacialmask_' \
            + resol + '_endgame.pp'
    else:
        gmpath = run['ancil_root'] + '/conservation/glacialmask_' \
            + resol + '.pp'

    # NEMO:
    # Jan 2015:  At the present version, Maverick does not handle ocean fields.
    #            This means that presently, we cannot include water
    #            conservation in the ocean sub-model. The part of the code that
    #            calculates fluxes in this function is now commented out.

    # TODO : Include ocean water conservation

    # pp-file with constant iceberg-calving flux:
    #  (instead of a NEMO diagnostic)

    # TODO : put mask in central directory

    # wfpath = run['ancil_root']+'/conservation/qrclim.icecalve_'+resol+'.pp'

    # mesh file (it is dependent on GO version), for GO5:
    # mesh = run['ancil_root']+'/conservation/NEMOGO5_ocean_mask.nc'
    # areas = run['ancil_root']+'/conservation/NEMOGO5_ocean_area.nc'

    # load land-fraction mask and calculate ocean fraction mask:
    # VPREDOI
    # land fraction file is needed
    # /home/users/valeriu/base_masks_autoassess
    # lfm = iris.load_cube(lfpath)
    # ofm = -1.0 * lfm + 1.0
    lfm = pptest
    ofm = pptest

    # load glacial mask and obtain corresponding fraction masks:
    # VPREDOI
    # glacial mask is needed
    # /home/users/valeriu/base_masks_autoassess
    # gm = iris.load_cube(gmpath)
    gm = pptest
    # lfm_is = gm * lfm
    # ofm_is = gm * ofm
    lfm_is = pptest * pptest
    ofm_is = pptest * pptest

    # calculate water fluxes into the different sub-models

    # Template strings for writing information to file
    # hdr_temp = "{0:>39s}\n"
    # val_temp = "{0:>39s}  {1:7.3f} \n"
    # eql_temp = "{0:>48s}\n"

    # Calculate global WATER CONSERVATION IN VARIOUS SUB-MODELS

    # Ocean:

    # print 'Calculating fluxes into the ocean ... '
    # print

    # directory where ocean netcdf files are (run dependent):
    # opath = run['ss_ocean']

    # fval = gwb.fluxes_ocean_submodel(expid, mesh, areas, opath, wfpath,
    #                                  run['from_annual'].year,
    #                                  run['to_annual'].year)

    # name = 'global net water flux into ocean and sea ice'
    # metrics[name] = float(fval[-1])

    # with open('glob_freshwater_fluxes.tmp', "a") as table:
    #     table.write(hdr_temp.format('Ocean + Sea ice')
    #     table.write(hdr_temp.format('***************')
    #     table.write(val_temp.format('Precipitation minus evaporation',
    #                                 fval[0]))
    #     table.write(val_temp.format('River discharge', fval[1]))
    #     table.write(val_temp.format('Iceberg calving', fval[2]))
    #     table.write(eql_temp.format('******'))
    #     table.write(val_temp.format('Net Flux', fval[3]))
    #     table.write('\n\n')

    # Calculate TRIP fluxes

    # flbl = ['surface runoff', 'sub-surface runoff','inland basin runoff',
    #         'river outflow']
    stash_f = ['m01s08i234', 'm01s08i235', 'm01s08i245', 'm01s26i004']
    f_mult = [lfm, lfm, -1.0, -1.0]
    fval = gwb.fluxes_submodel(run, stash_f, f_mult)
    name = 'global net water flux TRIP'
    metrics[name] = float(fval[-1])

    runoff = fval[0] + fval[1]
    ibrunoff = fval[2]
    rdis_trip = fval[3]

    table.append(['TRIP', ''])
    table.append(['Surf.+sub-surf. runoff', val_format.format(runoff)])
    table.append(['River discharge', val_format.format(rdis_trip)])
    table.append(['Inland-basin runoff', val_format.format(ibrunoff)])
    table.append(['Net Flux', val_format.format(fval[4])])
    table.append(['', ''])

    # Calculate Ice sheet fluxes

    # flbl = ['l. scale snowfall', 'conv. snowfall', 'sublimation',
    #         'sublim. sea-ice', 'snow melt']
    stash_f = [
        'm01s04i204', 'm01s05i206', 'm01s03i298', 'm01s03i353', 'm01s08i231'
    ]
    f_mult = [lfm_is, lfm_is, -1.0 * gm, ofm_is, -1.0 * lfm_is]

    fval = gwb.fluxes_submodel(run, stash_f, f_mult)

    # There's no conservation metric for ice sheets.

    sfall_is = fval[0] + fval[1]
    sublim_is = fval[2] + fval[3]
    smelt_is = fval[4]
    net_is = fval[5]

    table.append(['Ice sheets', ''])
    table.append(['Snowfall', val_format.format(sfall_is)])
    table.append(['Sublimation', val_format.format(sublim_is)])
    table.append(['Snow melt', val_format.format(smelt_is)])
    table.append(['Net Flux', val_format.format(fval[5])])
    table.append(['', ''])

    # Calculate land-snow fluxes
    stash_f = [
        'm01s04i204', 'm01s05i206', 'm01s03i298', 'm01s03i353', 'm01s08i231'
    ]
    f_mult = [lfm, lfm, -1.0, ofm, -1.0 * lfm]

    fval = gwb.fluxes_submodel(run, stash_f, f_mult)
    name = 'global net water flux land snow'
    metrics[name] = float(fval[-1] - net_is)

    table.append(['Land Snow', ''])
    table.append(['Snowfall', val_format.format(fval[0] + fval[1] - sfall_is)])
    table.append(
        ['Sublimation',
         val_format.format(fval[2] + fval[3] - sublim_is)])
    table.append(['Snow melt', val_format.format(fval[4] - smelt_is)])
    table.append(['Net Flux', val_format.format(fval[5] - net_is)])
    table.append(['', ''])

    smelt_l = fval[4]

    # Calculate soil moisture fluxes
    # calculate rainfall minus evaporation and obtain the remaining fluxes from
    # results from sub-models previously calculated:
    stash_f = [
        'm01s04i203', 'm01s05i205', 'm01s03i298', 'm01s03i223', 'm01s03i232'
    ]
    f_mult = [lfm, lfm, 1.0, -1.0, ofm]
    fval = gwb.fluxes_submodel(run, stash_f, f_mult)

    name = 'global net water flux soil moisture'
    metrics[name] = float(fval[5] - smelt_l - runoff - ibrunoff)

    table.append(['Soil moisture', ''])
    table.append(['Rainfall minus evaporation', val_format.format(fval[5])])
    table.append(['Snow melt', val_format.format(-1.0 * smelt_l)])
    table.append(
        ['Surf. + sub-surv. runoff',
         val_format.format(-1.0 * runoff)])
    table.append(['Inland basin runoff', val_format.format(-1.0 * ibrunoff)])
    table.append(
        ['Net Flux',
         val_format.format(fval[5] - smelt_l - runoff - ibrunoff)])
    table.append(['', ''])

    # Calculate atmospheric fluxes
    stash_f = ['m01s03i223', 'm01s05i216']
    f_mult = [1.0, -1.0]

    fval = gwb.fluxes_submodel(run, stash_f, f_mult)
    name = 'global net water flux into atmosphere'
    metrics[name] = float(fval[-1])

    table.append(['Atmosphere', ''])
    table.append(['Evaporation', val_format.format(fval[0])])
    table.append(['Precipitation', val_format.format(fval[1])])
    table.append(['Net Flux', val_format.format(fval[2])])

    _fig = render_mpl_table(
        table,
        header_rows=2,
        header_columns=0,
        col_width=5,
        highlight_cells=[(2, 0), (6, 1), (8, 0), (12, 1), (14, 0), (18, 1),
                         (20, 0), (25, 1), (27, 0), (30, 1)])
    plt.savefig(fluxes_table + '.png', format='png')

    return metrics
