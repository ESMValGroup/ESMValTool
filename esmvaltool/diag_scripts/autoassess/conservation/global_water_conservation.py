"""
Autoassess Conservation GC Water Conservation Module.

Module with routines to estimate conservation of water in all sub-models
in a GC configuration. Presently it calculates conservation as long-term
water fluxes across various sub-models, using annual mean data. It is
expected that in the future it will include exact calculations.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from . import global_water_budget as gwb
# from esmvaltool.diag_scripts.autoassess.loaddata import load_run_ss
from .matplotlib_table import render_mpl_table


def resolution(cube):
    """
    Get data resolution.

    Get resolution for cube and whether it is on ENDGame grid
    Assume on full p grid
    Different algorithm required if longitude or latitude missing
    """
    resol = "n" + str(cube.coord('longitude').points.size / 2)
    endgame = (cube.coord('latitude').points.size % 2 == 0)
    return resol, endgame


def global_freshwater_fluxes(run):
    """
    Calculate fluxes.

    Function to calculate long-term water fluxes in varios sub-models
    The conservation will be measured using as units 1e9 Kg/m2s~ Sv
    It is assumed that data comes from a model configuration of GA5
    and GO5 or higher.
    """
    metrics = dict()

    # Preliminaries to produce table with fluxes:
    expid = run['runid']
    fluxes_table = expid + '_global_freshwater_fluxes_table'
    table = []
    val_format = '{:.3f}'
    table.append([expid, ''])
    table.append(['GLOBAL FRESHWATER FLUXES', '(1e9 Kg/s ~ Sv)'])

    # determine atmospheric horizontal model resolution:
    # pptest = load_run_ss(run, 'seasonal', 'precipitation_flux')
    # m01s05i216: total precipitation rate; arbitrary cube
    # resol, endgame = resolution(pptest)

    # land fraction mask:
    mask_dir = os.path.join(
        os.path.dirname(os.path.dirname(__file__)), 'autoassess_source')
    # lfpath = os.path.join(mask_dir, 'qrparm.landfrac_' + resol + '.pp')
    # gmpath = os.path.join(mask_dir, 'glacialmask_' + resol + '.pp')

    # use generic landsea.nc file instead of dedicated masks
    # TODO: replace with dedicated mask files when available
    mask_file = os.path.join(mask_dir, 'landsea.nc')
    import iris
    masks = iris.load_cube(mask_file)
    # lfm = iris.load_cube(lfpath)
    # gm = iris.load_cube(gmpath)
    lfm = masks
    gm = masks
    lfm.data = np.ma.masked_array(lfm.data, mask=(lfm.data != 1.))
    gm.data = np.ma.masked_array(gm.data, mask=(gm.data != 1.))
    ofm = -1.0 * lfm + 1.0
    lfm_is = gm * lfm
    ofm_is = gm * ofm

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

    header_rows = 2
    header_columns = 0
    col_width = 5
    highlight_cells = [(2, 0), (6, 1), (8, 0), (12, 1), (14, 0), (18, 1),
                       (20, 0), (25, 1), (27, 0), (30, 1)]
    render_mpl_table(table, header_rows, header_columns,
                     col_width, highlight_cells)
    plt.savefig(fluxes_table + '.png', format='png')

    return metrics
