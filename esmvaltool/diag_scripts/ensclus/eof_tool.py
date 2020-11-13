"""Computing EOFs and PCs."""

import datetime

import matplotlib.pyplot as plt
import numpy as np

import cartopy.crs as ccrs
from eofs.standard import Eof


def eof_computation(var, lat):
    """Computing the EOFs and PCs.

    EOF analysis of a data array with spatial dimensions that
    represent latitude and longitude with weighting. In this example
    the data array is dimensioned (ntime, nlat, nlon), and in order
    for the latitude weights to be broadcastable to this shape, an
    extra length-1 dimension is added to the end
    """
    print('_________________________________________________________')
    print('Computing the EOFs and PCs')
    weights_array = np.sqrt(np.cos(np.deg2rad(lat)))[:, np.newaxis]

    start = datetime.datetime.now()
    solver = Eof(var, weights=weights_array)
    end = datetime.datetime.now()
    print('EOF computation took me %s seconds' % (end - start))

    # ALL VARIANCE FRACTIONS
    varfrac = solver.varianceFraction()
    # acc = np.cumsum(varfrac * 100)

    # ---------------------------------------PCs unscaled  (case 0 of scaling)
    pcs_unscal0 = solver.pcs()
    # ---------------------------------------EOFs unscaled  (case 0 of scaling)
    eofs_unscal0 = solver.eofs()

    # ---------------------------------------PCs scaled  (case 1 of scaling)
    pcs_scal1 = solver.pcs(pcscaling=1)

    # ---------------------------------------EOFs scaled (case 2 of scaling)
    eofs_scal2 = solver.eofs(eofscaling=2)

    return solver, pcs_scal1, eofs_scal2, pcs_unscal0, eofs_unscal0, varfrac


def eof_plots(neof, pcs_scal1, eofs_scal2, var, varunits, lat, lon,
              tit, numens, varfrac):
    """Plot of the nth the EOFs and PCs.

    Plot the PC scaled (divided by the square-root of their eigenvalues)
    in the selected domain
    """
    print('_________________________________________________________')
    print('Plotting the EOFs and PCs')
    print('Variable: {0} Units: {1}'.format(var, varunits))
    print('Ensemble members: {0}'.format(numens))

    # ------------------------------------------PCs scaled  (case 1 of scaling)
    figpc_scal1 = plt.figure(figsize=(24, 14))
    axes = figpc_scal1.gca()
    plt.plot(pcs_scal1[:, neof])
    plt.axhline(y=0, color='k', linestyle='--')
    tt_pc = '{0}   PC{1}: explained variance {2}%\n'\
        .format(tit, neof + 1, "%.2f" % (varfrac[neof] * 100))
    plt.title(tt_pc, fontsize=34, fontweight='bold')
    plt.grid(True)
    for tickx in axes.xaxis.get_major_ticks():
        tickx.label.set_fontsize(28)
    for ticky in axes.yaxis.get_major_ticks():
        ticky.label.set_fontsize(28)
    plt.ylabel('PC{0} {1}'.format(neof, varunits), fontsize=28)
    plt.xlabel('ensemble members', fontsize=28)

    # Plot the EOF scaled (multiplied by the square-root of their eigenvalues)
    # in the selected domain

    # ------------------------------------------EOFs scaled (case 2 of scaling)

    # rangecolorbar=np.arange(-180, 200, 20)
    figeof_scal2 = plt.figure(figsize=(14, 14))
    # ax = figeof_scal2.gca()
    proj = ccrs.PlateCarree()
    axes = plt.axes(projection=proj)
    axes.set_global()
    axes.coastlines()
    axes.gridlines()
    fill2 = axes.contourf(lon, lat, eofs_scal2[neof, ...], cmap=plt.cm.RdBu_r,
                          transform=ccrs.PlateCarree())
    cbar = plt.colorbar(fill2, orientation='horizontal')
    # cb.ax.set_position([0.9, 0.1, 0.001, 0.7])#([0.9, 0.1, 0.02, 0.8])
    cbar.set_label(varunits, rotation=0, fontsize=20)
    cbar.ax.tick_params(labelsize=20)
    tt_eof = '{0}\nEOF{1}: explained variance {2}%\n'\
        .format(tit, neof + 1, "%.2f" % (varfrac[neof] * 100))
    plt.title(tt_eof, fontsize=34, fontweight='bold')
    plt.tight_layout()

    return figpc_scal1, figeof_scal2
