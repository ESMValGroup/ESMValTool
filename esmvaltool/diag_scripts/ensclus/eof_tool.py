"""Computing EOFs and PCs."""

# Standard packages
import datetime
import warnings
import numpy as np
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
from matplotlib.cbook import MatplotlibDeprecationWarning
from eofs.standard import Eof
warnings.simplefilter('ignore', MatplotlibDeprecationWarning)


def eof_computation(var, varunits, lat, lon):
    """Computing the EOFs and PCs.

    EOF analysis of a data array with spatial dimensions that
    represent latitude and longitude with weighting. In this example
    the data array is dimensioned (ntime, nlat, nlon), and in order
    for the latitude weights to be broadcastable to this shape, an
    extra length-1 dimension is added to the end
    """
    print('_________________________________________________________'
          '___________________________________________________________')
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
    print('___________________________________________________'
          '_________________________________________________________________')
    print('Plotting the EOFs and PCs')

    # ------------------------------------------PCs scaled  (case 1 of scaling)
    figPC_scal1 = plt.figure(figsize=(24, 14))
    ax = figPC_scal1.gca()
    plt.plot(pcs_scal1[:, neof])
    plt.axhline(y=0, color='k', linestyle='--')
    ttPC = '{0}   PC{1}: explained variance {2}%\n'\
        .format(tit, neof + 1, "%.2f" % (varfrac[neof] * 100))
    plt.title(ttPC, fontsize=34, fontweight='bold')
    plt.grid(True)
    for tickx in ax.xaxis.get_major_ticks():
        tickx.label.set_fontsize(28)
    for ticky in ax.yaxis.get_major_ticks():
        ticky.label.set_fontsize(28)
    plt.ylabel('PC{0} {1}'.format(neof, varunits), fontsize=28)
    plt.xlabel('ensemble members', fontsize=28)

    # Plot the EOF scaled (multiplied by the square-root of their eigenvalues)
    # in the selected domain

    # ------------------------------------------EOFs scaled (case 2 of scaling)

    # rangecolorbar=np.arange(-180, 200, 20)
    figEOF_scal2 = plt.figure(figsize=(14, 14))
    # ax = figEOF_scal2.gca()
    proj = ccrs.PlateCarree()
    ax = plt.axes(projection=proj)
    ax.set_global()
    ax.coastlines()
    ax.gridlines()
    fill2 = ax.contourf(lon, lat, eofs_scal2[neof, ...], cmap=plt.cm.RdBu_r,
                        transform=ccrs.PlateCarree())
    cb = plt.colorbar(fill2, orientation='horizontal')
    # cb.ax.set_position([0.9, 0.1, 0.001, 0.7])#([0.9, 0.1, 0.02, 0.8])
    cb.set_label(varunits, rotation=0, fontsize=20)
    cb.ax.tick_params(labelsize=20)
    ttEOF = '{0}\nEOF{1}: explained variance {2}%\n'\
        .format(tit, neof + 1, "%.2f" % (varfrac[neof] * 100))
    plt.title(ttEOF, fontsize=34, fontweight='bold')
    plt.tight_layout()

    return figPC_scal1, figEOF_scal2
