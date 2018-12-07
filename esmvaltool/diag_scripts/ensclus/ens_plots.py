"""Plot the chosen field for each ensemble."""

import math
import os

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

# User-defined libraries
from read_netcdf import read_n_2d_fields


def ens_plots(dir_output, dir_plot, name_outputs, numclus,
              field_to_plot, plot_type):
    """Plot the chosen field for each ensemble."""
    tit = field_to_plot
    print('Number of clusters: {0}'.format(numclus))

    print(name_outputs)
    varname = name_outputs.split("_")[0]
    kind = name_outputs.split("_")[-2]
    exp = name_outputs.split("_")[-1]

    # Reading the netCDF file of N 2Dfields of anomalies, saved by ens_anom.py
    ifile = os.path.join(dir_output, 'ens_anomalies_{0}.nc'
                         .format(name_outputs))
    vartoplot, varunits, lat, lon = read_n_2d_fields(ifile)
    print('vartoplot dim: (numens x lat x lon)={0}'.format(vartoplot.shape))
    numens = vartoplot.shape[0]

    # ____________Load labels
    namef = os.path.join(dir_output, 'labels_{0}.txt'.format(name_outputs))
    labels = np.loadtxt(namef, dtype=int)
    print(labels)

    vmi = vartoplot.min()
    vma = vartoplot.max()

    if field_to_plot == 'anomalies':
        # compute range colorbar for anomalies
        delta = 0.05
        if abs(math.floor(vmi * 100) / 100) < math.ceil(vma * 100) / 100:
            rangecbarmin = -math.ceil(vma * 100) / 100
            rangecbarmax = math.ceil(vma * 100) / 100 + delta
        else:
            rangecbarmin = math.floor(vmi * 100) / 100
            rangecbarmax = abs(math.floor(vmi * 100) / 100) + delta
    else:
        # compute range colorbar for climatologies
        delta = 0.2
        rangecbarmin = math.floor(vmi)
        rangecbarmax = math.ceil(vma) + delta

    clevels = np.arange(rangecbarmin, rangecbarmax, delta)
    # clevels=np.arange(2,44,delta)
    # clevels=np.arange(-0.7,0.75,delta)
    # clevels=np.arange(0,6.2,delta)

    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'DarkOrange', 'grey']

    proj = ccrs.PlateCarree()

    x = int(np.ceil(np.sqrt(numens * 1.6)))
    y = int(np.ceil(numens / x))
    print(x, y)
    fig = plt.figure(figsize=(24, 14))
    for nens in range(numens):
        # print('//////////ENSEMBLE MEMBER {0}'.format(nens))
        axes = plt.subplot(x, y, nens + 1, projection=ccrs.PlateCarree())
        axes.coastlines("110m")
        # ax.set_extent([-10, 60, -30,90],ccrs.PlateCarree())

        # Plot Data
        if field_to_plot == 'anomalies':
            map_plot = plt.contourf(lon, lat, vartoplot[nens], clevels,
                                    cmap=plt.cm.RdBu_r,
                                    transform=proj)
        else:
            map_plot = plt.contourf(lon, lat, vartoplot[nens], clevels,
                                    transform=proj)
        # print('min={0}'.format(vartoplot[nens].min()))
        # print('max={0}\n'.format(vartoplot[nens].max()))

        # Add Title
        subtit = nens
        title_obj = plt.title(subtit, fontsize=32, fontweight='bold')
        for nclus in range(numclus):
            if nens in np.where(labels == nclus)[0]:
                title_obj.set_backgroundcolor(colors[nclus])

    cax = plt.axes([0.1, 0.03, 0.8, 0.03])  # horizontal
    cbar = plt.colorbar(map_plot, cax=cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=18)

    plt.suptitle(exp + ' ' + kind + ' ' + varname + ' ' + tit + ' (' +
                 varunits + ')', fontsize=45, fontweight='bold')

    plt.subplots_adjust(top=0.85)
    top = 0.89     # the top of the subplots of the figure
    bottom = 0.12  # the bottom of the subplots of the figure
    left = 0.02    # the left side of the subplots of the figure
    right = 0.98   # the right side of the subplots of the figure
    hspace = 0.36  # amount of height reserved for white space between subplots
    wspace = 0.14  # amount of width reserved for blank space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                        wspace=wspace, hspace=hspace)

    # plot the selected fields
    namef = os.path.join(dir_plot, ('{0}_{1}.' + plot_type)
                         .format(field_to_plot, name_outputs))
    fig.savefig(namef)  # bbox_inches='tight')
    print('A ', plot_type, ' figure for the selected fields is saved in {0}'
          .format(dir_plot))
    print('___________________________________________________'
          '_________________________________________________________________')

    return
