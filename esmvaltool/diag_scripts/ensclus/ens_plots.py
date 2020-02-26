"""Plot the chosen field for each ensemble."""

import math
import os

import matplotlib.pyplot as plt
import numpy as np
import cartopy.crs as ccrs

# User-defined libraries
from read_netcdf import read_n_2d_fields


def ens_plots(dir_output, dir_plot, name_outputs, numclus,
              field_to_plot, plot_type, season, area, extreme, numensmax):
    """Plot the chosen field for each ensemble."""
    print('Number of clusters: {0}'.format(numclus))

    varname = name_outputs.split("_")[0]
    kind = name_outputs.split("_")[-2]
    exp = name_outputs.split("_")[-1]

    # Reading the netCDF file of N 2Dfields of anomalies, saved by ens_anom.py
    namef = os.path.join(dir_output, 'ens_anomalies_{0}.nc'
                         .format(name_outputs))
    vartoplot, varunits, lat, lon = read_n_2d_fields(namef)
    print('vartoplot dim: (numens x lat x lon)={0}'.format(vartoplot.shape))
    numens = vartoplot.shape[0]

    # ____________Load labels
    namef = os.path.join(dir_output, 'labels_{0}.txt'.format(name_outputs))
    labels = np.loadtxt(namef, dtype=int)

    namef = os.path.join(dir_output, 'repr_ens_{0}.txt'.format(name_outputs))
    reprens = np.loadtxt(namef, dtype=int)

    namef = os.path.join(dir_output, 'legend_{0}.txt'.format(name_outputs))
    legends = np.loadtxt(namef, dtype=str)
    years = np.loadtxt(namef, dtype=str)[:, 6]
    legends = np.loadtxt(namef, dtype=str)[:, 1]

    vmi = round_down(np.nanpercentile(vartoplot, 0.1))
    vma = round_up(np.nanpercentile(vartoplot, 99.9))

    if field_to_plot == 'anomalies':
        # compute range colorbar for anomalies
        if abs(vmi) < abs(vma):
            rangecbarmin = -abs(vma)
            rangecbarmax = abs(vma)
        else:
            rangecbarmin = -abs(vmi)
            rangecbarmax = abs(vmi)
    else:
        # compute range colorbar for climatologies
        rangecbarmin = vmi
        rangecbarmax = vma

    delta = round_down((rangecbarmax - rangecbarmin) / 100)
    clevels = np.arange(rangecbarmin, rangecbarmax + delta, delta)
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'DarkOrange', 'grey']
    proj = ccrs.PlateCarree()
    if min(lon) < 180. < max(lon):
        clon = 180.
    else:
        clon = 0.

    numens0 = min(numens, numensmax)
    nfigs = int(np.ceil(numens / numens0))

    ypos = int(np.ceil(np.sqrt(numens0 / 2.0)))
    xpos = int(np.ceil(numens0 / ypos))
    fsize = int(min(max(4. / ypos * 15, 10), 32))

    namef_list = []
    for ifig in range(nfigs):
        fig = plt.figure(figsize=(24, 14))
        for iens in range(ifig * numens0, min((ifig + 1) * numens0, numens)):
            axes = plt.subplot(
                ypos, xpos, iens + 1 - ifig * numens0,
                projection=ccrs.PlateCarree(central_longitude=clon)
            )
            axes.set_extent([min(lon), max(lon), min(lat), max(lat)],
                            crs=ccrs.PlateCarree())
            axes.coastlines("110m")

            # Plot Data
            if field_to_plot == 'anomalies':
                map_plot = plt.contourf(lon, lat, vartoplot[iens], clevels,
                                        cmap=plt.cm.RdBu_r,
                                        transform=proj, extend='both')
            else:
                map_plot = plt.contourf(lon, lat, vartoplot[iens], clevels,
                                        transform=proj, extend='both')

            if iens in reprens:
                rect = plt.Rectangle((-0.01, -0.01), 1.02, 1.02, fill=False,
                                     transform=axes.transAxes, clip_on=False,
                                     zorder=10)
                rect.set_edgecolor(colors[labels[iens]])
                rect.set_linewidth(6.0)
                axes.add_artist(rect)

            # Add Title
            title_obj = plt.title(iens, fontsize=int(fsize * 1.8),
                                  fontweight='bold', loc='left')
            title_obj.set_backgroundcolor(colors[labels[iens]])
            title_obj = plt.title(legends[iens], fontsize=fsize, loc='right')

        cax = plt.axes([0.1, 0.03, 0.8, 0.03])  # horizontal
        cbar = plt.colorbar(map_plot, cax=cax, orientation='horizontal')
        cbar.ax.tick_params(labelsize=24)
        cbar.set_ticks(np.arange(rangecbarmin,
                                 rangecbarmax + delta, delta * 20))
        cbar.ax.set_ylabel(varname + '\n[' + varunits + ']', fontsize=24,
                           fontweight='bold', rotation='horizontal',
                           verticalalignment='center')
        cbar.ax.yaxis.set_label_position('right')
        cbar.ax.yaxis.set_label_coords(1.05, 1.4)

        plt.suptitle(field_to_plot.capitalize() + ' ' + varname + ' ' +
                     extreme + ' ' + area + ' ' + season + ' ' +
                     kind + ' ' + exp + ' ' + years[0],
                     fontsize=40, fontweight='bold')

        top = 0.89     # the top of the subplots of the figure
        bottom = 0.12  # the bottom of the subplots of the figure
        left = 0.02    # the left side of the subplots of the figure
        right = 0.98   # the right side of the subplots of the figure
        hspace = 0.36  # height reserved for white space between subplots
        wspace = 0.14  # width reserved for blank space between subplots
        plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top,
                            wspace=wspace, hspace=hspace)

        # plot the selected fields
        namef = os.path.join(dir_plot, ('{0}_{1}.fig{2}.' + plot_type)
                             .format(field_to_plot, name_outputs, ifig + 1))
        fig.savefig(namef)  # bbox_inches='tight')
        print('A ', plot_type, ' figure for the selected fields saved in {0}'
              .format(dir_plot))
        namef_list.append(namef)
    return namef_list


def round_up(x, sig=2):
    """Round up to a given number of significant digits."""
    dig = pow(10., sig - int(math.floor(math.log10(abs(x)))) - 1)
    return math.ceil(x * dig) / dig


def round_down(x, sig=2):
    """Round down to a given number of significant digits."""
    dig = pow(10., sig - int(math.floor(math.log10(abs(x)))) - 1)
    return math.floor(x * dig) / dig
