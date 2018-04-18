import numpy as np
import matplotlib as mpl
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches
import matplotlib.colors as colors
import iris
import iris.quickplot as qplt
import iris.plot as iplt
import os, sys
from mycolormaps import getcolors
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
mpl.rc('font', family='Times New Roman')

def MapPlot(var2plot, title=None, clevs=None, pltname='iris_test_plot.ps', colorReverse=False):
    # Plot Lat-Lon maps
    fig = plt.figure(figsize=(8, 3), dpi=100)
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    if colorReverse:
        cf = iplt.contourf(var2plot, clevs, cmap=getcolors('WhBlGrYeRe', colorReverse=True), extend='both')
    else:
        cf = iplt.contourf(var2plot, clevs, cmap=getcolors('WhBlGrYeRe'), extend='both')

    #iplt.pcolormesh(var2plot, cmap=getcolors('WhBlGrYeRe') )
    ax.coastlines()
    ax.gridlines()
    #ax.set_global()
    plt.title(title)

    ax.set_xticks([0, 60, 120, 180, 240, 300, 360], crs=ccrs.PlateCarree())
    ax.set_yticks([-30, 0, 30], crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(zero_direction_label=True)
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)

    labels = ax.get_xticklabels()
    plt.setp(labels, rotation=0, fontsize=10)
    labels = ax.get_yticklabels()
    plt.setp(labels, rotation=0, fontsize=10)

    plt_ax = plt.gca()
    left, bottom, width, height = plt_ax.get_position().bounds
    first_plot_left = plt_ax.get_position().bounds[0]

    # the width of the colorbar should now be simple
    width = left - first_plot_left + width * 0.9

    # Add axes to the figure, to place the colour bar
    colorbar_axes = fig.add_axes([first_plot_left + 0.0375, bottom + 0.1, width, 0.05])
    # Add the colour bar
    cbar = plt.colorbar(cf, colorbar_axes, orientation='horizontal')
    cbar.ax.tick_params(labelsize=10)

    plt.savefig(pltname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def HovTimeLon(var, time, lons, levels=None, title=None, figname=None):
    L, T = np.meshgrid(lons, time)

    cmap = getcolors('ncl_default')
    norm = colors.BoundaryNorm(levels, len(cmap.colors))

    fig = plt.figure()
    ax = fig.add_subplot(111)
    CS = plt.contourf(L, T, var, levels=levels,
                      cmap=cmap, norm=norm)
    plt.colorbar(CS, orientation='horizontal')
    CS = plt.contour(L, T, var, levels=levels, colors='k')
    #ax.set_aspect('equal')
    # set axes range
    plt.xlim(min(lons), max(lons))
    plt.ylim(min(time), max(time))
    plt.title(title)#
    plt.xlabel('Longitude (degrees east)')
    plt.ylabel('Lag/lead (days)')
    plt.grid()
    plt.savefig(figname)#, bbox_inches='tight')
    #plt.show()
    plt.close()

def PlotSpecRaw(spec, freq, wave, levels=np.linspace(-2, 2, 21), \
             title='', figname='spec_test.ps'):
    plt.clf()
    minfrq4plt = 0.
    maxfrq4plt = 0.8
    minwav4plt = -15
    maxwav4plt = 15

    minfrq = minfrq4plt
    maxfrq = min([maxfrq4plt, max(freq)])
    F, W = np.meshgrid(wave, freq)

    cmap = getcolors('ncl_default')
    norm = colors.BoundaryNorm(levels, len(cmap.colors))
    norm = colors.Normalize(clip=False)
    CS = plt.contourf(F, W, spec, levels=levels,
                      cmap=cmap, norm=norm, extend='both')
    bar = plt.colorbar(CS)

    tick_locs = levels
    tick_labels = levels
    bar.locator = ticker.FixedLocator(tick_locs)
    bar.formatter = ticker.FixedFormatter(tick_labels)
    bar.update_ticks()

    # set axes range
    plt.xlim(minwav4plt, maxwav4plt)
    plt.ylim(minfrq, maxfrq)

    # Lines
    plt.plot([0, 0], [0, 0.5], 'k--', lw=0.5)

    # Line markers of periods
    frqs = [80, 30, 6, 3]
    for frq in frqs:
        plt.plot([-15, 15], [1. / frq, 1. / frq], 'k--', lw=0.5)
        plt.text(-14.7, 1. / frq, str(frq) + ' days', {'color' : 'k'})

    plt.title(title)#
    plt.xlabel('Westward     Zonal Wave Number     Eastward')
    plt.ylabel('Frequency (CPD)')

    plt.savefig(figname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def plotAntiSymmetric(spec, freq, wave, Apzwn, \
         Afreq, levels=np.linspace(-2, 2, 21), \
         title='', figname='specAntiSym_test.ps'):
    plt.clf()
    minfrq4plt = 0.
    maxfrq4plt = 0.8
    minwav4plt = -15
    maxwav4plt = 15

    minfrq = minfrq4plt
    maxfrq = min([maxfrq4plt, max(freq)])
    F, W = np.meshgrid(wave, freq)

    cmap = getcolors('ncl_default')
    norm = colors.BoundaryNorm(levels, len(cmap.colors))

    CS = plt.contourf(F, W, spec, levels=levels,
                      cmap=cmap, norm=norm, extend='both')
    bar = plt.colorbar(CS)

    tick_locs = levels
    tick_labels = levels
    bar.locator = ticker.FixedLocator(tick_locs)
    bar.formatter = ticker.FixedFormatter(tick_labels)
    bar.update_ticks()

    # set axes range
    plt.xlim(minwav4plt, maxwav4plt)
    plt.ylim(minfrq, maxfrq)

    # Lines
    plt.plot([0, 0], [0, 0.5], 'k--', lw=0.5)

    # Line markers of periods
    frqs = [80, 30, 6, 3]
    for frq in frqs:
        plt.plot([-15, 15], [1. / frq, 1. / frq], 'k--', lw=0.5)
        plt.text(-14.7, 1. / frq, str(frq) + ' days', {'color' : 'k'})

    plt.title(title)#
    plt.xlabel('Westward     Zonal Wave Number     Eastward')
    plt.ylabel('Frequency (CPD)')

    # Symmetric waves
    # Equatorial Rossby
    for i in range(3):
        for j in range(3):
            plt.plot(Apzwn[i, j, :], Afreq[i, j, :], 'k', lw=0.5)

    plt.text(-10., 0.15, "MRG", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-3., 0.58, "n=2 IG", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(6., 0.4, "n=0 EIG", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-3., 0.475, "h=12", {'color' : 'k', 'backgroundcolor':'w'})

    plt.savefig(figname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def plotSymmetric(spec, freq, wave, Apzwn, \
         Afreq, levels=np.linspace(-2, 2, 21), \
         title='', figname='specSym_test.ps'):

    plt.clf()
    minfrq4plt = 0.
    maxfrq4plt = 0.8
    minwav4plt = -15
    maxwav4plt = 15

    minfrq = minfrq4plt
    maxfrq = min([maxfrq4plt, max(freq)])
    F, W = np.meshgrid(wave, freq)

    cmap = getcolors('ncl_default')
    norm = colors.BoundaryNorm(levels, len(cmap.colors))

    CS = plt.contourf(F, W, spec, levels=levels,
                      cmap=cmap, norm=norm, extend='both')
    bar = plt.colorbar(CS)

    tick_locs = levels
    tick_labels = levels
    bar.locator = ticker.FixedLocator(tick_locs)
    bar.formatter = ticker.FixedFormatter(tick_labels)
    bar.update_ticks()

    # set axes range
    plt.xlim(minwav4plt, maxwav4plt)
    plt.ylim(minfrq, maxfrq)

    # Lines
    plt.plot([0, 0], [0, 0.5], 'k--', lw=0.5)

    # Line markers of periods
    frqs = [80, 30, 6, 3]
    for frq in frqs:
        plt.plot([-15, 15], [1. / frq, 1. / frq], 'k--', lw=0.5)
        plt.text(-14.7, 1. / frq, str(frq) + ' days', {'color' : 'k'})

    plt.title(title)#
    plt.xlabel('Westward     Zonal Wave Number     Eastward')
    plt.ylabel('Frequency (CPD)')

    # Symmetric waves
    # Equatorial Rossby
    for i in range(3, 6):
        for j in range(3):
            plt.plot(Apzwn[i, j, :], Afreq[i, j, :], 'k', lw=0.5)

    plt.text(11.5, 0.4, "Kelvin", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-10.7, 0.07, "n=1 ER", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-3., 0.45, "n=1 IG", {'color' : 'k', 'backgroundcolor':'w'})
    plt.text(-14., 0.46, "h=12", {'color' : 'k', 'backgroundcolor':'w'})

    plt.savefig(figname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def mjo_wavenum_freq_season_plot(pow_cube, levels=np.arange(0, 3, 0.2), \
             title='', figname='wavenum_freq_season_plot.ps'):
    NW = 6
    fMin = -0.05
    fMax = 0.05

    constraint = iris.Constraint(frequency=lambda cell: fMin <= cell <= fMax, \
                                 wavenumber=lambda cell: 0 <= cell <= NW)
    pow_cube = pow_cube.extract(constraint)

    freq = pow_cube.coord('frequency').points
    wave = pow_cube.coord('wavenumber').points

    W, F = np.meshgrid(freq, wave)

    # set zeroth frequency to minimum value
    izero = np.where(freq == 0)[0][0]
    pow_cube.data[:,izero] = np.min(pow_cube.data) # 0th freq

    CS = plt.contourf(W, F, pow_cube.data, levels=levels,
                      cmap=getcolors('WhBlGrYeRe'), extend="both")
    plt.colorbar(CS)

    # set axes range
    plt.ylim(0, NW)
    plt.xlim(fMin, fMax)

    # Line markers of periods
    frqs = [80, 30]
    for frq in frqs:
        plt.plot([1. / frq, 1. / frq], [-0, 15], 'k--', lw=0.5)
        plt.text(1. / frq, 5.5, str(frq) + 'd', {'color' : 'k'})

    plt.plot([0, 0], [-0, 15], 'k:', lw=0.25)
    plt.title(title)  #
    plt.xlabel('Westward     Frequency     Eastward')
    plt.ylabel('Zonal wavenumber')
    plt.savefig(figname)#,bbox_inches='tight')
    #plt.show()
    plt.close()

def rmm_plot(rmmfile, figname):

    x = 4
    y = 0.707107
    linewidth = 0.25
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal', axisbg='lightgrey')
    plt.plot(np.arange(10, 11))
    plt.plot([-x, -y], [-x, -y], 'k', lw=linewidth)
    plt.plot([y, x], [y, x], 'k', lw=linewidth)
    plt.plot([-x, -y], [x, y], 'k', lw=linewidth)
    plt.plot([y, x], [-y, -x], 'k', lw=linewidth)
    plt.plot([-x, -1], [0, 0], 'k', lw=linewidth)
    plt.plot([1, x], [0, 0], 'k', lw=linewidth)
    plt.plot([0, 0], [-x, -1], 'k', lw=linewidth)
    plt.plot([0, 0], [1, x], 'k', lw=linewidth)

    c = mpatches.Circle((0, 0), 1, fc="lightgrey", ec="k", lw=linewidth)
    ax.add_patch(c)

    plt.text(-3.75, -2, 'Phase 1', va='center', rotation=90, fontsize='smaller')
    plt.text(-2, -3.75, 'Phase 2', ha='center', fontsize='smaller')
    plt.text(2, -3.75, 'Phase 3', ha='center', fontsize='smaller')
    plt.text(3.75, -2, 'Phase 4', ha='center', va='center', rotation=90, fontsize='smaller')
    plt.text(3.75, 2, 'Phase 5', ha='center', va='center', rotation=90, fontsize='smaller')
    plt.text(2, 3.75, 'Phase 6', ha='center', va='center', fontsize='smaller')
    plt.text(-2, 3.75, 'Phase 7', ha='center', va='center', fontsize='smaller')
    plt.text(-3.75, 2, 'Phase 8', va='center', rotation=90, fontsize='smaller')

    plt.text(-4.75, 0, 'Western Hem, Africa', va='center', rotation=90)
    plt.text(0, -4.75, 'Indian Ocean', ha='center', rotation=0)
    plt.text(4.5, 0, 'Maritime continent', ha='center', va='center', rotation=90)
    plt.text(0, 4.5, 'Western Pacific', va='center', ha='center', rotation=0)

    plt.ylim(-x, x)
    plt.xlim(-x, x)

    if os.path.isfile(rmmfile):
        C = np.loadtxt(rmmfile)
        year = C[:, 0]
        pc1 = C[:, 3]
        pc2 = C[:, 4]
        pha = C[:, 6]

        inds = np.where(pc1 != -999.)[0]
        plt.plot(pc1[inds], pc2[inds], color='blue', alpha=0.15, marker='o',
                 markersize=3)


    plt.savefig(figname)#,bbox_inches='tight')
    #plt.title('MJO Phase space')#
    #plt.xlabel('RMM1')
    #plt.ylabel('RMM2')
    #plt.show()
    plt.close()


def phasefreq_plot(rmmfile, figname):
    if os.path.isfile(rmmfile):
        C = np.loadtxt(rmmfile)
        year = C[:, 0]
        month = C[:, 1]
        pc1 = C[:, 3]
        pc2 = C[:, 4]
        pha = C[:, 5]
        amp = C[:, 6]

        freq = np.zeros((2, 9)) # winter/summer, 0 - all phases, 1:9 - 8 phases
        for phase in range(1, 9):
            # Summer
            months = [5, 6, 7, 8, 9, 10]
            freq[0, phase] = len(
                [i for i in range(len(pha))
                    if amp[i] > 1.0 and pha[i] == phase and month[i] in months]
            )
            freq[0, 0] = len(
                [i for i in range(len(pha))
                    if amp[i] <= 1.0 and month[i] in months]
            )

            # Winter
            months = [11, 12, 1, 2, 3, 4]
            freq[1, phase] = len([i for i in range(len(pha)) if amp[i] > 1.0 and
                pha[i] == phase and month[i] in months])
            freq[1, 0] = len(
                [i for i in range(len(pha))
                    if amp[i] <= 1.0 and month[i] in months]
            )

        freq = 100. * freq / freq.sum()
        phases = np.arange(0, 9, 1)
        months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep',
                  'Oct', 'Nov', 'Dec']

        fig = plt.figure(1)
        ax = fig.add_subplot(111, aspect=0.9, axisbg='lightgrey')
        plot1 = plt.bar(phases[1:] - 0.35, freq[0, 1:], color='w', width=0.35)
        plot2 = plt.bar(phases[1:] + 0.0, freq[1, 1:], color='g', width=0.35)
        plt.text(1.5, 4.25, 'No MJO summer : ' + str(round(freq[0, 0])))
        plt.text(1.5, 4., 'No MJO winter : ' + str(round(freq[1, 0])))
        ax.yaxis.grid(color='lightgray', linestyle='-')
        plt.ylim([0, 5])
        plt.legend((plot1[0], plot2[0]), ('Summer', 'Winter'))
        plt.xlabel('MJO Phase')
        plt.ylabel('% frequency')
        plt.title('MJO Phase Frequency')
        plt.savefig(figname, bbox_inches='tight')
        plt.close()
