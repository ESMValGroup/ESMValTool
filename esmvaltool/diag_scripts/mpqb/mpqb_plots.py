
import iris
import matplotlib
import matplotlib.pyplot as plt
import copy
import datetime
import yaml
import os
from esmvaltool.diag_scripts.shared._base import get_diagnostic_filename
import numpy as np


dataset_plotnames = {
  'ERA-Interim-Land' : 'ERA-Interim/Land',
  'CDS-SATELLITE-SOIL-MOISTURE' : 'ESA-CCI',
  'cds-era5-land-monthly' : 'ERA5-Land',
  'cds-era5-monthly' : 'ERA5',
  'MERRA2' : 'MERRA-2',
  'cds-satellite-lai-fapar' : 'SPOT-VGT',
  'CDS-SATELLITE-ALBEDO' : 'SPOT-VGT',
}



def get_ecv_plot_config(ecv_name):
    cfg_filename = os.path.join(os.path.split(__file__)[0], 'mpqb_cfg.yml')
    with open(cfg_filename,'r') as handle:
        plot_config = yaml.safe_load(handle)
    ecv_plot_config = plot_config[ecv_name]
    return ecv_plot_config

def get_plottitle_timeperiod(cube):
    time_coord = cube.coord('time')
    time_asdatetime = time_coord.units.num2date(time_coord.bounds).flatten()
    starttime = min(time_asdatetime) + datetime.timedelta(days=1)
    endtime = max(time_asdatetime) - datetime.timedelta(days=1)
    start, end = starttime.year, endtime.year
    return f"{start} - {end}"



def mpqb_mapplot(cube,filename,**plotkwargs):
    from matplotlib.ticker import MaxNLocator
    from matplotlib.colors import BoundaryNorm

    plottitle = dataset_plotnames[plotkwargs.pop('title')]

    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(projection=iris.plot.default_projection(cube))

    # replace the cmap key with the cmap object, and add grey shading for masked values
    cmapname = plotkwargs.pop('cmap')
    cmap = matplotlib.cm.get_cmap(cmapname)
    cmap.set_bad("grey", 0.1)

    from matplotlib.colors import from_levels_and_colors

    if np.abs(plotkwargs['vmin'])==np.abs(plotkwargs['vmax']):
        # Diverging colorbar centred around zero
        discrete=True
        N=7 # number of colors
        if discrete:
            color_list = cmap(np.linspace(0, 1, N))
            cmap_name = cmap.name + str(N)
            cmap = cmap.from_list(cmap_name, color_list, N)

        levels = MaxNLocator(nbins=cmap.N+1, symmetric=True).tick_values(plotkwargs['vmin'],plotkwargs['vmax'])
        # Remove zero from levels
        levels = np.delete(levels, len(levels)/2)

        cmap, norm = from_levels_and_colors(levels, color_list)
        cmap.set_bad("grey", 0.1)
        plotkwargs['norm'] = norm

    plotkwargs['rasterized'] = True
    plotkwargs['cmap'] = cmap

    pcols = iris.plot.pcolormesh(cube, **plotkwargs)
    # Take out small grid lines like this
    pcols.set_edgecolor('face')
    plt.gca().coastlines()

    # Colorbar
    cb = plt.colorbar(pcols, orientation='horizontal')
    cb.set_label(cube.units)
    cb.ax.tick_params(labelsize=8)

    # Add timeperiod to plot title
    timeperiod = get_plottitle_timeperiod(cube)
    plottitle = f"{plottitle} {timeperiod}"
    plt.title(plottitle)
    fig.savefig(filename, bbox_inches='tight')
    plt.close(fig)
