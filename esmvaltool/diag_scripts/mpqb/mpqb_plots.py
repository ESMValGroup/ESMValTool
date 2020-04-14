
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
  'missing': 'missing',
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
    starttime = min(time_asdatetime) + datetime.timedelta(seconds=1)
    endtime = max(time_asdatetime) - datetime.timedelta(seconds=1)
    start, end = starttime.year, endtime.year
    return f"{start} - {end}"

def mpqb_mapplot(cube,filename,**plotkwargs):
    plottitle = dataset_plotnames.pop(plotkwargs.pop('title','missing'),'missing')
    
    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(projection=iris.plot.default_projection(cube))
    # replace the cmap key with the cmap object, and add grey shading for masked values
    cmapname = plotkwargs.pop('cmap')
    cmap = matplotlib.cm.get_cmap(cmapname)
    cmap.set_bad("grey", 0.1)
    plotkwargs['cmap'] = cmap
    plotkwargs['rasterized'] = True
    # Take out small grid lines like this
    pcols = iris.quickplot.pcolormesh(cube,**plotkwargs)
    pcols.set_edgecolor('face')
    plt.gca().coastlines()
    # Add timeperiod to plot title
    timeperiod = get_plottitle_timeperiod(cube)
    plottitle = f"{plottitle} {timeperiod}"
    plt.title(plottitle)
    fig.savefig(filename, bbox_inches='tight')
    plt.close(fig)
