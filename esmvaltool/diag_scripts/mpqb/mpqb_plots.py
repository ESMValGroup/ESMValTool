
import iris
import matplotlib
import matplotlib.pyplot as plt
import copy
import yaml
import os
from esmvaltool.diag_scripts.shared._base import get_diagnostic_filename


def get_ecv_plot_config(ecv_name):
    cfg_filename = os.path.join(os.path.split(__file__)[0],'mpqb_cfg.yml')
    with open(cfg_filename,'r') as handle:
        plot_config = yaml.safe_load(handle)
    ecv_plot_config = plot_config[ecv_name]
    return ecv_plot_config


def mpqb_mapplot(cube,filename,**plotkwargs):
    plottitle = plotkwargs.pop('title')
    fig = plt.figure(dpi=200)
    ax = fig.add_subplot(projection=iris.plot.default_projection(cube))
    # replace the cmap key with the cmap object, and add grey shading for masked values
    cmapname = plotkwargs.pop('cmap')
    cmap = matplotlib.cm.get_cmap(cmapname)
    cmap.set_bad("grey", 0.4)
    plotkwargs['cmap'] = cmap
    plotkwargs['antialiased'] = True
    # Take out small grid lines like this
    plotkwargs['linewidth'] = 0
    iris.quickplot.pcolormesh(cube,**plotkwargs)
    plt.gca().coastlines()
    plt.title(plottitle)
    fig.savefig(filename)
    plt.close(fig)
