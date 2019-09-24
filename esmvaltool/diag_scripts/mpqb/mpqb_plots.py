
import iris
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
    fig = plt.figure()
    ax = fig.add_subplot(projection=iris.plot.default_projection(cube))
    iris.quickplot.pcolormesh(cube,**plotkwargs)
    plt.gca().coastlines()
    plt.title(plottitle)
    fig.savefig(filename)
    plt.close(fig)
