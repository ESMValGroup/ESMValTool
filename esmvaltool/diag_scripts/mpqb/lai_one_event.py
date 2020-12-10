
import logging
import os

from esmvalcore import preprocessor as pp
import cf_units
import iris
import numpy as np
import yaml
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import from_levels_and_colors


from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic
from esmvaltool.diag_scripts.shared._base import get_plot_filename
from esmvaltool.diag_scripts.shared.trend_mpqb_common.diag1d import (
    mannkendall1d, theilslopes1d)
from esmvaltool.diag_scripts.shared.trend_mpqb_common.sharedutils import \
    parallel_apply_along_axis
from mpqb_plots import mpqb_mapplot
from mpqb_utils import get_mpqb_cfg

logger = logging.getLogger(os.path.basename(__file__))


def global_mapplot(cube, dataset_cfg, cfg):
    fig = plt.figure(dpi=200)
    fig.add_subplot(projection=iris.plot.default_projection(cube))

    levels = [-2.5,-1.5,-.5,-.2,.2,.5,1.5,2.5]
    nbins = len(levels)
    cmap = matplotlib.cm.get_cmap('BrBG')
    color_list = cmap(np.linspace(0, 1, nbins+1))
    cmap_name = cmap.name + str(nbins)
    cmap = cmap.from_list(cmap_name, color_list, len(levels))
    cmap, norm = from_levels_and_colors(levels, color_list, extend='both')
    cmap.set_bad('grey', 0.1)

    pcols = iris.plot.pcolormesh(cube, cmap=cmap, norm=norm)
    # Take out small grid lines like this
    pcols.set_edgecolor('face')
    plt.gca().coastlines()

    # Colorbar
    colorbar = plt.colorbar(pcols, orientation='horizontal', extend='both')
    colorbar.set_label(cube.units)
    colorbar.ax.tick_params(labelsize=8)


    plt.gca().margins(0)

    plt.title(f"JJA 2003")
    baseplotname = f"lai_event2003_" \
               f"_{dataset_cfg['variable_group']}" \
               f"_{dataset_cfg['start_year']}-" \
               f"{dataset_cfg['end_year']}"



    filename = get_plot_filename(baseplotname, cfg)
    fig.savefig(filename, bbox_inches='tight')
    plt.close(fig)


def main(cfg):

    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    grouped_input_data = group_metadata(input_data, 'alias', sort='alias')

    # Loop through all datasets
    for alias in grouped_input_data.keys():
        dataset_cfg = grouped_input_data[alias][0]
        dataset = grouped_input_data[alias][0]['dataset']

        logger.info("Opening dataset: %s", dataset)
        # Opening the pair
        cube = iris.load_cube(dataset_cfg['filename'])

        # Cut down to 2003 event
        cube = pp.extract_time(cube, 2003, 6, 10, 2003, 8, 20)

        # Calc mean over event.
        cube = cube.collapsed('time', iris.analysis.MIN)

        # Plot the results (if configured to plot)
        if cfg['write_plots']:
            global_mapplot(cube, dataset_cfg, cfg)
        else:
            logger.warning("This script wants to create a plot, but isn't allowed to.")
    logger.info("Finished!")




if __name__ == '__main__':
    with run_diagnostic() as global_cfg:
        main(global_cfg)