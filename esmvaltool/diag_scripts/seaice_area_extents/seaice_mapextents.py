"""Diagnostic script to plot extent and differences.

based on code from Anton Steketee's COSIMA recipes notebook
https://cosima-recipes.readthedocs.io/en/latest/Examples/Sea_Ice_Area_Concentration_Volume_with_Obs.html
"""

import calendar
import logging
import os

import iris
import iris.plot as iplt
import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
from cartopy.crs import SouthPolarStereo
from esmvalcore.preprocessor import extract_month

from esmvaltool.diag_scripts.shared import (
    group_metadata,
    select_metadata,
    run_diagnostic,
    save_figure,
)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def map_fig_diff(model_dict, obs_si, months):
    """Create map with model dictionary: labels, cubes."""
    # figure set up, width 9 for 2 models
    figure = plt.figure(figsize=(9, len(months) * 3.5))
    j = 0  # to iterate through positions on figure

    for mon in months:  # eg.[2,9]
        i = 1
        obs_cube = extract_month(obs_si, mon)
        for mod_label, mod_si in model_dict.items():
            mod_cube = extract_month(mod_si, mon)

            out = mod_cube.copy()
            diff = mod_cube.data - obs_cube.data
            out.data = diff

            plt.subplot(len(months), 3, i + j * 3,
                        projection=SouthPolarStereo(true_scale_latitude=-70))

            diffmap = iplt.contourf(out, levels=np.arange(-90, 91, 20),
                                    cmap='RdBu')

            iplt.contour(obs_cube, levels=[15], colors=['yellow'])
            iplt.contour(mod_cube, levels=[15],
                         linewidths=1.0, colors=['black'])

            plt.title(calendar.month_abbr[mon] + ' ' + mod_label)

            i += 1
        j += 1

    line_cdr = mlines.Line2D([], [], color='yellow', label="Observed Extent")
    line_mod = mlines.Line2D([], [], color='black', label="Modelled Extent")

    plt.legend(handles=[line_cdr, line_mod], loc='center left',
               bbox_to_anchor=(1.2, 0.5))
    cax = plt.axes([0.7, 0.55, 0.04, 0.3])
    _ = plt.colorbar(diffmap, cax=cax,
                     label='Difference in \nSea Ice Concentration')

    plt.subplots_adjust(left=0.05, bottom=0.05,
                        right=0.95, top=0.95,
                        wspace=0.05, hspace=0.05)

    return figure


def main(cfg):
    """Compute."""
    # Get the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()

    groups = group_metadata(input_data, 'variable_group')
    # assign obs_si - select
    selection = select_metadata(input_data, project='OBS6')
    obs_si = iris.load_cube(selection[0]['filename'])

    for group_name in groups.keys():
        mod_dict = {}
        ancestor_filels = []
        logger.info("Processing variable %s", group_name)

        for attr in groups[group_name]:
            if attr['project'].startswith('OBS'):
                logger.info("Load OBS dataset %s", attr['dataset'])
                obs_si = iris.load_cube(attr['filename'])
            else:
                logger.info("load model dataset %s", attr['dataset'])
                mod_dict[attr['dataset']] = iris.load_cube(attr['filename'])
            ancestor_filels.append(attr['filename'])

        logger.info("creating map differences")
        mapfig = map_fig_diff(mod_dict, obs_si, cfg['months'])

        provenance_record = get_provenance_record(ancestor_filels)
        save_figure(group_name, provenance_record, cfg, figure=mapfig)


def get_provenance_record(ancestor_files):
    """Build provenance dictionary."""
    record = {
        'ancestors': ancestor_files,
        'authors': ['chun_felicity', 'steketee_anton'],
        'caption': 'siconc observations difference from model',
        'domains': ['shpolar'],
        'plot_types': ['polar'],
        'references': [],
        'statistics': ['diff'],
    }
    return record


if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
