"""Python example diagnostic."""
import inspect
import logging
import os
import sys

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import yaml

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def determine_profiles_str(cube):
    options = ['latitude', 'longitude']
    for o in options:
        coord = cube.coord(o)
        if len(coord.points) > 1:
            continue
        value = coord.points.mean()
        if o == 'latitude':
            return str(value) + ' N'
        if o == 'longitude':
            if value > 180.:
                return str(value - 360.) + ' W'
            return str(value) + ' E'
    return ''


def profilesPlots(
        cfg,
        md,
        fn,
):
    """
        This function makes a simple plot for an indivudual model.
        The cfg is the opened global config,
        md is the metadata dictionairy
        fn is the preprocessing model file.
        """
    # Load cube and set up units
    cube = iris.load_cube(fn)
    cube = diagtools.sensibleUnits(cube, md['short_name'])

    # Make annual means from:
    cube = cube.aggregated_by('year', iris.analysis.MEAN)

    # Is this data is a multi-model dataset?
    multiModel = md['model'].find('MultiModel') > -1

    #
    times = cube.coord('time')
    timefloats = diagtools.timecoord_to_float(times)

    cmap = plt.cm.get_cmap('jet')

    plotDetails = {}
    for t, time in enumerate(timefloats):

        c = cmap((time - timefloats[0]) / (timefloats[-1] - timefloats[0]))

        qplt.plot(cube[t, :], cube[t, :].coord('depth'), c=c)
        plotDetails[t] = {'c': c, 'ls': '-', 'lw': 1, 'label': str(int(time))}

    # Add title to plot
    title = ' '.join([
        md['model'],
        md['long_name'],
    ])
    plt.title(title)

    # Add Legend outside right.
    diagtools.add_legend_outside_right(plotDetails, plt.gca())

    # Determine png filename:
    if multiModel:
        path = diagtools.folder(
            cfg['plot_dir']) + os.path.basename(fn).replace(
                '.nc', '_profile.png')
    else:
        path = diagtools.get_image_path(
            cfg,
            md,
            suffix='profile',
            image_extention='png',
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def main(cfg):
    #####
    input_files = diagtools.get_input_files(cfg)

    print('cfg:\tContents:')
    for k in cfg.keys():
        print('CFG:\t', k, '\t', cfg[k])

    for i, metadatafilename in enumerate(cfg['input_files']):
        print(
            '\nmetadata filename:',
            metadatafilename,
        )

        metadata = diagtools.get_input_files(cfg, index=i)
        for fn in sorted(metadata.keys()):

            print('-----------------')
            print(
                'model filenames:\t',
                fn,
            )

            ######
            # Time series of individual model
            profilesPlots(cfg, metadata[fn], fn)

    logger.debug("\n\nThis works\n\n")
    print('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
