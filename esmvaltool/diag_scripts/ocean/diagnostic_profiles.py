import logging
import os
import sys

import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def determine_profiles_str(cube):
    """
        Determines a string from the cube, to describe the profile.
        """
    options = ['latitude', 'longitude']
    for option in options:
        coord = cube.coord(option)
        if len(coord.points) > 1:
            continue
        value = coord.points.mean()
        if option == 'latitude':
            return str(value) + ' N'
        if option == 'longitude':
            if value > 180.:
                return str(value - 360.) + ' W'
            return str(value) + ' E'
    return ''


def profiles_plots(
        cfg,
        md,
        fn,
):
    """
        This function makes a simple profile plot for an individual model.

        The cfg is the opened global config,
        md is the metadata dictionairy
        fn is the preprocessing model file.
        """
    # Load cube and set up units
    cube = iris.load_cube(fn)
    cube = diagtools.bgc_units(cube, md['short_name'])

    # Make annual means from:
    cube = cube.aggregated_by('year', iris.analysis.MEAN)

    # Is this data is a multi-model dataset?
    multi_model = md['model'].find('MultiModel') > -1

    #
    times = cube.coord('time')
    time_floats = diagtools.timecoord_to_float(times)

    cmap = plt.cm.get_cmap('jet')

    plot_details = {}
    for time_index, time in enumerate(time_floats):

        c = cmap((time - time_floats[0]) / (time_floats[-1] - time_floats[0]))

        qplt.plot(cube[time_index, :], cube[time_index, :].coord('depth'), c=c)
        plot_details[time_index] = {'c': c, 'ls': '-', 'lw': 1,
                                    'label': str(int(time))}

    # Add title to plot
    title = ' '.join([
        md['model'],
        md['long_name'],
    ])
    plt.title(title)

    # Add Legend outside right.
    diagtools.add_legend_outside_right(plot_details, plt.gca())

    # Determine png filename:
    if multi_model:
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
    """
        Main function to load the config file, and send it to the plot maker.

        The cfg is the opened global config.
        """
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
            profiles_plots(cfg, metadata[fn], fn)

    logger.debug("\n\nThis works\n\n")
    print('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
