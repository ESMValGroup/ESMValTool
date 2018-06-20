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


def make_profiles_plots(
        cfg,
        metadata,
        filename,
):
    """
    This function makes a simple profile plot for an individual model.

    The cfg is the opened global config,
    metadata is the metadata dictionairy
    filename is the preprocessing model file.
    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Make annual means from:
    cube = cube.aggregated_by('year', iris.analysis.MEAN)

    # Is this data is a multi-model dataset?
    multi_model = metadata['model'].find('MultiModel') > -1

    #
    times = cube.coord('time')
    times_float = diagtools.timecoord_to_float(times)
    time_0 = times_float[0]

    cmap = plt.cm.get_cmap('jet')

    plot_details = {}
    for time_index, time in enumerate(times_float):

        color = cmap((time - time_0) / (times_float[-1] - time_0))

        qplt.plot(cube[time_index, :], cube[time_index, :].coord('depth'),
                  c=color)

        plot_details[time_index] = {'c': color, 'ls': '-', 'lw': 1,
                                    'label': str(int(time))}

    # Add title to plot
    title = ' '.join([
        metadata['model'],
        metadata['long_name'],
    ])
    plt.title(title)

    # Add Legend outside right.
    diagtools.add_legend_outside_right(plot_details, plt.gca())

    # Determine png filename:
    if multi_model:
        path = diagtools.folder(
            cfg['plot_dir']) + os.path.basename(filename).replace(
                '.nc', '_profile.png')
    else:
        path = diagtools.get_image_path(
            cfg,
            metadata,
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
    for index, metadata_filename in enumerate(cfg['input_files']):
        print(
            '\nmetadata filename:',
            metadata_filename,
        )

        metadatas = diagtools.get_input_files(cfg, index=index)
        for filename in sorted(metadatas.keys()):

            print('-----------------')
            print(
                'model filenames:\t',
                filename,
            )

            ######
            # Time series of individual model
            make_profiles_plots(cfg, metadatas[filename], filename)

    logger.debug("\n\nThis works\n\n")
    print('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
