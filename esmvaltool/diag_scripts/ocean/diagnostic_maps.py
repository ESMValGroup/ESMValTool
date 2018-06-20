"""Python example diagnostic."""
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


def make_map_plots(
        cfg,
        metadata,
        filename,
):
    """
    This function makes a simple map plot for an indivudual model.

    The cfg is the opened global config,
    metadata is the metadata dictionairy
    filename is the preprocessing model file.
    """
    # Load cube and set up units
    cube = iris.load_cube(filename)
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = metadata['model'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Making plots for each layer
    for layer_index, (layer, cube_layer) in enumerate(cubes.items()):
        layer = str(layer)

        qplt.contourf(cube_layer, 25)

        try:
            plt.gca().coastlines()
        except AttributeError:
            print('Not able to add coastlines')

        # Add title to plot
        title = ' '.join([metadata['model'], metadata['long_name']])
        if layer:
            title = ' '.join(
                [title, '(', layer,
                 str(cube_layer.coords('depth')[0].units), ')'])
        plt.title(title)

        # Determine png filename:
        if multi_model:
            path = diagtools.folder(
                cfg['plot_dir']) + os.path.basename(filename).replace(
                    '.nc', '_map_' + str(layer_index) + '.png')
        else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='map_' + str(layer_index),
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
    ####
    for k in cfg.keys():
        print('CFG:\t', k, '\t', cfg[k])

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
            make_map_plots(cfg, metadatas[filename], filename)

    logger.debug("\n\nThis works\n\n")
    print('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
