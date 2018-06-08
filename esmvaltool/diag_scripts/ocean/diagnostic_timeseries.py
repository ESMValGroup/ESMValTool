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


def TimeSeriesPlots(
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

    # Is this data is a multi-model dataset?
    multiModel = md['model'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Making plots for each layer
    for l, (layer, c) in enumerate(cubes.items()):
        layer = str(layer)

        if multiModel:
            qplt.plot(c, label=md['model'], ls=':')
        else:
            qplt.plot(c, label=md['model'])

        # Add title, legend to plots
        title = ' '.join([md['model'], md['long_name']])
        if layer:
            title = ' '.join(
                [title, '(', layer,
                 str(c.coords('depth')[0].units), ')'])
        plt.title(title)
        plt.legend(loc='best')

        # Determine png filename:
        if multiModel:
            path = diagtools.folder(
                cfg['plot_dir']) + os.path.basename(fn).replace(
                    '.nc', '_timeseries_' + str(l) + '.png')
        else:
            path = diagtools.get_image_path(
                cfg,
                md,
                suffix='timeseries_' + str(l),
                image_extention='png',
            )

        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def multiModelTimeSeries(
        cfg,
        metadata,
    ):
    """
        This method makes a time series plot showing several models.
        This function makes a simple plot for an indivudual model.
        The cfg is the opened global config,
        metadata is the metadata dictionairy.
        """

    # Load the data for each layer as a separate cube
    model_cubes = {}
    layers = {}
    for fn in sorted(metadata.keys()):
        cube = iris.load_cube(fn)
        cube = diagtools.sensibleUnits(cube, metadata[fn]['short_name'])
        cubes = diagtools.make_cube_layer_dict(cube)
        model_cubes[fn] = cubes
        for l, c in cubes.items():
            layers[l] = True

    # Make a plot for each layer
    for layer in layers:

        long_name = ''
        z_units = ''
        # Plot each file in the group
        for fn in sorted(metadata.keys()):

            multiModel = metadata[fn]['model'].find('MultiModel') > -1

            if multiModel:
                qplt.plot(
                    model_cubes[fn][layer],
                    label=metadata[fn]['model'],
                    ls=':')
            else:
                qplt.plot(model_cubes[fn][layer], label=metadata[fn]['model'])
            long_name = metadata[fn]['long_name']
            if layer != '':
                z_units = model_cubes[fn][layer].coords('depth')[0].units

        # Add title, legend to plots
        title = long_name
        if layer:
            title = ' '.join([title, '(', str(layer), str(z_units), ')'])
        plt.title(title)
        plt.legend(loc='best')

        # Saving files:
        if cfg['write_plots']:
            path = diagtools.get_image_path(
                cfg,
                metadata[fn],
                prefix='MultiModel',
                suffix='_'.join(['timeseries', str(layer)]),
                image_extention='png',
                basenamelist=['field', 'short_name', 'start_year', 'end_year'],
            )
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

        #######
        # Multi model time series
        multiModelTimeSeries(
            cfg,
            metadata,
        )

        for fn in sorted(metadata.keys()):

            print('-----------------')
            print(
                'model filenames:\t',
                fn,
            )

            ######
            # Time series of individual model
            TimeSeriesPlots(cfg, metadata[fn], fn)

    logger.debug("\n\nThis works\n\n")
    print('Success')


# metadata = get_input_files(cfg)

if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
