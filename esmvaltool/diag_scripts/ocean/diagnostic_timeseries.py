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


def time_series_plots(
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
    cube = diagtools.bgc_units(cube, md['short_name'])

    # Is this data is a multi-model dataset?
    multi_model = md['model'].find('MultiModel') > -1

    # Make a dict of cubes for each layer.
    cubes = diagtools.make_cube_layer_dict(cube)

    # Making plots for each layer
    for la, (layer, c) in enumerate(cubes.items()):
        layer = str(layer)

        if multi_model:
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
        if multi_model:
            # path = diagtools.folder(
            #     cfg['plot_dir']) + os.path.basename(fn).replace(
            #         '.nc', '_timeseries_' + str(l) + '.png')
            path = diagtools.get_image_path(
                cfg,
                md,
                prefix='MultiModel',
                suffix='_'.join(['timeseries', str(layer)]),
                image_extention='png',
                basenamelist=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        else:
            path = diagtools.get_image_path(
                cfg,
                md,
                suffix='timeseries_' + str(la),
                image_extention='png',
            )

        # Saving files:
        if cfg['write_plots']:

            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def multi_model_time_series(
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
        cube = diagtools.bgc_units(cube, metadata[fn]['short_name'])
        cubes = diagtools.make_cube_layer_dict(cube)
        model_cubes[fn] = cubes
        for l, c in cubes.items():
            layers[l] = True

    # Make a plot for each layer
    for layer in layers:

        long_name = ''
        z_units = ''
        plot_details = {}
        cmap = plt.cm.get_cmap('viridis')

        # Plot each file in the group
        for i, fn in enumerate(sorted(metadata.keys())):
            c = cmap(
                (float(i) / (len(metadata.keys()) - 1.)
                 ))  # - timefloats[0]) / (timefloats[-1] - timefloats[0]))

            multi_model = metadata[fn]['model'].find('MultiModel') > -1

            if multi_model:
                qplt.plot(
                    model_cubes[fn][layer],
                    c=c,
                    # label=metadata[fn]['model'],
                    ls=':',
                    lw=2.,
                )
                plot_details[fn] = {
                    'c': c,
                    'ls': ':',
                    'lw': 2.,
                    'label': metadata[fn]['model']
                }
            else:
                qplt.plot(
                    model_cubes[fn][layer],
                    c=c,
                    # label=metadata[fn]['model'])
                    ls='-',
                    lw=2.,
                )
                plot_details[fn] = {
                    'c': c,
                    'ls': '-',
                    'lw': 2.,
                    'label': metadata[fn]['model']
                }

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
                prefix='MultipleModels_',
                suffix='_'.join(['timeseries', str(layer)]),
                image_extention='png',
                basenamelist=[
                    'field', 'short_name', 'preprocessor', 'diagnostic',
                    'start_year', 'end_year'
                ],
            )

        # Resize and add legend outside thew axes.
        plt.gcf().set_size_inches(9., 6.)
        diagtools.add_legend_outside_right(
            plot_details, plt.gca(), column_width=0.15)

        logger.info('Saving plots to %s', path)
        plt.savefig(path)
        plt.close()


def main(cfg):
    """
        Main function to load the config file, and send it to the plot maker.

        The cfg is the opened global config.
        """
    for i, metadata_filename in enumerate(cfg['input_files']):
        print(
            '\nmetadata filename:',
            metadata_filename,
        )

        metadata = diagtools.get_input_files(cfg, index=i)

        #######
        # Multi model time series
        multi_model_time_series(
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
            time_series_plots(cfg, metadata[fn], fn)

    logger.debug("\n\nThis works\n\n")
    print('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
