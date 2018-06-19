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


def determine_transect_str(cube):
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


def transectsPlots(
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

    qplt.contourf(cube, 25)
    plt.axes().set_yscale('log')

    # Add title to plot
    title = ' '.join(
        [md['model'], md['long_name'],
         determine_transect_str(cube)])
    plt.title(title)

    # Determine png filename:
    if multiModel:
        path = diagtools.folder(
            cfg['plot_dir']) + os.path.basename(fn).replace(
                '.nc', '_transect.png')
    else:
        path = diagtools.get_image_path(
            cfg,
            md,
            suffix='transect',
            image_extention='png',
        )

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def main(cfg):
    #####
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
            transectsPlots(cfg, metadata[fn], fn)

    logger.debug("\n\nThis works\n\n")
    print('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
