"""
Maps diagnostics
================

Diagnostic to produce images of a map with coastlines from a cube.
These plost show latitude vs longitude and the cube value is used as the colour
scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 3D+time field would be::

  preprocessors:
    prep_map:
      extract_levels:
        levels:  [100., ]
         scheme: linear_extrap
      climate_statistics:
        operator: mean


Note that this recipe may not function on machines with no access to the
internet, as cartopy may try to download the shapefiles. The solution to
this issue is the put the relevant cartopy shapefiles on a disk visible to your
machine, then link that path to ESMValTool via the `auxiliary_data_dir`
variable. The cartopy masking files can be downloaded from::

  https://www.naturalearthdata.com/downloads/

Here, cartopy uses the 1:10, physical coastlines and land files::

      110m_coastline.dbf  110m_coastline.shp  110m_coastline.shx
      110m_land.dbf  110m_land.shp  110m_land.shx

This tool is part of the ocean diagnostic tools package in the ESMValTool.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
"""
import logging
import os
import sys
from itertools import product
import matplotlib.pyplot as plt

import iris
import iris.quickplot as qplt
import cartopy

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def make_map_plots(
        cfg,
        metadata,
        cube,
        key,
        cmap='YlOrRd'
):
    """
    Make a simple map plot for an individual model.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.
    metadata: dict
        the metadata dictionary
    filename: str
        the preprocessed model file.

    """
    # Is this data is a multi-model dataset?
    multi_model = metadata['dataset'].find('MultiModel') > -1
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    if multi_model:
        path = diagtools.folder(cfg['plot_dir'])
        path +=os.path.basename(filename).replace(
                    '.nc', '_map_'+ key.replace(' ', ''), image_extention
                )
    else:
            path = diagtools.get_image_path(
                cfg,
                metadata,
                suffix='map_' + key.replace(' ', '') + image_extention,
            )

    # Load cube and set up units
    cube = diagtools.bgc_units(cube, metadata['short_name'])

    # Making plots for each layer
    qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)

    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # Add title to plot
    title = ' '.join([metadata['dataset'], key])
    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()


def split_variable_groups(variable_group):
    """
    Split variable group into variable and experiment.
    """
    variable, exp, threshold = variable_group.split('_')
    if variable == 'tas':
        variable = 'Surface Temperature'
    exp = exp.upper()
    exp = ' '.join([exp[:3], exp[3], exp[4:]])
    if threshold == '15':
        threshold = '1.5'
    threshold += u'\N{DEGREE SIGN}'
    return variable, exp, threshold


def make_ensemble_map_plots(
        cfg,
        cube,
        variable_group,
        cmap='YlOrRd'
):
    """
    Make a simple map plot for an individual model.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Determine image filename:
    path = diagtools.folder([cfg['plot_dir'], 'ensemble_plots'])+'Ensemblemean_'+variable_group+image_extention

    # Making plots for each layer
    qplt.contourf(cube, 12, linewidth=0, rasterized=True, cmap=cmap)

    try:
        plt.gca().coastlines()
    except AttributeError:
        logger.warning('Not able to add coastlines')

    # tas_ssp585_6
    variable, exp, threshold = split_variable_groups(variable_group)

    # Add title to plot
    title = ' '.join([variable, '- ensemble mean of', exp, 'after', threshold, 'warming'])
    plt.title(title)

    # Saving files:
    if cfg['write_plots']:
        logger.info('Saving plots to %s', path)
        plt.savefig(path)

    plt.close()





def make_gwt_map_plots(cfg):
    """
    Make plots

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    metadatas = diagtools.get_input_files(cfg)
    #print('\n', cfg.keys())

    files_dict = {}
    short_names = set()
    ensembles = set()
    exps = set()
    variable_groups = set()

    for fn, details in metadatas.items():
        #print(fn, details.keys())
        short_names.add(details['short_name'])
        ensembles.add(details['ensemble'])
        exps.add(details['exp'])
        variable_groups.add(details['variable_group'])

        unique_key = (details['variable_group'], details['ensemble'])
        try:
            files_dict[unique_key].append(fn)
        except:
            files_dict[unique_key] = [fn, ]

    print(files_dict.keys())

    # lets look at  minus the historical
    anomaly_cubes = {variable_group:{} for variable_group in variable_groups}

    # Calculate the anomaly for each ensemble/threshold combination
    for ensemble in ensembles:
        for variable_group in variable_groups:
            if variable_group == 'tas_hist':
                continue
            print('Plotting:', ensemble, variable_group)

            if (variable_group, ensemble) not in files_dict:
                continue
            fn = files_dict[(variable_group, ensemble)][0]
            fn_hist = files_dict[('tas_hist', ensemble)][0]

            details = metadatas[fn]
            cube = iris.load_cube( fn)
            cube = diagtools.bgc_units(cube, details['short_name'])

            cube_hist =  iris.load_cube( fn_hist)
            cube_hist = diagtools.bgc_units(cube_hist, details['short_name'])

            cube.data = cube.data - cube_hist.data
            anomaly_cubes[variable_group][ensemble] = cube
            key = variable_group.replace('_',' ') + ' '+ensemble

            # Produce a plot of the anomaly.
            make_map_plots(cfg, details, cube, key)

    for variable_group in variable_groups:
        if variable_group == 'tas_hist':
            continue
        cube_list = []
        for vari, cube in anomaly_cubes[variable_group].items():
            print(variable_group, vari )
            cube_list.append(cube)

        if cube_list == []: continue
        if len(cube_list)<=1:
            ensemble_mean = cube_list[0]
        else:
            cube_data = cube_list[0].data
            for c in cube_list[1:]:
                cube_data += c.data
            cube_data = cube_data/float(len(cube_list))
            ensemble_mean = cube_list[0]
            ensemble_mean.data = cube_data

        make_ensemble_map_plots(cfg, ensemble_mean, variable_group)





    # thresholds = ['15', '2']
    # for short_name in short_names:
    #     for exp in exps:
    #         if exp == 'historical': continue
    #         exp = exp[11:]
    #         for ensemble in ensembles:
    #             for threshold in thresholds:
    #                 #print(short_name, exp, ensemble)
    #                 unique_key = short_name+'_'+exp+'_'+ensemble
    #                 print(unique_key, unique_key in files_dict)




def main(cfg):
    """
    Load the config file, and send it to the plot makers.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionary, passed by ESMValTool.

    """
    cartopy.config['data_dir'] = cfg['auxiliary_data_dir']

    make_gwt_map_plots(cfg)


    # for index, metadata_filename in enumerate(cfg['input_files']):
    #     logger.info(
    #         'metadata filename:\t%s',
    #         metadata_filename,
    #     )
    #
    #     metadatas = diagtools.get_input_files(cfg, index=index)
    #     #thresholds = diagtools.load_thresholds(cfg, metadatas)
    #
    #
    #     # for filename in sorted(metadatas.keys()):
    #     #
    #     #     logger.info('-----------------')
    #     #     logger.info(
    #     #         'model filenames:\t%s',
    #     #         filename,
    #     #     )
    #     #
    #     #     ######
    #     #     # Contour maps of individual model
    #     #     if thresholds:
    #     #         make_map_contour(cfg, metadatas[filename], filename)
    #     #
    #     #     ######
        #     # Maps of individual model
        #     make_map_plots(cfg, metadatas[filename], filename)

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
