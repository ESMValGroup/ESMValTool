"""
Model 1 vs Model 2 vs Observations diagnostics.
===============================================

Diagnostic to produce an image showing four maps, based on a comparison of two
different models results against an observational dataset. This process is
often used to compare a new iteration of a model under development against
a previous version of the same model. The four map plots are:

* Top left: model 1
* Top right: model 1 minus model 2
* Bottom left: model 2 minus obs
* Bottom right: model 1 minus obs

All four plots show latitude vs longitude and the cube value is used as the
colour scale.

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An appropriate preprocessor for a 3D+time field would be::
   preprocessors:
     prep_map:
       extract_levels:
         levels:  [100., ]
         scheme: linear_extrap
       climate_statistics:
         operator: mean

This diagnostic also requires the ``exper_model``, ``exper_model`` and
``observational_dataset`` keys in the recipe::

  diagnostics:
     diag_name:
       ...
       scripts:
         Global_Ocean_map:
           script: ocean/diagnostic_maps_quad.py
           exper_model:  {Model 1 dataset details}
           control_model: {Model 2 dataset details}
           observational_dataset: {Observational dataset details}

This tool is part of the ocean diagnostic tools package in the ESMValTool,
and was based on the plots produced by the Ocean Assess/Marine Assess toolkit.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk

"""
import logging
import os
import sys
from pprint import pformat

import cartopy.crs as ccrs
import iris
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import numpy as np

from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import (  # ProvenanceLogger,
    group_metadata, run_diagnostic, save_figure,
)

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def add_map_subplot(subplot, cube, nspace, title='', cmap=''):
    """Add a map subplot to the current pyplot figure.

    Parameters
    ----------
    subplot: int
        The matplotlib.pyplot subplot number. (ie 221)
    cube: iris.cube.Cube
        the iris cube to be plotted.
    nspace: numpy.array
        An array of the ticks of the colour part.
    title: str
        A string to set as the subplot title.
    cmap: str
        A string to describe the matplotlib colour map.
    """
    plt.subplot(subplot)
    new_cube, extent = iris.analysis.cartography.project(cube,
                                                         ccrs.PlateCarree(),
                                                         nx=400,
                                                         ny=200)
    qplot = qplt.contourf(new_cube,
                          nspace,
                          linewidth=0,
                          cmap=plt.cm.get_cmap(cmap))
    qplot.colorbar.set_ticks(
        [nspace.min(), (nspace.max() + nspace.min()) / 2.,
         nspace.max()])

    plt.gca().coastlines()
    plt.title(title)


def multi_model_maps(cfg):
    """Make the four pane model vs model vs obs comparison plot.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.
    input_files: dict
        the metadata dictionary
    """
    from pprint import pprint
    filenames = {}
    ctl_key = 'control_model'
    exp_key = 'exper_model'
    obs_key = 'observational_dataset'
    pprint(dict(cfg))
    input_data = cfg['input_data'].values()
    grouped_input_data = group_metadata(input_data, 'dataset', sort='dataset')
    logger.info(
        "Example of how to group and sort input data by variable groups from "
        "the recipe:\n%s", pformat(grouped_input_data))

    input_files = diagtools.get_input_files(cfg)

    model_types = [ctl_key, exp_key, obs_key]
    for model_type in model_types:
        logger.debug(model_type, cfg[model_type])
        filenames[model_type] = diagtools.match_model_to_key(
            model_type, cfg[model_type], input_files)

    # Load the data for each cube
    model_type_cube = {}
    cube_depth = {}
    for model_type, input_file in filenames.items():
        nemo_cube = iris.load_cube(input_file)
        remove_extra_time_axis(nemo_cube)
        nemo_cube = diagtools.bgc_units(nemo_cube,
                                        input_files[input_file]['short_name'])
        #        cube_depth[model_type] = False
        cube_depth[model_type] = [0]
        for coord in nemo_cube.coords():
            if coord.standard_name == 'depth' and coord.shape[0] > 1:
                cube_depth[model_type] = coord.points
        model_type_cube[model_type] = nemo_cube

    # Workaround to allow us to subtract cubes...
    for model_type in [exp_key, ctl_key]:
        model_type_cube[model_type].coord(
            "latitude").points[:] = model_type_cube[obs_key].coord(
                "latitude").points
        model_type_cube[model_type].coord(
            "longitude").points[:] = model_type_cube[obs_key].coord(
                "longitude").points

    # ####
    # load names:
    exper = input_files[filenames[exp_key]]['dataset']
    control = input_files[filenames[ctl_key]]['dataset']
    obs = input_files[filenames[obs_key]]['dataset']
    long_name = model_type_cube[exp_key].long_name

    # Load image format extension
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    if len(cube_depth[exp_key]) == 1:
        indices = [0]
    else:
        indices = [2, 5, 10, 15]

    for index in indices:
        fig = plt.figure()
        fig.set_size_inches(9, 6)

        # Take slices of the cubes for this index
        cube_slice = {}
        for model_type in model_types:
            if len(cube_depth[model_type]) == 1:
                # 2D cube
                cube_slice[model_type] = iris.util.squeeze(
                    model_type_cube[model_type])
            else:
                # 3D cube - select relevant layer
                slices = [
                    slice(None) for ii in model_type_cube[model_type].shape
                ]
                print("model_type, model_type_cube[model_type].coords() : ",
                      model_type, model_type_cube[model_type].coords())
                coord_dim = model_type_cube[model_type].coord_dims('depth')[0]
                slices[coord_dim] = index
                cube_slice[model_type] = iris.util.squeeze(
                    model_type_cube[model_type][tuple(slices)])

        cube221 = cube_slice[exp_key]
        cube222 = cube_slice[exp_key] - cube_slice[ctl_key]
        cube223 = cube_slice[ctl_key] - cube_slice[obs_key]
        cube224 = cube_slice[exp_key] - cube_slice[obs_key]

        # create the z axis for plots 2, 3, 4.
        zrange1 = diagtools.get_cube_range([
            cube221,
        ])
        if model_type_cube[exp_key].long_name == "Sea Surface Temperature":
            zrange2 = [-5.0, 5.0]
        else:
            zrange2 = [-2.0, 2.0]

        linspace1 = np.linspace(zrange1[0], zrange1[1], 12, endpoint=True)
        linspace2 = np.linspace(zrange2[0], zrange2[1], 12, endpoint=True)

        # Add the sub plots to the figure.
        add_map_subplot(221, cube221, linspace1, cmap='viridis', title=exper)
        add_map_subplot(222,
                        cube222,
                        linspace2,
                        cmap='bwr',
                        title=' '.join([exper, 'minus', control]))
        add_map_subplot(223,
                        cube223,
                        linspace2,
                        cmap='bwr',
                        title=' '.join([control, 'minus', obs]))
        add_map_subplot(224,
                        cube224,
                        linspace2,
                        cmap='bwr',
                        title=' '.join([exper, 'minus', obs]))

        # Add overall title


#        formatted_layer = f"{layer:.3f}"
#        fig.suptitle(long_name, fontsize=14)

# Add overall title
    layers = cube_depth["exper_model"]
    for index in indices:
        layer = layers[index]
        print('layers', layers)
        print('cube depth : ', cube_depth.items())
        value = model_type_cube.values()
        print('value :', value)
        keys = model_type_cube.keys()
        print('key:', keys)
        key_length = len(keys)
        print('key_length', key_length)

        if len(layers) == 1:
            fig.suptitle(long_name, fontsize=14)
        else:
            try:
                print("PRINTED CORRECTLY")
                #                layer_value = float(layer)
                formatted_layer = f"{layer:.3f}"
                fig.suptitle(long_name + ' at ' + formatted_layer + 'm',
                             fontsize=14)
            except ValueError:
                print(f"Error: Could not convert layer '{layer}' to float.")

        # Determine image filename:
        fn_list = [long_name, exper, control, obs, str(index), model_type]
        path = diagtools.folder(
            cfg['plot_dir']) + '_'.join(fn_list) + str(layer)
        path = path.replace(' ', '') + image_extention

        # Saving files:
        logger.info('Saving plots to %s', path)
        provenance_record = diagtools.prepare_provenance_record(
            cfg,
            caption=f'Quadmap models comparison against {obs}{layer}',
            statistics=[
                'mean',
                'diff',
            ],
            domain=['global'],
            plot_type=['map'],
            ancestors=list(input_files.keys()),
        )
        save_figure('_'.join(fn_list), provenance_record, cfg, fig, close=True)


def remove_extra_time_axis(cube):
    """Remove the extra time axis from the |input variable| provided by the
    ``cube`` parameter.

    If two or more time coordinates with the same standard name exist
    in the |input variable|, the auxiliary time coordinates are
    removed. This situation occurs when neither the `time_counter` nor
    `time_centered` coordinates in |model output files| from NEMO is a
    record.

    Otherwise, remove the `time_counter` coordinate. This situation
    occurs when the `time_counter` coordinate in |model output files|
    from NEMO is a record; the points are sequential whole numbers
    corresponding to the number of time slices in the |input variable|,
    which prevents concatenation (see
    https://groups.google.com/forum/#!topic/scitools-iris/BHhKs_DQSUA).
    In this case, the `time_counter` coordinate is a dimension
    coordinate, so the remaining auxiliary time coordinate is promoted
    to a dimension coordinate.

    :param cube: the |input variable|
    :type cube: :class:`iris.cube.Cube`
    """
    count_time_axes = 0
    time_axes_names = []
    for coord in cube.coords():
        if iris.util.guess_coord_axis(coord) == 'T':
            count_time_axes = count_time_axes + 1
            time_axes_names.append(coord.standard_name)

    if count_time_axes >= 2 and len(set(time_axes_names)) == 1:
        for aux_coord in cube.coords(dim_coords=False):
            if iris.util.guess_coord_axis(aux_coord) == 'T':
                cube.remove_coord(aux_coord)
    else:
        for coord in cube.coords():
            if coord.var_name == 'time_counter':
                cube.remove_coord(coord)


def main(cfg):
    """Load the config file, and send it to the plot maker.

    Parameters
    ----------
    cfg: dict
        the opened global config dictionairy, passed by ESMValTool.



    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )
        input_files = diagtools.get_input_files(cfg, index=index)
        # #####
        # Multi model time series
    """

    multi_model_maps(cfg)
    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
