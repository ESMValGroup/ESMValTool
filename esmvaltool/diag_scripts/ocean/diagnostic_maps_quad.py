"""
Diagnostic Maps quad:

Diagnostic to produce an image showing four maps.
These plost show latitude vs longitude and the cube value is used as the colour
scale.
        model1              model 1 minus model2
        model2 minus obs    model1 minus obs

Note that this diagnostic assumes that the preprocessors do the bulk of the
hard work, and that the cube received by this diagnostic (via the settings.yml
and metadata.yml files) has no time component, a small number of depth layers,
and a latitude and longitude coordinates.

An approproate preprocessor for a 3D+time field would be:
preprocessors:
  prep_map:
    extract_levels:
      levels:  [100., ]
      scheme: linear_extrap
    time_average:

This tool is part of the ocean diagnostic tools package in the ESMValTool,
and was based on the plots produced by the Ocean Assess/Marine Assess toolkit.

Author: Lee de Mora (PML)
        ledm@pml.ac.uk
        """
import logging
import os
import sys
import matplotlib
matplotlib.use('Agg')  # noqa
import iris
import cartopy

import matplotlib.pyplot as plt
import iris.quickplot as qplt
import numpy as np

import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import run_diagnostic

# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def makemapplot(fig, ax, lons, lats, data, title, zrange=[-100, 100],
                lon0=0., drawCbar=True, cbarlabel='', doLog=False, ):

    if len(lons) == 0:
        return fig, ax
    try:
        if len(lons.compressed()) == 0:
            return False, False
    except AttributeError:
        logger.warning('makemapplot: latitude and longitude not provided.')

    lons = np.array(lons)
    lats = np.array(lats)
    data = np.ma.array(data)

    if doLog and zrange[0] * zrange[1] <= 0.:
        data = np.ma.masked_less_equal(ma.array(data), 0.)

    if data.ndim == 1:
        if doLog:
            im = ax.scatter(lons, lats, c=data, lw=0, marker='s',
                            transform=cartopy.crs.PlateCarree(),
                            norm=LogNorm(),
                            vmin=zrange[0],
                            vmax=zrange[1])
        else:
            im = ax.scatter(lons, lats, c=data, lw=0, marker='s',
                            transform=cartopy.crs.PlateCarree(),
                            vmin=zrange[0],
                            vmax=zrange[1])
    else:
        crojp2, data, newLon, newLat = regrid(data, lats, lons)

        if doLog:
            im = ax.pcolormesh(newLon, newLat, data,
                               transform=cartopy.crs.PlateCarree(),
                               norm=LogNorm(vmin=zrange[0],
                                            vmax=zrange[1]), )
        else:
            im = ax.pcolormesh(newLon, newLat, data,
                               transform=cartopy.crs.PlateCarree(),
                               vmin=zrange[0], vmax=zrange[1])

    ax.add_feature(cartopy.feature.LAND,  facecolor='0.85')

    if drawCbar:
        c1 = fig.colorbar(im, pad=0.05, shrink=0.75)
        if len(cbarlabel) > 0:
            c1.set_label(cbarlabel)

    plt.title(title)

    ax.set_axis_off()
    plt.axis('off')
    ax.axis('off')

    return fig, ax


def match_moddel_to_key(m, cfg_dict, input_files_dict, ):
    """Match up the three models and observations dataset from the configs."""
    if not isinstance(cfg_dict, dict):
        print('problem!', m)
        assert 0
    print(cfg_dict.keys())
    print(input_files_dict.keys())

    for input_file, intput_dict in input_files_dict.items():
        match = True
        for key, value in cfg_dict.items():
            if key not in intput_dict:
                match = False
                continue
            if value != intput_dict[key]:
                match = False
                continue
            print("found a match:\t", m,   key, ':', value, '==',
                  intput_dict[key])
        if match:
            print("found a matching file:", m,
                  os.path.basename(input_file), cfg_dict)
            return input_file
    return ''


def get_cube_range(cubes):
    """Determinue the minimum and maximum values of an array of cubes."""
    mi = []
    ma = []
    for cube in cubes:
        mi.append(cube.data.min())
        ma.append(cube.data.max())
    return [np.min(mi), np.max(ma), ]


def get_cube_range_diff(cubes):
    """Determinue the largest deviation from zero in an array of cubes."""
    ma = []
    for cube in cubes:
        ma.append(np.abs(cube.data.min()))
        ma.append(np.abs(cube.data.max()))
    return [-1. * np.max(ma), np.max(ma)]


def multi_model_maps(
        cfg,
        input_files,
):
    """
    Make a simple map plot for an individual model.

    The cfg is the opened global config,
    input_files is the input files dictionairy
    filename is the preprocessing model file.
    """
    filenames = {}
    ctl_key = 'control_model'
    exp_key = 'exper_model'
    obs_key = 'observational_dataset'
    model_types = [ctl_key, exp_key, obs_key]
    for m in model_types:
        print(m, cfg[m])
        filenames[m] = match_moddel_to_key(m, cfg[m], input_files, )

    # ####
    # Load the data for each layer as a separate cube
    layers = {}
    cubes = {}
    for m, input_file in filenames.items():
        print(input_file)
        cube = iris.load_cube(input_file)
        cube = diagtools.bgc_units(cube, input_files[input_file]['short_name'])

        cubes[m] = diagtools.make_cube_layer_dict(cube)
        for layer in cubes[m]:
            layers[layer] = True

    print('layers:', layers)
    print('cubes:', cubes.keys())

    # ####
    # load names:
    exper = input_files[filenames[exp_key]]['dataset']
    control = input_files[filenames[ctl_key]]['dataset']
    obs = input_files[filenames[obs_key]]['dataset']
    long_name = cubes[exp_key][list(layers.keys())[0]].long_name

    # Load image format extention
    image_extention = diagtools.get_image_format(cfg)

    # Make a plot for each layer
    for layer in layers:
        fig = plt.figure()
        fig.set_size_inches(9, 6)

        cube221 = cubes[exp_key][layer]
        cube222 = cubes[exp_key][layer] - cubes[ctl_key][layer]
        cube223 = cubes[ctl_key][layer] - cubes[obs_key][layer]
        cube224 = cubes[exp_key][layer] - cubes[obs_key][layer]

        zrange = get_cube_range_diff([cube222, cube223, cube224])
        n_points = 15
        linspace = np.linspace(zrange[0], zrange[1], n_points, endpoint=True)

        ax = plt.subplot(221)
        qplt.contourf(cube221, n_points, linewidth=0, )
        plt.gca().coastlines()
        plt.title(exper)

        ax = plt.subplot(222, )
        qplt.contourf(cube222, linspace, cmap=plt.cm.get_cmap('bwr'))
        plt.gca().coastlines()
        plt.title(' '.join([exper, 'minus', control]))

        ax = plt.subplot(223, )
        qplt.contourf(cube223, linspace, cmap=plt.cm.get_cmap('bwr'))
        plt.gca().coastlines()
        plt.title(' '.join([control, 'minus', obs]))

        ax = plt.subplot(224, )
        qplt.contourf(cube224, linspace, cmap=plt.cm.get_cmap('bwr'))
        plt.gca().coastlines()
        plt.title(' '.join([exper, 'minus', obs]))

        fig.suptitle(long_name, fontsize=14)

        # Determine image filename:
        fn_list = [long_name, exper, control, obs, str(layer)]
        path = diagtools.folder(cfg['plot_dir']) + '_'.join(fn_list)
        path = path.replace(' ', '') + image_extention

        # Saving files:
        if cfg['write_plots']:
            logger.info('Saving plots to %s', path)
            plt.savefig(path)

        plt.close()


def main(cfg):
    """
    Load the config file, and send it to the plot maker.

    The cfg is the opened global config.
    """
    for index, metadata_filename in enumerate(cfg['input_files']):
        logger.info(
            'metadata filename:\t%s',
            metadata_filename,
        )
        input_files = diagtools.get_input_files(cfg, index=index)
        # #####
        # Multi model time series
        multi_model_maps(
            cfg,
            input_files,
        )

    logger.info('Success')


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
