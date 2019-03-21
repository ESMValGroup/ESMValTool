#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
###############################################################################
ipcc_fig_9_15.py
Author: Klaus Zimmermann (SMHI, Sweden)
CRESCENDO project
###############################################################################

Description
    Calculate and plot August bottom potential temperature and August and
    February sea ice extent in the southern ocean. For the reference as
    specified in the namelist, the actual temperature is plotted; for the rest
    of the specified datasets the error agains the reference is plotted.
    Additionally, the RMSE is given in the respective label.
    This has been modelled after IPCC AR5 WG1 Ch. 9, Fig. 9.15.

Required diag_script_info attributes (diagnostics specific)

Optional diag_script_info attributes (diagnostic specific)
    [main_plot]
        filename : basename of the plot file

    [ref_line_style]
        linestyle = -
        linewidth = 4

Required variable_info attributes (variable specific)
    none

Required variable attributes (defined in namelist)
    none

Caveats

Modification history
    20180417-A_zimm_kl: written

###############################################################################
"""

from collections import OrderedDict
from functools import partial
import logging
import os

import cartopy.crs as ccrs
import iris
from iris.analysis.cartography import area_weights
from iris import Constraint
import iris.plot as iplt
import matplotlib.pyplot as plt
import numpy as np

# ESMValTool python packages
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger, extract_variables, get_diagnostic_filename,
    get_plot_filename, group_metadata, io, plot, run_diagnostic,
    variables_available)

logger = logging.getLogger(os.path.basename(__file__))


def cm_to_inch(cm):
    CM_PER_INCH = 2.54
    return cm/CM_PER_INCH

def inch_to_cm(inch):
    CM_PER_INCH = 2.54
    return inch*CM_PER_INCH

def calc_error(data, reference=None):
    if reference is None:
        return None
    error = data - reference
    error.metadata = data.metadata
    error.long_name += ' error'
    return error


def setup_figure(no_plots, auto_arrange=False):
    fig = plt.figure(
        figsize=(cm_to_inch(20.), cm_to_inch(22.)),
        tight_layout={'rect': (0., 0.06, 1., 1.)},
        dpi=300,
    )
    if auto_arrange:
        no_columns = min(np.floor(np.sqrt(no_plots)), 4.)
        no_rows = int(np.ceil(no_plots/no_columns))
        no_columns = int(no_columns)
    else:
        no_rows = 4
        no_columns = 4
    proj = ccrs.NearsidePerspective(central_latitude=-90.0,
                                    satellite_height=10000000)
    ax = [fig.add_subplot(no_rows, no_columns, i, projection=proj)
          for i in range(1, no_plots+1)]
    return fig, ax


def index_to_label(i):
    if i > 25:
        raise RuntimeError('Can produce only 26 labels')
    return chr(ord('a')+i)


def plot_field(ax, data, index=None, show_rms=True, **kwargs):
    model, field = data
    field = field.extract(Constraint(latitude=lambda cell: cell.point < -40))
    mesh = iplt.pcolormesh(field, axes=ax, **kwargs)
    mesh.set_rasterized(True)
    ax.coastlines()
    # Calculate fontsize. Empirically determined as 1/12 of axes height to
    # optimally use antarctica space.
    POINTS_PER_INCH = 72.
    scale_transform = ax.figure.dpi_scale_trans.inverted()
    ax_height = ax.get_window_extent().transformed(scale_transform).height
    fontsize = int(ax_height/12.*POINTS_PER_INCH)
    ax.text(
        0.54, 0.485,
        model,
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transAxes,
        fontsize=fontsize,
    )
    if show_rms:
        weights = area_weights(field)
        rms = field.collapsed(['latitude', 'longitude'],
                              iris.analysis.RMS, weights=weights)
        ax.text(
            0.71, 0.38,
            '{:.2f}'.format(rms.data),
            horizontalalignment='right',
            verticalalignment='center',
            transform=ax.transAxes,
            fontsize=fontsize,
        )
    if index is not None:
        ax.text(
            0., 1.,
            '{})'.format(index_to_label(index)),
            horizontalalignment='left',
            verticalalignment='top',
            transform=ax.transAxes,
        )
    return mesh


def plot_sea_ice_extents(ax, sic, **kwargs):
    lat = sic.coord('latitude').points
    lon = sic.coord('longitude').points
    if lat.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    ax.contour(lon, lat, sic.data, [15.],
               transform=ccrs.PlateCarree(), **kwargs)


def plot_antarctic(ax, data, i_plot, **kwargs):
    name = data[0]
    data_dict = data[1]
    bottom_layer = data_dict.get('tob', None)
    mesh = plot_field(ax, (name, bottom_layer), i_plot, **kwargs)
    sic_feb = data_dict.get('sic_feb', None)
    sic_aug = data_dict.get('sic_aug', None)
    if sic_feb is not None:
        plot_sea_ice_extents(ax, sic_feb,
                             colors='k', linestyles='solid', linewidths=1.)
    if sic_aug is not None:
        plot_sea_ice_extents(ax, sic_aug,
                             colors='k', linestyles='dashed', linewidths=1.)
    return mesh


def add_colorbars(fig, ref_mesh, model_mesh):
    ref_cbaxes = fig.add_axes((0.05, .03, .40, .04))
    ref_cb = fig.colorbar(ref_mesh,
                          cax=ref_cbaxes, orientation='horizontal')
    ref_cb.solids.set_rasterized(True)
    model_cbaxes = fig.add_axes((.55, .03, .40, .04))
    model_cb = fig.colorbar(model_mesh,
                            cax=model_cbaxes, orientation='horizontal')
    model_cb.solids.set_rasterized(True)
    fig.text(.5, .05, u'(\u00b0C)', {'fontsize': 12},
             horizontalalignment='center', verticalalignment='center')


def produce_plots(cfg, data):
    no_models = len(data)
    output_path = get_plot_filename('ch09_fig09_15', cfg)
    fig, ax = setup_figure(no_models)
    reference_data = data.pop('REFERENCE')
    ref_mesh = plot_antarctic(ax[0], ('WOCE', reference_data), 0,
                              show_rms=False, cmap='viridis', vmin=-2., vmax=3.)
    model_data = [(dataset_id, data[dataset_id]) for dataset_id in data.keys()]
    meshes = [plot_antarctic(a, d, i+1, cmap='coolwarm', vmin=-3., vmax=3.)
              for i, (a, d) in enumerate(zip(ax[1:], model_data))]
    add_colorbars(fig, ref_mesh, meshes[0])
    fig.savefig(output_path)
    # tags = [
    #     "DM_polar",
    #     "PT_other",
    #     "ST_diff", "ST_rmsd", "ST_clim",
    # ]
    # tag_output_file(output_path, E, "A_zimm_kl",
    #                 "Southern ocean bottom layer potential temperature error. "
    #                 "Similar to Flato et al. 2013, fig. 9.15.",
    #                 tags)


def load_data(config):
    for key in config['input_data'].keys():
        fn = config['input_data'][key]['filename']
        config['input_data'][key]['cube'] = iris.load_cube(fn)


def prepare_reference(config, groups):
    thetao_obs = groups.pop(config['thetao_reference'])
    assert len(thetao_obs) == 1
    thetao_obs = thetao_obs[0]['cube']
    sic_obs = groups.pop(config['sic_reference'])
    sic_obs_feb = [entry['cube'] for entry in sic_obs
                   if entry['variable_group'] == 'sic_obs_feb'][0]
    sic_obs_aug = [entry['cube'] for entry in sic_obs
                   if entry['variable_group'] == 'sic_obs_aug'][0]
    data = {
        'tob': thetao_obs,
        'sic_feb': sic_obs_feb,
        'sic_aug': sic_obs_aug,
    }
    return data


def prepare_dataset(name, dataset, reference):
    thetao = [entry['cube'] for entry in dataset
              if entry['variable_group'].startswith('thetao')][0]
    sic_feb = [entry['cube'] for entry in dataset
               if (entry['variable_group'].startswith('sic')
                   and entry['variable_group'].endswith('feb'))]
    sic_aug = [entry['cube'] for entry in dataset
               if (entry['variable_group'].startswith('sic')
                   and entry['variable_group'].endswith('aug'))]
    error = calc_error(thetao, reference['tob'])
    data = {
        name: {
            'tob': error,
        }
    }
    if len(sic_feb) == 1:
        data[name]['sic_feb'] = sic_feb[0]
    if len(sic_aug) == 1:
        data[name]['sic_aug'] = sic_aug[0]
    return data


def prepare_data(config):
    groups = group_metadata(config['input_data'].values(), 'dataset')
    reference = prepare_reference(config, groups)
    models = [prepare_dataset(name, dataset, reference)
              for name, dataset in groups.items()]
    data = {
        'REFERENCE': reference,
    }
    for model in models:
        data.update(model)
    return data


def main(config):
    """
    Arguments
        project_info : Dictionary containing project information

    Description
        This is the main routine of the diagnostic.
    """

    load_data(config)
    data = prepare_data(config)
    produce_plots(config, data)
    # E, modelconfig = prepare_project_info(project_info,
    #                                       authors=["A_zimm_kl"],
    #                                       diagnostics=["D_flato13ipcc"],
    #                                       projects=["P_crescendo"])
    # data = prepare_data(E, modelconfig)
    # if E.get_write_plots():
    #     produce_plots(E, modelconfig, data)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
