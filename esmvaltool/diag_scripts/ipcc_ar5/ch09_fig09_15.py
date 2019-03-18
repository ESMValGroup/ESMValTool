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

# from auxiliary import info, warning, error  # noqa: F401
# from smhi import regrid_esmpy
# from smhi.pipeline import (build_var_constraint,
#                            load,
#                            LoadingStep, ProcessingStep)
# from smhi.tools import (
#     calc_error,
#     ensure_dir_exists,
#     get_plot_options,
#     prepare_project_info,
#     tag_output_file,
# )


def calc_error(data, reference=None):
    if reference is None:
        return None
    error = data[1] - reference[1]
    error.metadata = data[1].metadata
    error.long_name += ' error'
    return data[0], error


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
    ax.text(
        0.54, 0.49,
        model,
        horizontalalignment='center',
        verticalalignment='center',
        transform=ax.transAxes,
    )
    if show_rms:
        weights = area_weights(field)
        rms = field.collapsed(['latitude', 'longitude'],
                              iris.analysis.RMS, weights=weights)
        ax.text(
            0.72, 0.37,
            '{:.2f}'.format(rms.data),
            horizontalalignment='right',
            verticalalignment='center',
            transform=ax.transAxes,
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
    bottom_layer = (data[0], data[1]['tob'])
    mesh = plot_field(ax, bottom_layer, i_plot, **kwargs)
    plot_sea_ice_extents(ax, data[1]['sic_feb'],
                         colors='k', linestyles='solid', linewidths=1.)
    plot_sea_ice_extents(ax, data[1]['sic_aug'],
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
    fig.text(.5, .05, u'(\u00b0C)', {'fontsize': 22},
             horizontalalignment='center', verticalalignment='center')


def produce_plots(cfg, data):
    no_models = len(data)
    output_path = get_plot_filename('ch09_fig09_15', cfg)
    fig, ax = setup_figure(no_models)
    reference_data = data.pop('REFERENCE')
    ref_mesh = plot_antarctic(ax[0], ('WOCE', reference_data), 0,
                              show_rms=False, cmap='jet', vmin=-2., vmax=3.)
    model_data = [(dataset_id, data[dataset_id]) for dataset_id in data.keys()]
    meshes = [plot_antarctic(a, d, i+1, cmap='bwr', vmin=-2., vmax=2.)
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


def get_dataset(entry, dataset_id=None):
    cube = entry['cube']
    if dataset_id == None:
        # try CMIP6 style source identification
        dataset_id = cube.attributes.get('source_id', None)
    if dataset_id is None:
        # try CMIP5 style source identification
        dataset_id = cube.attributes.get('model_id', None)
    if dataset_id is None:
        # take name from recipe
        dataset_id = entry['dataset']
    return dataset_id, cube


def prepare_tob(obs, models):
    assert len(obs) == 1
    obs = get_dataset(obs[0], 'REFERENCE')
    models = map(get_dataset, models)
    errors = map(partial(calc_error, reference=obs), models)
    res = {
        obs[0]: {
            'tob': obs[1],
        }
    }
    res.update({
        model[0]: {
            'tob': model[1],
        } for model in errors
    })
    return res


def prepare_sic(obs_feb, obs_aug, models_feb, models_aug):
    assert len(obs_feb) == 1
    assert len(obs_aug) == 1
    obs_feb = get_dataset(obs_feb[0], 'REFERENCE')
    obs_aug = get_dataset(obs_aug[0], 'REFERENCE')
    assert obs_feb[0] == obs_aug[0]
    models_feb = map(get_dataset, models_feb)
    models_aug = map(get_dataset, models_aug)
    res = {
        obs_feb[0]: {
            'sic_feb': obs_feb[1],
            'sic_aug': obs_feb[1],
        }
    }
    res.update({
        model[0][0]: {
            'sic_feb': model[0][1],
            'sic_aug': model[1][1],
        } for model in zip(models_feb, models_aug)
    })
    return res


def prepare_data(config):
    groups = group_metadata(config['input_data'].values(), 'variable_group')
    tob = prepare_tob(groups['thetao_obs'],
                      groups['thetao_models'])
    sic = prepare_sic(groups['sic_obs_feb'],
                      groups['sic_obs_aug'],
                      groups['sic_models_feb'],
                      groups['sic_models_aug'])
    data = {}
    for key in tob.keys():
        data[key] = {**tob[key], **sic[key]}
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
