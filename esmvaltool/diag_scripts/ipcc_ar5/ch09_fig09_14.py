# -*- coding: utf-8 -*-

"""Diagnostic script to plot figure 9.42a of IPCC AR5 chapter 9.

Description
-----------
Calculate and plot the following quantities with regards to sea
surface temperature: zonal mean error, equatorial mean error,
equatorial mean.  The errors are calculated agains the reference given
in the namelist.  Equatorial here means between 5 degrees north and 5
degrees south.  This has been modelled after IPCC AR5 WG1 Ch. 9,
Fig. 9.14.

Author
------
Klaus Zimmermann (SMHI, Sweden)

Project
-------
CRESCENDO
"""

import logging
import os

import iris
from iris.experimental.equalise_cubes import equalise_attributes
from iris.exceptions import CoordinateNotFoundError
import iris.plot as iplt
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import matplotlib.ticker as mticker
import numpy as np

from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger, get_diagnostic_filename,
    get_plot_filename, group_metadata, run_diagnostic)


matplotlib.rcParams.update({'font.size': 9})

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    record = {
        'caption':
        ('(a) Zonally averaged sea surface temperature (SST) error in CMIP5 '
         'models. (b) Equatorial SST error in CMIP5 models. (c) Zonally '
         'averaged multi-model mean SST error for CMIP5 (red line) together '
         'with inter-model standard deviation (shading). (d) Equatorial '
         'multi-model mean SST in CMIP5(red line) together with inter-model '
         'standard deviation (shading) and observations (black).  Model '
         'climatologies are derived from the 1979-1999 mean of the historical '
         'simulations. The Hadley Centre Sea Ice and Sea Surface Temperature '
         '(HadISST)(Rayner et al., 2003) observational climatology for '
         '1979-1999 is used as reference for the error calculation (a), (b), '
         'and (c); and for observations in (d).'),
        'statistics': ['anomaly', 'mean', 'stddev', 'clim'],
        'domains': ['eq', 'global'],
        'plot_types': ['geo', 'sect', 'zonal'],
        'authors': ['zimm_kl'],
        'projects': ['crescendo'],
        'references': ['flato13ipcc', 'hadisst'],
        'realms': ['ocean'],
        'themes': ['phys'],
        'ancestors':
        ancestor_files,
    }
    return record


DEGREE_SYMBOL = u'\u00B0'


def _fix_lons(lons):
    """
    Fix the given longitudes into the range ``[-180, 180]``.
    """
    lons = np.array(lons, copy=False, ndmin=1)
    fixed_lons = ((lons + 180) % 360) - 180
    # Make the positive 180s positive again.
    fixed_lons[(fixed_lons == -180) & (lons > 0)] *= -1
    return fixed_lons


def _lon_heimisphere(longitude):
    """Return the hemisphere (E, W or '' for 0) for the given longitude."""
    longitude = _fix_lons(longitude)
    if longitude == 0 or longitude == 180:
        hemisphere = ''
    elif longitude > 0:
        hemisphere = 'E'
    elif longitude < 0:
        hemisphere = 'W'
    else:
        hemisphere = ''
    return hemisphere


def _lat_heimisphere(latitude):
    """Return the hemisphere (N, S or '' for 0) for the given latitude."""
    if latitude > 0:
        hemisphere = 'N'
    elif latitude < 0:
        hemisphere = 'S'
    else:
        hemisphere = ''
    return hemisphere


def _east_west_formatted(longitude, num_format='g'):
    fmt_string = u'{longitude:{num_format}}{degree}{hemisphere}'
    longitude = _fix_lons(longitude)[0]
    return fmt_string.format(longitude=abs(longitude), num_format=num_format,
                             hemisphere=_lon_heimisphere(longitude),
                             degree=DEGREE_SYMBOL)


def _north_south_formatted(latitude, num_format='g'):
    fmt_string = u'{latitude:{num_format}}{degree}{hemisphere}'
    return fmt_string.format(latitude=abs(latitude), num_format=num_format,
                             hemisphere=_lat_heimisphere(latitude),
                             degree=DEGREE_SYMBOL)


#: A formatter which turns longitude values into nice longitudes such as 110W
LONGITUDE_FORMATTER = mticker.FuncFormatter(lambda v, pos:
                                            _east_west_formatted(v))
#: A formatter which turns longitude values into nice longitudes such as 45S
LATITUDE_FORMATTER = mticker.FuncFormatter(lambda v, pos:
                                           _north_south_formatted(v))


CM_PER_INCH = 2.54


def cm_to_inch(cm):
    return cm/CM_PER_INCH


def inch_to_cm(inch):
    return inch*CM_PER_INCH


def calc_error(data, reference=None):
    if reference is None:
        return None
    error = data - reference
    error.metadata = data.metadata
    name = data.long_name
    if name is None:
        name = data.name()
    error.long_name = '{} error'.format(name)
    return error


def multi_model_merge(cubes):
    def promote_model_name(cube):
        nc = cube.copy()
        model_name = nc.attributes['model_id']
        coord = iris.coords.AuxCoord(np.array([model_name]),
                                     standard_name=None,
                                     units='no_unit',
                                     long_name=u'model',
                                     var_name='model')
        nc.add_aux_coord(coord)
        return nc
    cl = iris.cube.CubeList([promote_model_name(m) for m in cubes])
    equalise_attributes(cl)
    for c in cl:
        c.cell_methods = tuple()
        for co in ['day_of_month', 'month_number', 'year']:
            try:
                c.remove_coord(co)
            except CoordinateNotFoundError:
                pass
    return cl.merge_cube()


def load_data(config):
    for key in config['input_data'].keys():
        fn = config['input_data'][key]['filename']
        config['input_data'][key]['cube'] = iris.load_cube(fn)


def prepare_reference(group):
    ref_name = group[0]['reference_dataset']
    reference_candidates = [ds for ds in group if ds['dataset'] == ref_name]
    assert len(reference_candidates) == 1
    reference = reference_candidates[0]
    group.remove(reference)
    return reference


def mask_equatorial(equ):
    lon = equ.coord('longitude').points
    equ.data.mask[np.logical_and(98. <= lon, lon <= 121.)] = True
    return equ


def prepare_data(config):
    groups = group_metadata(config['input_data'].values(), 'variable_group')
    zm_g = groups["tos_zm"]
    zm_ref = prepare_reference(zm_g)['cube']
    zm_errors = [calc_error(dataset['cube'], zm_ref) for dataset in zm_g]
    eq_g = groups["tos_eq"]
    eq_ref = mask_equatorial(prepare_reference(eq_g)['cube'])
    eqs = [mask_equatorial(ds['cube']) for ds in eq_g]
    eq_errors = [calc_error(eq, eq_ref) for eq in eqs]
    return (zm_errors, eqs, eq_ref, eq_errors)


def setup_figure():
    fig = plt.figure(figsize=(cm_to_inch(18), cm_to_inch(15)))
    ax = np.array(
        [[fig.add_axes([0.10, 0.56, 0.30, 0.35]),
          fig.add_axes([0.50, 0.56, 0.30, 0.35])],
         [fig.add_axes([0.10, 0.10, 0.30, 0.35]),
          fig.add_axes([0.50, 0.10, 0.30, 0.35])]]
    )
    return fig, ax


def plot_zonal_mean_errors_ensemble(ax, zonal_mean_errors, ref_line_style):
    ax.set_title('(a) Zonal mean SST error CMIP5')
    ax.yaxis.set_label_text(u'SST error (°C)')
    ax.yaxis.set_minor_locator(MultipleLocator(.5))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.set_ylim(-5., 5.)
    ax.set_xlim(-90., 90.)
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_label_text(u'Latitude')
    ls = []
    labels = []
    cl = multi_model_merge(zonal_mean_errors)
    for e in zonal_mean_errors:
        ls.append(iplt.plot(e, axes=ax)[0])
        labels.append(e.attributes['model_id'])
    ensemble_mean = cl.collapsed('model', iris.analysis.MEAN)
    m = iplt.plot(ensemble_mean, axes=ax, **ref_line_style)[0]
    ls = [m] + ls
    labels = ['CMIP5 mean'] + labels
    return (ls, labels)


def plot_equatorial_errors(ax, equatorial_errors, ref_line_style):
    ax.set_title('(b) Equatorial SST error CMIP5')
    ax.yaxis.set_label_text(u'SST error (°C)')
    ax.yaxis.set_minor_locator(MultipleLocator(.5))
    ax.xaxis.set_minor_locator(MultipleLocator(30))
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.set_ylim(-5., 5.)
    ax.set_xlim(25., 360.)
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_label_text(u'Longitude')
    for e in equatorial_errors:
        iplt.plot(e, label=e.attributes['model_id'], axes=ax)
    cl = multi_model_merge(equatorial_errors)
    ensemble_mean = cl.collapsed('model', iris.analysis.MEAN)
    iplt.plot(ensemble_mean, label='CMIP5 mean', axes=ax, **ref_line_style)


def plot_zonal_mean_errors_comparison(ax, zonal_mean_errors, ref_line_style):
    ax.set_title('(c) Zonal mean SST error CMIP5')
    ax.yaxis.set_label_text(u'SST error (°C)')
    ax.yaxis.set_minor_locator(MultipleLocator(.5))
    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.xaxis.set_major_locator(MultipleLocator(30))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    ax.set_ylim(-5., 5.)
    ax.set_xlim(-90., 90.)
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_label_text(u'Latitude')
    lat = zonal_mean_errors[0].coord('latitude').points
    data = np.ma.vstack([m.data for m in zonal_mean_errors])
    std = data.std(axis=0)
    avg = data.mean(axis=0)
    ax.fill_between(lat, avg-std, avg+std, alpha=.5)
    ax.plot(lat, avg, **ref_line_style)


def plot_equatorials(ax, reference, equatorials, ref_line_style):
    ax.set_title('(d) Equatorial SST CMIP5')
    ax.yaxis.set_label_text(u'SST (°C)')
    ax.yaxis.set_minor_locator(MultipleLocator(.5))
    ax.xaxis.set_minor_locator(MultipleLocator(30))
    ax.xaxis.set_major_locator(MultipleLocator(60))
    ax.yaxis.set_major_locator(MultipleLocator(2))
    ax.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    ax.set_ylim(22., 31.)
    ax.set_xlim(25., 360.)
    ax.tick_params(which='both', direction='in')
    ax.xaxis.set_label_text(u'Longitude')
    lon = reference.coord('longitude').points
    data = np.ma.vstack([m.data for m in equatorials[1:]])
    std = data.std(axis=0)
    avg = data.mean(axis=0)
    ax.fill_between(lon, avg-std, avg+std, alpha=.5)
    ax.plot(lon, avg, **ref_line_style)
    ls = ax.plot(lon, reference.data, 'k', **ref_line_style)
    return (ls, ['HadISST'])


def draw_legend(fig, ls, labels):
    return fig.legend(ls, labels,
                      loc='upper left',
                      fontsize=7.,
                      bbox_to_anchor=(.81, .92))


def produce_plots(config, data):
    (zonal_mean_errors,
     equatorials,
     equatorial_ref,
     equatorial_errors) = data
    ref_line_style = {'linestyle': '-', 'linewidth': 2.}
    fig, ax = setup_figure()
    ls, labels = plot_zonal_mean_errors_ensemble(ax[0, 0],
                                                 zonal_mean_errors,
                                                 ref_line_style)
    plot_equatorial_errors(ax[0, 1], equatorial_errors, ref_line_style)
    plot_zonal_mean_errors_comparison(ax[1, 0],
                                      zonal_mean_errors, ref_line_style)
    ref_ls, ref_labels = plot_equatorials(ax[1, 1], equatorial_ref,
                                          equatorials, ref_line_style)
    ls = ref_ls + ls
    labels = ref_labels + labels
    legend = draw_legend(fig, ls, labels)
    path = get_plot_filename('ch09_fig09_14', config)
    fig.savefig(path, additional_artists=[legend], tight_layout=True)
    return path


def write_data(config, data):
    cubes = iris.cube.CubeList([data[2]]+data[0]+data[1]+data[3])
    path = get_diagnostic_filename('ch09_fig09_14', config)
    iris.save(cubes, path)
    return path


def main(config):
    """
    Arguments
        config: Dictionary containing project information

    Description
        This is the main routine of the diagnostic.
    """
    load_data(config)
    data = prepare_data(config)
    if config['write_plots']:
        plot_path = produce_plots(config, data)
    ancestor_files = list(config['input_data'].keys())
    provenance_record = get_provenance_record(ancestor_files)
    if plot_path is not None:
        provenance_record['plot_file'] = plot_path
    netcdf_path = write_data(config, data)
    with ProvenanceLogger(config) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
