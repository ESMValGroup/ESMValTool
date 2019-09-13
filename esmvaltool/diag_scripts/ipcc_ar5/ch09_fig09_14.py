# -*- coding: utf-8 -*-

"""Diagnostic script to plot figure 9.14 of IPCC AR5 chapter 9.

Description
-----------
Calculate and plot the following quantities with regards to sea
surface temperature: zonal mean error, equatorial mean error,
equatorial mean.  The errors are calculated against the reference given
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
        'authors': ['zimmermann_klaus'],
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
    """Fix the given longitudes into the range ``[-180, 180]``."""
    lons = np.array(lons, copy=False, ndmin=1)
    fixed_lons = ((lons + 180) % 360) - 180
    # Make the positive 180s positive again.
    fixed_lons[(fixed_lons == -180) & (lons > 0)] *= -1
    return fixed_lons


def _lon_heimisphere(longitude):
    """Return the hemisphere (E, W or '' for 0) for the given longitude."""
    longitude = _fix_lons(longitude)
    if longitude in (0, 180):
        hemisphere = ''
    elif longitude > 0:
        hemisphere = ' E'
    elif longitude < 0:
        hemisphere = ' W'
    else:
        hemisphere = ''
    return hemisphere


def _lat_heimisphere(latitude):
    """Return the hemisphere (N, S or '' for 0) for the given latitude."""
    if latitude > 0:
        hemisphere = ' N'
    elif latitude < 0:
        hemisphere = ' S'
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


def cm_to_inch(cms):
    """Convert cm to inch."""
    return cms / CM_PER_INCH


def calc_error(data, reference=None):
    """Calculate the error against a reference."""
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
    """
    Merge cubes of different models into one cube.

    This merges cubes from different models/datsets into one big cube
    by promoting the cmip model_id attribute to a scalar coordinate and then
    performing a merge along that coordinate. Conflicting attributes and
    coordinates are simply removed.
    """
    def promote_model_name(cube):
        """Promote model_id attribute to scalar variable."""
        new_cube = cube.copy()
        model_name = new_cube.attributes['model_id']
        coord = iris.coords.AuxCoord(np.array([model_name]),
                                     standard_name=None,
                                     units='no_unit',
                                     long_name=u'model',
                                     var_name='model')
        new_cube.add_aux_coord(coord)
        return new_cube
    cube_list = iris.cube.CubeList([promote_model_name(m) for m in cubes])
    equalise_attributes(cube_list)
    for cube in cube_list:
        cube.cell_methods = tuple()
        for coord in ['day_of_month', 'month_number', 'year']:
            try:
                cube.remove_coord(coord)
            except CoordinateNotFoundError:
                pass
    return cube_list.merge_cube()


def load_data(config):
    """Load cubes into config dict."""
    for key in config['input_data'].keys():
        filename = config['input_data'][key]['filename']
        config['input_data'][key]['cube'] = iris.load_cube(filename)


def prepare_reference(group):
    """Prepare reference cube and remove from the group."""
    ref_name = group[0]['reference_dataset']
    reference_candidates = [ds for ds in group if ds['dataset'] == ref_name]
    assert len(reference_candidates) == 1
    reference = reference_candidates[0]
    group.remove(reference)
    return reference


def mask_equatorial(equ):
    """Mask out Indonesian island area."""
    lon = equ.coord('longitude').points
    equ.data.mask = equ.data.mask | np.logical_and(lon >= 98., lon <= 121.)
    return equ


def prepare_data(config):
    """Perform data calculations."""
    groups = group_metadata(config['input_data'].values(), 'variable_group')
    zm_g = groups["tos_zm"]
    zm_ref = prepare_reference(zm_g)['cube']
    zm_errors = [calc_error(dataset['cube'], zm_ref) for dataset in zm_g]
    eq_g = groups["tos_eq"]
    eq_ref = mask_equatorial(prepare_reference(eq_g)['cube'])
    eqs = [mask_equatorial(ds['cube']) for ds in eq_g]
    eq_errors = [calc_error(eq, eq_ref) for eq in eqs]
    data = {
        'zonal_mean_errors': zm_errors,
        'equatorials': eqs,
        'equatorial_ref': eq_ref,
        'equatorial_errors': eq_errors,
    }
    return data


def setup_figure():
    """Setup basic figure."""
    fig = plt.figure(figsize=(cm_to_inch(18), cm_to_inch(15)))
    axes = np.array(
        [[fig.add_axes([0.10, 0.56, 0.30, 0.35]),
          fig.add_axes([0.50, 0.56, 0.30, 0.35])],
         [fig.add_axes([0.10, 0.10, 0.30, 0.35]),
          fig.add_axes([0.50, 0.10, 0.30, 0.35])]]
    )
    return fig, axes


def plot_zonal_mean_errors_ensemble(axes, zonal_mean_errors, ref_line_style):
    """Plot zonal mean error plot (subfigure a)."""
    axes.set_title('(a) Zonal mean SST error CMIP5')
    axes.yaxis.set_label_text(u'SST error (°C)')
    axes.yaxis.set_minor_locator(MultipleLocator(.5))
    axes.xaxis.set_minor_locator(MultipleLocator(10))
    axes.xaxis.set_major_locator(MultipleLocator(30))
    axes.yaxis.set_major_locator(MultipleLocator(2))
    axes.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    axes.set_ylim(-5., 5.)
    axes.set_xlim(-90., 90.)
    axes.tick_params(which='both', direction='in', top=True, right=True,
                     labelsize=7.)
    axes.xaxis.set_label_text(u'Latitude')
    lines = []
    labels = []
    cube_list = multi_model_merge(zonal_mean_errors)
    for error in zonal_mean_errors:
        lines.append(iplt.plot(error, axes=axes)[0])
        labels.append(error.attributes['model_id'])
    ensemble_mean = cube_list.collapsed('model', iris.analysis.MEAN)
    mean_line = iplt.plot(ensemble_mean,
                          axes=axes, color='#e61f25', **ref_line_style)[0]
    lines = [mean_line] + lines
    labels = ['CMIP5 mean'] + labels
    return (lines, labels)


def plot_equatorial_errors(axes, equatorial_errors, ref_line_style):
    """Plot equatorial errors (subfigure b)."""
    axes.set_title('(b) Equatorial SST error CMIP5')
    axes.yaxis.set_label_text(u'SST error (°C)')
    axes.yaxis.set_minor_locator(MultipleLocator(.5))
    axes.xaxis.set_minor_locator(MultipleLocator(30))
    axes.xaxis.set_major_locator(MultipleLocator(60))
    axes.yaxis.set_major_locator(MultipleLocator(2))
    axes.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    axes.set_ylim(-5., 5.)
    axes.set_xlim(25., 360.)
    axes.tick_params(which='both', direction='in', top=True, right=True,
                     labelsize=7.)
    axes.xaxis.set_label_text(u'Longitude')
    for error in equatorial_errors:
        iplt.plot(error, label=error.attributes['model_id'], axes=axes)
    cube_list = multi_model_merge(equatorial_errors)
    ensemble_mean = cube_list.collapsed('model', iris.analysis.MEAN)
    iplt.plot(ensemble_mean, label='CMIP5 mean',
              axes=axes, color='#e61f25', **ref_line_style)


def plot_zonal_mean_errors_project(axes, zonal_mean_errors, ref_line_style):
    """Plot zonal error multi model mean (subfigure c)."""
    axes.set_title('(c) Zonal mean SST error CMIP5')
    axes.yaxis.set_label_text(u'SST error (°C)')
    axes.yaxis.set_minor_locator(MultipleLocator(.5))
    axes.xaxis.set_minor_locator(MultipleLocator(10))
    axes.xaxis.set_major_locator(MultipleLocator(30))
    axes.yaxis.set_major_locator(MultipleLocator(2))
    axes.xaxis.set_major_formatter(LATITUDE_FORMATTER)
    axes.set_ylim(-5., 5.)
    axes.set_xlim(-90., 90.)
    axes.tick_params(which='both', direction='in', top=True, right=True,
                     labelsize=7.)
    axes.xaxis.set_label_text(u'Latitude')
    lat = zonal_mean_errors[0].coord('latitude').points
    data = np.ma.vstack([m.data for m in zonal_mean_errors])
    std = data.std(axis=0)
    avg = data.mean(axis=0)
    axes.fill_between(lat, avg - std, avg + std, facecolor='#e61f25', alpha=.5)
    axes.plot(lat, avg, color='#e61f25', **ref_line_style)


def plot_equatorials(axes, reference, equatorials, ref_line_style):
    """Plot equatorial multi model mean (subfigure d)."""
    axes.set_title('(d) Equatorial SST CMIP5')
    axes.yaxis.set_label_text(u'SST (°C)')
    axes.yaxis.set_minor_locator(MultipleLocator(.5))
    axes.xaxis.set_minor_locator(MultipleLocator(30))
    axes.xaxis.set_major_locator(MultipleLocator(60))
    axes.yaxis.set_major_locator(MultipleLocator(2))
    axes.xaxis.set_major_formatter(LONGITUDE_FORMATTER)
    axes.set_ylim(22., 31.)
    axes.set_xlim(25., 360.)
    axes.tick_params(which='both', direction='in', top=True, right=True,
                     labelsize=7.)
    axes.xaxis.set_label_text(u'Longitude')
    lon = reference.coord('longitude').points
    data = np.ma.vstack([m.data for m in equatorials[1:]])
    std = data.std(axis=0)
    avg = data.mean(axis=0)
    axes.fill_between(lon, avg - std, avg + std, facecolor='#e61f25', alpha=.5)
    axes.plot(lon, avg, color='#e61f25', **ref_line_style)
    lines = axes.plot(lon, reference.data, 'k', **ref_line_style)
    return (lines, ['HadISST'])


def draw_legend(fig, lines, labels):
    """Draw the legend."""
    return fig.legend(lines, labels,
                      loc='upper left',
                      fontsize=6.,
                      bbox_to_anchor=(.81, .92))


def produce_plots(config, data):
    """Produce all elements of the full plot."""
    ref_line_style = {'linestyle': '-', 'linewidth': 2.}
    fig, axes = setup_figure()
    lines, labels = plot_zonal_mean_errors_ensemble(
        axes[0, 0],
        data['zonal_mean_errors'],
        ref_line_style)
    plot_equatorial_errors(axes[0, 1],
                           data['equatorial_errors'],
                           ref_line_style)
    plot_zonal_mean_errors_project(axes[1, 0],
                                   data['zonal_mean_errors'],
                                   ref_line_style)
    ref_ls, ref_labels = plot_equatorials(axes[1, 1],
                                          data['equatorial_ref'],
                                          data['equatorials'],
                                          ref_line_style)
    all_lines = ref_ls + lines
    all_labels = ref_labels + labels
    legend = draw_legend(fig, all_lines, all_labels)
    path = get_plot_filename('ch09_fig09_14', config)
    fig.savefig(path, additional_artists=[legend], tight_layout=True)
    return path


def write_data(config, data):
    """Write all the calculated data to output file."""
    cubes = iris.cube.CubeList([data['equatorial_ref']]
                               + data['zonal_mean_errors']
                               + data['equatorials']
                               + data['equatorial_errors'])
    path = get_diagnostic_filename('ch09_fig09_14', config)
    iris.save(cubes, path)
    return path


def main(config):
    """
    Run sst zonal mean and equatorial errors diagnostic.

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
    if config['write_plots']:
        provenance_record['plot_file'] = plot_path
    netcdf_path = write_data(config, data)
    with ProvenanceLogger(config) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as cfg:
        main(cfg)
