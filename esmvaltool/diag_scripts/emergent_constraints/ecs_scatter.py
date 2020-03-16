#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to calculate various emergent constraints for ECS.

Description
-----------
Calculate the X-axis of various emergent constraint for the Equilibrium
Climate Sensitivity (ECS).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
diag : str
  Emergent constraint to calculate (must be one of ``'brient_shal'``, ``'su'``,
  ``'volodin'``, ``'zhai'``.
metric : str, optional (default: 'regression_slope')
    Metric to measure model error. Only relevant for Su et al. (2014)
    constraint. Must be one of ``'regression_slope'``,
    ``'correlation_coefficient'``.
output_attributes : dict, optional
    Write additional attributes to netcdf files.
pattern : str, optional
    Pattern matched against ancestor file names.
savefig_kwargs : dict, optional
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings : dict, optional
    Options for :func:`seaborn.set` (affects all plots), see
    <https://seaborn.pydata.org/generated/seaborn.set.html>.

"""

import logging
import os
from copy import deepcopy
from inspect import isfunction
from multiprocessing import Pool
from pprint import pformat

import dask.array as da
import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.interpolate import interp1d
from scipy.stats import linregress

import esmvaltool.diag_scripts.emergent_constraints as ec
import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvalcore.cmor._fixes.shared import (add_plev_from_altitude,
                                           add_sigma_factory)
from esmvaltool.diag_scripts.shared import (ProvenanceLogger,
                                            get_diagnostic_filename,
                                            get_plot_filename, group_metadata,
                                            io, run_diagnostic,
                                            select_metadata)

logger = logging.getLogger(os.path.basename(__file__))


def _check_variables(datasets, necessary_short_names):
    """Check if ``datasets`` contain necessary variables."""
    dataset_name = datasets[0]['dataset']
    necessary_short_names = set(necessary_short_names)
    short_names = set(group_metadata(datasets, 'short_name').keys())
    if short_names != necessary_short_names:
        raise ValueError(
            f"Expected variables {necessary_short_names} for dataset "
            f"'{dataset_name}', got {short_names}")


# TODO: Remove when bug in iris is fixed
def _fix_iris_bug_derived_coord(cube):
    """Fix iris bug concerning derived coords and aggregation."""
    air_pressure_coord = cube.coord('air_pressure')
    try:
        altitude_coord = cube.coord('altitude')
    except iris.exceptions.CoordinateNotFoundError:
        altitude_coord = None
    for aux_factory in cube.aux_factories:
        cube.remove_aux_factory(aux_factory)
    try:
        cube.coord('air_pressure')
    except iris.exceptions.CoordinateNotFoundError:
        cube.add_aux_coord(air_pressure_coord, np.arange(cube.ndim))
    if altitude_coord is not None:
        cube.add_aux_coord(altitude_coord, (1, 2, 3))


def _get_cube(datasets, short_name):
    """Get cube with specific ``'short_name'`` from datasets."""
    datasets = select_metadata(datasets, short_name=short_name)
    if len(datasets) != 1:
        raise ValueError(
            f"Expected exactly one dataset with short_name '%s', got %d:\n%s",
            short_name, len(datasets), datasets)
    return iris.load_cube(datasets[0]['filename'],
                          ih.var_name_constraint(short_name))


def _get_level_width(air_pressure_bounds, ref_lev, ref_zg):
    """Get level widths of a single grid cell."""
    ref_lev = ref_lev.compressed()
    ref_zg = ref_zg.compressed()
    if len(ref_lev) < 2:
        return np.full(air_pressure_bounds.shape[0], np.nan)
    func = interp1d(ref_lev, ref_zg, kind='cubic', fill_value='extrapolate')
    level_widths = []
    for bounds in air_pressure_bounds:
        level_width = abs(func(bounds[0]) - func(bounds[1]))
        level_widths.append(level_width)
    return np.array(level_widths)


def _get_level_widths(cube, zg_cube):
    """Get all level widths for whole :class:`iris.cube.Cube`."""
    logger.info("Calculating level widths from 'air_pressure' coordinate")

    # Get air_pressure bounds
    (_, z_idx) = _get_z_coord(cube)
    air_pressure_coord = cube.coord('air_pressure')
    if air_pressure_coord.bounds is None:
        raise ValueError(
            f"Derived coordiante 'air_pressure' of cube "
            f"{cube.summary(shorten=True)} does not have bounds")
    if air_pressure_coord.shape == cube.shape:
        air_pressure_bounds = air_pressure_coord.bounds
    else:
        air_pressure_bounds = np.expand_dims(air_pressure_coord.bounds, 0)
        air_pressure_bounds = np.broadcast_to(air_pressure_bounds,
                                              cube.shape + (2, ))
    air_pressure_bounds = np.moveaxis(air_pressure_bounds, z_idx, -2)
    air_pressure_shape = air_pressure_bounds.shape[:-1]
    air_pressure_bounds = air_pressure_bounds.reshape(-1, cube.shape[z_idx], 2)

    # Geopotential height (pressure level -> altitude)
    (z_coord_zg, z_idx_zg) = _get_z_coord(zg_cube)
    ref_zg = np.moveaxis(zg_cube.data, z_idx_zg,
                         -1).reshape(-1, zg_cube.shape[z_idx_zg])
    mask = np.ma.getmaskarray(ref_zg)
    ref_lev = np.expand_dims(z_coord_zg.points, 0)
    ref_lev = np.ma.array(np.broadcast_to(ref_lev, ref_zg.shape), mask=mask)

    # Check shapes
    if air_pressure_bounds.shape[0] != ref_zg.shape[0]:
        raise ValueError(f"Expected identical first dimensions for cubes "
                         f"{cube.summary(shorten=True)} and "
                         f"{zg_cube.summary(shorten=True)}, got shapes "
                         f"{air_pressure_bounds.shape} and {ref_zg.shape}")

    # Cacluate level widths in parallel
    pool = Pool(64)
    level_widths = np.ma.masked_invalid(
        pool.starmap(_get_level_width, zip(air_pressure_bounds, ref_lev,
                                           ref_zg)))
    level_widths = level_widths.reshape(air_pressure_shape)
    level_widths = np.moveaxis(level_widths, -1, z_idx)
    return level_widths


def _get_level_width_coord(cube, zg_cube):
    """Get auxiliary coordinate which describes vertical level widths [m]."""
    try:
        altitude_coord = cube.coord('altitude')
    except iris.exceptions.CoordinateNotFoundError:
        level_widths = _get_level_widths(cube, zg_cube)
    else:
        logger.info("Calculating level widths from 'altitude' coordinate")
        if altitude_coord.bounds is None:
            raise ValueError(
                f"Height coordiante 'altitude' of cube "
                f"{cube.summary(shorten=True)} does not have bounds")
        level_widths = np.abs(altitude_coord.bounds[..., 1] -
                              altitude_coord.bounds[..., 0])
        if level_widths.shape != cube.shape:
            level_widths = np.expand_dims(level_widths, 0)
            level_widths = np.broadcast_to(level_widths, cube.shape)

    # Create coordinate
    aux_coord = iris.coords.AuxCoord(level_widths,
                                     var_name='level_width',
                                     long_name='Width of vertical layer',
                                     units='m')
    return aux_coord


def _get_mean_over_subsidence(cube, wap_cube, lat_constraint=None):
    """Get mean over subsidence regions."""
    if lat_constraint is not None:
        cube = cube.intersection(latitude=lat_constraint,
                                 longitude=(0.0, 360.0),
                                 ignore_bounds=True)
        wap_cube = wap_cube.intersection(latitude=lat_constraint,
                                         longitude=(0.0, 360.0),
                                         ignore_bounds=True)
    else:
        cube = cube.copy()
        wap_cube = wap_cube.copy()

    # Get monthly climatologies
    iris.coord_categorisation.add_month_number(cube, 'time')
    iris.coord_categorisation.add_month_number(wap_cube, 'time')
    cube = cube.aggregated_by('month_number', iris.analysis.MEAN)
    wap_cube = wap_cube.aggregated_by('month_number', iris.analysis.MEAN)

    # Mask subsidence regions (positive wap at 500 hPa)
    mask = da.where(wap_cube.core_data() > 0, False, True)
    cube.data = da.ma.masked_array(cube.core_data(), mask=mask)
    area_weights = iris.analysis.cartography.area_weights(cube)
    cube = cube.collapsed(['latitude', 'longitude'],
                          iris.analysis.MEAN,
                          weights=area_weights)
    return cube


def _get_seasonal_mblc_fraction(cl_cube, wap_cube, lat_constraint):
    """Calculate MBLC fraction."""
    cl_cube = cl_cube.intersection(latitude=lat_constraint,
                                   longitude=(0.0, 360.0),
                                   ignore_bounds=True)
    wap_cube = wap_cube.intersection(latitude=lat_constraint,
                                     longitude=(0.0, 360.0),
                                     ignore_bounds=True)

    # Calculate total cloud area fraction below 700 hPa
    levs = cl_cube.coord('air_pressure').core_points()
    mask = np.where(levs >= 70000, False, True)
    if mask.shape != cl_cube.shape:
        mask = np.broadcast_to(np.expand_dims(mask, 0), cl_cube.shape)
    cl_cube.data = da.ma.masked_array(cl_cube.core_data(), mask=mask)
    inv_cl_cube = cl_cube.copy(data=1.0 - cl_cube.core_data() / 100.0)
    (z_coord, z_idx) = _get_z_coord(inv_cl_cube)
    total_cl = (1.0 - inv_cl_cube.core_data().prod(axis=z_idx)) * 100.0
    clt_cube = inv_cl_cube.collapsed(z_coord, iris.analysis.MEAN)  # dummy
    clt_cube.data = total_cl
    clt_cube.cell_methods = clt_cube.cell_methods[:-1]

    # TODO: Remove when bug in iris is fixed
    for aux_factory in clt_cube.aux_factories:
        clt_cube.remove_aux_factory(aux_factory)

    # Get mean over subsidence regions
    return _get_mean_over_subsidence(clt_cube, wap_cube)


def _get_su_cube_dict(grouped_data, var_name, reference_datasets):
    """Extract cubes for Su et al. (2014) constraint."""
    ref_data = None

    # Reference data
    ref_filenames = []
    for ref_dataset_name in reference_datasets.split('|'):
        if ref_dataset_name not in grouped_data:
            raise ValueError(
                f"Reference dataset '{ref_dataset_name}' not found for "
                f"variable '{var_name}'")
        cube = iris.load_cube(grouped_data[ref_dataset_name][0]['filename'])
        if ref_data is None:
            ref_data = np.ma.array(cube.data)
        else:
            ref_data = np.ma.where(np.ma.getmaskarray(ref_data),
                                   np.ma.array(cube.data), ref_data)
        ref_filenames.append(grouped_data[ref_dataset_name][0]['filename'])
    ref_cube = cube.copy(ref_data)
    ref_cube.attributes['dataset'] = reference_datasets
    ref_cube.attributes['ancestors'] = '|'.join(ref_filenames)
    ref_cube.coord('air_pressure').attributes['positive'] = 'down'

    # All other cubes
    cube_dict = {reference_datasets: ref_cube}
    for (dataset_name, datasets) in grouped_data.items():
        if dataset_name in reference_datasets:
            continue
        cube = iris.load_cube(datasets[0]['filename'])
        cube.attributes['ancestors'] = datasets[0]['filename']
        cube_dict[dataset_name] = cube

    return cube_dict


def _get_su_variable(grouped_data):
    """Get variable and reference datasets of Su et al. (2014) constraint."""
    var_name = None
    reference_datasets = None
    for (dataset_name, datasets) in grouped_data.items():
        if len(datasets) != 1:
            raise ValueError(
                f"Expected exactly one file for dataset '{dataset_name}', got "
                f"{len(datasets):d}")
        new_var_name = datasets[0]['short_name']
        new_reference_datasets = datasets[0].get('reference_dataset')
        if var_name is None:
            var_name = new_var_name
        else:
            if new_var_name != var_name:
                raise ValueError(
                    f"Expected identical 'short_name' for all datasets of Su "
                    f"et al. (2014) constraint, got '{var_name}' and "
                    f"'{new_var_name}'")
        if reference_datasets is None:
            reference_datasets = new_reference_datasets
        else:
            if new_reference_datasets != reference_datasets:
                raise ValueError(
                    f"Expected identical 'reference_dataset' for all datasets "
                    f"of Su et al. (2014) constraint, got "
                    f"'{reference_datasets}' and '{new_reference_datasets}'")
    if reference_datasets is None:
        raise ValueError(f"'reference_dataset' not given for variable "
                         f"'{var_name}'")
    logger.info(
        "Found variable '%s' for Su et al. (2014) constraint", var_name)
    logger.info("Found reference datasets '%s'", reference_datasets)
    return (var_name, reference_datasets)


def _get_weighted_cloud_fractions(cl_cube, zg_cube, level_limits):
    """Calculate mass-weighted cloud fraction."""
    level_width_coord = _get_level_width_coord(cl_cube, zg_cube)
    cl_cube.add_aux_coord(level_width_coord, np.arange(cl_cube.ndim))

    # Mask data appropriately
    levs = cl_cube.coord('air_pressure')
    cloud_fractions = []
    for limits in level_limits:
        clt_cube = cl_cube.copy()
        mask = np.where(levs.points <= limits[0], False, True)
        mask |= np.where(levs.points >= limits[1], False, True)
        if mask.shape != clt_cube.shape:
            mask = np.broadcast_to(np.expand_dims(mask, 0), clt_cube.shape)
        clt_cube.data = da.ma.masked_array(clt_cube.core_data(), mask=mask)
        (z_coord, _) = _get_z_coord(clt_cube)

        # (Mass-weighted) vertical averaging
        clt_cube = clt_cube.collapsed(
            z_coord,
            iris.analysis.MEAN,
            weights=clt_cube.coord(var_name='level_width').points)

        # Temporal averaging
        clt_cube = clt_cube.collapsed('time', iris.analysis.MEAN)

        # (Area-weighted) horizontal averaging
        area_weights = iris.analysis.cartography.area_weights(clt_cube)
        clt_cube = clt_cube.collapsed(['latitude', 'longitude'],
                                      iris.analysis.MEAN,
                                      weights=area_weights)
        cloud_fractions.append(clt_cube.data)
    return cloud_fractions


def _get_z_coord(cube):
    """Get index of Z coordinate."""
    for coord in cube.coords(dim_coords=True):
        if iris.util.guess_coord_axis(coord) == 'Z':
            z_coord = coord
            break
    else:
        raise ValueError(f"Cannot determine height axis (Z) of cube "
                         f"{cube.summary(shorten=True)}")
    return (z_coord, cube.coord_dims(z_coord)[0])


def _get_zhai_data_frame(datasets, lat_constraint):
    """Get :class:`pandas.DataFrame` including the data for ``zhai``."""
    cl_cube = _get_cube(datasets, 'cl')
    wap_cube = _get_cube(datasets, 'wap')
    tos_cube = _get_cube(datasets, 'tos')

    # Add air_pressure coordinate if necessary
    if not cl_cube.coords('air_pressure'):
        if cl_cube.coords('altitude'):
            add_plev_from_altitude(cl_cube)
        elif cl_cube.coords('atmosphere_sigma_coordinate'):
            add_sigma_factory(cl_cube)
        else:
            raise ValueError(
                f"No 'air_pressure' coord available in cube "
                f"{cl_cube.summary(shorten=True)}")

    # Apply common mask (only ocean)
    mask_2d = da.ma.getmaskarray(tos_cube.core_data())
    mask_3d = mask_2d[:, np.newaxis, ...]
    mask_3d = da.broadcast_to(mask_3d, cl_cube.shape)
    wap_cube.data = da.ma.masked_array(wap_cube.core_data(), mask=mask_2d)
    cl_cube.data = da.ma.masked_array(cl_cube.core_data(), mask=mask_3d)

    # Calculate SST mean and MBLC fraction
    tos_cube = _get_mean_over_subsidence(tos_cube, wap_cube, lat_constraint)
    mblc_cube = _get_seasonal_mblc_fraction(cl_cube, wap_cube, lat_constraint)
    return pd.DataFrame(
        {'tos': tos_cube.data, 'mblc_fraction': mblc_cube.data},
        index=pd.Index(np.arange(12) + 1, name='month'),
    )


def _pearson_correlation_coeff(x_data, y_data):
    """Similarity metric using Pearson correlatin coefficient."""
    reg = linregress(x_data, y_data)
    return reg.rvalue


def _regression_slope_metric(x_data, y_data):
    """Similarity metric using the slope of a linear regression."""
    reg = linregress(x_data, y_data)
    return reg.slope


def _similarity_metric(cube, ref_cube, metric):
    """Calculate similarity metric between two cubes."""
    if metric == 'regression_slope':
        metric_func = _regression_slope_metric
    elif metric == 'correlation_coefficient':
        metric_func = _pearson_correlation_coeff
    else:
        raise ValueError(
            f"Expected one of 'regression_slope', 'correlation_coefficient' "
            f"similarity metric for diagnostic 'su', got '{metric}'")
    new_data = np.ma.array(cube.data, copy=True).ravel()
    ref_data = np.ma.array(ref_cube.data, copy=True).ravel()
    mask = np.ma.getmaskarray(ref_data) | np.ma.getmaskarray(new_data)
    return metric_func(np.ma.array(new_data, mask=mask).compressed(),
                       np.ma.array(ref_data, mask=mask).compressed())


def brient_shal(grouped_data, _):
    """Brient et al. (2016) constraint."""
    diag_data = {}

    # Variable attributes
    var_attrs = {
        'short_name': 'gamma',
        'long_name': 'Fraction of tropical (30°S - 30°N) low clouds with tops '
                     'below 850 hPa whose tops are also below 950 hPa (over '
                     'oceanic weak subsidence regions)',
        'units': '%',
    }
    attrs = {
        'plot_xlabel': r'Cloud shallowness index $\gamma$ [%]',
        'plot_title': 'Brient et al. (2016) constraint',
        'provenance_authors': ['schlund_manuel'],
        'provenance_domains': ['trop'],
        'provenance_realms': ['atmos'],
        'provenance_references': ['brient16climdyn'],
        'provenance_statistics': ['mean'],
        'provenance_themes': ['EC'],

    }

    # Calculate constraint
    for (dataset_name, datasets) in grouped_data.items():
        logger.info("Processing dataset '%s'", dataset_name)
        _check_variables(datasets, {'cl', 'wap', 'zg'})

        # Load cubes
        cl_cube = _get_cube(datasets, 'cl')
        wap_cube = _get_cube(datasets, 'wap')
        zg_cube = _get_cube(datasets, 'zg')

        # Add air_pressure coordinate if necessary
        if not cl_cube.coords('air_pressure'):
            if cl_cube.coords('altitude'):
                add_plev_from_altitude(cl_cube)
            elif cl_cube.coords('atmosphere_sigma_coordinate'):
                add_sigma_factory(cl_cube)
            else:
                raise ValueError(
                    f"No 'air_pressure' coord available in cube "
                    f"{cl_cube.summary(shorten=True)}")

        # TODO: Remove when bug in iris is fixed
        _fix_iris_bug_derived_coord(cl_cube)

        # Calculate monthly climatologies
        iris.coord_categorisation.add_month_number(cl_cube, 'time')
        iris.coord_categorisation.add_month_number(wap_cube, 'time')
        iris.coord_categorisation.add_month_number(zg_cube, 'time')
        cl_cube = cl_cube.aggregated_by('month_number', iris.analysis.MEAN)
        wap_cube = wap_cube.aggregated_by('month_number', iris.analysis.MEAN)
        zg_cube = zg_cube.aggregated_by('month_number', iris.analysis.MEAN)

        # Mask weak subsidence regions
        wap_cube.convert_units('hPa day-1')
        mask = np.where(wap_cube.data >= 10.0, False, True)
        mask |= np.where(wap_cube.data <= 30.0, False, True)
        cl_mask = np.broadcast_to(np.expand_dims(mask, 1), cl_cube.shape)
        cl_cube.data = da.ma.masked_array(cl_cube.core_data(), mask=cl_mask)
        wap_cube.data = da.ma.masked_array(wap_cube.core_data(), mask=mask)
        zg_mask = np.broadcast_to(np.expand_dims(mask, 1), zg_cube.shape)
        zg_cube.data = da.ma.masked_array(zg_cube.core_data(), mask=zg_mask)

        # Get mass-weighted cloud fractions
        [cf_950,
         cf_850] = _get_weighted_cloud_fractions(cl_cube, zg_cube,
                                                 [(100000, 90000),
                                                  (90000, 80000)])
        diag_data[dataset_name] = 100.0 * cf_950 / (cf_950 + cf_850)

    return (diag_data, var_attrs, attrs)


def su(grouped_data, cfg):
    """Su et al. (2014) constraint."""
    metric = cfg['metric']
    logger.info("Found metric '%s' for Su et al. (2014) constraint", metric)

    # Extract cubes
    (var_name, reference_datasets) = _get_su_variable(grouped_data)
    cube_dict = _get_su_cube_dict(grouped_data, var_name, reference_datasets)
    diag_data = {}
    ref_cube = cube_dict[reference_datasets]

    # Variable attributes
    var_attrs = {
        'short_name': 'alpha' if metric == 'regression_slope' else 'rho',
        'long_name': f"Error in vertically-resolved tropospheric "
                     f"zonal-average {ref_cube.long_name} between 40°N and "
                     f"45°S expressed as {metric.replace('_', ' ')} between "
                     f"model data and observations",
        'units': '1',
    }
    attrs = {
        'plot_xlabel': f'Model performance in {ref_cube.long_name} [1]',
        'plot_title': 'Su et al. (2014) constraint',
        'provenance_authors': ['schlund_manuel'],
        'provenance_domains': ['trop', 'midlat'],
        'provenance_realms': ['atmos'],
        'provenance_references': ['su14jgr'],
        'provenance_statistics': ['corr'],
        'provenance_themes': ['EC'],
    }

    # Calculate constraint
    for (dataset_name, cube) in cube_dict.items():
        logger.info("Processing dataset '%s'", dataset_name)

        # Plot cube
        if cube.ndim == 2:
            iris.quickplot.contourf(cube)
            filename = f"su_{dataset_name.replace('|', '_')}"
            plot_path = get_plot_filename(filename, cfg)
            plt.savefig(plot_path, **cfg['savefig_kwargs'])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Provenance
            netcdf_path = get_diagnostic_filename(filename, cfg)
            io.iris_save(cube, netcdf_path)
            ancestors = cube.attributes['ancestors'].split('|')
            provenance_record = ec.get_provenance_record(
                {'su': attrs}, ['su'],
                caption=f'{cube.long_name} for {dataset_name}.',
                plot_type='zonal', plot_file=plot_path, ancestors=ancestors)
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(netcdf_path, provenance_record)

        # Similarity metric
        diag_data[dataset_name] = _similarity_metric(cube, ref_cube, metric)

    return (diag_data, var_attrs, attrs)


def volodin(grouped_data, _):
    """Volodin (2008) constraint."""
    diag_data = {}

    # Variable attributes
    var_attrs = {
        'short_name': 'clt_diff',
        'long_name': 'Difference in total cloud fraction between tropics '
                     '(28°S - 28°N) and Southern midlatitudes (56°S - 36°S)',
        'units': '%',
    }
    attrs = {
        'plot_xlabel': 'Difference in tropical and midlatitude cloud fraction '
                       '[%]',
        'plot_title': 'Volodin (2008) constraint',
        'provenance_authors': ['schlund_manuel'],
        'provenance_domains': ['trop', 'shmidlat'],
        'provenance_realms': ['atmos'],
        'provenance_references': ['volodin08izvestiya'],
        'provenance_statistics': ['diff', 'mean'],
        'provenance_themes': ['EC'],
    }

    # Calculate constraint
    for (dataset_name, datasets) in grouped_data.items():
        logger.info("Processing dataset '%s'", dataset_name)
        _check_variables(datasets, {'clt'})

        # Check if tropical and midlatitudes clt is present
        if len(datasets) != 2:
            raise ValueError(
                f"Expected exactly two 'clt' datasets for dataset "
                f"'{dataset_name}', got {len(datasets):d}")
        for dataset in datasets:
            if ('trop' in dataset['variable_group']
                    or 'trop' in dataset['preprocessor']):
                trop = iris.load_cube(dataset['filename'])
                if trop.shape != ():
                    raise ValueError(
                        f"Expected scalar data for tropical 'clt' of dataset "
                        f"'{dataset_name}', got shape {trop.shape}")
                break
        else:
            raise ValueError(
                f"Expected exactly one dataset for tropical 'clt' (defined "
                f"by the string 'trop' in the variable group name or the "
                f"preprocessor name) for dataset '{dataset_name}', got none")
        for dataset in datasets:
            if ('midlat' in dataset['variable_group']
                    or 'midlat' in dataset['preprocessor']):
                midlat = iris.load_cube(dataset['filename'])
                if midlat.shape != ():
                    raise ValueError(
                        f"Expected scalar data for Southern midlatitudes "
                        f"'clt' of dataset '{dataset_name}', got shape "
                        f"{midlat.shape}")
                break
        else:
            raise ValueError(
                f"Expected exactly one dataset for Southern midlatitudes "
                f"'clt' (defined by the string 'midlat' in the variable group "
                f"name or the preprocessor name) for dataset "
                f"'{dataset_name}', got none")

        # Cloud fraction difference
        diag_data[dataset_name] = trop.data - midlat.data

    return (diag_data, var_attrs, attrs)


def zhai(grouped_data, cfg):
    """Zhai et al. (2015) constraint."""
    diag_data = {}

    # Variable attributes
    var_attrs = {
        'short_name': 'mblc_sst_response',
        'long_name': 'Response of seasonal Marine Boundary Layer Cloud (MBLC) '
                     'fraction to change in Sea Surface Temperature (SST) ',
        'units': '% K-1',
    }
    attrs = {
        'plot_xlabel': r'Response of MBLC fraction to SST changes '
                       r'[% K$^{-1}$]',
        'plot_title': 'Zhai et al. (2015) constraint',
        'provenance_authors': ['schlund_manuel'],
        'provenance_domains': ['trop', 'shmidlat'],
        'provenance_realms': ['atmos'],
        'provenance_references': ['zhai15grl'],
        'provenance_statistics': ['mean'],
        'provenance_themes': ['EC'],
    }

    # Calculate constraint
    for (dataset_name, datasets) in grouped_data.items():
        diag_data[dataset_name] = []
        logger.info("Processing dataset '%s'", dataset_name)
        _check_variables(datasets, {'cl', 'wap', 'tos'})

        # Consider both hemispheres seperately
        n_h = (20.0, 40.0)
        s_h = (-40.0, -20.0)
        for lat_constraint in (n_h, s_h):
            data_frame = _get_zhai_data_frame(datasets, lat_constraint)

            # MBLC fraction response to SST changes
            reg = linregress(data_frame['tos'].values,
                             data_frame['mblc_fraction'].values)
            diag_data[dataset_name].append(reg.slope)

            # Plot regression
            axes = sns.regplot(x='tos', y='mblc_fraction', data=data_frame)
            axes.text(
                0.05,
                0.05,
                rf"$\alpha={reg.slope:.3f}$ %/K ($R^2={reg.rvalue**2:.2f}$, "
                rf"$p={reg.pvalue:.4f}$)",
                transform=axes.transAxes)
            if lat_constraint == n_h:
                hem = 'Northern hemisphere'
                filename = f'zhai_{dataset_name}_nh'
            else:
                hem = 'Southern hemisphere'
                filename = f'zhai_{dataset_name}_sh'
            plot_path = get_plot_filename(filename, cfg)
            axes.set_title(f'{dataset_name} ({hem})')
            plt.savefig(plot_path, **cfg['savefig_kwargs'])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Provenance
            netcdf_path = get_diagnostic_filename(filename, cfg)
            cubes = iris.cube.CubeList([
                ec.pandas_object_to_cube(
                    data_frame['tos'], var_name='tos',
                    standard_name='sea_surface_temperature', units='K',
                    attributes={'region': hem}),
                ec.pandas_object_to_cube(
                    data_frame['mblc_fraction'], var_name='mblc_fraction',
                    long_name='Marine Boundary Layer Cloud fraction',
                    units='%', attributes={'region': hem}),
            ])
            io.iris_save(cubes, netcdf_path)
            provenance_record = ec.get_provenance_record(
                {'zhai': attrs}, ['zhai'],
                caption=f"Regression plot of 'mblc_fraction' vs 'tos' ({hem})",
                plot_type='scatter', plot_file=plot_path,
                ancestors=[d['filename'] for d in datasets])
            with ProvenanceLogger(cfg) as provenance_logger:
                provenance_logger.log(netcdf_path, provenance_record)

        # Mean over both hemispheres
        diag_data[dataset_name] = np.mean(diag_data[dataset_name])

    return (diag_data, var_attrs, attrs)


def check_cfg(cfg):
    """Check configuration :obj:`dict`."""
    necessary_options = ['diag']
    for opt in necessary_options:
        if opt not in cfg:
            raise KeyError(f"Necessary option '{opt}' not given")
    diag = cfg['diag']
    if diag not in globals() or not isfunction(globals()[diag]):
        raise ValueError(
            f"Selected diagnostic '{diag}' not available, it needs to be "
            f"implemented as a function of this diagnostic script")
    logger.info("Calculating constraint '%s'", diag)
    return diag


def check_input_data(input_data):
    """Check input data."""
    if not input_data:
        raise ValueError("No input data found")


def get_default_settings(cfg):
    """Get default configuration settings."""
    cfg.setdefault('metric', 'regression_slope')
    cfg.setdefault('savefig_kwargs', {
        'bbox_inches': 'tight',
        'dpi': 600,
        'orientation': 'landscape',
    })
    return cfg


def get_global_attributes(input_data, cfg):
    """Get attributes for psi cube for all datasets."""
    datasets = "|".join({str(d['dataset']) for d in input_data})
    projects = "|".join({str(d['project']) for d in input_data})
    ref = "|".join({str(d.get('reference_dataset')) for d in input_data})
    attrs = {
        'dataset': datasets,
        'project': projects,
        'reference_dataset': ref,
    }
    attrs.update(cfg.get('output_attributes', {}))
    return attrs


def main(cfg):
    """Run the diagnostic."""
    cfg = get_default_settings(cfg)
    diag = check_cfg(cfg)
    sns.set(**cfg.get('seaborn_settings', {}))

    # Get input data
    input_data = list(cfg['input_data'].values())
    input_data.extend(io.netcdf_to_metadata(cfg, pattern=cfg.get('pattern')))
    input_data = deepcopy(input_data)
    check_input_data(input_data)
    grouped_data = group_metadata(input_data, 'dataset')

    # Calculate X-axis of emergent constraint
    diag_func = globals()[diag]
    (diag_data, var_attrs, attrs) = diag_func(grouped_data, cfg)
    attrs.update(get_global_attributes(input_data, cfg))

    # Save data
    netcdf_path = get_diagnostic_filename(diag, cfg)
    io.save_scalar_data(diag_data, netcdf_path, var_attrs, attributes=attrs)
    logger.info("Found data:\n%s", pformat(diag_data))

    # Provenance
    provenance_record = ec.get_provenance_record(
        {diag: attrs}, [diag], caption=attrs['plot_xlabel'],
        ancestors=[d['filename'] for d in input_data])
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(netcdf_path, provenance_record)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
