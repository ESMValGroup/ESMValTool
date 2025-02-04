#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Diagnostic script to plot figure 10 of IPCC AR5 SPM.

Description
-----------
Calculate and plot the global mean surface air temperature increase as a
function of cumulative total global CO2 emissions (see IPCC AR5 WG1 fig.
SPM10).

Author
------
Manuel Schlund (DLR, Germany)

Project
-------
CRESCENDO

Configuration options in recipe
-------------------------------
save : dict, optional
    Keyword arguments for the `fig.saveplot()` function.
axes_functions : dict, optional
    Keyword arguments for the plot appearance functions.
dataset_style : str, optional
    Dataset style file (located in
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python`).
matplotlib_style : str, optional
    Dataset style file (located in
    :mod:`esmvaltool.diag_scripts.shared.plot.styles_python.matplotlib`).

"""

import logging
import os

import iris
import iris.coord_categorisation as iris_cat
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from cf_units import Unit
from scipy.spatial import ConvexHull

import esmvaltool.diag_scripts.shared.iris_helpers as ih
from esmvaltool.diag_scripts.shared import (
    get_plot_filename,
    group_metadata,
    plot,
    run_diagnostic,
    select_metadata,
    variables_available,
)

logger = logging.getLogger(os.path.basename(__file__))
plt.style.use(plot.get_path_to_mpl_style('small_font'))


def main(cfg):
    """Run the diagnostic."""
    input_data = cfg['input_data'].values()
    grouped_data = group_metadata(input_data, 'short_name')

    # Check if necessary variables are available
    necessary_vars = ['co2mass', 'tas', 'nbp_grid', 'fgco2_grid']
    if not variables_available(cfg, necessary_vars):
        raise ValueError(
            f"This diagnostic needs the variables {necessary_vars}")

    # Get name of historical experiment
    hist_exps = [
        name for name in group_metadata(input_data, 'exp').keys()
        if 'historical' in name
    ]
    if not hist_exps:
        raise ValueError("This diagnostic needs a historical experiment")
    hist_exp = hist_exps[0]
    if len(hist_exps) > 1:
        logger.warning("Found multiple historical experiments, using '%s'",
                       hist_exp)

    # Preprocess data
    max_decades = {}
    for (var, datasets) in grouped_data.items():
        logger.info("Processing %s", var)
        for dataset in datasets:
            cube = iris.load_cube(dataset['filename'])

            # co2mass: convert units from kgCO2 to kgC
            if var == 'co2mass':
                cube *= 12.0 / 44.0

            # tas: weighted mean
            elif var == 'tas':
                if dataset.get('fx_files', {}).get('areacella'):
                    area_cube = iris.load_cube(
                        dataset['fx_files']['areacella'])
                    area_weights = np.broadcast_to(area_cube.data, cube.shape)
                else:
                    for coord_name in ('latitude', 'longitude'):
                        coord = cube.coord(coord_name)
                        if not coord.has_bounds():
                            coord.guess_bounds()
                    area_weights = iris.analysis.cartography.area_weights(cube)
                    logger.warning(
                        "'areacella' file for experiment '%s' of '%s' not "
                        "found, estimated cell areas by bounds",
                        dataset['exp'], dataset['dataset'])
                cube = cube.collapsed(['latitude', 'longitude'],
                                      iris.analysis.MEAN,
                                      weights=area_weights)

            # nbp_grid and fgco2_grid: simple sum
            elif var in ('nbp_grid', 'fgco2_grid'):
                cube = cube.collapsed(['latitude', 'longitude'],
                                      iris.analysis.SUM)

            # Other variables not supported
            else:
                logger.warning("Variable '%s' not supported, skipping", var)
                continue

            def time_to_decade(coord, value):
                """Convert time into decade, e.g. 2050 means 2040-2049."""
                dec = 10
                year = coord.units.num2date(value).year
                return year + dec - year % dec

            # Add aux coordinate for for decade
            iris_cat.add_categorised_coord(
                cube, 'decade', 'time', time_to_decade, units='yr')

            # Decadal aggregation
            if var == 'co2mass':
                initial_co2mass = cube.data[0]

                def end_of_decade_value(data, axis):
                    """Return last element of decade."""
                    return np.take(data, -1, axis=axis)

                aggregator = iris.analysis.Aggregator('decade: point',
                                                      end_of_decade_value)
                cube = cube.aggregated_by('decade', aggregator)
                cube.data -= np.hstack((initial_co2mass, cube.data[:-1]))
            else:
                cube = cube.aggregated_by('decade', iris.analysis.MEAN)

            # Replace time coordinate
            time_idx = cube.coord_dims('time')
            new_coord = iris.coords.DimCoord.from_coord(cube.coord('decade'))
            cube.remove_coord('time')
            cube.remove_coord('decade')
            cube.add_dim_coord(new_coord, time_idx)

            # Add dataset coordinate and experiment
            cube.attributes = {}
            cube.cell_methods = ()
            for coord in cube.coords(dim_coords=False):
                cube.remove_coord(coord)
            dataset_coord = iris.coords.AuxCoord(
                dataset['dataset'], var_name='dataset', long_name='dataset')
            cube.add_aux_coord(dataset_coord, [])

            # nbp_grid and fgco2_grid: temporal integration
            if var in ('nbp_grid', 'fgco2_grid'):
                years = Unit('year')
                seconds = Unit('s')
                seconds_per_year = years.convert(10.0, seconds)
                cube *= seconds_per_year
                cube.units *= seconds

            # Save cube
            dataset['cube'] = cube

            # Save maximum number of decades
            exp = dataset['exp']
            if (max_decades.get(exp, np.array([])).size <
                    cube.coord('decade').shape[0]):
                max_decades[exp] = cube.coord('decade').points

    # Adapt time coordinate of 1pctCO2 run to historical run
    for dataset in select_metadata(input_data, exp='1pctCO2'):
        coord = dataset['cube'].coord('decade')
        start_year = max_decades[hist_exp][0]
        end_year = start_year + (coord.shape[0] - 1) * 10
        coord.points = np.linspace(start_year, end_year, coord.shape[0])

    # Unify decade coordinates
    for (exp, datasets) in group_metadata(input_data, 'exp').items():
        cubes = [d['cube'] for d in datasets]
        cubes = ih.unify_1d_cubes(cubes, 'decade')
        for (idx, dataset) in enumerate(datasets):
            dataset['cube'] = cubes[idx]

    # Calculate temperature anomalies for 1pctCO2 and historical
    # (relative to 1861-1880)
    tas = {}
    for (exp, datasets) in group_metadata(grouped_data['tas'], 'exp').items():
        tas[exp] = {}
        for dataset in datasets:
            cube = dataset['cube'].copy()

            # PiControl: relative to itself
            if exp == '1pctCO2':
                base_tas = np.ma.mean(cube.data[:2])

            # Others: relative to historical run
            else:
                hist_data = select_metadata(
                    grouped_data['tas'],
                    exp=hist_exp,
                    dataset=dataset['dataset'])
                if not hist_data:
                    logger.warning(
                        "Temperature data for historical experiment of '%s' "
                        "not available, skipping", dataset['dataset'])
                    continue
                base_tas = np.ma.mean(hist_data[0]['cube'].data[:2])
            cube -= base_tas

            # Drop first element for historical and 1pctCO2 run
            if exp in (hist_exp, '1pctCO2'):
                cube = cube[1:]

            # For future runs, include one earlier decade to match hist
            else:
                decade_coord = iris.coords.DimCoord(
                    np.insert(
                        cube.coord('decade').points, 0,
                        max_decades[hist_exp][-1]),
                    var_name='decade',
                    units='yr')
                dataset_coord = iris.coords.AuxCoord(
                    dataset['dataset'],
                    var_name='dataset',
                    long_name='dataset')
                metadata = cube.metadata
                cube = iris.cube.Cube(
                    np.ma.hstack([0.0, cube.data]),
                    dim_coords_and_dims=[(decade_coord, 0)],
                    aux_coords_and_dims=[(dataset_coord, [])])
                cube.metadata = metadata
            cube.attributes['exp'] = exp
            tas[exp][dataset['dataset']] = cube

    # Calculate cumultive CO2 emissions for every experiment
    # dC_E = dC_A + dC_L + dC_O = co2mass + nbp + fgco2
    co2_emiss = {}
    for (exp, datasets) in group_metadata(input_data, 'exp').items():
        co2_emiss[exp] = {}

        # CO2 mass is equal for every dataset
        co2_mass_datasets = select_metadata(datasets, short_name='co2mass')
        if co2_mass_datasets:
            co2_mass = co2_mass_datasets[0]['cube']
        else:
            logger.warning("No 'co2mass' data found for exp '%s', skipping",
                           exp)
            continue

        # nbp_grid and fgco2_grid
        for (name, name_datasets) in group_metadata(datasets,
                                                    'dataset').items():
            data_avail = True
            co2_emiss[exp][name] = co2_mass.copy()
            for var in ('nbp_grid', 'fgco2_grid'):
                var_data = select_metadata(name_datasets, short_name=var)
                if not var_data:
                    logger.warning(
                        "No '%s' data available for experiment '%s' of '%s', "
                        "skipping", var, exp, name)
                    data_avail = False
                else:
                    co2_emiss[exp][name] += var_data[0]['cube'].data
            if not data_avail:
                co2_emiss[exp].pop(name)
            else:
                cube = co2_emiss[exp][name]
                cube.data = cube.data.cumsum()

                # For historical run and 1pctCO2: relative to 1870
                if exp in (hist_exp, '1pctCO2'):
                    cube -= cube[0]
                    cube = cube[1:]

                # For future runs, include one earlier decade to match hist
                else:
                    decade_coord = iris.coords.DimCoord(
                        np.insert(
                            cube.coord('decade').points, 0,
                            max_decades[hist_exp][-1]),
                        var_name='decade',
                        units='yr')
                    dataset_coord = iris.coords.AuxCoord(
                        name, var_name='dataset', long_name='dataset')
                    metadata = cube.metadata
                    cube = iris.cube.Cube(
                        np.ma.hstack([0.0, cube.data]),
                        dim_coords_and_dims=[(decade_coord, 0)],
                        aux_coords_and_dims=[(dataset_coord, [])])
                    cube.metadata = metadata
                cube.coord('dataset').points = name
                cube.standard_name = None
                cube.long_name = 'Anthropogenic CO2 emissions'
                cube.convert_units('Gt')
                cube.attributes['exp'] = exp

                # Save cube and experiment
                cube.attributes['exp'] = exp
                co2_emiss[exp][name] = cube

    # Calculate multi-model means
    multi_model_sizes = {}
    for dict_ in (tas, co2_emiss):
        for exp in dict_:
            cubes = iris.cube.CubeList(dict_[exp].values())

            # Check length of ensembles
            if exp not in multi_model_sizes:
                multi_model_sizes[exp] = len(cubes)
            else:
                if multi_model_sizes[exp] != len(cubes):
                    logger.warning(
                        "Got different number of climate models for "
                        "multi-model averaging for temperature and CO2 "
                        "emissions of experiment '%s', expected %i, got %i",
                        exp, multi_model_sizes[exp], len(cubes))
            mmm_cube = cubes.merge_cube()
            mmm_cube = mmm_cube.collapsed('dataset', iris.analysis.MEAN)
            dict_[exp]['MultiModelMean'] = mmm_cube

    # Start future runs at temperature of historical run
    for exp in tas:
        if not exp.startswith('rcp'):
            continue
        for dataset in tas[exp]:
            if dataset not in tas[hist_exp]:
                logger.warning(
                    "Temperature for historical experiment of '%s' not "
                    "available, skipping", dataset)
                tas[exp].pop(dataset)
                continue
            tas[exp][dataset].data[0] = tas[hist_exp][dataset].data[-1]

    # Start future runs at emission level of historical run
    for exp in co2_emiss:
        if not exp.startswith('rcp'):
            continue
        for dataset in co2_emiss[exp]:
            if dataset not in co2_emiss[hist_exp]:
                logger.warning(
                    "CO2 emissions for historical experiment of '%s' not "
                    "available, skipping", dataset)
                co2_emiss[exp].pop(dataset)
                continue
            co2_emiss[exp][dataset].data += (
                co2_emiss[hist_exp][dataset].data[-1])

    # Calculate ranges (with given model range)
    model_range = cfg.get('range', 1.0)

    # Calculate historical and RCP range
    hist_rcp_tas = []
    hist_rcp_co2_emiss = []
    for exp in tas:
        if exp == '1pctCO2':
            continue
        valid_models_tas = []
        valid_models_co2_emiss = []
        for dataset in tas[exp]:
            if dataset not in co2_emiss[exp]:
                logger.warning(
                    "For temperature data of '%s' for experiment '%s': no "
                    "emission data available, skipping it for range "
                    "calculation (it is included in multi-model mean, though)",
                    dataset, exp)
                continue
            valid_models_tas.append(tas[exp][dataset].data)
            valid_models_co2_emiss.append(co2_emiss[exp][dataset].data)
        valid_models_tas = np.ma.array(valid_models_tas)
        valid_models_co2_emiss = np.ma.array(valid_models_co2_emiss)

        # Get all points in desired range for every dataset
        err = model_range * (np.ma.amax(valid_models_tas, axis=0) - np.ma.amin(
            valid_models_tas, axis=0)) * 0.5
        mean_tas = np.ma.mean(valid_models_tas, axis=0)
        mean_co2_emiss = np.ma.mean(valid_models_co2_emiss, axis=0)
        hist_rcp_tas.append(mean_tas - err)
        hist_rcp_tas.append(mean_tas + err)
        hist_rcp_co2_emiss.append(mean_co2_emiss)
        hist_rcp_co2_emiss.append(mean_co2_emiss)

    # Get all points on the range boundary
    hist_rcp_tas = np.ma.hstack(hist_rcp_tas).ravel()
    hist_rcp_co2_emiss = np.ma.hstack(hist_rcp_co2_emiss).ravel()
    hist_rcp_points = np.ma.vstack([hist_rcp_co2_emiss, hist_rcp_tas])
    hist_rcp_points = np.ma.swapaxes(hist_rcp_points, 0, 1)
    hist_rcp_points = np.ma.compress_rows(hist_rcp_points)
    hist_rcp_range = ConvexHull(hist_rcp_points)

    # Calculate 1pctCO2 range
    exp = '1pctCO2'
    valid_models_tas = []
    valid_models_co2_emiss = []
    for dataset in tas[exp]:
        if dataset not in co2_emiss[exp]:
            logger.warning(
                "For temperature data of '%s' for experiment '%s': no "
                "emission data available, skipping it for range "
                "calculation (it is included in multi-model mean, though)",
                dataset, exp)
            continue
        valid_models_tas.append(tas[exp][dataset].data)
        valid_models_co2_emiss.append(co2_emiss[exp][dataset].data)
    valid_models_tas = np.ma.array(valid_models_tas)
    valid_models_co2_emiss = np.ma.array(valid_models_co2_emiss)

    # Get all points in desired range for every dataset
    err = model_range * (np.ma.amax(valid_models_tas, axis=0) - np.ma.amin(
        valid_models_tas, axis=0)) * 0.5
    mean_tas = np.ma.mean(valid_models_tas, axis=0)
    onepct_mean_co2_emiss = np.ma.mean(valid_models_co2_emiss, axis=0)
    onepct_mean_tas_minus_err = mean_tas - err
    onepct_mean_tas_plus_err = mean_tas + err

    # Plot
    if cfg['write_plots']:
        (fig, axes) = plt.subplots()
        handles = []

        # Plot 1pctco2 range
        alph_1pctco2_range = 0.3
        c_1pctco2_range = 'black'
        axes.fill_between(
            onepct_mean_co2_emiss,
            onepct_mean_tas_minus_err,
            onepct_mean_tas_plus_err,
            color=c_1pctco2_range,
            alpha=alph_1pctco2_range,
            linewidth=0.0)

        # Plot historical and rcp range
        alph_hist_rcp_range = 0.3
        c_hist_rcp_range = 'red'
        axes.fill(
            hist_rcp_points[:, 0][hist_rcp_range.vertices],
            hist_rcp_points[:, 1][hist_rcp_range.vertices],
            color=c_hist_rcp_range,
            alpha=alph_hist_rcp_range,
            linewidth=0.0)

        # Years which should be plotted
        text_years = {
            '1pctCO2': [],
            hist_exp: [1890, 1950, 1980, 2000, 2010],
            'rcp26': [2030, 2050, 2100],
            'rcp45': [2030, 2050, 2100],
            'rcp60': [2050, 2100],
            'rcp85': [2050, 2100],
        }
        text_shifts = {
            '1pctCO2': [],
            hist_exp: [(50, -0.15), (-63, 0.14), (20, -0.12), (-160, 0.1),
                       (-160, 0.1)],
            'rcp26': [(-170, 0.1), (-150, 0.1), (0, -0.18)],
            'rcp45': [(20, -0.1), (30, -0.08), (20, -0.1)],
            'rcp60': [(20, -0.1), (40, 0)],
            'rcp85': [(-170, 0.1), (40, 0)],
        }

        # Plot all the selected models for all experiments
        for exp in tas:
            color = cfg.get(f'{exp}_color', 'red')
            linewidth = cfg.get(f'{exp}_linewidth', 2)
            marker = cfg.get(f'{exp}_marker', 'o')
            label = cfg.get(f'{exp}_label', exp)
            axes.plot(
                co2_emiss[exp]['MultiModelMean'].data,
                tas[exp]['MultiModelMean'].data,
                color=color,
                linewidth=linewidth,
                marker=marker,
                markerfacecolor=color,
                markeredgecolor=color)

            # Append historical and RCPs to legend
            if exp != '1pctCO2':
                handles.append(
                    mlines.Line2D([], [],
                                  color=color,
                                  label=label,
                                  linewidth=linewidth))

            # Add years
            for idx in range(len(text_years[exp])):
                txt_year = text_years[exp][idx]
                try:
                    x_idx = np.where(co2_emiss[exp]['MultiModelMean'].coord(
                        'decade').points == txt_year)[0]
                    y_idx = np.where(tas[exp]['MultiModelMean'].coord(
                        'decade').points == txt_year)[0]
                    x_pos = (co2_emiss[exp]['MultiModelMean'].data[x_idx][0] +
                             text_shifts[exp][idx][0])
                    y_pos = (tas[exp]['MultiModelMean'].data[y_idx][0] +
                             text_shifts[exp][idx][1])
                    axes.text(
                        x_pos,
                        y_pos,
                        '{:d}'.format(txt_year),
                        ha='left',
                        va='center',
                        color=color,
                        fontsize='smaller')
                except IndexError:
                    pass

        # Append historical and RCP range to legend
        handles.append(
            mpatches.Patch(
                color=c_hist_rcp_range,
                alpha=alph_hist_rcp_range,
                linewidth=0.0,
                label=cfg.get('historical_future_label', 'RCP range')))

        # Append 1pctCO2 to legend
        color = cfg.get('1pctCO2_color', 'red')
        linewidth = cfg.get('1pctCO2_linewidth', 2)
        marker = cfg.get('1pctCO2_marker', 'o')
        label = cfg.get('1pctCO2_label', '1pctCO2')
        handles.append(
            mlines.Line2D([], [],
                          color=color,
                          label=label,
                          linewidth=linewidth))
        handles.append(
            mpatches.Patch(
                color=c_1pctco2_range,
                alpha=alph_1pctco2_range,
                linewidth=0.0,
                label=cfg.get('1pctCO2_range_label', '1pctCO2 range')))

        # Create legend
        legend = axes.legend(
            handles=handles, loc='lower right', fontsize='smaller', ncol=2)

        # General plot appearance
        xlim_left = cfg.get('xlim_left', 0.0)
        xlim_right = cfg.get('xlim_right', 2500.0)
        ylim_top = cfg.get('ylim_top', 5.0)
        axes.set_xlim(left=xlim_left, right=xlim_right)
        axes.set_ylim(top=ylim_top)
        axes.set_xlabel(r'Cumulative total anthropogenic CO$_2$ emissions '
                        'from 1870 [GtC]')
        axes.set_ylabel(
            r'Temperature anomaly relative to 1861-1880 [$^\circ$C]')
        if cfg.get('plot_title'):
            axes.set_title(cfg['plot_title'])

        # Save file
        path = get_plot_filename('TCRE', cfg)
        logger.info("Wrote %s", path)
        fig.savefig(
            path,
            additional_artists=[legend],
            bbox_inches='tight',
            orientation='landscape')
        plt.close()


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
