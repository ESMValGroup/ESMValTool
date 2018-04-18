#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
   radiation assessment area
"""
# use Agg backend for running without X-server
import matplotlib as mpl
mpl.use('Agg')

import argparse
import csv
import datetime as dt
import os
import os.path
import pickle

import iris
import matplotlib.pyplot as plt
import numpy as np
from iris.experimental.equalise_cubes import equalise_attributes


AUTOASSESS_PARSER = argparse.ArgumentParser(description='Description')
AUTOASSESS_PARSER.add_argument('--ref-suite-id', required=True,
                               help='Reference suite ID.')
AUTOASSESS_PARSER.add_argument('--out-dir', required=True,
                               help='Output directory')
AUTOASSESS_PARSER.add_argument('--data-dir', required=True,
                               help='Directory containing model data.')
AUTOASSESS_PARSER.add_argument('--start-date', required=False,
                               help='Start date for assessment.')
AUTOASSESS_PARSER.add_argument('--end-date', required=False,
                               help='End date for assessment.')
AUTOASSESS_PARSER.add_argument('--obs-dir', required=False)
AUTOASSESS_PARSER.add_argument('--tmp-dir', required=False)
AUTOASSESS_PARSER.add_argument('--ancil-dir', required=False)
AUTOASSESS_PARSER.add_argument('--ncpu', default=1,
                               help='Number of available processors.')


def parse_args(parser):
    """Parse command line arguments, and check all paths are absolute."""
    args = parser.parse_args()

    for arg, val in vars(args).iteritems():
        if '_dir' in arg and not val is None and not os.path.isabs(val):
            msg = ' '.join(
                ['Cli argument ', str(arg), ' not absolute path:', str(val)])
            raise NotAbsolutePath(msg)
    return args


def unpickle_cubes(path):
    """Load cube list from path."""
    with open(path, 'r') as fh:
        cubes = pickle.load(fh)
    return cubes


class NotAbsolutePath(Exception):
    pass


def cube_time_mean(cubes):
    '''  Do time mean for cubes and coordinates bounds. '''

    time_meaned_cubes = iris.cube.CubeList()

    for cube in cubes:
        cube = cube.collapsed('time', iris.analysis.MEAN)
        time_meaned_cubes.append(cube)

    coord_sys = iris.coord_systems.GeogCS(iris.fileformats.pp.EARTH_RADIUS)
    for cube in time_meaned_cubes:
        if cube.coord_system() is None:
            cube.coord('latitude').coord_system = coord_sys
            cube.coord('longitude').coord_system = coord_sys

        cube.coord('latitude').guess_bounds()
        cube.coord('longitude').guess_bounds()

        lat_bounds = cube.coord('latitude').bounds.copy()
        lat_bounds[0, 0] = np.sign(lat_bounds[0, 0]) * 90.
        lat_bounds[-1, 1] = np.sign(lat_bounds[-1, 1]) * 90.
        cube.coord('latitude').bounds = lat_bounds

    return time_meaned_cubes


def area_weighted_mean(cube, mask=None):
    ''' Return a cube averaged over the total area.'''

    _area_weights = iris.analysis.cartography.area_weights(cube)
    if mask is not None:
        _area_weights = _area_weights * mask.data
    out = cube.collapsed(
        ['latitude', 'longitude'], iris.analysis.MEAN, weights=_area_weights)
    out.units = cube.units
    return out


def compute_energy_budget(cube_list, mask=None):
    '''
    Calculate budget elements from the supplied cube_list. Budget elements are

     0:radiation_net_toa,
     1:incoming_sw_rad_toa,          #cf-name : toa_incoming_shortwave_flux
     2:outgoing_sw_rad_toa           #cf-name : toa_outgoing_shortwave_flux
     3:total_downward_sw_surface,    #cf-name : surface_downwelling_shortwave_flux_in_air
     4:net_downward_sw_surface,      #cf-name : surface_net_downward_shortwave_flux
     5:upward_sw_reflected_surface
     6:sw_reflected_clouds,
     7:sw_absorbed_atm,
     8:outgoing_lw_rad_toa           #cf-name : toa_outgoing_longwave_flux
     9:downward_lw_absorbed_surface, #cf-name : surface_downwelling_longwave_flux_in_air
    10:net_downward_lw_surface,      #cf-name : surface_net_downward_longwave_flux/
    11:upward_lw_emitted_surface
    12:net_surface_radiation,
    13:sensible_heat,                #cf-name: surface_upward_sensible_heat_flux
    14:latent_heat                   #cf-name: surface_upward_latent_heat_flux
    15:radiation_adsorbed_surface,

    16:sw_cloud_forcing,
    17:lw_cloud_forcing,
    18:clear_sky_sw_rad_toa,         #cf-name: toa_outgoing_shortwave_flux_assuming_clear_sky
    19:clear_sky_lw_rad_toa,         #cf-name: toa_outgoing_longwave_flux_assuming_clear_sky

    All ratio are calculated in W m-2.

    '''
    incoming_sw_rad_toa_field = cube_list.extract_strict('toa_incoming_shortwave_flux')
    outgoing_sw_rad_toa_field = cube_list.extract_strict('toa_outgoing_shortwave_flux')
    total_downward_sw_surface_field = cube_list.extract_strict('surface_downwelling_shortwave_flux_in_air')
    clear_sky_sw_rad_toa_field = cube_list.extract_strict('toa_outgoing_shortwave_flux_assuming_clear_sky')
    clear_sky_lw_rad_toa_field = cube_list.extract_strict('toa_outgoing_longwave_flux_assuming_clear_sky')

    net_downward_sw_surface_field = cube_list.extract_strict('surface_net_downward_shortwave_flux')
    olr_field = cube_list.extract_strict('toa_outgoing_longwave_flux')
    # TODO standard name is 'surface_downwelling_longwave_flux_in_air'
    ilr_field = cube_list.extract_strict('surface_downwelling_longwave_flux')

    net_downward_lw_surface_field = cube_list.extract_strict('surface_net_downward_longwave_flux')
    sensible_heat_field = cube_list.extract_strict('surface_upward_sensible_heat_flux')
    latent_heat_field = cube_list.extract_strict('surface_upward_latent_heat_flux')

    # shortwave
    incoming_sw_rad_toa = area_weighted_mean(incoming_sw_rad_toa_field, mask)
    outgoing_sw_rad_toa = area_weighted_mean(outgoing_sw_rad_toa_field, mask)
    total_downward_sw_surface = area_weighted_mean(total_downward_sw_surface_field, mask)

    clear_sky_sw_rad_toa = area_weighted_mean(clear_sky_sw_rad_toa_field, mask)
    total_sw_cloud_forcing = outgoing_sw_rad_toa - clear_sky_sw_rad_toa

    net_downward_sw_surface = area_weighted_mean(net_downward_sw_surface_field, mask)
    upward_sw_reflected_surface_field = total_downward_sw_surface_field - net_downward_sw_surface_field
    upward_sw_reflected_surface = area_weighted_mean(upward_sw_reflected_surface_field, mask)
    sw_reflected_clouds_field = outgoing_sw_rad_toa_field - upward_sw_reflected_surface_field
    sw_reflected_clouds = area_weighted_mean(sw_reflected_clouds_field, mask)
    sw_absorbed_atm_field = incoming_sw_rad_toa_field - sw_reflected_clouds_field - total_downward_sw_surface_field
    sw_absorbed_atm = area_weighted_mean(sw_absorbed_atm_field, mask)

    # longwave
    outgoing_lw_rad_toa = area_weighted_mean(olr_field, mask)
    downward_lw_absorbed_surface = area_weighted_mean(ilr_field, mask)
    net_downward_lw_surface = area_weighted_mean(net_downward_lw_surface_field, mask)
    upward_lw_emitted_surface_field = ilr_field - net_downward_lw_surface_field
    upward_lw_emitted_surface = area_weighted_mean(upward_lw_emitted_surface_field, mask)

    clear_sky_lw_rad_toa = area_weighted_mean(clear_sky_lw_rad_toa_field, mask)
    total_lw_cloud_forcing = clear_sky_lw_rad_toa - outgoing_lw_rad_toa

    # net radiation
    net_surface_radiation_field = net_downward_sw_surface_field + net_downward_lw_surface_field
    net_surface_radiation = area_weighted_mean(net_surface_radiation_field, mask)

    # surface fluxes
    sensible_heat = area_weighted_mean(sensible_heat_field, mask)
    latent_heat = area_weighted_mean(latent_heat_field, mask)
    radiation_adsorbed_surface = net_surface_radiation - sensible_heat - latent_heat
    radiation_net_toa = incoming_sw_rad_toa - outgoing_sw_rad_toa - outgoing_lw_rad_toa

    # set long_name to a descriptive name (no CF standard name exists)
    radiation_net_toa.long_name = "radiation_net_toa"
    upward_sw_reflected_surface.long_name = "upward_sw_reflected_surface"
    sw_reflected_clouds.long_name = "sw_reflected_clouds"
    sw_absorbed_atm.long_name = "sw_absorbed_atm"
    upward_lw_emitted_surface.long_name = "upward_lw_emitted_surface"
    net_surface_radiation.long_name = "net_surface_radiation"
    radiation_adsorbed_surface.long_name = "radiation_adsorbed_surface"

    total_sw_cloud_forcing.long_name = "total_sw_cloud_forcing"
    total_lw_cloud_forcing.long_name = "total_lw_cloud_forcing"

    # set long_name to CF standard name
    incoming_sw_rad_toa.long_name = "toa_incoming_shortwave_flux"
    outgoing_sw_rad_toa.long_name = "toa_outgoing_shortwave_flux"
    total_downward_sw_surface.long_name = "surface_downwelling_shortwave_flux_in_air"
    net_downward_sw_surface.long_name = "surface_net_downward_shortwave_flux"
    outgoing_lw_rad_toa.long_name = "toa_outgoing_longwave_flux"
    downward_lw_absorbed_surface.long_name = "surface_downwelling_longwave_flux_in_air"
    net_downward_lw_surface.long_name = "surface_net_downward_longwave_flux"
    sensible_heat.long_name = "surface_upward_sensible_heat_flux"
    latent_heat.long_name = "surface_upward_latent_heat_flux"
    clear_sky_sw_rad_toa.long_name = "toa_outgoing_shortwave_flux_assuming_clear_sky"
    clear_sky_lw_rad_toa.long_name = "toa_outgoing_longwave_flux_assuming_clear_sky"

    return [
        radiation_net_toa, incoming_sw_rad_toa, outgoing_sw_rad_toa,
        clear_sky_sw_rad_toa, total_sw_cloud_forcing,
        total_downward_sw_surface, net_downward_sw_surface,
        upward_sw_reflected_surface, sw_reflected_clouds, sw_absorbed_atm,
        outgoing_lw_rad_toa, clear_sky_lw_rad_toa, total_lw_cloud_forcing,
        downward_lw_absorbed_surface, net_downward_lw_surface,
        upward_lw_emitted_surface, net_surface_radiation, sensible_heat,
        latent_heat, radiation_adsorbed_surface
    ]


def compute_clim(yearly_budget):
    '''Compute time average.'''

    yearly_budget = np.array(yearly_budget)
    clim_budget = []

    for j in range(yearly_budget.shape[1]):
        cube_list = iris.cube.CubeList(yearly_budget[:, j])
        equalise_attributes(cube_list)
        clim_budget.append(cube_list)

    return clim_budget


def climatological_energy_budget(cubes, outp_dir, MODEL, RESOL, YEARs, YEARe,
                                 season):
    '''Write the radiation budget data from model '''

    yearly_total_energy_budget = []

    fields = cube_time_mean(cubes)
    values = compute_energy_budget(fields)
    yearly_total_energy_budget.append(values)

    if not os.path.exists(outp_dir + '/output_' + MODEL):
        os.mkdir((outp_dir + '/output_' + MODEL), 0755)
    else:
        pass

    clim_total_energy_budget = compute_clim(yearly_total_energy_budget)
    model_output_file_name = str(
        outp_dir + '/output_' + MODEL + '/Global_Energy_Budget_' + MODEL + '_'
        + RESOL + '_' + str(YEARs) + '-' + str(YEARe) + '_' + season + '.txt')

    # "clim_total_energy_budget" are cubelist.
    with open(model_output_file_name, 'w') as fh:
        format_ = "%46s %7s %15.1f \n"
        for i in range(len(clim_total_energy_budget)):
            row_heading = clim_total_energy_budget[i][0].long_name
            units = clim_total_energy_budget[i][0].units
            budgetdata = clim_total_energy_budget[i][0].data
            fh.write(format_ %
                     (str(row_heading) + '', str(units) + ',', budgetdata))

    indir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'obs_data')
    stephens_data_tmp = []
    stephens_name = []
    stephens_data_err = []
    stephens_obs_file = 'Stephens_et_al_2012_obs_Energy_Budget.txt'
    with open(os.path.join(indir, stephens_obs_file), 'r') as stephens_f:
        for cols in [line.split() for line in stephens_f]:
            stephens_name.append(cols[0])
            stephens_data_tmp.append(cols[3])
            stephens_data_err.append(cols[4])

    demory_data_tmp = []
    demory_name = []
    demory_obs_file = 'Demory_et_al_2014_obs_Energy_Budget.txt'
    with open(os.path.join(indir, demory_obs_file), 'r') as demory_f:
        for cols in [line.split() for line in demory_f]:
            demory_name.append(cols[0])
            demory_data_tmp.append(cols[3])

    # Extract the CERES EBAF data
    obs_names = ['toa_incoming_shortwave_flux',
                 'toa_outgoing_shortwave_flux',
                 'toa_outgoing_shortwave_flux_assuming_clear_sky',
                 'surface_downwelling_shortwave_flux_in_air',
                 'toa_outgoing_longwave_flux',
                 'toa_outgoing_longwave_flux_assuming_clear_sky',
                 'surface_downwelling_longwave_flux_in_air']

    ceres_ebaf_list_tmp = []
    ceres_ebaf_data_dump = []

    # TODO local absolute path
    ceres_dir = '/project/earthobs/CERES/EBAF_CMIP5/supermeansnc/'
    ceres_filename = 'CERES-EBAF_L4_Ed2-6_200003-201002.' + season + '.nc'
    ceres_root = os.path.abspath(ceres_dir + ceres_filename)

    for obs_name in obs_names:
        ceres_ebaf_raw_cube = iris.load(ceres_root, obs_name)

        grid_areas = iris.analysis.cartography.area_weights(
            ceres_ebaf_raw_cube[0])
        avg_exp = ceres_ebaf_raw_cube[0].collapsed(
            ['latitude', 'longitude'], iris.analysis.MEAN, weights=grid_areas)
        avg_exp.long_name = obs_name
        if avg_exp.long_name == 'toa_outgoing_shortwave_flux_assuming_clear_sky':
            toa_out_sw_cs = avg_exp
        if avg_exp.long_name == 'toa_outgoing_shortwave_flux':
            toa_out_sw = avg_exp
        if avg_exp.long_name == 'toa_outgoing_longwave_flux_assuming_clear_sky':
            toa_out_lw_cs = avg_exp
        if avg_exp.long_name == 'toa_outgoing_longwave_flux':
            toa_out_lw = avg_exp

        ceres_ebaf_list_tmp.append(obs_name)
        ceres_ebaf_data_dump.append("%0.2f" % avg_exp.data)

    ebaf_sw_cloud_forcing = toa_out_sw - toa_out_sw_cs
    ebaf_lw_cloud_forcing = toa_out_lw_cs - toa_out_lw
    ebaf_sw_cloud_forcing.long_name = 'total_sw_cloud_forcing'
    ebaf_lw_cloud_forcing.long_name = 'total_lw_cloud_forcing'

    ceres_ebaf_list_tmp.append(ebaf_sw_cloud_forcing.long_name)
    ceres_ebaf_data_dump.append("%0.2f" % ebaf_sw_cloud_forcing.data)
    ceres_ebaf_list_tmp.append(ebaf_lw_cloud_forcing.long_name)
    ceres_ebaf_data_dump.append("%0.2f" % ebaf_lw_cloud_forcing.data)

    model_dump_tmp = []
    model_list_tmp = []
    with open(model_output_file_name, 'r') as model_f:
        for cols in [line.split() for line in model_f]:
            model_list_tmp.append(cols[0])
            model_dump_tmp.append(cols[3])

    ceres_ebaf_list = stephens_name[:]
    ceres_ebaf_data_tmp = stephens_data_tmp[:]

    model_list = stephens_name[:]
    model_data_tmp = stephens_data_tmp[:]

    for i in range(len(stephens_name)):
        ceres_ebaf_data_tmp[i] = float('NaN')
        model_data_tmp[i] = float('NaN')
        for j in range(len(ceres_ebaf_list_tmp)):
            if stephens_name[i] == ceres_ebaf_list_tmp[j]:
                ceres_ebaf_data_tmp[i] = ceres_ebaf_data_dump[j]
            else:
                print "There's no match list. Skip the CERES EBAF."  # TODO
        for k in range(len(model_list_tmp)):
            if stephens_name[i] == model_list_tmp[k]:
                model_data_tmp[i] = model_dump_tmp[k]
            else:
                print "There's no match list. Skip the MODEL."  # TODO

    model_data = [float(model_data) for model_data in model_data_tmp]
    stephens_data = [
        float(stephens_data) for stephens_data in stephens_data_tmp
    ]
    demory_data = [float(demory_data) for demory_data in demory_data_tmp]
    ceres_ebaf_data = [
        float(ceres_ebaf_data) for ceres_ebaf_data in ceres_ebaf_data_tmp
    ]
    stephens_error = [
        float(stephens_error) for stephens_error in stephens_data_err
    ]

    model_minus_demory = list(np.array(model_data) - np.array(demory_data))
    model_minus_ebf = list(np.array(model_data) - np.array(ceres_ebaf_data))
    model_minus_stephens = list(np.array(model_data) - np.array(stephens_data))

    n_groups = 20
    fig, ax = plt.subplots(figsize=(12, 8))
    index = np.arange(n_groups) * 2.0

    bar_width, opacity = 0.5, 0.4
    plt.ylim(-20, 20)
    rects1 = plt.bar(
        index + 0.2,
        model_minus_stephens,
        bar_width,
        alpha=opacity + 0.2,
        color='cornflowerblue',
        label=MODEL + '(' + season + ')' + ' - Stephens(2012)',
        yerr=stephens_error)
    rects2 = plt.bar(
        index + 0.2 + bar_width,
        model_minus_ebf,
        bar_width,
        alpha=opacity + 0.2,
        color='orange',
        label=MODEL + '(' + season + ')' + ' - CERES' + '(' + season + ')')
    rects3 = plt.bar(
        index + 0.2 + bar_width * 2,
        model_minus_demory,
        bar_width,
        alpha=opacity - 0.2,
        color='black',
        label=MODEL + '(' + season + ')' + ' - Demory(2014)')
    ax.spines['bottom'].set_position(('data', 0))
    ax.spines['top'].set_position(('data', 0))

    plt.xticks(
        index + bar_width + 0.5,
        model_list,
        ha='center',
        rotation=90,
        fontsize=10)
    plt.legend(frameon=False, fontsize=10, loc='upper left')

    # TODO: use static filename, and put additional information into figure
    # title
    plt.savefig(outp_dir + '/output_' + MODEL + '/Ebudget_' + MODEL + '_' +
                RESOL + '_' + str(YEARs) + '_' + str(
                    YEARe) + '_' + season + '.png')
    plt.close()


def main():
    "AutoAssess Radiation Budget assessment function."
    AREA_NAME = 'radiation_budget'
    args = parse_args(AUTOASSESS_PARSER)

    model_run_dirs = [
        os.path.join(args.data_dir, dir_) for dir_ in os.walk(args.data_dir).next()[1]
        ]

    cubes_paths = [
        os.path.join(dir_, AREA_NAME, 'cubes.pickle') for dir_ in model_run_dirs
        ]

    # load pickled cubes of all model runs
    cube_lists = []
    for cubes_path in cubes_paths:
        cube_lists.append(unpickle_cubes(cubes_path))

    # create output directories for single run assessments
    suite_ids = [cube_list[0].attributes['MODEL_RUN_ID'] for cube_list in cube_lists]

    # reference run is last
    suite_ids.remove(args.ref_suite_id)
    suite_ids.append(args.ref_suite_id)
    AREA_OUT_DIR = os.path.join(args.out_dir, '_vs_'.join(suite_ids), AREA_NAME)

    for suite_id in suite_ids:
        out_dir = os.path.join(AREA_OUT_DIR, suite_id)
        os.makedirs(out_dir)

    ### Assessment
    model = 'AMIP'
    resol = ''  # only used for output file names
    season = 'ann'
    date_format = '%Y/%m/%d'
    start_year = dt.datetime.strptime(args.start_date, date_format).year
    end_year = dt.datetime.strptime(args.end_date, date_format).year

    ### assess single runs
    # energy budget
    for cube_list in cube_lists:
        cube_list = work_around_radiation_forecast_period(cube_list)
        suite_id = cube_list[0].attributes['MODEL_RUN_ID']
        out_dir = os.path.join(AREA_OUT_DIR, suite_id)
        climatological_energy_budget(
            cube_list,
            out_dir,
            model,
            resol,
            start_year,
            end_year,
            season,
            )

    # radiation metrics
    import radiation_metrics as rm
    for suite_id in suite_ids:
        metrics = rm.global_rd(args.data_dir, suite_id, AREA_NAME)
        out_dir = os.path.join(AREA_OUT_DIR, suite_id)

        # write metrics to csv file for each suite_id
        with open(os.path.join(out_dir, 'metrics.csv'), 'w') as fh:
            writer = csv.writer(fh)
            for metric in metrics.items():
                writer.writerow(metric)


def work_around_radiation_forecast_period(cube_list):
    """
    forecast period is offset by 1h for some radiation cubes
    see UM ticket: https://code.metoffice.gov.uk/trac/um/ticket/2737
    Workaround: Delete forecast period coordinate
    """
    for cube in cube_list:
        cube.remove_coord('forecast_period')
    return cube_list


if __name__ == "__main__":
    main()
