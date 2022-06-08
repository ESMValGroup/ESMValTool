"""Diagnostic script to build climate patterns from CMIP6 models.

Description
-----------
Builds patterns, anomaly and climatology cubes from CMIP6 models.
This diagnostic needs preprocessed mean monthly cubes, with no
gridding requirements. Default re-grid specification exists to
decrease CPU-load and run-time.

Author
------
Gregory Munday (Met Office, UK)


Configuration options in recipe
-------------------------------
grid: str, optional (default: constrained)
    options: constrained, full
    def: removes Antarctica from grid
imogen_mode: bool, optional (default: off)
    options: on, off
    def: outputs extra data (anomaly, climatology) per variable
         to drive JULES-IMOGEN configuration
output_r2_scores: bool, optional (default: off)
    options: on, off
    def: outputs determinant values per variable to measure pattern robustness
parallelise: bool, optional (default: off)
    options: on, off
    def: parallelises code to run N models at once
parallel_threads: int, optional (default: null)
    options: any int, up to the amount of CPU cores accessible by user
    def: number of threads/cores to parallelise over. If 'parallelise: on' and
         'parallel_threads': null, diagnostic will automatically parallelise
         over N-1 accessible threads
"""

import logging
from pathlib import Path

import iris
import iris.coord_categorisation
import iris.cube
import numpy as np
import sklearn.linear_model
import sub_functions as sf
from plotting import plot_cp_timeseries, plot_patterns, plot_scores
from rename_variables import (
    rename_anom_variables,
    rename_clim_variables,
    rename_regression_variables,
    rename_variables_base,
)

from esmvaltool.diag_scripts.shared import run_diagnostic

logger = logging.getLogger(Path(__file__).stem)


def compute_diagnostic(filename):
    """Loads cube, removes any dimensions of length: 1.

    Parameters
    ----------
    filename (path): path to load cube file

    Returns
    -------
    cube (cube): a cube
    """
    logger.debug("Loading %s", filename)
    cube = iris.load_cube(filename)

    logger.debug("Running example computation")
    cube = iris.util.squeeze(cube)

    return cube


def climatology(cube):
    """Handles aggregation to make climatology.

    Parameters
    ----------
    cube (cube): cube loaded from config dictionary

    Returns
    -------
    cube_aggregated (cube): 40 year climatology cube from 1850-1889
    """
    cube_40yr = cube.extract(
        iris.Constraint(time=lambda t: 1850 <= t.point.year <= 1889,
                        month_number=lambda t: 1 <= t.point <= 12))
    cube_aggregated = make_monthly_climatology(cube_40yr)

    return cube_aggregated


def constrain_latitude(cube):
    """Constrains latitude to decrease run-time when output fed to IMOGEN.

    Parameters
    ----------
    cube (cube): cube loaded from config dictionary

    Returns
    -------
    cube_clipped (cube): cube with latitudes set from -55 to 82.5 degrees
    """
    cube_clipped = cube.extract(
        iris.Constraint(latitude=lambda cell: 82.5 >= cell >= -55))

    return cube_clipped


def make_monthly_climatology(cube):
    """Generates a climatology by month_number.

    Parameters
    ----------
    cube (cube): cube loaded from config dictionary

    Returns
    -------
    cube_month_climatol (cube): cube aggregated by month_number
    """
    if not cube.coords("month_number"):
        iris.coord_categorisation.add_month_number(cube, "time",
                                                   "month_number")
    cube_month_climatol = cube.aggregated_by("month_number",
                                             iris.analysis.MEAN)

    return cube_month_climatol


def diurnal_temp_range(cubelist):
    """Calculates diurnal range from daily max and min temperatures.

    Parameters
    ----------
    cubelist (cubelist): cubelist of tasmin and tasmax

    Returns
    -------
    range_cube (cube): cube of calculated diurnal range
    """
    range_cube = cubelist[0] - cubelist[1]
    range_cube.rename("Diurnal Range")
    range_cube.var_name = ("range_tl1")

    return range_cube


def calculate_diurnal_range(clim_list, ts_list):
    """Facilitates diurnal range calculation and appending.

    Parameters
    ----------
    clim_list (cubelist): cubelist of climatology cubes
    ts_list (cubelist): cubelist of standard timeseries cubes

    Returns
    -------
    clim_list_final (cubelist): cubelist of climatology cubes
        including diurnal range
    ts_list_final (cubelist): cubelist of standard timeseries cubes
        including diurnal range
    """
    temp_range_list_clim = iris.load([])
    temp_range_list_ts = iris.load([])
    comb_list = [clim_list, ts_list]

    for i in range(len(comb_list)):
        for cube in comb_list[i]:
            if (cube.var_name == "tasmax"
                    or cube.var_name == "tasmin") and cube in clim_list:
                temp_range_list_clim.append(cube)
            elif (cube.var_name == "tasmax"
                  or cube.var_name == "tasmin") and cube in ts_list:
                temp_range_list_ts.append(cube)
            else:
                pass

    derived_diurnal_clim = diurnal_temp_range(temp_range_list_clim)
    derived_diurnal_ts = diurnal_temp_range(temp_range_list_ts)

    # append diurnal range to lists
    clim_list_final, ts_list_final = append_diurnal_range(
        derived_diurnal_clim, derived_diurnal_ts, clim_list, ts_list)

    return clim_list_final, ts_list_final


def append_diurnal_range(derived_diurnal_clim, derived_diurnal_ts, clim_list,
                         ts_list):
    """Appends diurnal range to cubelists.

    Parameters
    ----------
    derived_diurnal_clim (cube): derived diurnal climatology cube
    derived_diurnal_ts (cube): derived diurnal timeseries cube
    clim_list (cubelist): existing climatology cubelist, no range
    ts_list (cubelist): existing timeseries cubelist, no range

    Returns
    -------
    clim_list_final (cubelist): cubelist of climatology cubes
        including diurnal range
    ts_list_final (cubelist): cubelist of standard timeseries cubes
        including diurnal range
    """
    # creating cube list without tasmax or tasmin
    # (since we just wanted the diurnal range)
    clim_list_final = iris.load([])
    ts_list_final = iris.load([])

    for cube in clim_list:
        if not (cube.var_name == "tasmax" or cube.var_name == "tasmin"):
            clim_list_final.append(cube)

    for cube in ts_list:
        if not (cube.var_name == "tasmax" or cube.var_name == "tasmin"):
            ts_list_final.append(cube)

    clim_list_final.append(derived_diurnal_clim)
    ts_list_final.append(derived_diurnal_ts)

    return clim_list_final, ts_list_final


def calculate_anomaly(clim_list, ts_list):
    """Calculates variables as anomalies, and adds diurnal range as variable.

    Parameters
    ----------
    clim_list (cubelist): cubelist of climatology variables
    ts_list (cubelist): cubelist of standard variable timeseries

    Returns
    -------
    clim_list_final (cubelist): cubelist of clim. vars, inc. diurnal range
    anom_list_final (cubelist): cubelist of anomaly vars, inc. diurnal range
    """
    # calculate diurnal temperature range cube
    clim_list_final, ts_list_final = calculate_diurnal_range(
        clim_list, ts_list)

    anom_list_final = ts_list_final.copy()

    # calc the anom by subtracting the monthly climatology from
    # the time series
    for i in range(len(ts_list_final)):
        i_months = anom_list_final[i].coord(
            'month_number').points - 1  # -1 because months are numbered 1..12
        anom_list_final[i].data -= clim_list_final[i][i_months].data

    return clim_list_final, anom_list_final


def regression(tas, cube_data):
    """Calculates the coefficients of regression between global surface temp
    and a variable.

    Parameters
    ----------
    tas (cube): near-surface air temperature
    cube_data (arr): cube.data array of a variable

    Returns
    -------
    slope_array (arr): array of grid cells with same shape as initial cube,
        containing the regression slope
    score_array (arr): array of grid cells with same shape as initial cube,
        containing the regression score as a measure of robustness
    """
    slope_array = np.full(tas.data.shape[1:], np.nan)
    score_array = np.full(tas.data.shape[1:], np.nan)

    # calculate global average
    tas_data = sf.area_avg(tas, return_cube=False)

    for i in range(tas.data.shape[1]):
        for j in range(tas.data.shape[2]):
            if tas.data[0, i, j] is not np.ma.masked:
                model = sklearn.linear_model.LinearRegression(
                    fit_intercept=False, copy_X=True)

                x = tas_data.reshape(-1, 1)
                y = cube_data[:, i, j]

                model.fit(x, y)
                slope_array[i, j] = model.coef_
                score_array[i, j] = model.score(x, y)

    return slope_array, score_array


def regression_units(tas, cube):
    """Calculates regression coefficient units.

    Parameters
    ----------
    tas (cube): near-surface air temperature
    cube (cube): cube of a given variable

    Returns
    -------
    units (str): string of calculated regression units
    """
    print('Cube Units: ', cube.units)
    units = cube.units / tas.units
    print('Regression Units: ', units)

    return units


def calculate_regressions(anom_list):
    """Facilitates the calculation of regression coefficients (climate
    patterns) and the creation of a new cube of patterns per variable.

    Parameters
    ----------
    anom_list (cubelist): cube list of variables as anomalies

    Returns
    -------
    regr_var_list (cubelist): cube list of newly created regression slope
        value cubes, for each variable
    score_list (cubelist): cube list of newly created regression score cubes,
        for each variable
    """
    regr_var_list = iris.cube.CubeList([])
    score_list = iris.cube.CubeList([])

    for cube in anom_list:
        if cube.var_name == 'tl1_anom':

            # getting last 85 years from full timeseries
            tas = cube[-1020:]

    for cube in anom_list:
        cube_ssp = cube[-1020:]
        month_list = iris.cube.CubeList([])
        score_month_list = iris.cube.CubeList([])

        # exctracting months, regressing, and merging
        for i in range(1, 13):
            month_constraint = iris.Constraint(imogen_drive=i)
            month_cube_ssp = cube_ssp.extract(month_constraint)
            month_tas = tas.extract(month_constraint)

            regr_array, score_array = regression(month_tas,
                                                 month_cube_ssp.data)

            # re-creating cube
            if (cube.var_name == 'swdown_anom'
                    or cube.var_name == 'lwdown_anom'):
                units = 'W m-2 K-1'
            else:
                units = regression_units(tas, cube_ssp)

            # assigning dim_coords
            coord1 = tas.coord(contains_dimension=1)
            coord2 = tas.coord(contains_dimension=2)
            dim_coords_and_dims = [(coord1, 0), (coord2, 1)]

            # assigning aux_coord
            coord_month = iris.coords.AuxCoord(i, var_name='imogen_drive')
            aux_coords_and_dims = [(coord_month, ())]

            cube = rename_regression_variables(cube)

            # creating cube of regression values
            regr_cube = iris.cube.Cube(regr_array,
                                       units=units,
                                       dim_coords_and_dims=dim_coords_and_dims,
                                       aux_coords_and_dims=aux_coords_and_dims,
                                       var_name=cube.var_name,
                                       standard_name=cube.standard_name)

            # calculating cube of r2 scores
            score_cube = iris.cube.Cube(
                score_array,
                units='R2',
                dim_coords_and_dims=dim_coords_and_dims,
                aux_coords_and_dims=aux_coords_and_dims,
                var_name=cube.var_name,
                standard_name=cube.standard_name)

            month_list.append(regr_cube)
            score_month_list.append(score_cube)

        conc_cube = month_list.merge_cube()
        regr_var_list.append(conc_cube)

        conc_score_cube = score_month_list.merge_cube()
        score_list.append(conc_score_cube)

    return regr_var_list, score_list


def write_scores(scores, work_path):
    """Saves the global average regression scores per variable in a text file.

    Parameters
    ----------
    scores (cubelist): cube list of regression score cubes, for each variable
    work_path (path): path to work_dir, to save scores

    Returns
    -------
    None
    """
    for cube in scores:
        score = sf.area_avg(cube, return_cube=False)
        mean_score = np.mean(score)

        # saving scores
        file = open(work_path + 'scores', 'a')
        data = '{0:10.3f}'.format(mean_score)
        name = cube.var_name
        file.write(name + ': ' + data + '\n')
        file.close()


def cube_saver(list_of_cubelists, work_path, name_list, mode):
    """Saves desired cubelists to work_dir, depending on switch settings.

    Parameters
    ----------
    list_of_cubelists (list): list containing desired cubelists
    work_path (path): path to work_dir, to save cubelists
    name_list (list): list of filename strings for saving
    mode (str): switch option passed through by ESMValTool config dict

    Returns
    -------
    None
    """
    if mode == 'imogen_scores':
        for i in range(0, 4):
            iris.save(list_of_cubelists[i], work_path + name_list[i])

    if mode == 'imogen':
        for i in range(0, 3):
            iris.save(list_of_cubelists[i], work_path + name_list[i])

    if mode == 'scores':
        for i in range(2, 4):
            for cube in list_of_cubelists[i]:
                rename_variables_base(cube)
            iris.save(list_of_cubelists[i], work_path + name_list[i])

    if mode == 'base':
        for cube in list_of_cubelists[2]:
            rename_variables_base(cube)
        iris.save(list_of_cubelists[2], work_path + name_list[2])


def saving_outputs(clim_list_final, anom_list_final, regressions, scores,
                   imogen_mode, r2_scores, plot_path, work_path):
    """Saves data and plots to relevant directories.

    Parameters
    ----------
    clim_list_final (cubelist): cube list of all variable climatologies
    anom_list_final (cubelist): cube list of all variable anomalies
    regressions (cubelist): cube list of all variable regression slopes
    scores (cubelist): cube list of all variable regression scores

    Returns
    -------
    None
    """
    list_of_cubelists = [clim_list_final, anom_list_final, regressions, scores]
    name_list = [
        "climatology_variables.nc", "anomaly_variables.nc", "patterns.nc",
        "scores.nc"
    ]

    # saving data + plotting
    if imogen_mode is True:
        if r2_scores is True:
            plot_scores(list_of_cubelists[3], plot_path)
            write_scores(scores, work_path)
            plot_cp_timeseries(list_of_cubelists, plot_path)
            cube_saver(list_of_cubelists,
                       work_path,
                       name_list,
                       mode='imogen_scores')

        elif r2_scores is False:
            plot_cp_timeseries(list_of_cubelists, plot_path)
            cube_saver(list_of_cubelists, work_path, name_list, mode='imogen')

    elif imogen_mode is False:
        if r2_scores is True:
            plot_scores(list_of_cubelists[3], plot_path)
            write_scores(scores, work_path)
            plot_patterns(list_of_cubelists[2], plot_path)
            cube_saver(list_of_cubelists, work_path, name_list, mode='scores')

        elif r2_scores is False:
            plot_patterns(list_of_cubelists[2], plot_path)
            cube_saver(list_of_cubelists, work_path, name_list, mode='base')


def patterns(model):
    """Driving function for script, taking in model data and saving parameters.

    Parameters
    ----------
    model (str): model name

    Returns
    -------
    None
    """
    with run_diagnostic() as cfg:
        input_data = cfg["input_data"].values()
        grid_spec = cfg["grid"]
        imogen_mode = cfg["imogen_mode"]
        r2_scores = cfg["output_r2_scores"]
        work_path = cfg["work_dir"] + "/"
        plot_path = cfg["plot_dir"] + "/"

    clim_list = iris.load([])
    ts_list = iris.load([])

    for dataset in input_data:
        if dataset['dataset'] == model:
            input_file = dataset["filename"]

            # preparing single cube
            cube_initial = compute_diagnostic(input_file)

            if grid_spec == "constrained":
                cube = constrain_latitude(cube_initial)
            elif grid_spec == "full":
                cube = cube_initial

            # appending to timeseries list
            ts_list.append(cube)

            # making climatology
            clim_cube = climatology(cube)
            clim_list.append(clim_cube)

    # calculate anomaly over historical + ssp timeseries
    clim_list_final, anom_list_final = calculate_anomaly(clim_list, ts_list)

    for i in range(len(clim_list_final)):
        rename_clim_variables(clim_list_final[i])
        rename_anom_variables(anom_list_final[i])

    regressions, scores = calculate_regressions(anom_list_final.copy())

    model_work_dir, model_plot_dir = sf.make_model_dirs(
        cube_initial, work_path, plot_path)

    saving_outputs(clim_list_final, anom_list_final, regressions, scores,
                   imogen_mode, r2_scores, model_plot_dir, model_work_dir)


def main(cfg):
    """Takes in driving data with parallelisation options.

    Parameters
    ----------
    cfg (dict): the global config dictionary, passed by ESMValTool.

    Returns
    -------
    None
    """
    input_data = cfg["input_data"].values()
    parallelise = cfg["parallelise"]
    threads = cfg["parallel_threads"]

    models = []
    for x in input_data:
        model = x['dataset']
        if model not in models:
            models.append(model)

    if parallelise is True:
        sf.parallelise(patterns, threads)(models)
    else:
        for model in models:
            patterns(model)


if __name__ == "__main__":

    with run_diagnostic() as config:
        main(config)
