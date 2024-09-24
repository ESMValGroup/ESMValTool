# weathertyping diagnostic for esmvaltool

# operating system manipulations (e.g. path constructions)
import json
import logging
import os
import warnings

# plotting imports
import cartopy.crs as ccrs
import cartopy.feature as cfeature
# to manipulate iris cubes
import iris
import iris.analysis.cartography
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
# general imports
import numpy as np
import pandas as pd
import seaborn as sns
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.colors import ListedColormap

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import ProvenanceLogger, group_metadata

iris.FUTURE.datum_support = True

logger = logging.getLogger(os.path.basename(__file__))

# Ignoring a warning that is produced when selecting timesteps of a weathertype
warnings.filterwarnings('ignore', '.*Collapsing a non-contiguous coordinate*')


def generate_grayscale_hex_values(x):
    """Generate grayscale values for plotting seasonal occurrence.

    Args:
        x (int): number of weathertypes

    Returns:
        np.list: array with grescale values as hex
    """
    grayscale_values = np.linspace(0, 1, x)
    grayscale_hex = [
        f'#{int(value * 255):02x}{int(value * 255):02x}{int(value * 255):02x}'
        for value in grayscale_values]

    return grayscale_hex


def plot_seasonal_occurrence(cfg: dict, wt_cubes: iris.cube.Cube,
                             dataset_name: str, ensemble: str = ""):
    """Plot relative monthly/seasonal occurrence of weathertypes.

    Args:
        cfg (dict): Configuration dict from recipe
        wt_cubes (iris.cube.Cube): Cube with weathertypes
        dataset_name (str): name of dataset
    """
    output_path = f'{cfg["plot_dir"]}/seasonal_occurrence'

    if not os.path.exists(f'{output_path}'):
        os.makedirs(f'{output_path}')

    month_list = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun',
                  'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']

    relative_occurrences = {}
    # {wt_string: {month: {wt1: rel_occurrence, wt2: rel_occurrence, ....}}}
    # first do absolute occurrence, then relative occurrence
    for cube in wt_cubes:
        month_dict = {}
        for month in range(1, 13):
            month_constraint = iris.Constraint(
                time=iris.time.PartialDateTime(month=month))
            array = cube.extract(month_constraint).data
            unique, counts = np.unique(array, return_counts=True)
            count_dict = dict(zip(unique, counts/sum(counts)))
            month_dict[month] = count_dict
        relative_occurrences[cube.long_name] = month_dict

    x = month_list

    for wt_string in relative_occurrences:
        wt_numbers = max(len(value) for value in
                         relative_occurrences.get(wt_string).values())
        colors = generate_grayscale_hex_values(wt_numbers)
        wt_stack = np.zeros((wt_numbers, 12))
        for month, month_value in relative_occurrences.get(wt_string).items():
            print(month_value)
            for wt in month_value.keys():
                print(month_value.get(wt))
                wt_stack[np.int8(wt-1), month-1] = month_value.get(wt)

        y = np.vstack(list(wt_stack))

        # plot
        _, ax = plt.subplots(figsize=(10, 10))

        ax.set_title(f'{dataset_name}')

        ax.stackplot(x, y, colors=colors)

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),
                  fancybox=True, shadow=True, ncol=9,
                  labels=tuple(f'WT {i+1}' for i in range(0, wt_numbers)))

        ax.set(xlim=(0, 11), xticks=np.arange(0, 12),
               ylim=(0, 1), yticks=np.arange(0, 1.1, 0.1))
        ax.set_xlabel('Month')
        ax.set_ylabel('Cumulative Relative occurrence')

        plt.savefig(f'{output_path}/{dataset_name}_{ensemble}_'
                    f'{wt_string}_rel_occurrence.png')
        plt.savefig(f'{output_path}/{dataset_name}_{ensemble}_'
                    f'{wt_string}_rel_occurrence.pdf')
        plt.close()


def get_cfg_vars(cfg: dict):
    """Get list of vars from configuration dict.

    Args:
        cfg (dict): Configuration dict from recipe.

    Returns:
        tuple: cfg vars
    """
    preproc_variables_dict = group_metadata(
        cfg.get('input_data').values(), 'dataset')

    correlation_threshold = cfg.get('correlation_threshold')
    rmse_threshold = cfg.get('rmse_threshold')
    work_dir = cfg.get('work_dir')
    plotting = cfg.get('plotting', False)
    automatic_slwt = cfg.get('automatic_slwt')
    predefined_slwt = cfg.get('predefined_slwt', False)

    return (preproc_variables_dict, correlation_threshold, rmse_threshold,
            work_dir, plotting, automatic_slwt, predefined_slwt)


def load_wt_preprocessors(dataset: str, preproc_variables_dict: dict):
    """Load preprocessor cubes for calculating Lamb and simplified
    weathertypes.

    Args:
        dataset (str): dataset name
        preproc_variables_dict (dict): dictionary of preprocessor variables

    Returns:
        _type_: list of preprocessor vars for weathertyping
    """
    wt_preproc = iris.load_cube(
        preproc_variables_dict.get(dataset)[0].get('filename'))
    wt_preproc_prcp = iris.load_cube(
        preproc_variables_dict.get(dataset)[1].get('filename'))
    wt_preproc_prcp_eobs = iris.load_cube(
        preproc_variables_dict.get('E-OBS')[0].get('filename'))

    return wt_preproc, wt_preproc_prcp, wt_preproc_prcp_eobs


def get_ancestors_era5_eobs(dataset: str, preproc_variables_dict: dict):
    """Get ancestors for ERA5/E-OBS.

    Args:
        dataset (str): dataset name
        preproc_variables_dict (dict): dictionary of preprocessor variables

    Returns:
        _type_: lists of ERA5/E-OBS ancestors
    """
    era5_ancestors = [
        preproc_variables_dict.get(dataset)[0].get('filename'),
        preproc_variables_dict.get(dataset)[1].get('filename')
    ]
    eobs_ancestors = [
        preproc_variables_dict.get(dataset)[0].get('filename'),
        preproc_variables_dict.get('E-OBS')[0].get('filename')
    ]
    return era5_ancestors, eobs_ancestors


def get_model_output_filepath(dataset: str, value: list):
    """Generate output filepath for model data.

    Args:
        dataset (str): Model name
        value (list): Model variables

    Returns:
        _type_: Output filepath and preprocessor path for
        future referencing.
    """
    model_name = dataset
    timerange = value[0].get('timerange').replace('/', '-')
    experiment = value[0].get('exp')
    ensemble = value[0].get('ensemble')

    output_file_path = f'{model_name}/{experiment}/{ensemble}/{timerange}'
    preproc_path = value[0].get('filename')

    return output_file_path, preproc_path


def get_preproc_lists(preproc_vars: list):
    """Put preprocessor variables and paths into list for further use.

    Args:
        preproc_vars (list): List of variables for specific
        dataset.

    Returns:
        _type_: List of preprocessor cubes for mean calculations
        as well as path to files.
    """
    preproc_path_psl = preproc_vars[-3].get('filename')
    preproc_path_prcp = preproc_vars[-2].get('filename')
    preproc_path_tas = preproc_vars[-1].get('filename')

    mean_preproc_psl = iris.load_cube(preproc_vars[-3].get('filename'))
    mean_preproc_prcp = iris.load_cube(preproc_vars[-2].get('filename'))
    mean_preproc_tas = iris.load_cube(preproc_vars[-1].get('filename'))

    preproc_list = [mean_preproc_psl, mean_preproc_prcp, mean_preproc_tas]
    preproc_path_list = [preproc_path_psl, preproc_path_prcp, preproc_path_tas]

    return preproc_list, preproc_path_list


def get_looping_dict(preproc_vars: list):
    """Put cubes into dictionary for further use.

    Args:
        preproc_vars (list): Values of dataset dictionary

    Returns:
        _type_: Dictionary of the form {'var': [preprov_var, preproc_path]}
    """
    preproc, preproc_path = get_preproc_lists(preproc_vars)
    dict_ = {
        'psl': [preproc[0], preproc_path[0]],
        'pr': [preproc[1], preproc_path[1]],
        'tas': [preproc[2], preproc_path[2]]
    }
    return dict_


def load_wt_files(path: str, only_lwt=False):
    """Load *.nc files of weathertype data. If only_lwt is true, only Lamb
    weathertypes will be loaded. (useful for automatic_slwt = False)

    Args:
        path (str): Path to weathertype data.
        only_lwt (bool, optional): If True,
        only Lamb weathertypes will be loaded. Defaults to False.
        (useful for automatic_slwt = False)

    Returns:
        _type_: List of weathertype cubes.
    """
    if not only_lwt:
        lwt_cube = iris.load_cube(path, 'lwt')
        slwt_era5_cube = iris.load_cube(path, 'slwt_era5')
        slwt_eobs_cube = iris.load_cube(path, 'slwt_eobs')
        wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]
    else:
        lwt_cube = iris.load_cube(path, 'lwt')
        wt_cubes = [lwt_cube]

    return wt_cubes


def plot_means(cfg: dict,
               preproc_var: np.array,
               wt_cubes: iris.cube.Cube,
               data_info: dict,
               only_lwt=False):
    """Wrapper function to plot various means/std/anomalies.

    Args:
        cfg (dict): cfg dictionary provided by recipe
        preproc_var (np.array): variable to be plotted
        wt_cubes (iris.cube.Cube): list of wt cubes
        data_info (dict): dictionary with info to dataset
        only_lwt (bool, optional): If True,
        only Lamb weathertypes will be loaded. Defaults to False.
        (useful for automatic_slwt = False)
    """

    if not only_lwt:
        data_info['wt_string'] = 'lwt'
        calc_wt_means(
            cfg,
            preproc_var,
            wt_cubes,
            data_info
        )
        data_info['wt_string'] = 'slwt_ERA5'
        calc_wt_means(
            cfg,
            preproc_var,
            wt_cubes,
            data_info
        )
        data_info['wt_string'] = 'slwt_EOBS'
        calc_wt_means(
            cfg,
            preproc_var,
            wt_cubes,
            data_info
        )
    else:
        data_info['wt_string'] = 'lwt'
        calc_wt_means(
            cfg,
            preproc_var,
            wt_cubes,
            data_info
        )


def get_provenance_record(caption: str, ancestors: list, long_names: list,
                          plot_types: bool | list,
                          statistics: bool | list) -> dict:
    """Get provenance record.

    Args:
        caption (str): Caption for plots
        ancestors (list): List of ancestors
        long_names (list): Variable long names
        plot_types (bool | list): Type of plot. Can be false
        if output is not a plot.
        statistics (bool | list): Types of statistics used.

    Returns:
        dict: Provenance dictionary.
    """
    record = {
        'caption': caption,
        'domains': ['reg'],
        'authors': ['jury_martin', 'kroissenbrunner_thomas'],
        'references': ['maraun21jgr', 'jones93ijc'],
        'projects': ['preval'],
        'long_names': long_names,
        'ancestors': ancestors,
    }
    if plot_types is not False:
        record['plot_types'] = plot_types
    if statistics is not False:
        record['statistics'] = statistics
    return record


def log_provenance(filename: str, cfg: dict, provenance_record: dict):
    """Log provenance. Produces xml file provenance info.

    Args:
        caption (str): Caption of plots.
        filename (str): Output name of provenance.
        cfg (dict): Configuration dictionary provided by recipe.
        ancestors (list): List of ancestors.
        long_names (list): Variable long_names
        plot_types (bool | list): Plot types. Can be false
        if output is not a plot.
        statistics (bool | list): Statistics used.
    """

    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    logger.info('Output stored as %s', filename)


def turn_list_to_mapping_dict(list_: list) -> dict:
    """Turns a list of combined WT to a dictionary for further processing.

    Args:
        list_ (list): List where entries are lists with related WT

    Returns:
        dict: Mapping dicitonary keys are simplified WT, values are Lamb WT
    """
    result_dict = {}

    for i, s in enumerate(list_):
        for elem in s:
            if elem not in result_dict:
                result_dict[elem] = i + 1
            else:
                result_dict[elem].append(i)

    return result_dict


def get_mapping_dict(selected_pairs: list) -> dict:
    """Get mapping dictionary from list of selected pairs.

    Args:
        selected_pairs (list): Selected pairs of WTs based on
        precipitation patterns over specified area and
        correlation and RSME thresholds defined in recipe.S

    Returns:
        dict: Mapping dicitonary keys are simplified WT, values are Lamb WT
    """
    mapping_array = []

    for i, elem in enumerate(selected_pairs):
        mapping_array.append(elem[0])

    s = [set(i) for i in mapping_array if i]

    def find_intersection(m_list: list) -> list:
        for i, v in enumerate(m_list):
            for j, k in enumerate(m_list[i + 1:], i + 1):
                if v & k:
                    s[i] = v.union(m_list.pop(j))
                    return find_intersection(m_list)
        return m_list

    merged_tuples = find_intersection(s)
    mapping_dict = turn_list_to_mapping_dict(merged_tuples)

    return mapping_dict


def write_mapping_dict(work_dir: str, dataset: str, mapping_dict: dict):
    mapping_dict_reformat = convert_dict(mapping_dict)

    with open(f'{work_dir}/wt_mapping_dict_{dataset}.json',
              'w',
              encoding='utf-8') as file:
        json.dump(mapping_dict_reformat, file)


def calc_slwt_obs(cfg: dict, lwt: np.array, cube: iris.cube.Cube,
                  dataset: str, ancestors: list) -> np.array:
    """Calculate simplified weathertypes for observation datasets based on
    precipitation patterns over specified area.

    Args:
        cfg (dict): Configuration dictionary from recipe
        lwt (np.array): Array of Lamb WT
        cube (iris.cube.Cube): preprocessor cube to keep time coordinate
        dataset (str): Name of dataset
        correlation_thresold (float): correlation_threshold
        rmse_threshold (float): rsme_threshold
        ancestors (list): list of ancestors

    Returns:
        np.array: _description_
    """

    logger.info('Calculating simplified Lamb Weathertypes for %s', dataset)

    work_dir = cfg.get('work_dir')
    correlation_threshold = cfg.get('correlation_threshold')
    rmse_threshold = cfg.get('rmse_threshold')
    tcoord = cube.coord('time')

    wt_data_prcp = []
    for wt in range(1, 28):
        target_indices = np.where(lwt == wt)
        if len(target_indices[0]) < 1:
            logger.info(
                'calc_slwt_obs - CAUTION: Skipped wt %s \
                for dataset %s!', wt, dataset)
            continue
        dates = [
            tcoord.units.num2date(tcoord.points[i]) for i in target_indices
        ]
        if dataset == 'E-OBS':
            extracted_cube = cube[target_indices]
        else:
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0]))
        wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)
        wt_data_prcp.append(wt_cube_mean.data.compressed())
    selected_pairs = process_prcp_mean(cfg, wt_data_prcp,
                                       correlation_threshold,
                                       rmse_threshold, dataset)

    with open(f'{work_dir}/wt_selected_pairs_{dataset}.json',
              'w',
              encoding='utf-8') as file:
        json.dump(selected_pairs, file)

    mapping_dict = get_mapping_dict(selected_pairs)

    write_mapping_dict(work_dir, dataset, mapping_dict)

    provenance_record = get_provenance_record('Lamb Weathertypes',
                                              ancestors,
                                              ['Lamb Weathertypes'],
                                              False, False)

    log_provenance(f'{dataset}_wt_prov', cfg, provenance_record)

    return map_lwt_to_slwt(lwt, mapping_dict)


def calc_const():
    '''Calculate constants for weathertyping algorithm.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Returns:
        tuple: The four constants needed for WT calculation.
    '''

    const1 = 1 / np.cos(45 * np.pi / 180)
    const2 = np.sin(45 * np.pi / 180) / np.sin(40 * np.pi / 180)
    const3 = np.sin(45 * np.pi / 180) / np.sin(50 * np.pi / 180)
    const4 = 1 / (2 * np.cos(45 * np.pi / 180)**2)

    return const1, const2, const3, const4


def calc_westerly_flow(cube: iris.cube.Cube) -> np.array:
    '''Calculate the westerly flow over area.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): Cube of psl data.

    Returns:
        np.array: westerly flow
    '''

    return 1 / 2 * (cube.data[:, 1, 2] + cube.data[:, 1, 4]) - 1 / 2 * (
        cube.data[:, 3, 2] + cube.data[:, 3, 4])


def calc_southerly_flow(cube: iris.cube.Cube, const1: float) -> np.array:
    '''Calculate the southerly flow over area.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): Cube of psl data.
        const1 (float): const1

    Returns:
        np.array: southerly flow
    '''

    return const1 * (
        1 / 4 *
        (cube.data[:, 3, 4] + 2 * cube.data[:, 2, 4] + cube.data[:, 1, 4]) -
        1 / 4 *
        (cube.data[:, 3, 2] + 2 * cube.data[:, 2, 2] + cube.data[:, 1, 2]))


def calc_resultant_flow(w: np.array, s: np.array) -> np.array:
    '''Calculate the resultant flow.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        w (np.array): westerly flow.
        s (np.array): southerly flow

    Returns:
        np.array: resultant flow
    '''
    return (s**2 + w**2)**(1 / 2)


def calc_westerly_shear_velocity(cube: iris.cube.Cube, const2: float,
                                 const3: float) -> np.array:
    '''Calculate westerly shear velocity.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): cube of psl data
        const2 (float): const2
        const3 (float): const3

    Returns:
        np.array: westerly shear velocity
    '''
    return const2 * (1 / 2 *
                     (cube.data[:, 0, 2] + cube.data[:, 0, 4]) - 1 / 2 *
                     (cube.data[:, 2, 2] + cube.data[:, 2, 4])) - const3 * (
                         1 / 2 *
                         (cube.data[:, 2, 2] + cube.data[:, 2, 4]) - 1 / 2 *
                         (cube.data[:, 4, 2] + cube.data[:, 4, 4]))


def calc_southerly_shear_velocity(cube: iris.cube.Cube,
                                  const4: float) -> np.array:
    '''Calculate southerly shear velocity.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): cube of psl data
        const4 (float): const4

    Returns:
        np.array: southerly shear velocity
    '''
    return const4 * (
        1 / 4 *
        (cube.data[:, 3, 6] + 2 * cube.data[:, 2, 6] + cube.data[:, 1, 6]) -
        1 / 4 *
        (cube.data[:, 3, 4] + 2 * cube.data[:, 2, 4] + cube.data[:, 1, 4]) -
        1 / 4 *
        (cube.data[:, 3, 2] + 2 * cube.data[:, 2, 2] + cube.data[:, 1, 2]) +
        1 / 4 *
        (cube.data[:, 3, 0] + 2 * cube.data[:, 2, 0] + cube.data[:, 1, 0]))


def calc_total_shear_velocity(zw: np.array, zs: np.array) -> np.array:
    '''Calculate total shear velocity.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        zw (np.array): westerly shear velocity
        zs (np.array): southerly shear velocity

    Returns:
        np.array: total shear velocity
    '''
    return zw + zs


def wt_algorithm(cube: iris.cube.Cube, dataset: str) -> np.array:
    '''Algorithm to calculate Lamb weathertypes.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993),
    A comparison of Lamb circulation types with an objective classification
    scheme.
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): PSL field of dataset
        dataset (str): Name of dataset

    Returns:
        np.array: Array of Lamb WT for each day
    '''

    # lats and lons corresponding to datapoints
    # 55, 5 -> 1
    # 55, 15 -> 2
    # 50, -5 -> 3
    # 50, 5 -> 4
    # 50, 15 -> 5
    # 50, 25 -> 6
    # 45, -5 -> 7
    # 45, 5 -> 8
    # 45, 15 -> 9
    # 45, 25 -> 10
    # 40, -5 -> 11
    # 40, 5 -> 12
    # 40, 15 -> 13
    # 40, 25 -> 14
    # 35, 5 -> 15
    # 35, 15 -> 16

    # lons: -5, 0, 5, 10, 15, 20, 25
    # lats: 35, 40, 45, 50, 55

    logger.info('Calculating Lamb Weathertypes for %s', dataset)

    const1, const2, const3, const4 = calc_const()

    # westerly flow
    w = calc_westerly_flow(cube)
    # southerly flow
    s = calc_southerly_flow(cube, const1)
    # resultant flow
    f = calc_resultant_flow(w, s)
    # westerly shear vorticity
    zw = calc_westerly_shear_velocity(cube, const2, const3)
    # southerly shear vorticity
    zs = calc_southerly_shear_velocity(cube, const4)
    # total shear vorticity
    z = calc_total_shear_velocity(zw, zs)

    wt = np.zeros(len(z))

    for i, z_i in enumerate(z):

        direction = np.arctan(w[i] / s[i]) * 180 / np.pi  # deg
        if s[i] >= 0:
            direction += 180  # deg

        if direction < 0:
            direction += 360  # deg

        # Lamb pure directional type
        if abs(z_i) < f[i]:
            if 337.5 <= direction or direction < 22.5:
                wt[i] = 1
            elif 22.5 <= direction < 67.5:
                wt[i] = 2
            elif 67.5 <= direction < 112.5:
                wt[i] = 3
            elif 112.5 <= direction < 157.5:
                wt[i] = 4
            elif 157.5 <= direction < 202.5:
                wt[i] = 5
            elif 202.5 <= direction < 247.5:
                wt[i] = 6
            elif 247.5 <= direction < 292.5:
                wt[i] = 7
            elif 292.5 <= direction < 337.5:
                wt[i] = 8
        # Lamb’s pure cyclonic and anticyclonic type
        elif (2 * f[i]) < abs(z_i):
            if z_i > 0:
                wt[i] = 9

            elif z_i < 0:
                wt[i] = 10
        # Lambs’s synoptic/direction hybrid types
        elif f[i] < abs(z_i) < (2 * f[i]):
            if z_i > 0:
                if 337.5 <= direction or direction < 22.5:
                    wt[i] = 11
                elif 22.5 <= direction < 67.5:
                    wt[i] = 12
                elif 67.5 <= direction < 112.5:
                    wt[i] = 13
                elif 112.5 <= direction < 157.5:
                    wt[i] = 14
                elif 157.5 <= direction < 202.5:
                    wt[i] = 15
                elif 202.5 <= direction < 247.5:
                    wt[i] = 16
                elif 247.5 <= direction < 292.5:
                    wt[i] = 17
                elif 292.5 <= direction < 337.5:
                    wt[i] = 18

            elif z_i < 0:
                if 337.5 <= direction or direction < 22.5:
                    wt[i] = 19
                elif 22.5 <= direction < 67.5:
                    wt[i] = 20
                elif 67.5 <= direction < 112.5:
                    wt[i] = 21
                elif 112.5 <= direction < 157.5:
                    wt[i] = 22
                elif 157.5 <= direction < 202.5:
                    wt[i] = 23
                elif 202.5 <= direction < 247.5:
                    wt[i] = 24
                elif 247.5 <= direction < 292.5:
                    wt[i] = 25
                elif 292.5 <= direction < 337.5:
                    wt[i] = 26
        # light indeterminate flow, corresponding to Lamb’s unclassified type U
        elif abs(z_i) < 6 and f[i] < 6:
            wt[i] = 27

    return wt


def calc_lwt_slwt_model(cfg: dict, cube: iris.cube.Cube,
                        data_info: dict,
                        predefined_slwt: bool | dict):
    """Calculate Lamb WT and simplified WT for model data.

    Args:
        cfg (dict): Configuration dicitonary from recipe
        cube (iris.cube.Cube): PSL field of dataset
        dataset (str): Name of dataset
        preproc_path (str): Path of ancestors
        output_file_path (str): Path to write output file
    """

    work_dir = cfg.get('work_dir')
    dataset = data_info.get('dataset')
    preproc_path = data_info.get('preproc_path')
    output_file_path = data_info.get('output_file_path')

    if not os.path.exists(f'{work_dir}/{output_file_path}'):
        os.makedirs(f'{work_dir}/{output_file_path}')

    lwt = wt_algorithm(cube, dataset)

    tcoord = cube.coord('time')
    time_points = tcoord.units.num2date(tcoord.points)

    logger.info('Calculating simplified Lamb Weathertypes for %s', dataset)

    if not predefined_slwt:
        with open(f'{work_dir}/wt_mapping_dict_ERA5.json',
                  'r',
                  encoding='utf-8') as file:
            mapping_dict_era5_f = json.load(file)

        with open(f'{work_dir}/wt_mapping_dict_E-OBS.json',
                  'r',
                  encoding='utf-8') as file:
            mapping_dict_eobs_f = json.load(file)
        mapping_dict_era5 = reverse_convert_dict(mapping_dict_era5_f)
        mapping_dict_eobs = reverse_convert_dict(mapping_dict_eobs_f)

        slwt_era5 = map_lwt_to_slwt(lwt, mapping_dict_era5)
        slwt_eobs = map_lwt_to_slwt(lwt, mapping_dict_eobs)
    else:
        predefined_slwt = check_mapping_dict_format(predefined_slwt)
        write_mapping_dict(work_dir, 'ERA5', predefined_slwt)
        write_mapping_dict(work_dir, 'E-OBS', predefined_slwt)
        slwt_era5 = map_lwt_to_slwt(lwt, predefined_slwt)
        slwt_eobs = map_lwt_to_slwt(lwt, predefined_slwt)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(lwt, long_name='lwt'))
    wt_cube.append(iris.cube.Cube(slwt_era5, long_name='slwt_era5'))
    wt_cube.append(iris.cube.Cube(slwt_eobs, long_name='slwt_eobs'))

    wt_cube[0].add_dim_coord(tcoord, 0)
    wt_cube[1].add_dim_coord(tcoord, 0)
    wt_cube[2].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f'{work_dir}/{output_file_path}/{dataset}.nc')

    # write to csv file
    d = {
        'date': time_points[:],
        'lwt': np.int8(lwt),
        'slwt_ERA5': np.int8(slwt_era5),
        'slwt_EOBS': np.int8(slwt_eobs),
    }
    df = pd.DataFrame(data=d)
    df.to_csv(f'{work_dir}/{output_file_path}/{dataset}.csv', index=False)

    ancestors = [
        preproc_path, f'{work_dir}/wt_mapping_dict_ERA5.json',
        f'{work_dir}/wt_mapping_dict_E-OBS.json'
    ]
    provenance_record = get_provenance_record('Lamb Weathertypes',
                                              ancestors,
                                              ['Lamb Weathertypes'],
                                              False, False)

    log_provenance(f'{dataset}_wt_prov', cfg, provenance_record)


def get_colormap(colormap_string: str) -> ListedColormap:
    """Get colormaps based on string.

    Args:
        colormap_string (str): String to get Colormaps for either
                                psl, tas or precipitation.

    Returns:
        ListedColormap: Choosen Colormap
    """

    misc_seq_2_disc = [
        (230 / 255, 240 / 255, 240 / 255),
        (182 / 255, 217 / 255, 228 / 255),
        (142 / 255, 192 / 255, 226 / 255),
        (118 / 255, 163 / 255, 228 / 255),
        (116 / 255, 130 / 255, 222 / 255),
        (121 / 255, 97 / 255, 199 / 255),
        (118 / 255, 66 / 255, 164 / 255),
        (107 / 255, 40 / 255, 121 / 255),
        (86 / 255, 22 / 255, 75 / 255),
        (54 / 255, 14 / 255, 36 / 255),
    ]

    temp_seq_disc = [
        (254 / 255, 254 / 255, 203 / 255),
        (251 / 255, 235 / 255, 153 / 255),
        (244 / 255, 204 / 255, 104 / 255),
        (235 / 255, 167 / 255, 84 / 255),
        (228 / 255, 134 / 255, 80 / 255),
        (209 / 255, 98 / 255, 76 / 255),
        (164 / 255, 70 / 255, 66 / 255),
        (114 / 255, 55 / 255, 46 / 255),
        (66 / 255, 40 / 255, 24 / 255),
        (25 / 255, 25 / 255, 0 / 255),
    ]

    prec_seq_disc = [
        (255 / 255, 255 / 255, 229 / 255),
        (217 / 255, 235 / 255, 213 / 255),
        (180 / 255, 216 / 255, 197 / 255),
        (142 / 255, 197 / 255, 181 / 255),
        (105 / 255, 177 / 255, 165 / 255),
        (67 / 255, 158 / 255, 149 / 255),
        (44 / 255, 135 / 255, 127 / 255),
        (29 / 255, 110 / 255, 100 / 255),
        (14 / 255, 85 / 255, 74 / 255),
        (0 / 255, 60 / 255, 48 / 255),
    ]

    if colormap_string == 'psl':
        return ListedColormap(misc_seq_2_disc)
    if colormap_string == 'prcp':
        return ListedColormap(prec_seq_disc)
    if colormap_string == 'temp':
        return ListedColormap(temp_seq_disc)

    return None


def plot_maps(wt: np.array, cfg: dict, cube: iris.cube.Cube,
              data_info: dict, mode: str):
    """Plot maps.

    Args:
        wt (np.array): WT for which statistic was calculated
        cfg (dict): Configuration dicitonary from recipe
        cube (iris.cube.Cube): Data to be plotted
        dataset (str): Name of dataset
        var_name (str): Name of variable to be plotted
        wt_string (str): lwt, slwt_ERA5 or slwt_EOBS
                        slwt are calculated based on ERA5 or EOBS
                        precipitation data, respectively
        mode (str): Statistics that is used
    """

    dataset = data_info.get('dataset')
    var_name = data_info.get('var')
    wt_string = data_info.get('wt_string')
    ensemble = data_info.get('ensemble', '')

    logger.info('Plotting %s %s %s for %s %s', dataset, var_name, mode,
                wt_string, wt)

    local_path = f"{cfg.get('plot_dir')}/{mode}"

    if not os.path.exists(f'{local_path}'):
        os.makedirs(f'{local_path}')

    ax = plt.axes(projection=ccrs.PlateCarree())

    if var_name == 'psl':
        psl_cmap = get_colormap('psl')
        plt.title(f'{var_name} {mode}, wt: {wt}')
        unit = '[hPa]'
        im = iplt.contourf(cube / 100, cmap=psl_cmap)
    elif var_name == 'pr':
        prcp_cmap = get_colormap('prcp')
        if dataset == 'ERA5':
            unit = '[m]'
            plt.title(f'total {var_name} {mode}, wt: {wt}')
        else:
            unit = '[kg m-2 s-1]'
            plt.title(f'{var_name} flux {mode}, wt: {wt}')
        im = iplt.contourf(cube, cmap=prcp_cmap)
    elif var_name == 'tas':
        temp_cmap = get_colormap('temp')
        unit = '[K]'
        plt.title(f'1000 hPa {var_name} {mode}, wt: {wt}')
        im = iplt.contourf(cube, cmap=temp_cmap)

    cb = plt.colorbar(im)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(label=f'{var_name} {mode} {unit}')

    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=0.5,
        color='gray',
        alpha=0.5,
        linestyle='--',
    )
    gl.left_labels = True
    gl.bottom_labels = True
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = True
    gl.ylocator = mticker.FixedLocator(np.arange(20, 70, 5))
    gl.xlocator = mticker.FixedLocator([-10, -5, 0, 5, 10, 15])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'color': 'black', 'size': 8}

    ax.set_extent([-15, 20, 27.5, 62.5])

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    plt.savefig(f'{local_path}/{wt_string}_{wt}_{dataset}_{ensemble}'
                f'_{var_name}_{mode}.png')
    plt.savefig(f'{local_path}/{wt_string}_{wt}_{dataset}_{ensemble}_'
                f'{var_name}_{mode}.pdf')
    plt.close()


def rmse(subarray1: np.array, subarray2: np.array) -> np.array:
    """Calculate rsme.

    Args:
        subarray1 (np.array): array1
        subarray2 (np.array): array2

    Returns:
        np.array: rsme array
    """
    return np.sqrt(np.mean((subarray1 - subarray2)**2))


def convert_dict(dict_: dict) -> dict:
    """Convert mapping dictionary from {lwt: slwt, ...} format to {slwt: [lwt1,
    lwt2], ...}.

    Args:
        dict_ (dict): Dict in the {lwt: slwt, ...} format

    Returns:
        dict: Dict in the {slwt: [lwt1, lwt2], ...} format
    """

    new_dict = {}
    for dataset, value in dict_.items():
        if value not in new_dict:
            new_dict[value] = []
        new_dict[value].append(dataset)
    return new_dict


def reverse_convert_dict(dict_: dict) -> dict:
    """Convert mapping dictionary from {slwt: [lwt1, lwt2], ...} format to
    {lwt: slwt, ...}.

    Args:
        original_dict (dict): Dict in the {slwt: [lwt1, lwt2], ...}format

    Returns:
        dict: Dict in the  format {lwt: slwt, ...}
    """
    new_dict = {}
    for key, value_list in dict_.items():
        for original_key in value_list:
            new_dict[original_key] = key
    return new_dict


def plot_corr_rmse_heatmaps(cfg: dict, pattern_correlation_matrix: np.array,
                            rmse_matrix: np.array, dataset: str):
    """Plot heatmaps for correlation and rmse matrices.

    Args:
        cfg (dict): cfg dict from recipe
        pattern_correlation_matrix (np.array): pattern correlation matrix
        rmse_matrix (np.array): rmse matrix
        dataset (str): string of dataset
    """

    output_path = f'{cfg.get("plot_dir")}/heatmaps'

    if not os.path.exists(f'{output_path}'):
        os.makedirs(f'{output_path}')

    labels = np.arange(1, 28)

    mask = np.zeros_like(pattern_correlation_matrix)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style('white'):
        plt.figure(figsize=(10, 10))
        plt.title('Correlation Matrix')
        levels = np.linspace(np.min(pattern_correlation_matrix),
                             np.max(pattern_correlation_matrix), 9)
        ax = sns.heatmap(pattern_correlation_matrix,
                         mask=mask,
                         square=True,
                         annot=True,
                         annot_kws={'size': 6},
                         cmap='seismic',
                         xticklabels=labels,
                         yticklabels=labels,
                         cbar_kws={
                             'ticks': levels,
                             'shrink': 0.8,
                             'format': '%.2f'
                         })
        ax.set_xlabel('lwt', fontsize=8)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ax.set_ylabel('lwt', fontsize=8)
        plt.tight_layout()
        plt.savefig(f'{output_path}/correlation_matrix_{dataset}.png')
        plt.savefig(f'{output_path}/correlation_matrix_{dataset}.pdf')
        plt.close()

    mask = np.zeros_like(rmse_matrix)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style('white'):
        plt.figure(figsize=(10, 10))
        plt.title('RMSE Matrix')
        levels = np.linspace(np.min(rmse_matrix), np.max(rmse_matrix), 9)
        ax = sns.heatmap(rmse_matrix,
                         mask=mask,
                         square=True,
                         annot=True,
                         annot_kws={'size': 5},
                         cmap='seismic',
                         xticklabels=labels,
                         yticklabels=labels,
                         cbar_kws={
                             'ticks': levels,
                             'shrink': 0.8,
                             'format': '%.2f'
                         })
        ax.set_xlabel('lwt', fontsize=8)
        ax.set_xticklabels(ax.get_xticklabels(), rotation=0)
        ax.set_yticklabels(ax.get_yticklabels(), rotation=0)
        ax.set_ylabel('lwt', fontsize=8)
        plt.tight_layout()
        plt.savefig(f'{output_path}/rmse_matrix_{dataset}.png')
        plt.savefig(f'{output_path}/rmse_matrix_{dataset}.pdf')
        plt.close()


def write_corr_rmse_to_csv(cfg: dict, pattern_correlation_matrix: np.array,
                           rmse_matrix: np.array, dataset: str):
    """Write correlation and rsme matrix to csv files.

    Args:
        cfg (dict): Configuration dictionary from recipe
        pattern_correlation_matrix (np.array): Correlation matrix
        rmse_matrix (np.array): RSME matrix
        dataset (str): Name of dataset
    """

    logger.info('Writing corr and rsme matrices for %s', dataset)

    work_dir = cfg.get('work_dir')

    df_corr = pd.DataFrame(pattern_correlation_matrix)
    df_corr.index = range(1, len(df_corr) + 1)
    df_corr.columns = range(1, len(df_corr.columns) + 1)
    df_corr.to_csv(f'{work_dir}/correlation_matrix_{dataset}.csv',
                   index_label='Index')

    df_rmse = pd.DataFrame(rmse_matrix)
    df_rmse.index = range(1, len(df_rmse) + 1)
    df_rmse.columns = range(1, len(df_rmse.columns) + 1)
    df_rmse.to_csv(f'{work_dir}/rmse_matrix_{dataset}.csv',
                   index_label='Index')


def process_prcp_mean(cfg: dict, data: np.array, correlation_threshold: float,
                      rmse_threshold: float, dataset: str) -> list:
    """Process precipitation fields for specified area to get a list of
    selected pairs of weathertypes with the highest correlation (higher than
    correlation_threshold) and smallest RSME (smaller than rsme_threshold) for
    further processing and simplifying the WT.

    Args:
        cfg (dict): Configuration dictionary from recipe
        data (np.array): Precipitation data
        correlation_threshold (float): Correlation threshold
        rmse_threshold (float): RMSE threshold
        dataset (str): Name of dataset

    Returns:
        list: Selected pairs of WT. This is passed to get_mapping_dict
    """

    logger.info('Calculating corr and rsme matrices for %s', dataset)

    selected_pairs = []
    pattern_correlation_matrix = np.ma.corrcoef(data)
    n = len(data)
    rmse_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            rmse_matrix[i][j] = rmse(data[i], data[j])
            rmse_matrix[j][i] = rmse_matrix[i][j]
            if (pattern_correlation_matrix[i][j] >= correlation_threshold
                    and rmse_matrix[i][j] <= rmse_threshold):
                selected_pairs.append(
                    ((i + 1, j + 1), pattern_correlation_matrix[i][j],
                     rmse_matrix[i][j]))

    # write matrices to csv
    write_corr_rmse_to_csv(cfg, pattern_correlation_matrix, rmse_matrix,
                           dataset)
    # plot heatmaps for matrices
    plot_corr_rmse_heatmaps(cfg, pattern_correlation_matrix, rmse_matrix,
                            dataset)

    return selected_pairs


def calc_wt_means(cfg: dict, cube: iris.cube.Cube,
                  wt_cubes: iris.cube.CubeList, data_info: dict):
    """Calculate means of fields of each weathertype.

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube with variable data
        wt_cubes (iris.cube.CubeList): List of cubes of lwt, slwt_ERA5
                                        and slwt_EOBS
        dataset (str): Name of dataset
        var_name (str): Name of variable
        wt_string (str): lwt, slwt_ERA5 or slwt_EOBS
                        slwt are calculated based on ERA5 or EOBS
                        precipitation data, respectively
        preproc_path (str): Ancestor path
    """

    dataset = data_info.get('dataset')
    var_name = data_info.get('var')
    wt_string = data_info.get('wt_string')
    preproc_path = data_info.get('preproc_path')

    logger.info('Calculating %s %s means for %s', dataset, var_name, wt_string)

    work_dir = cfg.get('work_dir')

    if wt_string == 'slwt_ERA5':
        slwt_era5_cube = wt_cubes[1]
        tcoord = slwt_era5_cube.coord('time')
        slwt_era5 = slwt_era5_cube.data[:]
        num_slwt = len(np.unique(slwt_era5))
    elif wt_string == 'slwt_EOBS':
        slwt_eobs_cube = wt_cubes[2]
        tcoord = slwt_eobs_cube.coord('time')
        slwt_eobs = slwt_eobs_cube.data[:]
        num_slwt = len(np.unique(slwt_eobs))
    elif wt_string == 'lwt':
        lwt_cube = wt_cubes[0]
        tcoord = lwt_cube.coord('time')
        lwt = lwt_cube.data[:]

    # num_slwt_era5 = len(np.unique(slwt_era5))
    # num_slwt_eobs = len(np.unique(slwt_eobs))

    # if num_slwt_eobs != num_slwt_era5:
    #    logger.info('calc_wt_means - CAUTION: unequal number of \
    #                slwt_era5 (%s) \ and slwt_eobs (%s)!',
    #                num_slwt_era5, num_slwt_eobs)

    if 'slwt' in wt_string:
        for wt in range(1, num_slwt + 1):
            if wt_string == 'slwt_ERA5':
                target_indices = np.where(slwt_era5 == wt)
            elif wt_string == 'slwt_EOBS':
                target_indices = np.where(slwt_eobs == wt)
            else:
                logger.info('WT_STRING not supported!')
            if len(target_indices[0]) < 1:
                logger.info(
                    'calc_wt_means - CAUTION: Skipped %s %s \
                    for dataset %s!', wt_string, wt, dataset)
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)
            plot_maps(wt, cfg, wt_cube_mean, data_info, 'mean')
    elif wt_string == 'lwt':
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    'calc_wt_means - CAUTION: Skipped lwt %s \
                    for dataset %s!', wt, dataset)
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)
            plot_maps(wt, cfg, wt_cube_mean, data_info, 'mean')
    else:
        logger.info('WT_STRING NOT SUPPORTED.')

    ancestors = [f'{preproc_path}', f'{work_dir}/ERA5.nc']
    provenance_record = get_provenance_record(f'{var_name} means for \
                                              {wt_string}',
                                              ancestors,
                                              [var_name],
                                              ['map'], ['mean'])

    log_provenance(f'{dataset}_{var_name}_{wt_string}_means_prov',
                   cfg, provenance_record)


def calc_wt_anomalies(cfg: dict, cube: iris.cube.Cube,
                      wt_cubes: iris.cube.CubeList, data_info: dict):
    """Calculate anomalies of fields of each weathertype.

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube with variable data
        wt_cubes (iris.cube.CubeList): List of cubes of lwt, slwt_ERA5
                                        and slwt_EOBS
        dataset (str): Name of dataset
        var_name (str): Name of variable
        wt_string (str): lwt, slwt_ERA5 or slwt_EOBS
                        slwt are calculated based on ERA5 or EOBS
                        precipitation data, respectively
        preproc_path (str): Ancestor path
    """

    work_dir = cfg.get('work_dir')
    dataset = data_info.get('dataset')
    var_name = data_info.get('var_name')
    wt_string = data_info.get('wt_string')
    preproc_path = data_info.get('preproc_path')

    logger.info('Calculating %s %s anomalies for %s', dataset, var_name,
                wt_string)

    if wt_string == 'slwt_ERA5':
        slwt_era5_cube = wt_cubes[1]
        tcoord = slwt_era5_cube.coord('time')
        slwt_era5 = slwt_era5_cube.data[:]
    elif wt_string == 'slwt_EOBS':
        slwt_eobs_cube = wt_cubes[2]
        tcoord = slwt_eobs_cube.coord('time')
        slwt_eobs = slwt_eobs_cube.data[:]
    elif wt_string == 'lwt':
        lwt_cube = wt_cubes[0]
        tcoord = lwt_cube.coord('time')
        lwt = lwt_cube.data[:]

    num_slwt_era5 = len(np.unique(slwt_era5))
    num_slwt_eobs = len(np.unique(slwt_eobs))

    if num_slwt_eobs != num_slwt_era5:
        logger.info(
            'calc_wt_means - CAUTION: unequal number of \
                    slwt_era5 (%s) and slwt_eobs (%s)!', num_slwt_era5,
            num_slwt_eobs)

    if 'slwt' in wt_string:
        for wt in range(1, max(num_slwt_era5, num_slwt_eobs)):
            if wt_string == 'slwt_ERA5':
                target_indices = np.where(slwt_era5 == wt)
            elif wt_string == 'slwt_EOBS':
                target_indices = np.where(slwt_eobs == wt)
            else:
                logger.info('WT_STRING not supported!')
            if len(target_indices[0]) < 1:
                logger.info(
                    'calc_wt_anomalies - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset)
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)
            plot_maps(
                wt, cfg,
                cube.collapsed('time', iris.analysis.MEAN) - wt_cube_mean,
                data_info, 'anomaly')
    elif wt_string == 'lwt':
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    'calc_wt_anomalies - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset)
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)
            plot_maps(
                wt, cfg,
                cube.collapsed('time', iris.analysis.MEAN) - wt_cube_mean,
                data_info, 'anomaly')
    else:
        logger.info('WT_STRING NOT SUPPORTED.')

    ancestors = [f'{preproc_path}', f'{work_dir}/ERA5.nc']
    provenance_record = get_provenance_record(f'{var_name} anomaly for \
                                              {wt_string}',
                                              ancestors,
                                              [var_name],
                                              ['map'], ['anomaly'])

    log_provenance(f'{dataset}_{var_name}_{wt_string}_anomaly_prov',
                   cfg, provenance_record)


def calc_wt_std(cfg: dict, cube: iris.cube.Cube,
                wt_cubes: iris.cube.CubeList, data_info: dict):
    """Calculate standard deviation of fields of each weathertype.

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube with variable data
        wt_cubes (iris.cube.CubeList): List of cubes of lwt, slwt_ERA5
                                        and slwt_EOBS
        dataset (str): Name of dataset
        var_name (str): Name of variable
        wt_string (str): lwt, slwt_ERA5 or slwt_EOBS
                        slwt are calculated based on ERA5 or EOBS
                        precipitation data, respectively
        preproc_path (str): Ancestor path
    """

    work_dir = cfg.get('work_dir')
    dataset = data_info.get('dataset')
    var_name = data_info.get('var_name')
    wt_string = data_info.get('wt_string')
    preproc_path = data_info.get('preproc_path')

    logger.info('Calculating %s %s standard deviation for %s', dataset,
                var_name, wt_string)

    if wt_string == 'slwt_ERA5':
        slwt_era5_cube = wt_cubes[1]
        tcoord = slwt_era5_cube.coord('time')
        slwt_era5 = slwt_era5_cube.data[:]
    elif wt_string == 'slwt_EOBS':
        slwt_eobs_cube = wt_cubes[2]
        tcoord = slwt_eobs_cube.coord('time')
        slwt_eobs = slwt_eobs_cube.data[:]
    elif wt_string == 'lwt':
        lwt_cube = wt_cubes[0]
        tcoord = lwt_cube.coord('time')
        lwt = lwt_cube.data[:]

    num_slwt_era5 = len(np.unique(slwt_era5))
    num_slwt_eobs = len(np.unique(slwt_eobs))

    if num_slwt_eobs != num_slwt_era5:
        logger.info(
            'calc_wt_means - CAUTION: unequal number of \
                    slwt_era5 (%s) and slwt_eobs (%s)!', num_slwt_era5,
            num_slwt_eobs)

    if 'slwt' in wt_string:
        for wt in range(1, max(num_slwt_era5, num_slwt_eobs)):
            if wt_string == 'slwt_ERA5':
                target_indices = np.where(slwt_era5 == wt)
            elif wt_string == 'slwt_EOBS':
                target_indices = np.where(slwt_eobs == wt)
            else:
                logger.info('WT_STRING not supported!')
            if len(target_indices[0]) < 1:
                logger.info(
                    'calc_slwt_obs - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset)
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_std = extracted_cube.collapsed('time',
                                                   iris.analysis.STD_DEV)
            plot_maps(wt, cfg, wt_cube_std, data_info,
                      'standard deviation')
    elif wt_string == 'lwt':
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    'calc_wt_std - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset)
                continue
            dates = [
                tcoord.units.num2date(tcoord.points[i]) for i in target_indices
            ]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_std = extracted_cube.collapsed('time',
                                                   iris.analysis.STD_DEV)
            plot_maps(wt, cfg, wt_cube_std, data_info,
                      'standard deviation')
    else:
        logger.info('WT_STRING NOT SUPPORTED.')

    ancestors = [f'{preproc_path}', f'{work_dir}/ERA5.nc']
    provenance_record = get_provenance_record(f'{var_name} standard \
                                              deviation for \
                                              {wt_string}',
                                              ancestors,
                                              [var_name],
                                              ['map'], ['stddev'])

    log_provenance(f'{dataset}_{var_name}_{wt_string}_std_prov',
                   cfg, provenance_record)


def run_predefined_slwt(work_dir: str, dataset_name: str,
                        lwt: np.array, predefined_slwt: dict):
    predefined_slwt = check_mapping_dict_format(predefined_slwt)
    write_mapping_dict(work_dir, dataset_name, predefined_slwt)
    write_mapping_dict(work_dir, "E-OBS", predefined_slwt)
    slwt_era5 = map_lwt_to_slwt(lwt, predefined_slwt)
    slwt_eobs = map_lwt_to_slwt(lwt, predefined_slwt)

    return slwt_era5, slwt_eobs


def combine_wt_to_file(cfg: dict, wt_list: list,
                       cube: iris.cube.Cube,
                       file_name: str):
    """Combine lwt and slwt arrays to one file.

    Args:
        cfg (dict): Configuration dictionary from recipe
        lwt (np.array): lwt array
        slwt_era5 (np.array): slwt_era5 array
        slwt_eobs (np.array): slwt_eobs array
        cube (iris.cube.Cube): Cube of data to keep time coordinate
        file_name (str): Name of output file
    """

    lwt = wt_list[0]
    slwt_era5 = wt_list[1]
    slwt_eobs = wt_list[2]

    logger.info('Writing weathertypes to %s', file_name)

    tcoord = cube.coord('time')
    time_points = tcoord.units.num2date(tcoord.points)

    write_path = cfg.get('work_dir')

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(lwt, long_name='lwt'))
    wt_cube.append(iris.cube.Cube(slwt_era5, long_name='slwt_era5'))
    wt_cube.append(iris.cube.Cube(slwt_eobs, long_name='slwt_eobs'))

    wt_cube[0].add_dim_coord(tcoord, 0)
    wt_cube[1].add_dim_coord(tcoord, 0)
    wt_cube[2].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f'{write_path}/{file_name}.nc')

    # write to csv file
    d = {
        'date': time_points[:],
        'lwt': np.int8(lwt),
        'slwt_era5': np.int8(slwt_era5),
        'slwt_eobs': np.int8(slwt_eobs),
    }
    df = pd.DataFrame(data=d)
    df.to_csv(write_path + f'/{file_name}.csv', index=False)


def write_lwt_to_file(cfg: dict, lwt: np.array, cube: iris.cube.Cube,
                      file_name: str):
    """Write only lwt to file.

    Args:
        cfg (dict): Configuration dictionary from recipe
        lwt (np.array): lwt array
        cube (iris.cube.Cube): Cube to keep time coordinate
        file_name (str): Output filename
    """

    logger.info('Writing Lamb Weathertype to %s', file_name)

    tcoord = cube.coord('time')
    time_points = tcoord.units.num2date(tcoord.points)

    write_path = cfg.get('work_dir')

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(np.int8(lwt), long_name='lwt'))

    wt_cube[0].add_dim_coord(tcoord, 0)
    iris.save(wt_cube, f'{write_path}/{file_name}.nc')

    # write to csv file
    d = {'date': time_points[:], 'lwt': np.int8(lwt)}
    df = pd.DataFrame(data=d)
    df.to_csv(write_path + f'/{file_name}.csv', index=False)


def calc_lwt_model(cfg: dict, cube: iris.cube.Cube, dataset: str,
                   output_file_path: str):
    """Calculate lwt for model data.

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube to keep time coordinate
        dataset (str): Name of dataset
        output_file_path (str): Path to write output
    """

    work_dir = cfg.get('work_dir')

    if not os.path.exists(f'{work_dir}/{output_file_path}'):
        os.makedirs(f'{work_dir}/{output_file_path}')

    wt = wt_algorithm(cube, dataset)

    tcoord = cube.coord('time')
    time_points = tcoord.units.num2date(tcoord.points)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(wt, long_name='lwt'))

    logger.info(
        'Writing Lamb Weathertype for %s \
                to file %s.nc', dataset, dataset)

    wt_cube[0].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f'{work_dir}/{output_file_path}/{dataset}.nc')

    # write to csv file
    d = {'date': time_points[:], 'lwt': np.int8(wt)}
    df = pd.DataFrame(data=d)
    df.to_csv(f'{work_dir}/{output_file_path}/{dataset}.csv', index=False)


def map_lwt_to_slwt(lwt: np.array, mapping_dict: dict) -> np.array:
    """Map lwt array to slwt array.

    Args:
        lwt (np.array): array of lwt
        mapping_dict (dict): mapping dictionary in {lwt: slwt, ...} format

    Returns:
        np.array: array of slwt
    """

    return np.array([np.int8(mapping_dict.get(value, 0)) for value in lwt])


def check_mapping_dict_format(mapping_dict: dict) -> dict:
    """Check format of mapping dict and return in {lwt: slwt, ...} format.

    Args:
        mapping_dict (dict): mapping dict in any format

    Returns:
        dict: mapping dict in {lwt: slwt, ...} format
    """

    if isinstance(mapping_dict.get(list(mapping_dict.keys())[0]), list):
        dict_ = reverse_convert_dict(mapping_dict)
    else:
        dict_ = mapping_dict

    return dict_
