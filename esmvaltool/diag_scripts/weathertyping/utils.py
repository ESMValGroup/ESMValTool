# weathertyping diagnostic for esmvaltool

# operating system manipulations (e.g. path constructions)
import os

# to manipulate iris cubes
import iris

# plotting imports
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import iris.analysis.cartography
import iris.plot as iplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from matplotlib.colors import ListedColormap
import seaborn as sns

# general imports
import numpy as np
import pandas as pd
import json
import logging
import warnings
import os

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import (
    group_metadata,
    run_diagnostic,
    ProvenanceLogger
)

iris.FUTURE.datum_support = True
iris.FUTURE.save_split_attrs = True

logger = logging.getLogger(os.path.basename(__file__))
warnings.filterwarnings('ignore', '.*Collapsing a non-contiguous coordinate*')


def get_provenance_record(caption: str, ancestors: list,
                          long_names: list, plot_types: list) -> dict:
    '''Create a provenance record describing the diagnostic data and plots.'''
    record = {
        'caption':
        caption,
        'domains': ['reg'],
        'authors': [
            'jury_martin',
            'kroissenbrunner_thomas'
        ],
        'references': [
            'TBD'
        ],
        'projects': [
            'preval'
        ],
        'long_names':
        long_names,
        'plot_types':
        plot_types,
        'ancestors':
        ancestors,
    }
    return record


def log_provenance(caption: str, filename: str, cfg: dict,
                   ancestors: list, long_names: list, plot_types: list):
    '''Log provenance info.'''
    provenance_record = get_provenance_record(
        caption, ancestors, long_names, plot_types)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    logger.info('Output stored as %s', filename)


def turn_list_to_mapping_dict(list_: list) -> dict:
    '''Turns a list of combined WT to a dictionary for further
        processing.

    Args:
        list_ (list): List where entries are lists with related WT

    Returns:
        dict: Mapping dicitonary keys are simplified WT, values are Lamb WT
    '''
    result_dict = {}

    for i, s in enumerate(list_):
        for elem in s:
            if elem not in result_dict:
                result_dict[elem] = i + 1
            else:
                result_dict[elem].append(i)

    return result_dict


def get_mapping_dict(selected_pairs: list) -> dict:
    '''Get mapping dictionary from list of selected pairs.

    Args:
        selected_pairs (list): Selected pairs of WTs based on 
        precipitation patterns over specified area and
        correlation and RSME thresholds defined in recipe.S

    Returns:
        dict: Mapping dicitonary keys are simplified WT, values are Lamb WT 
    '''
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

    with open(
        f'{work_dir}/wt_mapping_dict_{dataset}.json', 'w', encoding='utf-8'
    ) as file:
        json.dump(mapping_dict_reformat, file)


def calculate_slwt_obs(cfg: dict, lwt: np.array, cube: iris.cube.Cube,
                       dataset: str, correlation_thresold: float,
                       rmse_threshold: float, ancestors: list) -> np.array:
    '''Calculate simplified weathertypes for observation datasets based on 
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
    '''

    logger.info('Calculating simplified Lamb Weathertypes for %s', dataset)

    work_dir = cfg.get('work_dir')
    tcoord = cube.coord('time')

    wt_data_prcp = []
    for wt in range(1, 28):
        target_indices = np.where(lwt == wt)
        if len(target_indices[0]) < 1:
            logger.info(
                'calculate_slwt_obs - CAUTION: Skipped wt %s \
                for dataset %s!', wt, dataset
            )
            continue
        dates = [tcoord.units.num2date(tcoord.points[i])
                 for i in target_indices]
        if dataset == 'E-OBS':
            extracted_cube = cube[target_indices]
        else:
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0])
            )
        wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)
        wt_data_prcp.append(wt_cube_mean.data.compressed())
    selected_pairs = process_prcp_mean(
        cfg, wt_data_prcp, correlation_thresold, rmse_threshold, dataset
    )

    with open(
        f'{work_dir}/wt_selected_pairs_{dataset}.json', 'w', encoding='utf-8'
    ) as file:
        json.dump(selected_pairs, file)

    mapping_dict = get_mapping_dict(selected_pairs)

    write_mapping_dict(work_dir, dataset, mapping_dict)

    log_provenance('Lamb Weathertypes', f'{dataset}_wt_prov', cfg,
                   ancestors,
                   ['Lamb Weathertypes'], [''])

    return map_lwt_to_slwt(lwt, mapping_dict)


def calc_const():
    '''Calculate constants for weathertyping algorithm.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993), 
    A comparison of Lamb circulation types with an objective classification scheme. 
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Returns:
        tuple: The four constants needed for WT calculation.
    '''

    const1 = 1 / np.cos(45 * np.pi / 180)
    const2 = np.sin(45 * np.pi / 180) / np.sin(40 * np.pi / 180)
    const3 = np.sin(45 * np.pi / 180) / np.sin(50 * np.pi / 180)
    const4 = 1 / (2 * np.cos(45 * np.pi / 180) ** 2)

    return const1, const2, const3, const4


def calc_westerly_flow(cube: iris.cube.Cube) -> np.array:
    '''Calculate the westerly flow over area.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993), 
    A comparison of Lamb circulation types with an objective classification scheme. 
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): Cube of psl data.

    Returns:
        np.array: westerly flow
    '''

    return 1 / 2 * (cube.data[:, 1, 2] + cube.data[:, 1, 4]) - 1 / 2 * (
        cube.data[:, 3, 2] + cube.data[:, 3, 4]
    )


def calc_southerly_flow(cube: iris.cube.Cube, const1: float) -> np.array:
    '''Calculate the southerly flow over area.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993), 
    A comparison of Lamb circulation types with an objective classification scheme. 
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): Cube of psl data.
        const1 (float): const1

    Returns:
        np.array: southerly flow
    '''

    return const1 * (
        1 / 4 * (cube.data[:, 3, 4] + 2 *
                 cube.data[:, 2, 4] + cube.data[:, 1, 4])
        - 1 / 4 * (cube.data[:, 3, 2] + 2 *
                   cube.data[:, 2, 2] + cube.data[:, 1, 2])
    )


def calc_resultant_flow(w: np.array, s: np.array) -> np.array:
    '''Calculate the resultant flow.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993), 
    A comparison of Lamb circulation types with an objective classification scheme. 
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        w (np.array): westerly flow.
        s (np.array): southerly flow

    Returns:
        np.array: resultant flow
    '''
    return (s**2 + w**2) ** (1 / 2)


def calc_westerly_shear_velocity(cube: iris.cube.Cube, const2: float, const3: float) -> np.array:
    '''Calculate westerly shear velocity.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993), 
    A comparison of Lamb circulation types with an objective classification scheme. 
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): cube of psl data
        const2 (float): const2
        const3 (float): const3

    Returns:
        np.array: westerly shear velocity
    '''
    return const2 * (
        1 / 2 * (cube.data[:, 0, 2] + cube.data[:, 0, 4])
        - 1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])
    ) - const3 * (
        1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])
        - 1 / 2 * (cube.data[:, 4, 2] + cube.data[:, 4, 4])
    )


def calc_southerly_shear_velocity(cube: iris.cube.Cube, const4: float) -> np.array:
    '''Calculate southerly shear velocity.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993), 
    A comparison of Lamb circulation types with an objective classification scheme. 
    Int. J. Climatol., 13: 655-663. https://doi.org/10.1002/joc.3370130606

    Args:
        cube (iris.cube.Cube): cube of psl data
        const4 (float): const4

    Returns:
        np.array: southerly shear velocity
    '''
    return const4 * (
        1 / 4 * (cube.data[:, 3, 6] + 2 *
                 cube.data[:, 2, 6] + cube.data[:, 1, 6])
        - 1 / 4 * (cube.data[:, 3, 4] + 2 *
                   cube.data[:, 2, 4] + cube.data[:, 1, 4])
        - 1 / 4 * (cube.data[:, 3, 2] + 2 *
                   cube.data[:, 2, 2] + cube.data[:, 1, 2])
        + 1 / 4 * (cube.data[:, 3, 0] + 2 *
                   cube.data[:, 2, 0] + cube.data[:, 1, 0])
    )


def calc_total_shear_velocity(zw: np.array, zs: np.array) -> np.array:
    '''Calculate total shear velocity.
    Eq. taken from: Jones, P.D., Hulme, M. and Briffa, K.R. (1993), 
    A comparison of Lamb circulation types with an objective classification scheme. 
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
    A comparison of Lamb circulation types with an objective classification scheme. 
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


def calculate_lwt_slwt_model(cfg: dict, cube: iris.cube.Cube, dataset: str,
                             preproc_path: str, output_file_path: str):
    '''Calculate Lamb WT and simplified WT for model data.

    Args:
        cfg (dict): Configuration dicitonary from recipe
        cube (iris.cube.Cube): PSL field of dataset
        dataset (str): Name of dataset
        preproc_path (str): Path of ancestors
        output_file_path (str): Path to write output file
    '''

    work_dir = cfg.get('work_dir')

    if not os.path.exists(f'{work_dir}/{output_file_path}'):
        os.makedirs(f'{work_dir}/{output_file_path}')

    wt = wt_algorithm(cube, dataset)

    tcoord = cube.coord('time')
    time_points = tcoord.units.num2date(tcoord.points)

    with open(f'{work_dir}/wt_mapping_dict_ERA5.json',
              'r', encoding='utf-8') as file:
        mapping_dict_era5_f = json.load(file)

    with open(f'{work_dir}/wt_mapping_dict_E-OBS.json',
              'r', encoding='utf-8') as file:
        mapping_dict_eobs_f = json.load(file)

    logger.info('Calculating simplified Lamb Weathertypes for %s', dataset)

    mapping_dict_era5 = reverse_convert_dict(mapping_dict_era5_f)
    mapping_dict_eobs = reverse_convert_dict(mapping_dict_eobs_f)

    slwt_era5 = map_lwt_to_slwt(wt, mapping_dict_era5)
    slwt_eobs = map_lwt_to_slwt(wt, mapping_dict_eobs)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(wt, long_name='lwt'))
    wt_cube.append(iris.cube.Cube(slwt_era5, long_name='slwt_era5'))
    wt_cube.append(iris.cube.Cube(slwt_eobs, long_name='slwt_eobs'))

    wt_cube[0].add_dim_coord(tcoord, 0)
    wt_cube[1].add_dim_coord(tcoord, 0)
    wt_cube[2].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f'{work_dir}/{output_file_path}/{dataset}.nc')

    # write to csv file
    d = {
        'date': time_points[:],
        'lwt': np.int8(wt),
        'slwt_ERA5': np.int8(slwt_era5),
        'slwt_EOBS': np.int8(slwt_eobs),
    }
    df = pd.DataFrame(data=d)
    df.to_csv(f'{work_dir}/{output_file_path}/{dataset}.csv', index=False)

    log_provenance('Lamb Weathertypes', f'{dataset}_wt_prov', cfg,
                   [preproc_path, f'{work_dir}/wt_mapping_dict_ERA5.json',
                    f'{work_dir}/wt_mapping_dict_E-OBS.json'],
                   ['Lamb Weathertypes'], [''])


def get_colormap(colormap_string: str) -> ListedColormap:
    '''Get colormaps based on string.

    Args:
        colormap_string (str): String to get Colormaps for either
                                psl, tas or precipitation.

    Returns:
        ListedColormap: Choosen Colormap
    '''

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


def plot_maps(wt: np.array, cfg: dict, cube: iris.cube.Cube, dataset: str,
              var_name: str, wt_string: str, mode: str):
    '''Plot maps.

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
    '''

    logger.info('Plotting %s %s %s for %s %s', dataset, var_name,
                mode, wt_string, wt)

    local_path = cfg.get('plot_dir')

    ax = plt.axes(projection=ccrs.PlateCarree())

    if var_name == 'psl':
        psl_cmap = get_colormap('psl')
        plt.title(f'{var_name} {mode}, wt: {wt}')
        unit = '[hPa]'
        im = iplt.contourf(cube / 100, cmap=psl_cmap)
    elif var_name == 'prcp':
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

    plt.savefig(
        os.path.join(
            local_path, f'{wt_string}_{wt}_{dataset}_{var_name}_{mode}.png')
    )
    plt.close()


def rmse(subarray1: np.array, subarray2: np.array) -> np.array:
    '''Calculate rsme.

    Args:
        subarray1 (np.array): array1
        subarray2 (np.array): array2

    Returns:
        np.array: rsme array
    '''
    return np.sqrt(np.mean((subarray1 - subarray2) ** 2))


def convert_dict(dict_: dict) -> dict:
    '''Convert mapping dictionary from {lwt: slwt, ...} format to {slwt: [lwt1, lwt2], ...}. 

    Args:
        dict_ (dict): Dict in the {lwt: slwt, ...} format

    Returns:
        dict: Dict in the {slwt: [lwt1, lwt2], ...} format
    '''

    new_dict = {}
    for key, value in dict_.items():
        if value not in new_dict:
            new_dict[value] = []
        new_dict[value].append(key)
    return new_dict


def reverse_convert_dict(dict_: dict) -> dict:
    '''Convert mapping dictionary from {slwt: [lwt1, lwt2], ...} format to {lwt: slwt, ...}. 

    Args:
        original_dict (dict): Dict in the {slwt: [lwt1, lwt2], ...}format

    Returns:
        dict: Dict in the  format {lwt: slwt, ...}
    '''
    new_dict = {}
    for key, value_list in dict_.items():
        for original_key in value_list:
            new_dict[original_key] = key
    return new_dict


def plot_corr_rmse_heatmaps(cfg: dict, pattern_correlation_matrix: np.array,
                            rmse_matrix: np.array, dataset: str):
    '''Plot heatmaps for correlation and rmse matrices

    Args:
        cfg (dict): cfg dict from recipe
        pattern_correlation_matrix (np.array): pattern correlation matrix
        rmse_matrix (np.array): rmse matrix
        dataset (str): string of dataset
    '''

    work_dir = cfg.get('work_dir')

    labels = np.arange(1, 28)

    mask = np.zeros_like(pattern_correlation_matrix)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style('white'):
        plt.title('Correlation Matrix')
        ax = sns.heatmap(pattern_correlation_matrix, mask=mask,
                         square=True, annot=True, annot_kws={'size': 5},
                         cmap='seismic', xticklabels=labels,
                         yticklabels=labels)
        ax.set_xlabel('lwt', fontsize=8, rotation=90)
        ax.set_ylabel('lwt', fontsize=8)
        plt.savefig(
            os.path.join(
                work_dir, f'correlation_matrix_{dataset}.png')
        )
        plt.close()

    mask = np.zeros_like(rmse_matrix)
    mask[np.triu_indices_from(mask)] = True
    with sns.axes_style('white'):
        plt.title('RMSE Matrix')
        ax = sns.heatmap(rmse_matrix, mask=mask,
                         square=True, annot=True, annot_kws={'size': 5},
                         cmap='seismic', xticklabels=labels,
                         yticklabels=labels)
        ax.set_xlabel('lwt', fontsize=8, rotation=90)
        ax.set_ylabel('lwt', fontsize=8)
        plt.savefig(
            os.path.join(
                work_dir, f'rmse_matrix_{dataset}.png')
        )
        plt.close()


def write_corr_rmse_to_csv(cfg: dict, pattern_correlation_matrix: np.array,
                           rmse_matrix: np.array, dataset: str):
    '''Write correlation and rsme matrix to csv files

    Args:
        cfg (dict): Configuration dictionary from recipe
        pattern_correlation_matrix (np.array): Correlation matrix
        rmse_matrix (np.array): RSME matrix
        dataset (str): Name of dataset
    '''

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
    '''Process precipitation fields for specified area to get
        a list of selected pairs of weathertypes with the 
        highest correlation (higher than correlation_threshold) 
        and smallest RSME (smaller than rsme_threshold) for further 
        processing and simplifying the WT.

    Args:
        cfg (dict): Configuration dictionary from recipe
        data (np.array): Precipitation data
        correlation_threshold (float): Correlation threshold
        rmse_threshold (float): RMSE threshold
        dataset (str): Name of dataset

    Returns:
        list: Selected pairs of WT. This is passed to get_mapping_dict
    '''

    logger.info('Calculating corr and rsme matrices for %s', dataset)

    selected_pairs = []
    pattern_correlation_matrix = np.ma.corrcoef(data)
    n = len(data)
    rmse_matrix = np.zeros((n, n))

    for i in range(n):
        for j in range(i + 1, n):
            rmse_matrix[i][j] = rmse(data[i], data[j])
            rmse_matrix[j][i] = rmse_matrix[i][j]
            if (
                pattern_correlation_matrix[i][j] >= correlation_threshold
                and rmse_matrix[i][j] <= rmse_threshold
            ):
                selected_pairs.append(
                    ((i + 1, j + 1),
                     pattern_correlation_matrix[i][j], rmse_matrix[i][j])
                )

    # write matrices to csv
    write_corr_rmse_to_csv(
        cfg, pattern_correlation_matrix, rmse_matrix, dataset)
    # plot heatmaps for matrices
    plot_corr_rmse_heatmaps(
        cfg, pattern_correlation_matrix, rmse_matrix, dataset)

    return selected_pairs


def calculate_wt_means(cfg: dict, cube: iris.cube.Cube,
                       wt_cubes: iris.cube.CubeList, dataset: str,
                       var_name: str, wt_string: str, preproc_path: str):
    '''Calculate means of fields of each weathertype

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
    '''

    logger.info('Calculating %s %s means for %s', dataset, var_name, wt_string)

    work_dir = cfg.get('work_dir')

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
        logger.info('calculate_wt_means - CAUTION: unequal number of \
                    slwt_era5 (%s) \ and slwt_eobs (%s)!',
                    num_slwt_era5, num_slwt_eobs)

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
                    'calculate_wt_means - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset
                )
                continue
            dates = [tcoord.units.num2date(tcoord.points[i])
                     for i in target_indices]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0])
            )
            wt_cube_mean = extracted_cube.collapsed('time',
                                                    iris.analysis.MEAN)
            plot_maps(wt, cfg, wt_cube_mean, dataset,
                      var_name, wt_string, 'mean')
    elif wt_string == 'lwt':
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    'calculate_wt_means - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset
                )
                continue
            dates = [tcoord.units.num2date(tcoord.points[i])
                     for i in target_indices]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0])
            )
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)
            plot_maps(wt, cfg, wt_cube_mean, dataset,
                      var_name, wt_string, 'mean')
    else:
        logger.info('WT_STRING NOT SUPPORTED.')

    log_provenance(f'{var_name} means for {wt_string}',
                   f'{dataset}_{var_name}_{wt_string}_means_prov', cfg,
                   [f'{preproc_path}', f'{work_dir}/ERA5.nc'],
                   [var_name], ['geo'])


def calculate_wt_anomalies(cfg: dict, cube: iris.cube.Cube,
                           wt_cubes: iris.cube.CubeList, dataset: str,
                           var_name: str, wt_string: str, preproc_path: str):
    '''Calculate anomalies of fields of each weathertype

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
    '''

    work_dir = cfg.get('work_dir')

    logger.info('Calculating %s %s anomalies for %s', dataset,
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
        logger.info('calculate_wt_means - CAUTION: unequal number of \
                    slwt_era5 (%s) \ and slwt_eobs (%s)!',
                    num_slwt_era5, num_slwt_eobs)

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
                    'calculate_wt_anomalies - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset
                )
                continue
            dates = [tcoord.units.num2date(tcoord.points[i])
                     for i in target_indices]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0])
            )
            wt_cube_mean = extracted_cube.collapsed('time',
                                                    iris.analysis.MEAN)
            plot_maps(wt, cfg, cube.collapsed('time', iris.analysis.MEAN) -
                      wt_cube_mean, dataset, var_name, wt_string, 'anomaly')
    elif wt_string == 'lwt':
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    'calculate_wt_anomalies - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset
                )
                continue
            dates = [tcoord.units.num2date(tcoord.points[i])
                     for i in target_indices]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0])
            )
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)
            plot_maps(wt, cfg, cube.collapsed('time', iris.analysis.MEAN) -
                      wt_cube_mean, dataset, var_name, wt_string, 'anomaly')
    else:
        logger.info('WT_STRING NOT SUPPORTED.')

    log_provenance(f'{var_name} anomaly for {wt_string}',
                   f'{dataset}_{var_name}_{wt_string}_anomalies_prov', cfg,
                   [f'{preproc_path}', f'{work_dir}/ERA5.nc'], [var_name], ['geo'])


def calculate_wt_std(cfg: dict, cube: iris.cube.Cube,
                     wt_cubes: iris.cube.CubeList, dataset: str,
                     var_name: str, wt_string: str, preproc_path: str):
    '''Calculate standard deviation of fields of each weathertype

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
    '''

    work_dir = cfg.get('work_dir')

    logger.info('Calculating %s %s standard deviation for %s',
                dataset, var_name, wt_string)

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
        logger.info('calculate_wt_means - CAUTION: unequal number of \
                    slwt_era5 (%s) \ and slwt_eobs (%s)!',
                    num_slwt_era5, num_slwt_eobs)

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
                    'calculate_slwt_obs - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset
                )
                continue
            dates = [tcoord.units.num2date(tcoord.points[i])
                     for i in target_indices]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0])
            )
            wt_cube_std = extracted_cube.collapsed('time',
                                                   iris.analysis.STD_DEV)
            plot_maps(wt, cfg, wt_cube_std, dataset, var_name,
                      wt_string, 'standard deviation')
    elif wt_string == 'lwt':
        for wt in range(1, 28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                logger.info(
                    'calculate_wt_std - CAUTION: Skipped wt %s \
                    for dataset %s!', wt, dataset
                )
                continue
            dates = [tcoord.units.num2date(tcoord.points[i])
                     for i in target_indices]
            extracted_cube = cube.extract(
                iris.Constraint(time=lambda t: t.point in dates[0])
            )
            wt_cube_std = extracted_cube.collapsed(
                'time', iris.analysis.STD_DEV)
            plot_maps(wt, cfg, wt_cube_std, dataset, var_name,
                      wt_string, 'standard deviation')
    else:
        logger.info('WT_STRING NOT SUPPORTED.')

    log_provenance(f'{var_name} standard deviation for {wt_string}',
                   f'{dataset}_{var_name}_{wt_string}_std_prov', cfg,
                   [f'{preproc_path}', f'{work_dir}/ERA5.nc'], [var_name], ['geo'])


def combine_wt_to_file(cfg: dict, lwt: np.array, slwt_era5: np.array,
                       slwt_eobs: np.array, cube: iris.cube.Cube,
                       file_name: str):
    '''Combine lwt and slwt arrays to one file

    Args:
        cfg (dict): Configuration dictionary from recipe
        lwt (np.array): lwt array
        slwt_era5 (np.array): slwt_era5 array
        slwt_eobs (np.array): slwt_eobs array
        cube (iris.cube.Cube): Cube of data to keep time coordinate
        file_name (str): Name of output file
    '''

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
    '''Write only lwt to file

    Args:
        cfg (dict): Configuration dictionary from recipe
        lwt (np.array): lwt array
        cube (iris.cube.Cube): Cube to keep time coordinate
        file_name (str): Output filename
    '''

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


def calculate_lwt_model(cfg: dict, cube: iris.cube.Cube, dataset: str, output_file_path: str):
    '''Calculate lwt for model data

    Args:
        cfg (dict): Configuration dictionary from recipe
        cube (iris.cube.Cube): Cube to keep time coordinate
        dataset (str): Name of dataset
        output_file_path (str): Path to write output
    '''

    work_dir = cfg.get('work_dir')

    if not os.path.exists(f'{work_dir}/{output_file_path}'):
        os.makedirs(f'{work_dir}/{output_file_path}')

    wt = wt_algorithm(cube, dataset)

    tcoord = cube.coord('time')
    time_points = tcoord.units.num2date(tcoord.points)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(wt, long_name='lwt'))

    logger.info('Writing Lamb Weathertype for %s \
                to file %s.nc', dataset, dataset)

    wt_cube[0].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f'{work_dir}/{output_file_path}/{dataset}.nc')

    # write to csv file
    d = {'date': time_points[:], 'lwt': np.int8(wt)}
    df = pd.DataFrame(data=d)
    df.to_csv(f'{work_dir}/{output_file_path}/{dataset}.csv', index=False)


def map_lwt_to_slwt(lwt: np.array, mapping_dict: dict) -> np.array:
    '''map lwt array to slwt array.

    Args:
        lwt (np.array): array of lwt
        mapping_dict (dict): mapping dictionary in {lwt: slwt, ...} format

    Returns:
        np.array: array of slwt
    '''

    return np.array([np.int8(mapping_dict.get(value, 0)) for value in lwt])


def check_mapping_dict_format(mapping_dict: dict) -> dict:
    '''Check format of mapping dict and return in {lwt: slwt, ...} format

    Args:
        mapping_dict (dict): mapping dict in any format

    Returns:
        dict: mapping dict in {lwt: slwt, ...} format
    '''

    if isinstance(mapping_dict.get(list(mapping_dict.keys())[0]), list):
        dict_ = reverse_convert_dict(mapping_dict)
    else:
        dict_ = mapping_dict

    return dict_
