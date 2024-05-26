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


logger = logging.getLogger(os.path.basename(__file__))
warnings.filterwarnings("ignore", ".*Collapsing a non-contiguous coordinate*")

def get_provenance_record(caption: str, ancestors: list,
                          long_names: list, plot_types: list) -> dict:
    """Create a provenance record describing the diagnostic data and plots."""
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
    """Log provenance info."""
    provenance_record = get_provenance_record(
        caption, ancestors, long_names, plot_types)
    with ProvenanceLogger(cfg) as provenance_logger:
        provenance_logger.log(filename, provenance_record)

    logger.info('Output stored as %s', filename)


def turn_set_to_mapping_dict(list_: list) -> dict:
    result_dict = {}

    for i, s in enumerate(list_):
        for elem in s:
            if elem not in result_dict:
                result_dict[elem] = i + 1
            else:
                result_dict[elem].append(i)

    return result_dict


def get_mapping_dict(selected_pairs: list) -> dict:
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
    mapping_dict = turn_set_to_mapping_dict(merged_tuples)

    return mapping_dict


def calculate_slwt_obs(cfg: dict, lwt: np.array, cube: iris.cube.Cube,
                       dataset: str, correlation_thresold: float,
                       rmse_threshold: float, ancestors: list) -> np.array:
    
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

    with open(
        f'{work_dir}/wt_mapping_dict_{dataset}.json', 'w', encoding='utf-8'
    ) as file:
        json.dump(mapping_dict, file)

    log_provenance('Lamb Weathertypes', f'{dataset}_wt_prov', cfg,
                   ancestors,
                   ['Lamb Weathertypes'], [''])

    return np.array([np.int8(mapping_dict.get(value, 0)) for value in lwt])


def wt_algorithm(cube: iris.cube.Cube, dataset: str) -> np.array:

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

    const1 = 1 / np.cos(45 * np.pi / 180)
    const2 = np.sin(45 * np.pi / 180) / np.sin(40 * np.pi / 180)
    const3 = np.sin(45 * np.pi / 180) / np.sin(50 * np.pi / 180)
    const4 = 1 / (2 * np.cos(45 * np.pi / 180) ** 2)

    # westerly flow
    w = 1 / 2 * (cube.data[:, 1, 2] + cube.data[:, 1, 4]) - 1 / 2 * (
        cube.data[:, 3, 2] + cube.data[:, 3, 4]
    )
    # southerly flow
    s = const1 * (
        1 / 4 * (cube.data[:, 3, 4] + 2 *
                 cube.data[:, 2, 4] + cube.data[:, 1, 4])
        - 1 / 4 * (cube.data[:, 3, 2] + 2 *
                   cube.data[:, 2, 2] + cube.data[:, 1, 2])
    )
    # resultant flow
    f = (s**2 + w**2) ** (1 / 2)
    # westerly shear vorticity
    zw = const2 * (
        1 / 2 * (cube.data[:, 0, 2] + cube.data[:, 0, 4])
        - 1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])
    ) - const3 * (
        1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])
        - 1 / 2 * (cube.data[:, 4, 2] + cube.data[:, 4, 4])
    )
    # southerly shear vorticity
    zs = const4 * (
        1 / 4 * (cube.data[:, 3, 6] + 2 *
                 cube.data[:, 2, 6] + cube.data[:, 1, 6])
        - 1 / 4 * (cube.data[:, 3, 4] + 2 *
                   cube.data[:, 2, 4] + cube.data[:, 1, 4])
        - 1 / 4 * (cube.data[:, 3, 2] + 2 *
                   cube.data[:, 2, 2] + cube.data[:, 1, 2])
        + 1 / 4 * (cube.data[:, 3, 0] + 2 *
                   cube.data[:, 2, 0] + cube.data[:, 1, 0])
    )
    # total shear vorticity
    z = zw + zs

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

    work_dir = cfg.get('work_dir')

    if not os.path.exists(f'{work_dir}/{output_file_path}'):
        os.makedirs(f'{work_dir}/{output_file_path}')

    wt = wt_algorithm(cube, dataset)

    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
    tcoord = cube.coord('time')
    time_points = tcoord.units.num2date(tcoord.points)

    with open(f'{work_dir}/wt_mapping_dict_ERA5.json',
              'r', encoding='utf-8') as file:
        mapping_dict_era5 = json.load(file)

    with open(f'{work_dir}/wt_mapping_dict_E-OBS.json',
              'r', encoding='utf-8') as file:
        mapping_dict_eobs = json.load(file)

    logger.info('Calculating simplified Lamb Weathertypes for %s', dataset)

    slwt_era5 = np.array(
        [np.int8(mapping_dict_era5.get(f'{value}', 0)) for value in np.int8(wt)]
    )
    slwt_eobs = np.array(
        [np.int8(mapping_dict_eobs.get(f'{value}', 0)) for value in np.int8(wt)]
    )

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


def rmse(subarray1: np.array, subarray2: np.array):
    '''
    Calculate root mean square difference between two arrays.
    '''
    return np.sqrt(np.mean((subarray1 - subarray2) ** 2))


def write_corr_rmse_to_csv(cfg: dict, pattern_correlation_matrix: np.array,
                           rmse_matrix: np.array, dataset: str):

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

    return selected_pairs


def calculate_wt_means(cfg: dict, cube: iris.cube.Cube,
                       wt_cubes: iris.cube.CubeList, dataset: str,
                       var_name: str, wt_string: str, preproc_path: str):

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

    if 'slwt' in wt_string:
        for wt in range(1, 10):
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

    if 'slwt' in wt_string:
        for wt in range(1, 10):
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

    if 'slwt' in wt_string:
        for wt in range(1, 10):
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

    logger.info('Writing weathertypes to %s', file_name)

    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
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

    logger.info('Writing Lamb Weathertype to %s', file_name)

    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
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

    work_dir = cfg.get('work_dir')

    if not os.path.exists(f'{work_dir}/{output_file_path}'):
        os.makedirs(f'{work_dir}/{output_file_path}')

    wt = wt_algorithm(cube, dataset)

    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
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


def run_my_diagnostic(cfg: dict):
    '''
    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    '''
    # assemble the data dictionary keyed by dataset name
    # this makes use of the handy group_metadata function that
    # orders the data by 'dataset'; the resulting dictionary is
    # keyed on datasets e.g. dict = {'MPI-ESM-LR': [var1, var2...]}
    # where var1, var2 are dicts holding all needed information per variable
    preproc_variables_dict = group_metadata(
        cfg.get('input_data').values(), 'dataset')

    correlation_threshold = cfg.get('correlation_threshold')
    rmse_threshold = cfg.get('rmse_threshold')
    work_dir = cfg.get('work_dir')

    # load cubes and run functions
    # key = dataset name, value is dataset
    if cfg.get('automatic_slwt'):
        for key, value in preproc_variables_dict.items():
            if key == 'ERA5':
                wt_preproc = iris.load_cube(value[0].get('filename'))
                wt_preproc_prcp = iris.load_cube(value[1].get('filename'))
                mean_preproc_psl = iris.load_cube(value[2].get('filename'))
                mean_preproc_prcp = iris.load_cube(value[3].get('filename'))
                mean_preproc_tas = iris.load_cube(value[4].get('filename'))
                wt_preproc_prcp_eobs = iris.load_cube(
                    preproc_variables_dict.get('E-OBS')[0].get('filename')
                )

                # calculate lwt
                lwt = wt_algorithm(wt_preproc, key)

                era5_ancestors = [value[0].get('filename'),
                                  value[1].get('filename')]
                eobs_ancestors = [value[0].get('filename'),
                preproc_variables_dict.get('E-OBS')[0].get('filename')]

                # calculate simplified lwt based on precipitation patterns
                slwt_era5 = calculate_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp,
                    key,
                    correlation_threshold,
                    rmse_threshold,
                    era5_ancestors,
                )
                slwt_eobs = calculate_slwt_obs(
                    cfg,
                    lwt,
                    wt_preproc_prcp_eobs,
                    'E-OBS',
                    correlation_threshold,
                    rmse_threshold,
                    eobs_ancestors,
                )


                # write to file
                combine_wt_to_file(cfg, lwt, slwt_era5,
                                   slwt_eobs, wt_preproc, key)

                # load weathertype files as cubes
                lwt_cube = iris.load_cube(f'{work_dir}/{key}.nc', 'lwt')
                slwt_era5_cube = iris.load_cube(
                    f'{work_dir}/{key}.nc', 'slwt_era5')
                slwt_eobs_cube = iris.load_cube(
                    f'{work_dir}/{key}.nc', 'slwt_eobs')
                wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]

                preproc_path_psl = value[2].get('filename')
                preproc_path_prcp = value[3].get('filename')
                preproc_path_tas = value[4].get('filename')

                # plot means
                calculate_wt_means(
                    cfg,
                    mean_preproc_psl,
                    wt_cubes,
                    key,
                    var_name='psl',
                    wt_string='lwt',
                    preproc_path=preproc_path_psl,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_prcp,
                    wt_cubes,
                    key,
                    var_name='prcp',
                    wt_string='lwt',
                    preproc_path=preproc_path_prcp,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_tas,
                    wt_cubes,
                    key,
                    var_name='tas',
                    wt_string='lwt',
                    preproc_path=preproc_path_tas,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_psl,
                    wt_cubes,
                    key,
                    var_name='psl',
                    wt_string='slwt_ERA5',
                    preproc_path=preproc_path_psl,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_prcp,
                    wt_cubes,
                    key,
                    var_name='prcp',
                    wt_string='slwt_ERA5',
                    preproc_path=preproc_path_prcp,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_tas,
                    wt_cubes,
                    key,
                    var_name='tas',
                    wt_string='slwt_ERA5',
                    preproc_path=preproc_path_tas,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_psl,
                    wt_cubes,
                    key,
                    var_name='psl',
                    wt_string='slwt_EOBS',
                    preproc_path=preproc_path_psl,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_prcp,
                    wt_cubes,
                    key,
                    var_name='prcp',
                    wt_string='slwt_EOBS',
                    preproc_path=preproc_path_prcp,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_tas,
                    wt_cubes,
                    key,
                    var_name='tas',
                    wt_string='slwt_EOBS',
                    preproc_path=preproc_path_tas,
                )
            else:
                if key == 'E-OBS':
                    continue
                wt_preproc = iris.load_cube(value[0].get('filename'))
                mean_preproc_psl = iris.load_cube(value[1].get('filename'))
                mean_preproc_prcp = iris.load_cube(value[2].get('filename'))
                mean_preproc_tas = iris.load_cube(value[3].get('filename'))

                model_name = key
                timerange = value[0].get('timerange').replace('/', '-')
                experiment = value[0].get('exp')
                ensemble = value[0].get('ensemble')

                output_file_path = f'{model_name}/{experiment}/{ensemble}/{timerange}'
                preproc_path = value[0].get('filename')

                # calculate weathertypes
                calculate_lwt_slwt_model(cfg, wt_preproc, key, preproc_path, output_file_path)

                # load wt files
                lwt_cube = iris.load_cube(
                    f'{work_dir}/{output_file_path}/{key}.nc', 'lwt')
                slwt_era5_cube = iris.load_cube(
                    f'{work_dir}/{output_file_path}/{key}.nc', 'slwt_era5')
                slwt_eobs_cube = iris.load_cube(
                    f'{work_dir}/{output_file_path}/{key}.nc', 'slwt_eobs')
                wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]


                preproc_path_psl = value[1].get('filename')
                preproc_path_prcp = value[2].get('filename')
                preproc_path_tas = value[3].get('filename')

                # plot means
                calculate_wt_means(
                    cfg,
                    mean_preproc_psl,
                    wt_cubes,
                    key,
                    var_name='psl',
                    wt_string='lwt',
                    preproc_path=preproc_path_psl,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_prcp,
                    wt_cubes,
                    key,
                    var_name='prcp',
                    wt_string='lwt',
                    preproc_path=preproc_path_prcp,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_tas,
                    wt_cubes,
                    key,
                    var_name='tas',
                    wt_string='lwt',
                    preproc_path=preproc_path_tas,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_psl,
                    wt_cubes,
                    key,
                    var_name='psl',
                    wt_string='slwt_ERA5',
                    preproc_path=preproc_path_psl,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_prcp,
                    wt_cubes,
                    key,
                    var_name='prcp',
                    wt_string='slwt_ERA5',
                    preproc_path=preproc_path_prcp,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_tas,
                    wt_cubes,
                    key,
                    var_name='tas',
                    wt_string='slwt_ERA5',
                    preproc_path=preproc_path_tas,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_psl,
                    wt_cubes,
                    key,
                    var_name='psl',
                    wt_string='slwt_EOBS',
                    preproc_path=preproc_path_psl,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_prcp,
                    wt_cubes,
                    key,
                    var_name='prcp',
                    wt_string='slwt_EOBS',
                    preproc_path=preproc_path_prcp,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_tas,
                    wt_cubes,
                    key,
                    var_name='tas',
                    wt_string='slwt_EOBS',
                    preproc_path=preproc_path_tas,
                )
    else:  # if automatic_slwt is false, just do the lwt
        for key, value in preproc_variables_dict.items():
            if key == 'ERA5':
                wt_preproc = iris.load_cube(value[0].get('filename'))
                wt_preproc_prcp = iris.load_cube(value[1].get('filename'))
                mean_preproc_psl = iris.load_cube(value[2].get('filename'))
                mean_preproc_prcp = iris.load_cube(value[3].get('filename'))
                mean_preproc_tas = iris.load_cube(value[4].get('filename'))
                wt_preproc_prcp_eobs = iris.load_cube(
                    preproc_variables_dict.get('E-OBS')[0].get('filename')
                )

                # calculate lwt
                lwt = wt_algorithm(wt_preproc, key)

                # write to file
                write_lwt_to_file(cfg, lwt, wt_preproc, key)

                # load weathertype files as cubes
                lwt_cube = iris.load_cube(f'{work_dir}/{key}.nc', 'lwt')
                wt_cubes = [lwt_cube]

                preproc_path_psl = value[2].get('filename')
                preproc_path_prcp = value[3].get('filename')
                preproc_path_tas = value[4].get('filename')

                # plot means
                calculate_wt_means(
                    cfg,
                    mean_preproc_psl,
                    wt_cubes,
                    key,
                    var_name='psl',
                    wt_string='lwt',
                    preproc_path=preproc_path_psl,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_prcp,
                    wt_cubes,
                    key,
                    var_name='prcp',
                    wt_string='lwt',
                    preproc_path=preproc_path_prcp,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_tas,
                    wt_cubes,
                    key,
                    var_name='tas',
                    wt_string='lwt',
                    preproc_path=preproc_path_tas,
                )
            else:
                if key == 'E-OBS':
                    continue
                wt_preproc = iris.load_cube(value[0].get('filename'))
                mean_preproc_psl = iris.load_cube(value[1].get('filename'))
                mean_preproc_prcp = iris.load_cube(value[2].get('filename'))
                mean_preproc_tas = iris.load_cube(value[3].get('filename'))

                model_name = key
                timerange = value[0].get('timerange').replace('/', '-')
                experiment = value[0].get('exp')
                ensemble = value[0].get('ensemble')

                output_file_path = f'{model_name}/{experiment}/{ensemble}/{timerange}'

                # calculate weathertypes
                calculate_lwt_model(cfg, wt_preproc, key, output_file_path)

                # load wt files
                lwt_cube = iris.load_cube(f'{work_dir}/{key}.nc', 'lwt')
                wt_cubes = [lwt_cube]

                # plot means
                preproc_path_psl = value[2].get('filename')
                preproc_path_prcp = value[3].get('filename')
                preproc_path_tas = value[4].get('filename')

                # plot means
                calculate_wt_means(
                    cfg,
                    mean_preproc_psl,
                    wt_cubes,
                    key,
                    var_name='psl',
                    wt_string='lwt',
                    preproc_path=preproc_path_psl,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_prcp,
                    wt_cubes,
                    key,
                    var_name='prcp',
                    wt_string='lwt',
                    preproc_path=preproc_path_prcp,
                )
                calculate_wt_means(
                    cfg,
                    mean_preproc_tas,
                    wt_cubes,
                    key,
                    var_name='tas',
                    wt_string='lwt',
                    preproc_path=preproc_path_tas,
                )


if __name__ == '__main__':

    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
