#weathertyping diagnostic for esmvaltool


# operating system manipulations (e.g. path constructions)
import os
import numpy as np
import pandas as pd
import json
# to manipulate iris cubes
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import cf_units
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.cm as cm
from matplotlib.colors import ListedColormap
import cartopy.util
from iris.analysis.cartography import wrap_lons

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import group_metadata, run_diagnostic

def turn_set_to_mapping_dict(array):

    result_dict = {}

    for i, s in enumerate(array):
        for elem in s:
            if elem not in result_dict:
                result_dict[elem] = i+1
            else:
                result_dict[elem].append(i)

    return result_dict

def get_mapping_dict(selected_pairs):

    mapping_dict_douglas_paper = {
            1: 1, 2: 1, 19: 1,
            3: 2, 4: 2, 22: 2, 21: 2,
            5: 3, 6: 3, 15: 3, 16: 3,
            7: 4, 8: 4, 11: 4, 18: 4,
            9: 5, 17: 5,
            10: 6, 20: 6,
            12: 7, 13: 7, 14: 7,
            23: 8, 24: 8,
            25: 9, 26: 9, 27:27, 0:0
        }

    mapping_array = []

    for i in range(0,len(selected_pairs)):
        mapping_array.append(selected_pairs[i][0])

    s=[set(i) for i in mapping_array if i]

    def find_intersection(m_list):
        for i,v in enumerate(m_list) : 
            for j,k in enumerate(m_list[i+1:],i+1):  
                if v&k:
                    s[i]=v.union(m_list.pop(j))
                    return find_intersection(m_list)
        return m_list

    merged_tuples = find_intersection(s)

    mapping_dict = turn_set_to_mapping_dict(merged_tuples)

    return mapping_dict

def calculate_slwt_obs(cfg, LWT, cube, dataset: str, correlation_thresold, rmse_threshold):
        
        print(dataset)
        
        work_dir = cfg['work_dir']

        tcoord = cube.coord("time")

        wt_data_prcp = []
        for wt in range(1,28):
            target_indices = np.where(LWT == wt)
            if len(target_indices[0]) < 1:
                print(f"calculate_slwt_obs - CAUTION: Skipped wt {wt} for dataset {dataset}!")
                continue
            dates = [tcoord.units.num2date(tcoord.points[i]) for i in target_indices]
            if dataset == "E-OBS":
                extracted_cube = cube[target_indices]#.extract(iris.Constraint(time=lambda t: t.point in dates[0]))
            else:
                extracted_cube = cube.extract(iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)#cube[target_indices[0], :, :].collapsed('time', iris.analysis.MEAN)
            #plot_mean(wt, cfg, wt_cube_mean, dataset, f"{var_name}_", wt_string)
            wt_data_prcp.append(wt_cube_mean.data.compressed())
        selected_pairs = process_prcp_mean(wt, cfg, wt_data_prcp, correlation_thresold, rmse_threshold, dataset)

        with open(f'{work_dir}/wt_selected_pairs_{dataset}.json', 'w') as file:
            json.dump(selected_pairs, file)

        mapping_dict = get_mapping_dict(selected_pairs)

        with open(f'{work_dir}/wt_mapping_dict_{dataset}.json', 'w') as file:
            json.dump(mapping_dict, file)

        return np.array([np.int8(mapping_dict.get(value, 0)) for value in np.int8(LWT)])

def wt_algorithm(cube):

    #lats and lons corresponding to datapoints
    #55, 5 -> 1
    #55, 15 -> 2
    #50, -5 -> 3
    #50, 5 -> 4
    #50, 15 -> 5
    #50, 25 -> 6
    #45, -5 -> 7
    #45, 5 -> 8
    #45, 15 -> 9
    #45, 25 -> 10
    #40, -5 -> 11
    #40, 5 -> 12
    #40, 15 -> 13
    #40, 25 -> 14
    #35, 5 -> 15
    #35, 15 -> 16

    #lons: -5, 0, 5, 10, 15, 20, 25
    #lats: 35, 40, 45, 50, 55
    
    const1 = 1 / np.cos(45 * np.pi / 180)
    const2 = np.sin(45 * np.pi / 180) / np.sin(40 * np.pi / 180)
    const3 = np.sin(45 * np.pi / 180) / np.sin(50 * np.pi / 180)
    const4 = 1 / (2 * np.cos(45 * np.pi / 180) ** 2)

    # westerly flow
    W = 1 / 2 * (cube.data[:, 1, 2] + cube.data[:, 1, 4]) - 1 / 2 * (cube.data[:, 3, 2] + cube.data[:, 3, 4])
    #southerly flow
    S = const1 * (1 / 4 * (cube.data[:, 3, 4] + 2 * cube.data[:, 2, 4] + cube.data[:, 1, 4]) - 1 / 4 * (cube.data[:, 3, 2] + 2 * cube.data[:, 2, 2] + cube.data[:, 1, 2]))
    #resultant flow
    F = (S ** 2 + W ** 2) ** (1 / 2)
    #westerly shear vorticity
    ZW = const2 * (1 / 2 * (cube.data[:, 0, 2] + cube.data[:, 0, 4]) - 1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4])) - const3 * (
            1 / 2 * (cube.data[:, 2, 2] + cube.data[:, 2, 4]) - 1 / 2 * (cube.data[:, 4, 2] + cube.data[:, 4, 4]))
    #southerly shear vorticity
    ZS = const4 * (1 / 4 * (cube.data[:, 3, 6] + 2 * cube.data[:, 2, 6] + cube.data[:, 1, 6]) - 1 / 4 * (cube.data[:, 3, 4] + 2 * cube.data[:, 2, 4] + cube.data[:, 1, 4]) - 1 / 4 * (
            cube.data[:, 3, 2] + 2 * cube.data[:, 2, 2] + cube.data[:, 1, 2]) + 1 / 4 * (cube.data[:, 3, 0] + 2 * cube.data[:, 2, 0] + cube.data[:, 1, 0]))
    #total shear vorticity
    Z = ZW + ZS

    WT = np.zeros(len(Z))

    for i in range(len(Z)):

        direction = np.arctan(W[i] / S[i]) * 180 / np.pi  # deg
        if S[i] >= 0:
            direction += 180  # deg

        if direction < 0:
            direction += 360  # deg

    # Lamb pure directional type
        if abs(Z[i]) < F[i]:
            if 337.5 <= direction or direction < 22.5:
                WT[i] = 1
            elif 22.5 <= direction < 67.5:
                WT[i] = 2
            elif 67.5 <= direction < 112.5:
                WT[i] = 3
            elif 112.5 <= direction < 157.5:
                WT[i] = 4
            elif 157.5 <= direction < 202.5:
                WT[i] = 5
            elif 202.5 <= direction < 247.5:
                WT[i] = 6
            elif 247.5 <= direction < 292.5:
                WT[i] = 7
            elif 292.5 <= direction < 337.5:
                WT[i] = 8
    # Lamb’s pure cyclonic and anticyclonic type
        elif (2 * F[i]) < abs(Z[i]):
            if Z[i] > 0:
                WT[i] = 9

            elif Z[i] < 0:
                WT[i] = 10
    #Lambs’s synoptic/direction hybrid types
        elif F[i] < abs(Z[i]) < (2 * F[i]):
            if Z[i] > 0:
                if 337.5 <= direction or direction < 22.5:
                    WT[i] = 11
                elif 22.5 <= direction < 67.5:
                    WT[i] = 12
                elif 67.5 <= direction < 112.5:
                    WT[i] = 13
                elif 112.5 <= direction < 157.5:
                    WT[i] = 14
                elif 157.5 <= direction < 202.5:
                    WT[i] = 15
                elif 202.5 <= direction < 247.5:
                    WT[i] = 16
                elif 247.5 <= direction < 292.5:
                    WT[i] = 17
                elif 292.5 <= direction < 337.5:
                    WT[i] = 18

            elif Z[i] < 0:
                if 337.5 <= direction or direction < 22.5:
                    WT[i] = 19
                elif 22.5 <= direction < 67.5:
                    WT[i] = 20
                elif 67.5 <= direction < 112.5:
                    WT[i] = 21
                elif 112.5 <= direction < 157.5:
                    WT[i] = 22
                elif 157.5 <= direction < 202.5:
                    WT[i] = 23
                elif 202.5 <= direction < 247.5:
                    WT[i] = 24
                elif 247.5 <= direction < 292.5:
                    WT[i] = 25
                elif 292.5 <= direction < 337.5:
                    WT[i] = 26
    # light indeterminate flow, corresponding to Lamb’s unclassified type U
        elif abs(Z[i]) < 6 and F[i] < 6:
            WT[i] = 27

    return WT

def calculate_lwt_slwt_model(cfg, cube, dataset):

    work_dir = cfg['work_dir']

    WT = wt_algorithm(cube)

    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    with open(f'{work_dir}/wt_mapping_dict_ERA5.json', 'r') as file:
        mapping_dict_era5 = json.load(file)
    
    with open(f'{work_dir}/wt_mapping_dict_E-OBS.json', 'r') as file:
        mapping_dict_eobs = json.load(file)

    slwt_era5 = np.array([np.int8(mapping_dict_era5.get(f"{value}", 0)) for value in np.int8(WT)])
    slwt_eobs = np.array([np.int8(mapping_dict_eobs.get(f"{value}", 0)) for value in np.int8(WT)])

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(WT, long_name="lwt"))
    wt_cube.append(iris.cube.Cube(slwt_era5, long_name="slwt_era5"))
    wt_cube.append(iris.cube.Cube(slwt_eobs, long_name="slwt_eobs"))

    wt_cube[0].add_dim_coord(tcoord, 0)
    wt_cube[1].add_dim_coord(tcoord, 0)
    wt_cube[2].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f"{work_dir}/{dataset}.nc")

    #write to csv file
    d = {'date': time_points[:], 'lwt': np.int8(WT), 'slwt_ERA5': np.int8(slwt_era5), 'slwt_EOBS': np.int8(slwt_eobs)}
    df = pd.DataFrame(data=d)
    df.to_csv(f'{work_dir}/{dataset}.csv', index=False)

    return f'weathertyping for {dataset} done.'

def get_colormap(colormap_string):

    misc_seq_2_disc = [
            (230/255, 240/255, 240/255),
            (182/255, 217/255, 228/255),
            (142/255, 192/255, 226/255),
            (118/255, 163/255, 228/255),
            (116/255, 130/255, 222/255),
            (121/255, 97/255, 199/255),
            (118/255, 66/255, 164/255),
            (107/255, 40/255, 121/255),
            (86/255, 22/255, 75/255),
            (54/255, 14/255, 36/255)]

    temp_seq_disc = [
            (254/255, 254/255, 203/255),
            (251/255, 235/255, 153/255),
            (244/255, 204/255, 104/255),
            (235/255, 167/255, 84/255),
            (228/255, 134/255, 80/255),
            (209/255, 98/255, 76/255),
            (164/255, 70/255, 66/255),
            (114/255, 55/255, 46/255),
            (66/255, 40/255, 24/255),
            (25/255, 25/255, 0/255)]
    
    prec_seq_disc = [
            (255/255, 255/255, 229/255),
            (217/255, 235/255, 213/255),
            (180/255, 216/255, 197/255),
            (142/255, 197/255, 181/255),
            (105/255, 177/255, 165/255),
            (67/255, 158/255, 149/255),
            (44/255, 135/255, 127/255),
            (29/255, 110/255, 100/255),
            (14/255, 85/255, 74/255),
            (0/255, 60/255, 48/255)]

    if colormap_string == "psl":
        return ListedColormap(misc_seq_2_disc)
    elif colormap_string == "prcp":
        return ListedColormap(prec_seq_disc)
    elif colormap_string == "temp":
        return ListedColormap(temp_seq_disc)
    else:
        print("Colormap string not supported!")

    return 

def plot_mean(wt, cfg, cube, dataset, var_name, wt_string):

    local_path = cfg['plot_dir']

    ax = plt.axes(projection=ccrs.PlateCarree())


    
    psl_cmap = get_colormap("psl")
    prcp_cmap = get_colormap("prcp")
    temp_cmap = get_colormap("temp")

    if var_name == "psl":
        plt.title(f"{var_name} mean, wt: {wt}")
        unit = "[hPa]"
        im = iplt.contourf(cube/100, cmap=psl_cmap)
    elif var_name == "prcp":
        if dataset == "ERA5":
            unit = "[m]"
            plt.title(f"total {var_name} mean, wt: {wt}")
        else:
            unit = "[kg m-2 s-1]"
            plt.title(f"{var_name} flux mean, wt: {wt}")
        im = iplt.contourf(cube, cmap=prcp_cmap)
    elif var_name == "tas":
        unit = "[K]"
        plt.title(f"1000 hPa {var_name} mean, wt: {wt}")
        im = iplt.contourf(cube, cmap=temp_cmap)

    cb = plt.colorbar(im)
    cb.ax.tick_params(labelsize=8)
    cb.set_label(label=f"{var_name} mean {unit}")

    gl = ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=True, linewidth=0.5, color='gray', alpha=0.5, linestyle='--')
    gl.left_labels = True
    gl.bottom_labels = True
    gl.top_labels = False
    gl.right_labels = False
    gl.xlines = True
    gl.ylocator = mticker.FixedLocator(np.arange(20,70,5))
    gl.xlocator = mticker.FixedLocator([-10, -5, 0, 5, 10, 15])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER
    gl.xlabel_style = {'size': 8, 'color': 'black'}
    gl.ylabel_style = {'color': 'black', 'size': 8}

    ax.set_extent([-15, 20, 27.5, 62.5])

    ax.coastlines()
    ax.add_feature(cfeature.BORDERS, linestyle=':')

    plt.savefig(os.path.join(local_path, f"{wt_string}_{wt}_{dataset}_{var_name}_mean.png"))
    plt.close()

    return

def rmse(subarray1, subarray2):
    """
    Calculate root mean square difference between two arrays.
    """
    return np.sqrt(np.mean((subarray1 - subarray2)**2))

def write_to_csv(cfg, pattern_correlation, rmse_matrix, dataset):

    work_dir = cfg['work_dir']

    df_corr = pd.DataFrame(pattern_correlation)
    df_corr.index = range(1, len(df_corr) + 1)
    df_corr.columns = range(1, len(df_corr.columns) + 1)
    df_corr.to_csv(f'{work_dir}/correlation_matrix_{dataset}.csv', index_label='Index')

    df_rmse = pd.DataFrame(rmse_matrix)
    df_rmse.index = range(1, len(df_rmse) + 1)
    df_rmse.columns = range(1, len(df_rmse.columns) + 1)
    df_rmse.to_csv(f'{work_dir}/rmse_matrix_{dataset}.csv', index_label='Index')

    return

def process_prcp_mean(wt, cfg, data, correlation_threshold, rmse_threshold, dataset):

    selected_pairs = []

    pattern_correlation = np.ma.corrcoef(data)

    n = len(data)
    rmse_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            rmse_matrix[i][j] = rmse(data[i], data[j])
            rmse_matrix[j][i] = rmse_matrix[i][j]
            if pattern_correlation[i][j] >= correlation_threshold and rmse_matrix[i][j] <= rmse_threshold:
                selected_pairs.append(((i+1, j+1), pattern_correlation[i][j], rmse_matrix[i][j]))

    print(selected_pairs)

    #write matrices to csv
    write_to_csv(cfg, pattern_correlation, rmse_matrix, dataset)

    return selected_pairs

def calculate_wt_means(cfg, cube, wt_cubes, dataset: str, var_name: str, wt_string: str):

    if wt_string == "slwt_ERA5":
        slwt_era5_cube = wt_cubes[1]#iris.load(f"{work_dir}/{file_name}.nc", "slwt_era5")
        tcoord = slwt_era5_cube.coord("time")
        slwt_era5 = slwt_era5_cube.data[:]
        print(slwt_era5)
    elif wt_string == "slwt_EOBS":
        slwt_eobs_cube = wt_cubes[2]#iris.load(f"{work_dir}/{file_name}.nc", "slwt_e-obs")
        tcoord = slwt_eobs_cube.coord("time")
        slwt_eobs = slwt_eobs_cube.data[:]
        print(slwt_eobs)
    elif wt_string == "lwt":    
        lwt_cube = wt_cubes[0]#iris.load(f"{work_dir}/{file_name}.nc", "lwt")
        tcoord = lwt_cube.coord("time")
        lwt = lwt_cube.data[:]
        print(lwt)

    if "slwt" in wt_string:
        for wt in range(1,10):
            if wt_string == "slwt_ERA5":
                target_indices = np.where(slwt_era5 == wt)
            elif wt_string == "slwt_EOBS":
                target_indices = np.where(slwt_eobs == wt)
            else:
                print("WT_STRING not supported!")
            if len(target_indices[0]) < 1:
                print(f"calculate_wt_mean - CAUTION: Skipped wt {wt} for dataset {dataset}!")
                continue
            dates = [tcoord.units.num2date(tcoord.points[i]) for i in target_indices]  
            extracted_cube = cube.extract(iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)#cube[target_indices[0], :, :].collapsed('time', iris.analysis.MEAN)
            plot_mean(wt, cfg, wt_cube_mean, dataset, var_name, wt_string)
    elif wt_string == "lwt":
        for wt in range(1,28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                print(f"calculate_wt_mean - CAUTION: Skipped wt {wt} for dataset {dataset}!")
                continue
            dates = [tcoord.units.num2date(tcoord.points[i]) for i in target_indices]
            extracted_cube = cube.extract(iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)#cube[target_indices[0], :, :].collapsed('time', iris.analysis.MEAN)
            plot_mean(wt, cfg, wt_cube_mean, dataset, var_name, wt_string)
    else:
        print("WT_STRING NOT SUPPORTED.")

    return

def combine_wt_to_file(cfg, lwt, slwt_era5, slwt_eobs, cube, file_name):

    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    write_path = cfg['work_dir']

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(lwt, long_name="lwt"))
    wt_cube.append(iris.cube.Cube(slwt_era5, long_name="slwt_era5"))
    wt_cube.append(iris.cube.Cube(slwt_eobs, long_name="slwt_eobs"))

    wt_cube[0].add_dim_coord(tcoord, 0)
    wt_cube[1].add_dim_coord(tcoord, 0)
    wt_cube[2].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f"{write_path}/{file_name}.nc")

    #write to csv file
    d = {'date': time_points[:], 'lwt': np.int8(lwt), 'slwt_era5': np.int8(slwt_era5), 'slwt_eobs': np.int8(slwt_eobs)}#, 'slwt': slwt_data}
    df = pd.DataFrame(data=d)
    df.to_csv(write_path + f'/{file_name}.csv', index=False)

    return 

def write_lwt_to_file(cfg, lwt, cube, file_name):
    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    write_path = cfg['work_dir']

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(np.int8(lwt), long_name="lwt"))

    wt_cube[0].add_dim_coord(tcoord, 0)
    iris.save(wt_cube, f"{write_path}/{file_name}.nc")

    #write to csv file
    d = {'date': time_points[:], 'lwt': np.int8(lwt)}
    df = pd.DataFrame(data=d)
    df.to_csv(write_path + f'/{file_name}.csv', index=False)

    return 

def calculate_lwt_model(cfg, cube, dataset):

    work_dir = cfg['work_dir']

    WT = wt_algorithm(cube)

    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(WT, long_name="lwt"))

    wt_cube[0].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f"{work_dir}/{dataset}.nc")

    #write to csv file
    d = {'date': time_points[:], 'lwt': np.int8(WT)}
    df = pd.DataFrame(data=d)
    df.to_csv(f'{work_dir}/{dataset}.csv', index=False)

    return f'weathertyping for {dataset} done.'

def run_my_diagnostic(cfg):
    """
    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    """
    # assemble the data dictionary keyed by dataset name
    # this makes use of the handy group_metadata function that
    # orders the data by 'dataset'; the resulting dictionary is
    # keyed on datasets e.g. dict = {'MPI-ESM-LR': [var1, var2...]}
    # where var1, var2 are dicts holding all needed information per variable
    preproc_variables_dict = group_metadata(cfg['input_data'].values(), 'dataset')

    correlation_threshold = cfg['correlation_threshold']
    rmse_threshold = cfg['rmse_threshold']
    work_dir = cfg['work_dir']

    #load cubes and run functions
    #key = dataset name, value is dataset
    if cfg['automatic_slwt']:
        for key, value in preproc_variables_dict.items():
            if key == "ERA5":
                wt_preproc = iris.load_cube(value[0]['filename'])
                wt_preproc_prcp = iris.load_cube(value[1]['filename'])
                mean_preproc_psl = iris.load_cube(value[2]['filename'])
                mean_preproc_prcp = iris.load_cube(value[3]['filename'])
                mean_preproc_tas = iris.load_cube(value[4]['filename'])
                wt_preproc_prcp_EOBS = iris.load_cube(preproc_variables_dict['E-OBS'][0]['filename'])

                #calculate lwt
                lwt = wt_algorithm(wt_preproc)

                #calculate simplified lwt based on precipitation patterns
                slwt_era5 = calculate_slwt_obs(cfg, lwt, wt_preproc_prcp, key, correlation_threshold, rmse_threshold)
                slwt_eobs = calculate_slwt_obs(cfg, lwt, wt_preproc_prcp_EOBS, "E-OBS", correlation_threshold, rmse_threshold)

                #write to file
                combine_wt_to_file(cfg, lwt, slwt_era5, slwt_eobs, wt_preproc, key)

                #load weathertype files as cubes
                lwt_cube = iris.load_cube(f"{work_dir}/{key}.nc", "lwt")
                slwt_era5_cube = iris.load_cube(f"{work_dir}/{key}.nc", "slwt_era5")
                slwt_eobs_cube = iris.load_cube(f"{work_dir}/{key}.nc", "slwt_eobs")
                wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]

                #plot means
                #calculate_wt_means(cfg, mean_preproc_psl, wt_cubes, key, var_name="psl", wt_string="lwt")
                #calculate_wt_means(cfg, mean_preproc_prcp, wt_cubes, key, var_name="prcp", wt_string="lwt")
                #calculate_wt_means(cfg, mean_preproc_tas, wt_cubes, key, var_name="tas", wt_string="lwt")
                #calculate_wt_means(cfg, mean_preproc_psl, wt_cubes, key, var_name="psl", wt_string="slwt_ERA5")
                #calculate_wt_means(cfg, mean_preproc_prcp, wt_cubes, key, var_name="prcp", wt_string="slwt_ERA5")
                #calculate_wt_means(cfg, mean_preproc_tas, wt_cubes, key, var_name="tas", wt_string="slwt_ERA5")
                #calculate_wt_means(cfg, mean_preproc_psl, wt_cubes, key, var_name="psl", wt_string="slwt_EOBS")
                #calculate_wt_means(cfg, mean_preproc_prcp, wt_cubes, key, var_name="prcp", wt_string="slwt_EOBS")
                #calculate_wt_means(cfg, mean_preproc_tas, wt_cubes, key, var_name="tas", wt_string="slwt_EOBS")
            else: 
                if key == "E-OBS":
                    continue
                wt_preproc = iris.load_cube(value[0]['filename'])
                mean_preproc_psl = iris.load_cube(value[1]['filename'])
                mean_preproc_prcp = iris.load_cube(value[2]['filename'])          
                mean_preproc_tas = iris.load_cube(value[3]['filename'])  

                #calculate weathertypes
                calculate_lwt_slwt_model(cfg, wt_preproc, key)  

                #load wt files
                lwt_cube = iris.load_cube(f"{work_dir}/{key}.nc", "lwt")
                slwt_era5_cube = iris.load_cube(f"{work_dir}/{key}.nc", "slwt_era5")
                slwt_eobs_cube = iris.load_cube(f"{work_dir}/{key}.nc", "slwt_eobs")
                wt_cubes = [lwt_cube, slwt_era5_cube, slwt_eobs_cube]   

                #plot means
                #calculate_wt_means(cfg, mean_preproc_psl, wt_cubes, key, var_name="psl", wt_string="lwt")
                #calculate_wt_means(cfg, mean_preproc_prcp, wt_cubes,  key, var_name="prcp", wt_string="lwt")
                #calculate_wt_means(cfg, mean_preproc_tas, wt_cubes, key, var_name="tas", wt_string="lwt")
                #calculate_wt_means(cfg, mean_preproc_psl, wt_cubes, key, var_name="psl", wt_string="slwt_ERA5")
                #calculate_wt_means(cfg, mean_preproc_prcp, wt_cubes, key, var_name="prcp", wt_string="slwt_ERA5")
                #calculate_wt_means(cfg, mean_preproc_tas, wt_cubes, key, var_name="tas", wt_string="slwt_ERA5")
                #calculate_wt_means(cfg, mean_preproc_psl, wt_cubes, key, var_name="psl", wt_string="slwt_EOBS")
                #calculate_wt_means(cfg, mean_preproc_prcp, wt_cubes, key, var_name="prcp", wt_string="slwt_EOBS")
                #calculate_wt_means(cfg, mean_preproc_tas, wt_cubes, key, var_name="tas", wt_string="slwt_EOBS")
    else: #if automatic_slwt is false, just do the lwt, DOES NOT WORK YET
        for key, value in preproc_variables_dict.items():
            if key == "ERA5":
                wt_preproc = iris.load_cube(value[0]['filename'])
                wt_preproc_prcp = iris.load_cube(value[1]['filename'])
                mean_preproc_psl = iris.load_cube(value[2]['filename'])
                mean_preproc_prcp = iris.load_cube(value[3]['filename'])
                mean_preproc_tas = iris.load_cube(value[4]['filename'])
                wt_preproc_prcp_EOBS = iris.load_cube(preproc_variables_dict['E-OBS'][0]['filename'])

                #calculate lwt
                lwt = wt_algorithm(wt_preproc)

                #write to file
                write_lwt_to_file(cfg, lwt, wt_preproc, key)

                #load weathertype files as cubes
                lwt_cube = iris.load_cube(f"{work_dir}/{key}.nc", "lwt")
                wt_cubes = [lwt_cube]

                #plot means
                calculate_wt_means(cfg, mean_preproc_psl, wt_cubes, key, var_name="psl", wt_string="lwt")
                calculate_wt_means(cfg, mean_preproc_prcp, wt_cubes, key, var_name="prcp", wt_string="lwt")
                calculate_wt_means(cfg, mean_preproc_tas, wt_cubes, key, var_name="tas", wt_string="lwt")
            else: 
                if key == "E-OBS":
                    continue
                wt_preproc = iris.load_cube(value[0]['filename'])
                mean_preproc_psl = iris.load_cube(value[1]['filename'])
                mean_preproc_prcp = iris.load_cube(value[2]['filename'])          
                mean_preproc_tas = iris.load_cube(value[3]['filename'])  

                #calculate weathertypes
                calculate_lwt_model(cfg, wt_preproc, key)  

                #load wt files
                lwt_cube = iris.load_cube(f"{work_dir}/{key}.nc", "lwt")
                wt_cubes = [lwt_cube]   

                #plot means
                calculate_wt_means(cfg, mean_preproc_psl, wt_cubes, key, var_name="psl", wt_string="lwt")
                calculate_wt_means(cfg, mean_preproc_prcp, wt_cubes,  key, var_name="prcp", wt_string="lwt")
                calculate_wt_means(cfg, mean_preproc_tas, wt_cubes, key, var_name="tas", wt_string="lwt")

    return

if __name__ == '__main__':

    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
        
