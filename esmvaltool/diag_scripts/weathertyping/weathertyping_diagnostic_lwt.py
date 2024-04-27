#weathertyping diagnostic for esmvaltool


# operating system manipulations (e.g. path constructions)
import os
import numpy as np
import pandas as pd
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

def calc_slwt(wt):
        mapping = {
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

        return np.array([mapping[value] for value in wt])

def calculate_weathertypes(cfg, cube, dataset):
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

    slwt_data = calc_slwt(np.int8(WT))

    iris.FUTURE.datum_support = True
    iris.FUTURE.save_split_attrs = True
    tcoord = cube.coord("time")
    time_points = tcoord.units.num2date(tcoord.points)

    write_path = cfg['work_dir']

    wt_cube = iris.cube.CubeList()
    wt_cube.append(iris.cube.Cube(WT, long_name="Lamb Weathertypes"))
    wt_cube.append(iris.cube.Cube(slwt_data, long_name="Simplified Lamb Weathertypes"))

    wt_cube[0].add_dim_coord(tcoord, 0)
    wt_cube[1].add_dim_coord(tcoord, 0)

    iris.save(wt_cube, f"{write_path}/{dataset}.nc")

    #write to csv file
    d = {'date': time_points[:], 'lwt': np.int8(WT), 'slwt': slwt_data}
    df = pd.DataFrame(data=d)
    df.to_csv(write_path + f'/{dataset}.csv', index=False)

    return f'weathertyping for {dataset} done.'

def plot_mean(wt, cfg, cube, dataset, var_name, wt_string):

    local_path = cfg['plot_dir']

    ax = plt.axes(projection=ccrs.PlateCarree())

    imola_rgb = [(25/255, 51/255, 178/255),
                (36/255, 70/255, 168/255),
                (45/255, 89/255, 159/255),
                (57/255, 106/255, 147/255),
                (73/255, 123/255, 132/255),
                (95/255, 146/255, 123/255),
                (122/255, 173/255, 116/255),
                (152/255, 203/255, 108/255),
                (195/255, 233/255, 102/255),
                (255/255, 254/255, 102/255)]
    
    imola_cmap = ListedColormap(imola_rgb)

    if var_name == "psl":
        plt.title(f"{var_name} mean, wt: {wt}")
        unit = "[hPa]"
        im = iplt.contourf(cube/100, cmap=imola_cmap)
    elif var_name == "prcp":
        if dataset == "ERA5":
            unit = "[m]"
            plt.title(f"total {var_name} mean, wt: {wt}")
        else:
            unit = "[kg m-2 s-1]"
            plt.title(f"{var_name} flux mean, wt: {wt}")
        im = iplt.contourf(cube, cmap=imola_cmap)
    elif var_name == "tas":
        unit = "[K]"
        plt.title(f"1000 hPa {var_name} mean, wt: {wt}")
        im = iplt.contourf(cube, cmap=imola_cmap)

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

def rmsd(subarray1, subarray2):
    """
    Calculate root mean square difference between two arrays.
    """
    return np.sqrt(np.mean((subarray1 - subarray2)**2))

def process_prcp_mean(wt, cfg, data, dataset, var_name, wt_string):

    correlation_threshold = 0.7
    rmsd_threshold = 0.00002
    selected_pairs = []

    pattern_correlation = np.ma.corrcoef(data)
    print(pattern_correlation)

    n = len(data)
    rmsd_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(i+1, n):
            rmsd_matrix[i][j] = rmsd(data[i], data[j])
            rmsd_matrix[j][i] = rmsd_matrix[i][j]
            if pattern_correlation[i][j] >= correlation_threshold and rmsd_matrix[i][j] <= rmsd_threshold:
                selected_pairs.append(((i+1, j+1), pattern_correlation[i][j], rmsd_matrix[i][j]))
            
    print(rmsd_matrix)

    print(selected_pairs)

    work_dir = cfg['work_dir']

    np.savetxt(f"{work_dir}/rmsd_matrix.txt", rmsd_matrix)
    np.savetxt(f"{work_dir}/corr_matrix.txt", pattern_correlation)

    return

def calculate_wt_means(cfg, cube, dataset: str, var_name: str, wt_string: str="slwt"):

    work_dir = cfg['work_dir']

    lwt_cube = iris.load_cube(f"{work_dir}/{dataset}.nc", "Lamb Weathertypes")
    #lwt_cube = iris.load_cube(f"{work_dir}/ERA5.nc", "Lamb Weathertypes") #for EOBS
    slwt_cube = iris.load_cube(f"{work_dir}/{dataset}.nc", "Simplified Lamb Weathertypes")
    # = iris.load_cube(f"{work_dir}/ERA5.nc", "Simplified Lamb Weathertypes") #for EOBS
    lwt = lwt_cube.data[:]
    slwt = slwt_cube.data[:]

    if wt_string == "slwt":
        tcoord = slwt_cube.coord("time")
    elif wt_string == "lwt":    
        tcoord = lwt_cube.coord("time")

    if wt_string == "slwt":
        for wt in range(1,10):
            target_indices = np.where(slwt == wt)
            dates = [tcoord.units.num2date(tcoord.points[i]) for i in target_indices]  
            extracted_cube = cube.extract(iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)#cube[target_indices[0], :, :].collapsed('time', iris.analysis.MEAN)
            plot_mean(wt, cfg, wt_cube_mean, dataset, var_name, wt_string)
    elif wt_string == "lwt":
        wt_data_prcp = []
        for wt in range(1,28):
            target_indices = np.where(lwt == wt)
            if len(target_indices[0]) < 1:
                continue
            dates = [tcoord.units.num2date(tcoord.points[i]) for i in target_indices]
            extracted_cube = cube.extract(iris.Constraint(time=lambda t: t.point in dates[0]))
            wt_cube_mean = extracted_cube.collapsed('time', iris.analysis.MEAN)#cube[target_indices[0], :, :].collapsed('time', iris.analysis.MEAN)
            plot_mean(wt, cfg, wt_cube_mean, dataset, var_name, wt_string)
            #wt_data_prcp.append(wt_cube_mean.data.compressed())
        #(wt, cfg, wt_data_prcp, dataset, var_name, wt_string)
    else:
        print("WT_STRING NOT SUPPORTED.")

    return

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
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')

    #load cube
    #keys are datasets: models, era5
    #values are variables
    for key, value in my_files_dict.items():
        if key == "ERA5" and len(value) > 1:
            wt_preproc = iris.load_cube(value[0]['filename'])
            mean_preproc_psl = iris.load_cube(value[1]['filename'])
            wt_preproc_prcp = iris.load_cube(value[2]['filename'])
            mean_preproc_tas = iris.load_cube(value[3]['filename'])
            calculate_weathertypes(cfg, wt_preproc, key)
            calculate_wt_means(cfg, wt_preproc_prcp, key, var_name="prcp", wt_string="lwt")
            calculate_wt_means(cfg, mean_preproc_psl, key, var_name="psl", wt_string="lwt")
            calculate_wt_means(cfg, mean_preproc_tas, key, var_name="tas", wt_string="lwt")
        elif key == "ERA5" and len(value) == 1:
            wt_preproc = iris.load_cube(value[0]['filename'])
            calculate_weathertypes(cfg, wt_preproc, key)
        else: 
            wt_preproc = iris.load_cube(value[0]['filename'])
            mean_preproc_psl = iris.load_cube(value[1]['filename'])
            mean_preproc_prcp = iris.load_cube(value[2]['filename'])          
            mean_preproc_temp = iris.load_cube(value[3]['filename'])  
            calculate_weathertypes(cfg, wt_preproc, key)        
            calculate_wt_means(cfg, mean_preproc_psl, key, var_name="psl", wt_string="lwt")
            calculate_wt_means(cfg, mean_preproc_prcp, key, var_name="prcp", wt_string="lwt")       
            calculate_wt_means(cfg, mean_preproc_temp, key, var_name="tas", wt_string="lwt")       


    return

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
