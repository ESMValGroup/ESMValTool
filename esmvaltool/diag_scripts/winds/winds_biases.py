import cartopy.crs as ccrs 
import esmvalcore.preprocessor as eprep
import iris
import iris.cube
import iris.plot as iplt
import logging
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import numpy as np
import glob

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))

def update_mpl():

    mpl.rcParams['font.size'] = 18
    mpl.rcParams['axes.spines.top'] = False
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['figure.titlesize'] = 'medium'

    return

def get_levels(aver_info: str):

    if aver_info == 'Daily':
        mean_levs = np.arange(-25, 1)
        std_levs =  np.arange(0, 36)
    elif int(aver_info[:-2])<25:
        mean_levs = np.arange(-15, 1)
        std_levs =  np.arange(0, 26)
    else:
        mean_levs = np.arange(-25, 1)
        std_levs =  np.arange(0, 36)
    
    return mean_levs, std_levs


def get_derived_wind(dataset_info : list):

    try:
        wind_from_h = select_metadata(dataset_info, variable_group='derived')[0]['filename']
        wind_from_h_cb = iris.load_cube(wind_from_h)
    except:
        uas = select_metadata(dataset_info, variable_group='uas')[0]['filename']
        uas_cb = iris.load_cube(uas)
        vas = select_metadata(dataset_info, variable_group='vas')[0]['filename']
        vas_cb = iris.load_cube(vas)
        wind_cblst = iris.cube.CubeList([uas_cb, vas_cb])
        wind_from_h_cb = eprep._derive.sfcwind.DerivedVariable().calculate(wind_cblst)

    return wind_from_h_cb


def get_diff(wind_from_h_cb : iris.cube.Cube, wind_h_cb : iris.cube.Cube):

    try:
        diff_cb = (wind_from_h_cb - wind_h_cb)*100/wind_h_cb
    except:
        diff_cb = (wind_from_h_cb - wind_h_cb.data)*100/wind_h_cb.data

    return diff_cb


def plot_fig(cube : iris.cube.Cube, levels : np.ndarray, year: int,
                                    dataset : str, cfg : dict, metric : str):
    
    update_mpl()
    
    cmap = 'Blues_r' if metric == 'mean' else 'Reds'
    img_form = diagtools.get_image_format(cfg)
    hours = cfg.get('hours')

    fig, ax = plt.subplots(ncols=1, nrows=1,
                                    subplot_kw={'projection': ccrs.Robinson()})
    fig.set_size_inches(8, 6)
    levs = iplt.contourf(cube, axes=ax, levels=levels, cmap=cmap, extend='both')
    ax.coastlines(linewidth=0.5)
    # fig.suptitle(f'{hours} {dataset} {metric} wind speed difference in {year}',
    #                                                               fontsize='xx-large')
    fig.subplots_adjust(left=0.01, bottom=0.1, right=0.99, top=0.9999)
    cax = fig.add_axes([0.1,0.11,0.8,0.02])
    fig.colorbar(levs, cax=cax, orientation='horizontal', label='Difference, %')
    # fig.savefig(os.path.join(cfg['plot_dir'], 
    #                          f'{hours}_{metric}_diff_{dataset}_{year}{img_form}'))
    fig.savefig(os.path.join(cfg['plot_dir'], 
                             f'{hours}_{metric}_diff_{dataset}_{year}.png'))

    return

def main(cfg):

    mean_levs, std_levs = get_levels(cfg['hours'])

    input_data = cfg['input_data']

    datasets = group_metadata(input_data.values(), 'dataset', sort=True)

    for dataset in datasets.keys():
        vars_data = datasets[dataset]
        # get real wind
        wind_h = select_metadata(vars_data, variable_group='real')[0]['filename']
        wind_h_cb = iris.load_cube(wind_h)
        # get derived wind
        wind_from_h_cb = get_derived_wind(vars_data)

        y = str(input_data[wind_h]['start_year'])
        
        diff_cb = get_diff(wind_from_h_cb, wind_h_cb)

        mean_diff = diff_cb.collapsed('time', iris.analysis.MEAN)
        std_diff = diff_cb.collapsed('time', iris.analysis.STD_DEV)

        plot_fig(mean_diff, mean_levs, y, dataset, cfg, 'mean')
        # plot_fig(std_diff, std_levs, y, dataset, cfg, 'StDev')

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
