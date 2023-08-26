import logging
import iris
import iris.plot as iplt
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs 
import esmvalcore.preprocessor as eprep
import scipy
import os

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_data, ProvenanceLogger
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def calculate_cube_stats(data_cube):

    data_dic = {}

    mean_cb = eprep.climate_statistics(data_cube, operator='mean')
    std_cb = eprep.climate_statistics(data_cube, operator='std_dev')
    p95_cb = data_cube.collapsed('time', iris.analysis.PERCENTILE, percent=95)
    skew_cb = mean_cb ; skew_cb.data = scipy.stats.skew(data_cube.data, axis=0)

    data_dic = { 'mean': mean_cb,
                 'std': std_cb,
                 'skewnes': skew_cb,
                 'p95': p95_cb}

    return data_dic

def plot_maps(data_dic, dataset, cfg, ensemble=None, bias=False):

    n_seas = len(data_dic.keys())
    if n_seas == 1: 
        ncols = 2 ; nrows = 2
    else:
        ncols = 4; nrows = n_seas        

    plot_dic = {'mean': {'levels': np.linspace(0,65,25), 'cbar': 'GnBu'}, 
                'std': {'levels': np.linspace(0,21,21), 'cbar': 'GnBu'},
                'skew': {'levels': np.linspace(-2.2,2.2,24), 'cbar': 'RdBu'},
                'p95': {'levels': np.linspace(0,70,20), 'cbar': 'GnBu'}}
    if bias: 
        plot_dic = {'mean': {'levels': np.linspace(-5,5,20), 'cbar': 'RdBu'}, 
                    'std': {'levels': np.linspace(-2,2,10), 'cbar': 'RdBu'},
                    'skew': {'levels': np.linspace(-1,1,20), 'cbar': 'RdBu'},
                    'p95': {'levels': np.linspace(-5,5,20), 'cbar': 'RdBu'}}        
        
    fig_maps, ax_maps = plt.subplots(ncols=ncols, nrows=nrows,
                                      subplot_kw={'projection': ccrs.Robinson()})
    ax_maps = ax_maps.flatten()
    
    for nseas, seas in enumerate(data_dic.keys()):
        n = nseas
        for nstat, stat in enumerate(data_dic[seas].keys()):
            n = n + nstat
            iplt.contourf(data_dic[seas][stat], axes=ax_maps[n], 
                          levels=plot_dic[stat]['levels'], extend='both', 
                          cmap=plot_dic[stat]['cbar'])


    fig_name = 'wv_'
    if bias: 
        fig_name = fig_name + 'bias_'
    fig_name = fig_name + dataset
    if ensemble != None:
        fig_name = fig_name+ '_' +ensemble
    
    fig_maps.savefig(os.path.join(cfg['plot_dir'], fig_name + diagtools.get_image_format(cfg)))

    return


def main(cfg):

    input_data = cfg['input_data']

    data_groups = list(group_metadata(input_data.values(), 'variable_group').keys())

    ref_dic = {}

    ref_dataset = list(group_metadata(input_data.values(), 'reference_dataset').keys())[0]
    for group in data_groups: 
        ref_f_name = select_metadata(input_data.values(), 
                                       dataset=ref_dataset, 
                                       variable_group=group)[0]['filename']
        ref_cube = iris.load_cube(ref_f_name)
        ref_dic[group] = calculate_cube_stats(ref_cube)

    datasets = group_metadata(input_data.values(), 'dataset', sort=True)

    for dataset in datasets.keys():
        reals = group_metadata(datasets[dataset], 'ensemble')
        for real in reals:
            for group in data_groups:
                if real == None:
                    dset_fname = select_metadata(input_data.values(), 
                                            dataset=dataset, 
                                            variable_group=group)[0]['filename']
                else: 
                    dset_fname = select_metadata(input_data.values(), 
                                            dataset=dataset, ensemble = real,
                                            variable_group=group)[0]['filename']
                data_cube = iris.load_cube(dset_fname)
                data_dic = calculate_cube_stats(data_cube)
            plot_maps(data_dic, dataset, cfg, real)
            if dataset != ref_dataset: 
                bias_dic = {}
                for group in data_groups: 
                    bias_dic[group] = data_dic[group] - ref_dic[group]
                plot_maps(data_dic, dataset, cfg, real, bias=True)
            
            logger.info("Success")


        




    logger.info('Success')

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)