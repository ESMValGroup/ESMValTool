import logging
import cf_units
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
    diffs = data_cube - mean_cb
    sqrt_f = iris.analysis.maths.IFunc(np.sqrt, lambda cube: cf_units.Unit('1'))
    sqrt_var_cb = sqrt_f(eprep.climate_statistics(diffs**2, operator='mean'))
    skew_cb = eprep.climate_statistics(diffs**3, operator='mean')/sqrt_var_cb**3
    skew_cb.standard_name=mean_cb.standard_name
    skew_cb.long_name=mean_cb.long_name
    skew_cb.var_name=mean_cb.var_name
    skew_cb.units=mean_cb.units
    skew_cb.attributes=mean_cb.attributes

    data_dic = { 'mean': mean_cb,
                 'std': std_cb,
                 'skew': skew_cb,
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
        plot_dic = {'mean': {'levels': np.linspace(-20,20,20), 'cbar': 'RdBu'}, 
                    'std': {'levels': np.linspace(-8,8,15), 'cbar': 'RdBu'},
                    'skew': {'levels': np.linspace(-2,2,10), 'cbar': 'RdBu'},
                    'p95': {'levels': np.linspace(-20,20,20), 'cbar': 'RdBu'}}        
        
    fig_maps, ax_maps = plt.subplots(ncols=ncols, nrows=nrows,
                                      subplot_kw={'projection': ccrs.Robinson()})
    ax_maps = ax_maps.flatten()
    
    n=0
    for nseas, seas in enumerate(data_dic.keys()):
        for nstat, stat in enumerate(data_dic[seas].keys()):
            iplt.contourf(data_dic[seas][stat], axes=ax_maps[n], 
                          levels=plot_dic[stat]['levels'], extend='both', 
                          cmap=plot_dic[stat]['cbar'])
            ax_maps[n].coastlines(linewidth=0.5)
            ax_maps[n].set_title(stat)
            n = n + 1


    fig_name = 'wv_' ; fig_title = 'Water vapour '
    if bias: 
        fig_name = fig_name + 'bias_'
        fig_title = fig_title + 'bias '
    fig_name = fig_name + dataset
    fig_title = fig_title + dataset
    if ensemble != None:
        fig_name = fig_name+ '_' +ensemble
        fig_title = fig_title + ' ' +ensemble
    fig_title = fig_title +'('+cfg['time_range']+')'
    fig_maps.suptitle(fig_title)
    
    fig_maps.savefig(os.path.join(cfg['plot_dir'], fig_name + diagtools.get_image_format(cfg)))
    plt.close(fig_maps)

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
            data_dic = {}
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
                data_dic [group] = calculate_cube_stats(data_cube)
            plot_maps(data_dic, dataset, cfg, real)
            if dataset != ref_dataset: 
                bias_dic = {}
                for group in data_groups: 
                    bias_dic[group] = {}
                    for stat in data_dic[group]:
                        bias_dic[group][stat] = data_dic[group][stat] - ref_dic[group][stat]
                plot_maps(data_dic, dataset, cfg, real, bias=True)


    logger.info('Success')

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)