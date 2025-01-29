import csv
import iris
import climextremes as cex
import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import genextreme as gev

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def calculate_stds(data_dic, obs_gev_data, cfg):   

    era_var = obs_gev_data['ano_obs_cb'].data.std()
    
    col_mod = (25/255, 14/255, 79/255)
    col_obs = 'indianred'

    fig_stds, ax_stds = plt.subplots(nrows=1, ncols=1)

    fig_stds.set_size_inches(8., 9.)
    fig_stds.set_dpi(200)

    y_ticks = np.arange(0, len(data_dic.keys()))
    y_labs = np.zeros(len(data_dic.keys()), dtype='<U30')

    for nm, model in enumerate(data_dic.keys()):
        if model != 'Multi-Model-Mean':
            for i in range(0, len(data_dic[model]['var_data'])):
                ax_stds.scatter(data_dic[model]['var_data'][i].data.std(), nm+1, s=50, marker='o', c=col_mod, zorder=2)
            y_labs[nm+1] = model
        else:
            wghts = list(); wdata = list()
            for i in range(0, len(data_dic[model]['var_data'])):
                wghts.append(np.full(data_dic[model]['var_data'][i].shape, data_dic[model]['var_data'][i].attributes['ensemble_weight']*data_dic[model]['var_data'][i].attributes['reverse_dtsts_n']))
                wdata.append(data_dic[model]['var_data'][i].data)
            wghts = np.asarray(wghts).flatten(); wdata = np.asarray(wdata).flatten()
            if cfg['model_weighting']:
                mmm = np.average(wdata, weights=wghts)
                mmm_std = np.sqrt(np.average((wdata - mmm)**2, weights=wghts))
            else: 
                mmm_std = wdata.std()
            ax_stds.scatter(mmm_std, 0, c = col_mod, s=70, marker='s', zorder=3)
            y_labs[0] = 'CMIP6'
    
    ax_stds.axvline(era_var, -1, len(data_dic.keys()) + 1, c=col_obs, zorder=2, lw=3)
    ax_stds.set_ylim(len(data_dic.keys()) -0.8, -0.2)

    ax_stds.set_yticks(y_ticks, labels=y_labs)
    ax_stds.grid(which='both', c='silver', zorder=1)

    ax_stds.set_xlabel('StD of '+cfg['var_label'] + ' anomalies, C')
    ax_stds.text(0.02, 0.97, cfg.get('litera')+' '+cfg['region'], fontsize='xx-large', transform=ax_stds.transAxes)

    # fig_stds.suptitle('Standard deviations (StD) of '+cfg['var_label']+' anomalies in '+cfg['region']+' relative to '+ str(cfg['reference_period'][0]) \
    #         + '-' + str(cfg['reference_period'][1]))

    plt.tight_layout()

    fig_stds.savefig(os.path.join(cfg['plot_dir'], 'stds_'+cfg['region'].lower()+'_'+cfg['var_label'].lower()+'.'+cfg['output_file_type']))           


    return

def create_timeseries(data_dic, mixns, obs_gev_data, cfg):

    b_max = np.ceil(mixns['max']*1.05)
    b_min = np.floor(mixns['min']*1.05)

    col_mod = (25/255, 14/255, 79/255)
    col_obs = 'indianred'

    t_s = np.arange(cfg['year_span'][0], cfg['year_span'][1]+1)

    for dataset in data_dic.keys(): 
        fig_ts, ax_ts = plt.subplots(nrows=1, ncols=1)
        fig_ts.set_size_inches(10., 8.)
        fig_ts.set_dpi(200)

        model_var_data = list()
        weights = list()

        for i in range(0,len(data_dic[dataset]['var_data'])):
            model_var_data.append(data_dic[dataset]['var_data'][i].data)
            weights.append(data_dic[dataset]['var_data'][i].attributes['ensemble_weight']*
                            data_dic[dataset]['var_data'][i].attributes['reverse_dtsts_n'])
        weights = np.array(weights)
        model_var_data = np.array(model_var_data)
        mean_var_arr = np.average(model_var_data, axis=0)
        min_var_arr = np.min(model_var_data, axis=0)
        max_var_arr = np.max(model_var_data, axis=0)

        ax_ts.text(0.47, 0.95, cfg['var_label'], transform=ax_ts.transAxes)
        ax_ts.plot(t_s, obs_gev_data['ano_obs_cb'].data, c=col_obs, zorder=5, label='ERA5')
        ax_ts.plot(t_s, mean_var_arr, c=col_mod, zorder=3, label=dataset)
        if len(data_dic[dataset]['var_data'])>1:
            # clean this noncense before submitting! 
            ax_ts.fill_between(t_s, min_var_arr, max_var_arr, color=col_mod, lw=0, alpha=0.25)
        ax_ts.plot([t_s-0.5, t_s+0.5], [0,0], c='grey')
        ax_ts.set_xlim(t_s[0]-0.5, t_s[-1]+0.5)
        ax_ts.set_ylim(b_min, b_max)
        ax_ts.legend(loc=0, ncols=2, fancybox=False, frameon=False)
        ax_ts.set_ylabel(cfg['var_label'] + ' anomaly, C')
        if dataset == 'Multi-Model-Mean':
            fig_ts.suptitle(cfg.get('litera')+' '+cfg['var_label']+' anomalies in '+cfg['region']+' relative to '+ str(cfg['reference_period'][0]) \
            + '-' + str(cfg['reference_period'][1])+ ' from '+str(len(data_dic.keys()) - 1) +' CMIP6 models' , fontsize = 'x-large')
        else:
            fig_ts.suptitle(cfg.get('litera')+' '+cfg['var_label']+' anomalies in '+cfg['region']+' relative to '+ str(cfg['reference_period'][0]) \
            + '-' + str(cfg['reference_period'][1])+ '  from '+dataset, fontsize = 'x-large')

        plt.tight_layout()

        fig_ts.savefig(os.path.join(cfg['plot_dir'], 'ts_'+cfg['region'].lower()+'_'+cfg['var_label'].lower()+'_'+dataset+'.'+cfg['output_file_type']))

    return


def main(cfg):

    input_data = cfg['input_data']

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)

    ano_obs_cb, obs_gev_data = obtain_obs_info(groups, cfg)

    # change later, so far to identify the group name
    group = [k for k in groups.keys() if 'txnx' in k][0]

    plotting_dic = {}
    mixns = {}

    mins = list(); maxs = list()
    group_data = groups[group]
    datasets = group_metadata(group_data, 'dataset')
    ens_var_cubelist = iris.cube.CubeList()
    for dataset in datasets.keys():
        filepaths = list(group_metadata(datasets[dataset], 'filename').keys())
        n_real = len(filepaths)
        mod_var_cubelist = iris.cube.CubeList()
        for filepath in filepaths:
            mod_cb = iris.load_cube(filepath)
            # adding weights to the data cubes 
            mod_cb.attributes['ensemble_weight'] = 1 / n_real
            mod_cb.attributes['reverse_dtsts_n'] = 1/ len(datasets)
            ens_var_cubelist.append(mod_cb); mod_var_cubelist.append(mod_cb)
            mins.append(mod_cb.collapsed('time', iris.analysis.MIN).data)
            maxs.append(mod_cb.collapsed('time', iris.analysis.MAX).data)
        plotting_dic[dataset] = {'var_data': mod_var_cubelist}
    mixns['max'] = np.asarray(maxs).max()
    mixns['min'] = np.asarray(mins).min()
    plotting_dic['Multi-Model-Mean'] = {'var_data' : ens_var_cubelist}

    st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(st_file)

    create_timeseries(plotting_dic, mixns, obs_gev_data, cfg)

    calculate_stds(plotting_dic, obs_gev_data, cfg)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
