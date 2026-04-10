import iris
import math
import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
import os

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, group_metadata, save_figure
import esmvaltool.diag_scripts.shared.plot as eplot
# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def obtain_reference(data_group: list):
    '''
    This function cretes a dictionary with reference data.

    Parameters:
    -----------
    data_group: 
        list with metadata for reference group
    
    Returns:
    --------
    reference_dic:
        dictionary with the reference data. Dataset name is a keyword
    '''

    reference_dic = {}

    ref_fnames = group_metadata(data_group, 'filename')

    for dataset_f in ref_fnames.keys():
        dataset_n = ref_fnames[dataset_f][0]['dataset']
        if len(group_metadata(data_group, 'dataset')[dataset_n])>1:
            key = ref_fnames[dataset_f][0]['alias']
        else: 
            key = dataset_n
        ref_cb = iris.load_cube(dataset_f)
        reference_dic[key] = ref_cb.data

    return reference_dic


def create_provenance(caption: str):
    '''Creates provenance dictionary'''

    provenance_dic = {'authors': ['malinina_elizaveta'], 
                      'caption': caption,
                      'references': ['malinina24']}

    return provenance_dic


def perkins_skill_score(model_data: np.ndarray, obs_dict: dict, bin_width: float| None):
    '''
    This function calculates the Perkins skill score for a given model and observation data.
    Parameters:
    -----------
    model_data:
        numpy array with the model data from all realizations
    obs_dict:
        dictionary with the observation data, keys are dataset names and values are numpy arrays
    cfg:
        config dictionary, comes from ESMValCore

    Returns:
    --------    
    pss_dic:
        dictionary with the Perkins skill score for each observation dataset, keys are observation datasets
    '''

    bin_width = 1 if bin_width is None else bin_width # TODO add nbins option

    mod_max, mod_min = np.nanmax(model_data), np.nanmin(model_data)

    pss_dic = {}

    for obs_dataset in obs_dict.keys():
        pss_dic[obs_dataset] = {}
        obs_data = obs_dict[obs_dataset]
        obs_max, obs_min = np.nanmax(obs_data), np.nanmin(obs_data) 
        bins = np.arange(np.floor(min(mod_min, obs_min)), np.ceil(max(mod_max, obs_max) + bin_width), bin_width)
        mod_hist, bin_edges = np.histogram(model_data, bins=bins, density=True)
        obs_hist, _ = np.histogram(obs_data, bins=bins, density=True)

        pss = float(np.sum(np.minimum(mod_hist, obs_hist))*bin_width)
        
        pss_dic[obs_dataset]['pss'] = pss
        pss_dic[obs_dataset]['bin_edges'] = bin_edges
        pss_dic[obs_dataset]['mod_hist'] = mod_hist
        pss_dic[obs_dataset]['obs_hist'] = obs_hist

    return pss_dic


def save_skills(data_dic, cfg):

    region = cfg.get('region', 'region')
    variable = cfg.get('var_label')

    skills_info = []

    for model in data_dic.keys():
        for obs_dataset in data_dic[model].keys():
            pss = data_dic[model][obs_dataset]['pss']
            dataset_dic = {'region': region, 'dataset': model, 'reference': obs_dataset, 'pss': pss}
            skills_info.append(dataset_dic)

    skills_df = pd.DataFrame(skills_info)

    skills_df.to_csv(os.path.join(cfg['work_dir'], 
                    f'perkins_skill_score_{variable}_{region}.csv'), index=False)
    
    return


def plot_skills(data_dic, obs_name_list, cfg):
    '''
    This function plots perkins skill score for each model and observation dataset. 

    Parameters:
    -----------
    data_dic:
        dictionary with the data to plot
    obs_name_list:
        list with the names of the observation datasets
    cfg:
        config dictionary, comes from ESMValCore
    '''

    mpl_st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(mpl_st_file)

    variable = cfg.get('var_label')
    exp_variable = variable.replace('_', ' ')
    region = cfg['region'] if cfg.get('region') else 'region'
    
    # create near square subplots
    n_models = len(data_dic.keys())
    ncols = math.ceil(math.sqrt(n_models))
    nrows = math.ceil(n_models / ncols)

    for obs_name in obs_name_list:
        fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 3*nrows), squeeze=False)
        for i, model in enumerate(data_dic.keys()):
            row, col = divmod(i, ncols)
            ax = axes[row, col]
            pss = data_dic[model][obs_name]['pss']
            bin_edges = data_dic[model][obs_name]['bin_edges']
            mod_hist = data_dic[model][obs_name]['mod_hist']
            obs_hist = data_dic[model][obs_name]['obs_hist']

            ax.bar(bin_edges[:-1], mod_hist, width=np.diff(bin_edges), 
                   align='edge', alpha=0.5, label='Model', color='#1A0F50')
            ax.bar(bin_edges[:-1], obs_hist, width=np.diff(bin_edges),
                    align='edge', alpha=0.5, label=obs_name, color='#cd5c5c')
            ax.set_title(f'{model} (PSS: {pss:.2f})')
            ax.set_xlabel(cfg.get('var_label'))
            ax.grid(False)
            if i == 0:
                ax.legend()

        # Remove empty panels
        for j in range(n_models, nrows * ncols):
            row, col = divmod(j, ncols)
            fig.delaxes(axes[row][col])

        caption = f'{exp_variable} Perkins Skill Score in {region} for {obs_name}'
        fig.suptitle(caption)

        prov_dic = create_provenance(caption)
        plt.tight_layout()

        fig_path = os.path.join(cfg['plot_dir'], f'pss_{variable}_{region}_{obs_name}')

        save_figure(fig_path, prov_dic, cfg, fig, close=True)        

    return


def main(cfg):

    input_data = cfg['input_data']

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)

    reference_dic = obtain_reference(groups.pop('reference'))

    remaining_metadata = []
    for k in groups.keys():
        remaining_metadata.extend(groups[k])

    data_dic = {}

    datasets = group_metadata(remaining_metadata, 'dataset')
    for dataset in datasets.keys():
        filepaths = list(group_metadata(datasets[dataset], 'filename').keys())
        mod_var_list = []
        for filepath in filepaths:
            mod_cb = iris.load_cube(filepath)
            mod_var_list.append(mod_cb.data)
        full_ens_data = np.asarray(mod_var_list).flatten()
        mod_pss_dic = perkins_skill_score(full_ens_data, reference_dic, cfg.get('bin_width'))
        data_dic[dataset] = mod_pss_dic

    plot_skills(data_dic, list(reference_dic.keys()), cfg)
    save_skills(data_dic, cfg)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
