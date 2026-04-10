import iris
import iris.cube
import logging
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, group_metadata, save_figure
import esmvaltool.diag_scripts.shared.plot as eplot
# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def bootstrap_stdev(data: np.ndarray, block_size: int, n_bootstrap: int = 2500):

    nens, nyears = data.shape
    n_blocks = int(np.ceil(nyears/block_size))

    boot_list = []

    # Create overlapping blocks
    blocks = np.array([data[:,i:i+block_size]
               for i in range(nyears - block_size + 1)])
    
    rng = np.random.default_rng(501)
    
    for _ in range(n_bootstrap):
        # Sample blocks with replacement
        sampled_blocks = blocks[rng.integers(0, len(blocks), size=n_blocks)]
        sampled_blocks = sampled_blocks.transpose(1, 0, 2).reshape(nens, n_blocks*block_size)
        # Flatten the sampled blocks and calculate the standard deviation
        boot_sample = sampled_blocks[:, :nyears]
        boot_list.append(boot_sample.std(axis=1, ddof=1).mean())

    return np.asarray(boot_list)

def obtain_reference(data_group: list, block_size: int = 10):
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
        ref_data = ref_cb.data.reshape(-1, ref_cb.data.shape[0])
        reference_dic[key] = {}
        reference_dic[key]['best'] = ref_data.std(ddof=1)
        reference_dic[key]['bootstrap'] = bootstrap_stdev(ref_data, block_size)

    return reference_dic


def create_provenance(caption: str):
    '''Creates provenance dictionary'''

    provenance_dic = {'authors': ['malinina_elizaveta'], 
                      'caption': caption,
                      'references': ['malinina24']}

    return provenance_dic


def save_ratios(data_dic, reference_dic, cfg):

    percentiles = cfg.get('percentiles', [5, 95])
    region = cfg.get('region', 'region')
    variable = cfg.get('var_label')

    ratios_info = []

    for ref_dataset in reference_dic.keys():
        ref_boot_values = reference_dic[ref_dataset]['bootstrap']
        self_ratios = ref_boot_values[:, None] / ref_boot_values[None, :]
        # remove self-division (i == j)
        mask = ~np.eye(len(ref_boot_values), dtype=bool)
        self_ratios = self_ratios[mask]
        ref_perc_vals = np.nanpercentile(self_ratios, percentiles).tolist()
        dataset_dic = {'region': region, 'dataset': ref_dataset, 'reference': ref_dataset, 
                       'range_start': ref_perc_vals[0], 'range_end': ref_perc_vals[1]}
        ratios_info.append(dataset_dic)
        for model in data_dic.keys():
            mod_vals = data_dic[model]['bootstrap']/reference_dic[ref_dataset]['bootstrap']
            perc_vals = np.nanpercentile(mod_vals, percentiles).tolist()
            dataset_dic = {'region': region,'dataset': model, 'reference': ref_dataset, 
                       'range_start': perc_vals[0], 'range_end': perc_vals[1]}
            ratios_info.append(dataset_dic)

    ratios_df = pd.DataFrame(ratios_info)

    ratios_df.to_csv(os.path.join(cfg['work_dir'], 
                    f'variability_ratio_bootstrap_{variable}_{region}.csv'), index=False)
   
    return


def plot_stdevs(data_dic, reference_dic, cfg):
    '''
    This function plots timeseries and the max/min values

    Parameters:
    -----------
    data_dic:
        dictionary with the data to plot
    reference_dic:
        dictionary with the reference data
    cfg:
        config dictionary, comes from ESMValCore
    '''

    mpl_st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(mpl_st_file)
    
    percentiles = cfg.get('percentiles', [5, 95])

    fig_stds, ax_stds = plt.subplots(nrows=1, ncols=1)

    fig_stds.set_size_inches(8., 9.)
    fig_stds.set_dpi(200)

    y_ticks = np.arange(0, len(data_dic.keys()))
    y_labs = np.zeros(len(data_dic.keys()), dtype='<U30')
    for nr, ref_dataset in enumerate(reference_dic.keys()):
        ref_color_st = eplot.get_dataset_style(ref_dataset, cfg.get('color_style'))
        for nm, model in enumerate(data_dic.keys()):
            mod_vals = data_dic[model]['bootstrap']/reference_dic[ref_dataset]['bootstrap']
            single_dot = ax_stds.scatter(data_dic[model]['best']/reference_dic[ref_dataset]['best'], 
                                         nm-0.1 + nr*0.2, marker='o', c=ref_color_st['color'], zorder=2, 
                                         label=f"model/{ref_dataset}")
            perc_vals = np.percentile(mod_vals, percentiles)
            ax_stds.hlines(nm-0.1 + nr*0.2, perc_vals[0], perc_vals[1], 
                           colors=ref_color_st['color'], zorder=2)
            y_labs[nm] = model
    
    legend_handles = [single_dot]

    for ref_dataset in reference_dic.keys():
        ref_color_st = eplot.get_dataset_style(ref_dataset, cfg.get('color_style'))
        ref_boot_values = reference_dic[ref_dataset]['bootstrap']
        self_ratios = ref_boot_values[:, None] / ref_boot_values[None, :]
        # remove self-division (i == j)
        mask = ~np.eye(len(ref_boot_values), dtype=bool)
        self_ratios = self_ratios[mask]
        # remove inf / nan (e.g., division by zero)
        self_ratios = self_ratios[np.isfinite(self_ratios)]
        ref_perc_vals = np.percentile(self_ratios, percentiles)
        ax_stds.fill_betweenx([len(data_dic.keys()) -0.5, -0.5], 
                                 ref_perc_vals[0], ref_perc_vals[1],
                                 alpha=0.15, color=ref_color_st['color'], lw=0)
   
    ax_stds.set_ylim(len(data_dic.keys()) -0.5, -0.5)

    ax_stds.set_yticks(y_ticks, labels=y_labs)
    ax_stds.grid(which='both', c='silver', zorder=1)

    ax_stds.set_xscale('log')

    variable = cfg.get('var_label')
    exp_variable = variable.replace('_', ' ')
    units = cfg.get('units')
    region = cfg.get('region', 'region')

    ax_stds.set_xlabel(f'StD of {exp_variable}, {units}')

    default_caption = f'{variable} variability ratio in {region}'

    caption = cfg['figure_caption'] if cfg.get('figure_caption') else default_caption
    fig_stds.suptitle(caption)

    prov_dic = create_provenance(caption)

    plt.legend(handles=legend_handles, bbox_to_anchor =(0.5,-0.1), 
               loc='lower center', ncols=4, fancybox=False, frameon=False)
    
    plt.tight_layout()

    fig_path = os.path.join(cfg['plot_dir'], f'variability_ratio_bootstrap_{variable}_{region}')

    save_figure(fig_path, prov_dic, cfg, fig_stds, close=True)        

    return


def main(cfg):

    input_data = cfg['input_data']

    block_size = cfg.get('block_size', 10)

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)

    reference_dic = obtain_reference(groups.pop('reference'))

    remaining_metadata = []
    for k in groups.keys():
        remaining_metadata.extend(groups[k])

    data_dic = {}

    datasets = group_metadata(remaining_metadata, 'dataset')
    ens_var_cubelist = iris.cube.CubeList()
    for dataset in datasets.keys():
        data_dic[dataset] ={}
        filepaths = list(group_metadata(datasets[dataset], 'filename').keys())
        mod_var_list = []
        for filepath in filepaths:
            mod_cb = iris.load_cube(filepath)
            ens_var_cubelist.append(mod_cb)
            mod_var_list.append(mod_cb.data)
        mod_data = np.asarray(mod_var_list)
        if mod_data.ndim == 1:
            mod_data = mod_data.reshape(-1, mod_data.shape[0])
        data_dic[dataset]['best'] = mod_data.std(axis=1, ddof=1).mean()
        data_dic[dataset]['bootstrap'] = bootstrap_stdev(mod_data, block_size)

    save_ratios(data_dic, reference_dic, cfg)
    plot_stdevs(data_dic, reference_dic, cfg)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
