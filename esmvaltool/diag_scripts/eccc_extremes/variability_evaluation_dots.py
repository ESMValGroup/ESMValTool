import iris
import iris.cube
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
    
    fig_stds, ax_stds = plt.subplots(nrows=1, ncols=1)

    fig_stds.set_size_inches(8., 9.)
    fig_stds.set_dpi(200)

    y_ticks = np.arange(0, len(data_dic.keys()))
    y_labs = np.zeros(len(data_dic.keys()), dtype='<U30')

    for nm, model in enumerate(data_dic.keys()):
        color_st = eplot.get_dataset_style(model, cfg.get('color_style'))
        for i in range(0, len(data_dic[model])):
            single_dot = ax_stds.scatter(data_dic[model][i], nm, marker='o',
                                c=color_st['color'], zorder=2, 
                                label='individual member')
        y_labs[nm] = model

    legend_handles = [single_dot]

    for ref_dataset in reference_dic.keys():
        ref_color_st = eplot.get_dataset_style(ref_dataset, cfg.get('color_style'))
        ref_line = ax_stds.axvline(reference_dic[ref_dataset], c=ref_color_st['color'], zorder=2, 
                                                        label=ref_dataset)
        legend_handles.append(ref_line)
    ax_stds.set_ylim(len(data_dic.keys()) -0.8, -0.2)

    ax_stds.set_yticks(y_ticks, labels=y_labs)
    ax_stds.grid(which='both', c='silver', zorder=1)

    variable = cfg.get('var_label')
    exp_variable = variable.replace('_', ' ')
    units = cfg.get('units')
    region = cfg['region'] if cfg.get('units') else 'region'

    ax_stds.set_xlabel(f'StD of {exp_variable}, {units}')

    default_caption = f'{variable} variability in {region}'

    caption = cfg['figure_caption'] if cfg.get('figure_caption') else default_caption
    fig_stds.suptitle(caption)

    prov_dic = create_provenance(caption)

    plt.legend(handles=legend_handles, bbox_to_anchor =(0.5,-0.1), 
               loc='lower center', ncols=4, fancybox=False, frameon=False)
    
    plt.tight_layout()

    fig_path = os.path.join(cfg['plot_dir'], f'variability_{variable}_{region}')

    save_figure(fig_path, prov_dic, cfg, fig_stds, close=True)        

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
    ens_var_cubelist = iris.cube.CubeList()
    for dataset in datasets.keys():
        filepaths = list(group_metadata(datasets[dataset], 'filename').keys())
        mod_var_list = []
        for filepath in filepaths:
            mod_cb = iris.load_cube(filepath)
            ens_var_cubelist.append(mod_cb)
            mod_var_list.append(mod_cb.data)
        data_dic[dataset] = mod_var_list

    plot_stdevs(data_dic, reference_dic, cfg)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
