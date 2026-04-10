import iris
import iris.cube
import yaml
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


def calculate_ensemble_stdev(ens_cubelist: iris.cube.CubeList, 
                                                    weight_flag: bool):
    '''
    This function calculates standard deviation of multi model ensemble

    Parameters:
    -----------
    ens_cubelist: 
        CubeList with the ensemble timeseries
    weight_flag:
        flag for weighting or not the data
    
    Returns:
    --------
    ens_std: float
        standard deviation for the ensemble
    '''

    wghts , wdata = [], []
    for i in range(0, len(ens_cubelist)):
        cb = ens_cubelist[i]
        cb_wght = cb.attributes['ensemble_weight']*cb.attributes['reverse_dtsts_n']
        wghts.append(np.full(cb.shape, cb_wght))
        wdata.append(ens_cubelist[i].data)

    wghts = np.asarray(wghts).flatten()
    wdata = np.asarray(wdata).flatten()

    if weight_flag:
        mmm = np.average(wdata, weights=wghts)
        ens_std = np.sqrt(np.average((wdata - mmm)**2, weights=wghts))
    else: 
        ens_std = wdata.std()

    return ens_std

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

    legend_handles = []

    data = []

    for nr, ref_dataset in enumerate(reference_dic.keys()):
        ref_color_st = eplot.get_dataset_style(ref_dataset, cfg.get('color_style'))
        for nm, model in enumerate(data_dic.keys()):
                for i in range(0, len(data_dic[model])):
                    data.append(reference_dic[ref_dataset]/data_dic[model][i])
                    single_dot = ax_stds.scatter(reference_dic[ref_dataset]/data_dic[model][i], nm-0.1 + nr*0.2, marker='o',
                                        c=ref_color_st['color'], zorder=2, 
                                        label=f"model/{ref_dataset}")
                y_labs[nm] = model
        legend_handles.append(single_dot)

    data = np.asarray(data)
    xmin, xmax = np.min(data), np.max(data)

    # Round limits to nearest 0.1
    start = np.floor(xmin * 10) / 10
    end   = np.ceil(xmax * 10) / 10

    # Generate ticks
    ticks = np.arange(start, end + 0.1, 0.1)

    ax_stds.set_xscale('log')

    # Apply ticks + labels
    ax_stds.set_xticks(ticks)
    ax_stds.set_xticklabels([f"{t:.1f}" for t in ticks])

    # Force plain formatting (avoid log-style labels)
    ax_stds.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax_stds.ticklabel_format(style='plain', axis='x')


    # ref_ratio_dic = {'comment': 'obs/mmm'}  # for the printout
    # for ref_dataset in reference_dic.keys():
    #     ref_color_st = eplot.get_dataset_style(ref_dataset, cfg.get('color_style'))
    #     ref_line = ax_stds.axvline(reference_dic[ref_dataset], c=ref_color_st['color'], zorder=2, 
    #                                                     label=ref_dataset)
    #     legend_handles.append(ref_line)
    #     # add the printout of the ratio
    #     ref_ratio_dic[ref_dataset] = float(reference_dic[ref_dataset] / data_dic['Multi-Model'])

    # with open(os.path.join(cfg['work_dir'], 'ratio_info.yml')) as ratio_f:
    #     yaml.safe_dump(ref_ratio_dic, ratio_f)
    #     # end of the printout


    ax_stds.set_ylim(len(data_dic.keys()) -0.5, -0.5)

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

    fig_path = os.path.join(cfg['plot_dir'], f'variability_ratio_{variable}_{region}')

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
