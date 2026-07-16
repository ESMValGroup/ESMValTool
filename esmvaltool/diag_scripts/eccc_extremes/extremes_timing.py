import xarray as xr
import pandas as pd 
from datetime import date
import logging
import numpy as np
import matplotlib.pyplot as plt
import os

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, group_metadata, save_figure
import esmvaltool.diag_scripts.shared.plot as eplot
# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))


def get_one_extreme(data_group: list, operator_str: str):

    ref_fnames = group_metadata(data_group, 'filename')

    times = []

    for f in ref_fnames.keys():
        short_name = ref_fnames[f][0]['short_name']
        da = xr.open_dataset(f)[short_name]
        time_arr = da.resample(time="YE").map(lambda x: getattr(x, f"idx{operator_str}")(dim="time")).values
        ts = pd.Series(time_arr)
        times.append(ts)

    one_ts = pd.concat(times)

    return one_ts


def extremes_by_dataset(metadata_dic: dict, operator_str: str):

    datasets = group_metadata(metadata_dic, 'dataset', sort=True)
    data = {}

    for dataset, metadata in datasets.items():
        ts = get_one_extreme(metadata, operator_str)
        data[dataset] = ts

    return data


def create_provenance(caption: str):
    '''Creates provenance dictionary'''

    provenance_dic = {'authors': ['malinina_elizaveta'], 
                      'caption': caption,
                      'references': ['malinina24']}

    return provenance_dic


def plot_distribution(data_dic: dict, reference_dic: dict, cfg: dict):
    '''
    This function plots distributions of max/min

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

    n_ref = len(reference_dic.keys())
    n_mod = len(data_dic.keys())
                
    n_panels = n_ref + n_mod

    bins, bins_lit = [], []
    for m in range(1,13):
        mon_lit = date(2000, m, 1).strftime("%b")
        bins.append(date(2000, m, 1).timetuple().tm_yday)
        bins.append(date(2000, m, 16).timetuple().tm_yday)
        bins_lit.extend([f"1hlf {mon_lit}", f"2hlf {mon_lit}"])

    bins.append(366)

    freq = []

    fig_distrs, ax_distrs = plt.subplots(nrows=n_panels, ncols=1, sharex=True, sharey=True)

    fig_distrs.set_size_inches(6., 2*n_panels)
    fig_distrs.set_dpi(200)


    for nr, ref_dataset in enumerate(reference_dic.keys()):
        color_st = eplot.get_dataset_style(ref_dataset, cfg.get('color_style'))
        ref_data = reference_dic[ref_dataset]
        doys = ref_data.dt.dayofyear.to_list()
        f, _, _ = ax_distrs[nr].hist(doys, bins, density=True, color=color_st['color'], alpha=0.8)
        ax_distrs[nr].set_title(ref_dataset)
        freq.append(f)

    for nm, model in enumerate(data_dic.keys()):
        color_st = eplot.get_dataset_style(model, cfg.get('color_style'))
        data = data_dic[model]
        try:
            doys = data.dt.dayofyear.to_list()
        except:
            doys = [p.dayofyr for p in data_dic[model]]
        f, _, _ = ax_distrs[n_ref+nm].hist(doys, bins, density=True, color=color_st['color'], alpha=0.8)
        ax_distrs[n_ref+nm].set_title(model)
        freq.append(f)


    ax_distrs[-1].set_xticks(bins[:-1] + np.diff(bins)/2, bins_lit, rotation=30)
    freq = np.asarray(freq)
    idxs = np.flatnonzero(~np.all(freq == 0, axis=0))
    # max in mind only
    ax_distrs[-1].set_xlim(bins[idxs[0]], bins[idxs[-1]+1])

    region = cfg.get('region', 'region')
    caption = f"Timing of {cfg['operator']} in {region}"

    fig_distrs.suptitle(caption)

    prov_dic = create_provenance(caption)

    plt.tight_layout()

    fig_path = os.path.join(cfg['plot_dir'], f"timing_{cfg['operator']}_{region}")

    save_figure(fig_path, prov_dic, cfg, fig_distrs, close=True)        

    return


def main(cfg):

    input_data = cfg['input_data']

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)

    operator_str = cfg.get('operator', 'max')

    reference = groups.pop('reference')

    ref_data = extremes_by_dataset(reference, operator_str)

    remaining_metadata = []
    for k in groups.keys():
        remaining_metadata.extend(groups[k])

    model_data = extremes_by_dataset(remaining_metadata, operator_str)

    plot_distribution(model_data, ref_data, cfg)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
