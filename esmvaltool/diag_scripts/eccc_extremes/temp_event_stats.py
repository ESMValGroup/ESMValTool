import xarray as xr
import calendar
import operator
import logging
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
from datetime import datetime, timedelta, date

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, group_metadata, save_figure
import esmvaltool.diag_scripts.shared.plot as eplot
# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def get_deseasonalized_values(da: xr.DataArray, smooth: int):

    da = da.sortby("time")
    rolling_da = da.rolling(time=smooth, center=True).mean()
    smoothed_clim = rolling_da.groupby("time.dayofyear").mean("time")
    deseasonal_da = da.groupby("time.dayofyear") - smoothed_clim

    return deseasonal_da


def compute_perc_values(da: xr.DataArray,  percentile: float,
                            half_window: int):

    # remove leap days from the data array
    da = da.sel(
        time=~((da.time.dt.month == 2) & (da.time.dt.day == 29))
    )
    
    percentiles = []

    for doy in range(1, 366):
        days = ((np.arange(doy-half_window, doy+half_window+1) - 1) % 365) + 1

        window_data = da.where(
            da.time.dt.dayofyear.isin(days),
            drop=True,
        )

        q = window_data.quantile(percentile/100, dim="time")
        percentiles.append(q)

    clim_percentiles = xr.concat(percentiles, dim="dayofyear")
    clim_percentiles = clim_percentiles.assign_coords(dayofyear=np.arange(1, 366))

    return clim_percentiles

def calculate_events(da: xr.DataArray, min_length: int):

    changes = da.astype(int).diff(dim='time').convert_calendar("standard", use_cftime=False)

    starts = changes.where(changes == 1).dropna(dim='time').time.dt.date.values
    ends = changes.where(changes == -1).dropna(dim='time').time.dt.date.values
    
    events = pd.DataFrame({'start': [datetime.strptime(str(s), '%Y-%m-%d') for s in starts], 
                           'end': [datetime.strptime(str(e), '%Y-%m-%d') - timedelta(days=1) for e in ends]})

    events['duration'] = (events['end'] - events['start']).dt.days + 1

    events = events[events['duration'] >= min_length]

    events['mean T'] = None
    events['peak T'] = None
    events['mean T exceedance'] = None
    events['peak T exceedance'] = None

    return events


def calculate_event_strengths(events: pd.DataFrame, da: xr.DataArray, deseasonal_da: xr.DataArray, opeartor_str: str):
    '''
    Calculate the temperature stats from the event dates
    '''

    extremes_operator = 'max' if opeartor_str == 'gt' else 'min'


    da = da.convert_calendar("standard", use_cftime=False)
    deseasonal_da = deseasonal_da.convert_calendar("standard", use_cftime=False)


    for index, row in events.iterrows():
        start_date = str(row['start'].date())
        end_date = str(row['end'].date())
        event_data = da.sel(time=slice(start_date, end_date))
        event_deseasonal_data = deseasonal_da.sel(time=slice(start_date, end_date))

        mean_T = event_data.mean().item()
        peak_T = getattr(event_data, extremes_operator)().item()
        mean_T_exceedance = event_deseasonal_data.mean().item()
        peak_T_exceedance = getattr(event_deseasonal_data, extremes_operator)().item()

        events.at[index, 'mean T'] = mean_T
        events.at[index, 'peak T'] = peak_T
        events.at[index, 'mean T exceedance'] = mean_T_exceedance
        events.at[index, 'peak T exceedance'] = peak_T_exceedance

    return events


def create_provenance(caption: str):
    '''Creates provenance dictionary'''

    provenance_dic = {'authors': ['malinina_elizaveta'], 
                      'caption': caption,
                      'references': ['malinina24']}

    return provenance_dic


def plot_stats(data_list, ref_list, cfg, label='full'):
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
    
    fig_means, ax_means = plt.subplots(nrows=6, ncols=1, sharex=True)

    fig_means.set_size_inches(6., 10.)
    fig_means.set_dpi(200)

    x_ticks = np.arange(0, len(ref_list)+len(data_list)+1)
    x_labs = np.zeros(len(ref_list)+len(data_list)+1, dtype='<U30')

    for nr, ref_data in enumerate(ref_list):
        ref_color_st = eplot.get_dataset_style(ref_data['dataset'], cfg.get('color_style'))
        x_labs[nr] = ref_data['alias']
        ax_means[0].scatter(nr, len(ref_data['events']), c=ref_color_st['color'], zorder=2)
        # add longest event
        ax_means[1].scatter(nr, ref_data['events']['duration'].mean().item(), c=ref_color_st['color'], zorder=2)
        ax_means[2].scatter(nr, ref_data['events']['mean T'].mean().item(), c=ref_color_st['color'], zorder=2)
        ax_means[3].scatter(nr, ref_data['events']['peak T'].mean().item(), c=ref_color_st['color'], zorder=2)
        ax_means[4].scatter(nr, ref_data['events']['mean T exceedance'].mean().item(), c=ref_color_st['color'], zorder=2)
        ax_means[5].scatter(nr, ref_data['events']['peak T exceedance'].mean().item(), c=ref_color_st['color'], zorder=2)

    num, mean_dur, mean_T, peak_T, mean_dev, peak_dev = [], [], [], [], [], [] 

    for nd, data in enumerate(data_list):
        x = len(ref_list)+nd
        color_st = eplot.get_dataset_style(data['dataset'], cfg.get('color_style'))
        x_labs[x] = data['alias']
        ax_means[0].scatter(x, len(data['events']), c=color_st['color'])
        num.append(len(data['events']))
        # add longest event
        ax_means[1].scatter(x, data['events']['duration'].mean().item(), c=color_st['color'], zorder=2)
        mean_dur.append(data['events']['duration'].mean().item())
        ax_means[2].scatter(x, data['events']['mean T'].mean().item(), c=color_st['color'], zorder=2)
        mean_T.append(data['events']['mean T'].mean().item())
        ax_means[3].scatter(x, data['events']['peak T'].mean().item(), c=color_st['color'], zorder=2)
        peak_T.append(data['events']['peak T'].mean().item())
        ax_means[4].scatter(x, data['events']['mean T exceedance'].mean().item(), c=color_st['color'], zorder=2)
        mean_dev.append(data['events']['mean T exceedance'].mean().item())
        ax_means[5].scatter(x, data['events']['peak T exceedance'].mean().item(), c=color_st['color'], zorder=2)
        peak_dev.append(data['events']['peak T exceedance'].mean().item())

    num = np.mean(num) ; mean_dur = np.mean(mean_dur); mean_T = np.mean(mean_T); 
    peak_T= np.mean(peak_T); mean_dev = np.mean(mean_dev); peak_dev = np.mean(peak_dev)

    ax_means[0].scatter(x+1, num, c=color_st['color'], s=100, zorder=2)
    ax_means[1].scatter(x+1, mean_dur, c=color_st['color'], s=100, zorder=2)
    ax_means[2].scatter(x+1, mean_T, c=color_st['color'], s=100, zorder=2)
    ax_means[3].scatter(x+1, peak_T, c=color_st['color'], s=100, zorder=2)
    ax_means[4].scatter(x+1, mean_dev, c=color_st['color'], s=100, zorder=2)
    ax_means[5].scatter(x+1, peak_dev, c=color_st['color'], s=100, zorder=2)

    x_labs[-1] = 'Model Mean'

    ax_means[0].set_title('Number of events')
    ax_means[1].set_title('Mean Duration')
    ax_means[2].set_title('Mean T')
    ax_means[3].set_title('Peak T')
    ax_means[4].set_title('Mean deviation')
    ax_means[5].set_title('Peak deviation')

    ax_means[-1].set_xticks(x_ticks, labels=x_labs, rotation=30)
    [a.grid(which='major', c='silver', zorder=0) for a in ax_means]

    region = cfg.get('region', 'region')

    default_caption = f"{cfg['event_type']} stats in {region} ({label})"

    caption = cfg['figure_caption'] if cfg.get('figure_caption') else default_caption
    fig_means.suptitle(caption)

    prov_dic = create_provenance(caption)
    
    plt.tight_layout()

    fig_path = os.path.join(cfg['plot_dir'], f"{cfg['event_type']}_stats_{label}_{region}")

    save_figure(fig_path, prov_dic, cfg, fig_means, close=True)        

    return


def main(cfg):

    input_data = cfg['input_data']

    smooth = cfg.get('smooth', 5)
    percentile = cfg.get('percentile', 90) 
    half_window = cfg.get('half_window', 15) 
    operator_str = cfg.get('operator', 'gt')
    func = getattr(operator, operator_str)
    min_length = cfg.get('min_length', 3)

    for f in input_data.keys():
        short_name = input_data[f]['short_name']
        dataset = input_data[f]['dataset']
        alias = input_data[f]['alias']
        da = xr.open_dataset(f)[short_name]
        da = da.assign_coords(dayofyear=da.time.dt.dayofyear)
        deseasonal_da = get_deseasonalized_values(da, smooth)
        percentile_values = compute_perc_values(deseasonal_da, percentile, half_window)
        result_mask = func(deseasonal_da.groupby('dayofyear'), percentile_values)
        events = calculate_events(result_mask, min_length)
        events = calculate_event_strengths(events, da, deseasonal_da, operator_str)
        input_data[f]['result_mask'] = result_mask
        input_data[f]['data'] = da
        input_data[f]['events'] = events
        alias = f'{alias}_{dataset}' if dataset not in alias else alias
        events.to_csv(os.path.join(cfg['work_dir'], f'events_{alias}.csv'))

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)

    reference_dic = groups.pop('reference')

    remaining_metadata = []
    for k in groups.keys():
        remaining_metadata.extend(groups[k])

    datasets = group_metadata(remaining_metadata, 'dataset')
    for dataset in datasets.keys():
        plot_stats(datasets[dataset], reference_dic, cfg)
        if len(cfg['months'])>1:
            short_label = f"{calendar.month_abbr[cfg['months'][0]]}-{calendar.month_abbr[cfg['months'][-1]]}"
        elif cfg['months'] == 1:
            short_label = f"{calendar.month_abbr[cfg['months'][0]]}-{calendar.month_abbr[cfg['months'][-1]]}"
        for data in datasets[dataset]:
            data['events'] = data['events'][
                data['events']['start'].dt.month.isin(cfg['months'])|
                data['events']['end'].dt.month.isin(cfg['months']) ]
        for ref_data in reference_dic:
            ref_data['events'] = ref_data['events'][
                ref_data['events']['start'].dt.month.isin(cfg['months'])|
                ref_data['events']['end'].dt.month.isin(cfg['months']) ]
        plot_stats(datasets[dataset], reference_dic, cfg, short_label)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
