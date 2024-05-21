"""diagnostic script to plot minima and maxima trends based on code
    from Anton's COSIMA cookbook notebook

"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
import os
import logging
from esmvaltool.diag_scripts.shared import run_diagnostic, save_figure
from esmvaltool.diag_scripts.shared._base import get_plot_filename


# This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def sea_ice_area(sic,area,coordls, range=[0.15,1]): 
    #sic percent is in 0-100. Mulitply by portion so divide by 100 ##
    sic = sic/100 
    return (sic*area).where((sic>=range[0])*(sic<=range[1])).sum(coordls)

def sea_ice_area_obs(ds):
    sic = ds.siconc #
    area_km2 = ds.areacello/1e6 #
    result = sea_ice_area(sic,area_km2,['x','y']).to_dataset(name='cdr_area') 

    #Theres a couple of data gaps which should be nan
    result.loc[{'time':'1988-01-01'}] = np.nan
    result.loc[{'time':'1987-12'}] = np.nan

    return result.sel(time=slice('1979','2018')) 

def sea_ice_area_model_sh(ds):
    sic = ds.siconc.where(ds.siconc.lat < -20, drop=True) #
    
    area_km2=ds.areacello/1e6 ##area convert to km2

    return sea_ice_area(sic,area_km2,['i','j']).to_dataset(name='si_area')

def min_and_max(ds):
    def min_and_max_year(da):
        result = xr.Dataset()
        result['min'] = da.min()
        result['max'] = da.max()
        return result
    annual_min_max_ds=ds.si_area.groupby('time.year').apply(min_and_max_year)
    return annual_min_max_ds

def plot_trend(model_min_max, obs_a, minmax):
    """
    model_min_max: dictionary of model label and xarray ds with min and max.
    obs_a: xarray of observations area
    """

    figure, _axes = plt.subplots()
    
    # both min and max # multiple models? min_max dict 
    # add note for years change?
    for mod_label,model_min_max_dt in model_min_max.items():
        model_min_max_dt[minmax].plot(label=mod_label)

    if minmax == 'max':
        obs_a.cdr_area.groupby('time.year').max().plot(label='Obs CDR')
        plt.title('Trends in Sea-Ice Maxima')
    elif minmax == 'min':
        obs_a.cdr_area.groupby('time.year').min().plot(label='Obs CDR')
        plt.title('Trends in Sea-Ice Minima')

    
    plt.ylabel('Sea-Ice Area (km2)')

    _ = plt.legend()
    return figure


def main(cfg):
    """Compute sea ice area for each input dataset."""

    input_data = cfg['input_data'].values()

    data = []

    for dataset in input_data:
        # Load the data
        input_file = [dataset['filename'], dataset['short_name'], dataset['dataset']]
        # key for different models      
        logger.info(f"dataset: {dataset['long_name']}")
        data.append(input_file)

    df = pd.DataFrame(data, columns=['filename','short_name','dataset']) #

    logger.info(df[['short_name', 'dataset']])

    min_max = {}
    # sort to ensure order of reading
    for fp, sn, dt in df.sort_values(['dataset','short_name']).itertuples(index=False):
        if dt == 'NSIDC-G02202-sh':    
            if sn == 'areacello':
                area_obs = xr.open_dataset(fp) 
            else:
                obs_si = xr.open_dataset(fp)
        else:  # other models 
            if sn == 'areacello':
                area_mod = xr.open_dataset(fp)
                dta = dt 
            else:
                mod_si = xr.open_dataset(fp)

                # make sure correct area with model?
                if dta == dt:
                    mod_si['areacello'] = area_mod['areacello'] #
                    model_area_dt = sea_ice_area_model_sh(mod_si)    
                    model_min_max_dt = min_and_max(model_area_dt)
                    model_min_max_dt['year'] = model_min_max_dt.year + 1652  # make years compariable
                    min_max[dt] = model_min_max_dt
                else:
                    logger.warning(f"..{dt} missing a variable?")
    
    obs_si['areacello'] = area_obs['areacello']
    obs_a = sea_ice_area_obs(obs_si)

    provenance_record = get_provenance_record(df['filename'].to_list())
    for m in ['min','max']:
        fig = plot_trend(min_max, obs_a, m)
        # Save output
        output_path = get_plot_filename(f'{m}_trend', cfg)
        # fig.savefig(output_path) # use esmvaltool convenience function
        save_figure(output_path, provenance_record, cfg, figure=fig)


def get_provenance_record(ancestor_files):
    record = {
        'ancestors': ancestor_files,
        'authors': [
            'chun_felicity',
        ],
        'caption': 'added 1652 years to model years for comparability',
        'domains': ['shpolar'],
        'plot_types': ['times'],
        'references': [],
        'statistics': ['mean'],
        }
    return record

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
