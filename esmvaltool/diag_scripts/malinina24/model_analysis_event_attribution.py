'''
This diagnostic calculates return periods and risk ratios for an
extreme event of interest from models, thus providing actual
attribution.

This diagnostic requires an ancestor task to determine the event and
its rarity from the observations. 

The preprocessed data should be provided in form of spatially averaged
timeseries of the annual maxima of a variable (e.g., tasmax). While
some seasonal data can be used at the preprocessor stage, to the 
diagnostic only one value per year shall be provided. Example: 
extract JJA season first, calculate annual mean using JJA data only.

This code was developped for the results used in Malinina&Gillett(2024)
and Gillet et al. (2022) (Weather and Climate Extremes). 

Author: Elizaveta Malinina (elizaveta.malinina-rieger@ec.gc.ca)
Initial development: 2021-2022
Last updated: August 2024
'''

import xarray as xr
import climextremes as cex
import logging
import numpy as np
import matplotlib.pyplot as plt
import os
from datetime import datetime, timedelta
from scipy.stats import genextreme as gev
import yaml

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_figure, io
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger

# import the classes from OBS diag 
from esmvaltool.diag_scripts.malinina24.observational_return_periods import StationaryRP, select_bins

logger = logging.getLogger(os.path.basename(__file__))

def calculate_anomaly(data: np.array, input_data: list, file_meta: dict,
                                                            anomaly_type: str):
    '''
    Parameters
    ----------
    data : 
        a numpy array for whihc anomalies to be calculated
    input_data : 
        an ESMValTool metadata dictionary with all the available data
    file_meta : 
        a dictionary with the data metadata
    anomaly_type: 
        a type of anomaly to calculate
    
    Returns
    -------
    anomaly : 
        numpy array with the anomaly data
    
    Raises
    ------
    NotImplementedError
        if the anomly type is not 'anomaly' or 'ratio' 
    '''
    ano_file = select_metadata(input_data, dataset=file_meta['dataset'],
                               ensemble=file_meta['ensemble'], 
                               variable_group='anomaly')[0]['filename']
    ano_data = xr.open_dataset(ano_file)[file_meta['short_name']].data
    if anomaly_type == 'anomaly':
        anomaly = data - ano_data
    elif anomaly_type == 'ratio':
        anomaly = data/ano_data
    else:
        raise NotImplementedError(f"The anomaly type {anomaly_type} isn't "
                                  f"impleneted. Possible entries: 'anomaly' "
                                  "and 'ratio'.")
    logger.info(f"The anomaly of {anomaly_type} for " +\
                file_meta['dataset'] + " has been calculated")
    
    return anomaly


def calculate_initial_cond(bootstrap: list[StationaryRP]):
    '''
    Parameters
    ----------
    bootstrap: 
        list with the bootstrap samples with the fitted GEV
    
    Returns
    -------
    initial_cond : dict
        dictionary with initial conditions
    '''        
    shape = float(np.mean([b.shape for b in bootstrap] ))
    scale = float(np.mean([b.scale for b in bootstrap]))
    loc = float(np.mean([b.loc for b in bootstrap]))
    initial_cond = {'location': loc, 'scale': scale, 'shape': shape}

    return initial_cond


def determine_weights(datasets: np.ndarray | list, shape: tuple):
    '''
    Parameters
    ----------
    datasets :
        array or list with the datasets names
    shape :
        tuple with the target shape of the weights array
    
    Returns
    -------
    weights :
        weights with the shape of shape
    '''
    # weighs are assigned in a way that each model in a sample gets
    # a weight of 1, e.g., model has 4 samples, it's weight is 0.25
    weights_1d = [1/np.count_nonzero(datasets == d) for d in datasets]
    weights = np.repeat(np.expand_dims(weights_1d, 1),shape[1]
                                                        ).reshape(shape)

    return weights


def bootstrap_gev(data: np.ndarray, event: float, datasets: np.ndarray, 
                                                            cfg: dict = {}):
    '''
    Parameters
    ----------
    data : 
        data which will be used for bootstrap
    event :
        event for which RPs to be calculated
    datasets : 
        array with the dataset names
    cfg :
        internal ESMValTool dict with parameters for the diagnostic
    '''
    seed = cfg['bootstrap_seed'] if cfg.get('bootstrap_seed') else 1
    yblock = cfg['yblock'] if cfg.get('yblock') else 1

    pool_iter_n = np.max([len(np.unique(datasets)) *100, 1000])

    # the number of the datasets to be included in the pool
    # in case the dataset number is less than 3, using 3 as
    # the minimum ensemble number of to draw conclusions on
    # 3 most likely to be used in the case of only one model
    n_dataset = np.max([len(np.unique(datasets)), 3])

    n_years = data.shape[1]

    pool_size = int(np.around(n_dataset*n_years/yblock))

    # defining the seed for the random sequence
    rng = np.random.default_rng(seed)

    bootstrap_samples = list()

    for i in range(pool_iter_n):
        sample_data = list() ; sample_datasets = list()
        mod_idxs = rng.integers(0, high=n_dataset, size=pool_size)
        for mod_idx in mod_idxs:
            y_idx = rng.integers(0, high=n_years-yblock)
            sample_data.append(data[mod_idx, y_idx:y_idx+yblock])
            sample_datasets.append(datasets[mod_idx])
        sample_data = np.asarray(sample_data)
        if cfg['model_weighting']:
            sample_weights = determine_weights(sample_datasets, 
                                            sample_data.shape).flatten()
        else:
            sample_weights = None
        SampleGEV = StationaryRP(sample_data.flatten(), event, 
                                                        weights=sample_weights)
        bootstrap_samples.append(SampleGEV)

    return bootstrap_samples

def get_basic_event_info(obs_info : dict, cfg: dict):
    '''
    Parameters
    ----------
    obs_info : 
        dictionary with the observational information
    cfg : 
        internal ESMValTool dictionary with the input from the recipe

    Returns
    -------
    event : float
        event strength from observations
    obs_event_rp : float
        best guess event return period from observation
    '''

    obs_ref_dataset = cfg['obs_ref_dataset']
    event = obs_info[obs_ref_dataset]['event']
    if obs_info.get('obs_gev'): 
        obs_gev_type = obs_info['obs_gev']
    else:
        obs_gev_type = 'stationary'
    obs_event_rp = obs_info[obs_ref_dataset][obs_gev_type+'_gev']['RP']

    return event, obs_event_rp


def calculate_intensity(GEV : StationaryRP, x_fine : np.ndarray, 
                                            event_rp: float):
    '''
    Calculate intensity of the event if the event_rp rarity

    Parameters
    ----------
    GEV : 
        class with the stationary GEV info
    x_fine : 
        array with the strengths
    event_rp :
        return period of the event for which inetensity is calculated
    
    Returns
    -------
    intensity : float 
        the strength of the event of the event_rp probability
    '''

    rp_func = 1/gev.sf(x_fine, GEV.shape, loc=GEV.loc, scale=GEV.scale)
    intensity = float(x_fine[np.argmin(np.abs(rp_func-event_rp))])

    return intensity


class Climate:
    '''
    Attributes
    ----------
    name
    start_year
    end_year
    data
    weights
    bootstrap
    BestGuessGEV
    RP_CI
    intensity
    intensity_CI 
    '''

    def __init__(self, input_data: list, group: str, anomaly_calculation: bool,
                                                    cfg: dict, event: float):
        
        group_info = select_metadata(input_data, variable_group=group)

        self.name = group
        self.start_year = list(group_metadata(group_info, 'start_year'
                                                                  ).keys())[0]
        self.end_year = list(group_metadata(group_info, 'end_year').keys())[0]

        # obtain data and datasets to calculate weights
        files = group_metadata(group_info, 'filename')
        group_data = list(); datasets = list()
        for file in files.keys():
            file_meta = select_metadata(input_data, filename=file)[0]
            short_name = file_meta['short_name']
            file_dataset = file_meta['dataset']
            data = xr.open_dataset(file)[short_name].data
            if anomaly_calculation:
                data = calculate_anomaly(data, input_data, file_meta, 
                                            cfg['anomaly_type'])
            group_data.append(data)
            datasets.append(file_dataset)
        group_data = np.asarray(group_data) ; datasets = np.asarray(datasets)

        self.data = group_data
        if cfg['model_weighting']:
            self.weights = determine_weights(datasets, group_data.shape).flatten()
        else:
            self.weights = None

        self.bootstrap = bootstrap_gev(self.data, event, datasets, cfg)
        if cfg['initial_conditions']: 
            initial_cond = calculate_initial_cond(self.bootstrap)
        else:
            initial_cond = None
        self.BestGuessGEV = StationaryRP(self.data.flatten(), event,
                                                weights=self.weights,
                                                initial=initial_cond)
        self.rp_ci = {}
        # calculate RPs for the quantiles
        cfg['CI_quantiles']

    def calculate_intensities(self, x_fine, event_rp, CI):

        self.intensity = calculate_intensity(self.BestGuessGEV, x_fine, 
                                                                event_rp)
        ci_dict = {} 
        self.intensity_CI = ci_dict
        # calculate CIs 

    def update_event(self, event):
        self.BestGuessGEV.rp = float(1/gev.sf(event, self.BestGuessGEV.shape, 
                                               loc=self.BestGuessGEV.loc,
                                               scale=self.BestGuessGEV.scale))
        for i in range(len(self.bootstrap)):
            self.bootstrap[i].rp = float(1/gev.sf(event, 
                                                self.bootstrap[i].shape,
                                                loc=self.bootstrap[i].loc,
                                                scale=self.bootstrap[i].scale))

    def reform_to_yml(self):

        out_dic = {self.name: {'shape': self.BestGuessGEV.shape,
                               'scale': self.BestGuessGEV.scale,
                               'loc': self.BestGuessGEV.loc,
                               'RP': self.BestGuessGEV.rp
                            #    add CIs
                            # add intensities and their CIs
                               }}

        return out_dic


def calculate_risk_ratios(clim_list : list[Climate], ci_percs : list[float|int]):
    '''
    Calculate risk ratios from the list of climates

    Parameters
    ----------
    clim_list : 
        list with Climate information
    ci_dic : 
        list with the percentiles
    
    Returns
    -------
    risk_ratio_dic : dict
        dictionary with risk ratios

    Raises
    ------
    ValueError
        if the amount of climates is insuffiecient to calculate risk
        ratios (the number of climates is less than two)
    '''

    if len(clim_list)<2: 
        raise ValueError("The amount of climates in insufficient to calculate "
                         "a risk ratio. Required number of climates is 2, but "
                         f"only {len(clim_list)} were provided.")

    # sorting the climates from earliest to latest 
    clim_indices = np.argsort([c.start_year for c in clim_list])

    risk_ratio_dic = {}

    # starting with the second, because the oldest is always the reference 
    for clim_idx in clim_indices[1:]:
    # for climate in clim_list[clim_indices[1:]]:
        current_idx = np.where(clim_indices==clim_idx)[0][0]
        for dev_clim_idx in range(current_idx):
            old_clim_idx = clim_indices[dev_clim_idx]
            key = f"{clim_list[clim_idx].name}/{clim_list[old_clim_idx].name}"
            risk_ratio_dic[key] = {
                'best_guess': float(clim_list[old_clim_idx].BestGuessGEV.rp/\
                    clim_list[clim_idx].BestGuessGEV.rp)}
            bootstrap_rrs = list()
            # the number of bootstrap samples is the same in each climate
            for i in range(len(clim_list[clim_idx].bootstrap)):
                bootstrap_rrs = clim_list[old_clim_idx].bootstrap[i].rp/\
                                    clim_list[clim_idx].bootstrap[i].rp
            risk_ratio_dic[key]['CI'] = {}
            for ci_perc in ci_percs: 
                risk_ratio_dic[key]['CI'][ci_perc] = float(np.nanpercentile(
                        bootstrap_rrs, ci_perc, method='closest_observation'
                                                                    ).round(2))

    return risk_ratio_dic


class Climates:
    '''
    Class with climates information for a dataset

    Attributes
    ----------
    name : str
        name of the dataset
    climates : list
        list with the Climates classes
    risk_ratios
    max_value
    min_value
    '''

    def __init__(self, input_data: list, dataset_name: str, cfg: dict, 
                                                            obs_info: dict):
        
        event, event_rp = get_basic_event_info(obs_info, cfg)

        self.name = dataset_name
        self.climates = []

        groups = list(group_metadata(input_data, 'variable_group').keys())
        if 'anomaly' in groups:
            groups.remove('anomaly'); anomaly_calculation = True

        for group in groups:
            GroupClimate = Climate(input_data, group, anomaly_calculation, 
                                                                    cfg, event)
            self.climates.append(GroupClimate)

        max_value = np.max([cl.data.max() for cl in self.climates])
        min_value = np.min([cl.data.min() for cl in self.climates])

        self.max_value = np.ceil(max_value*1.1 if max_value> 0
                                               else max_value*0.9)
        self.min_value = np.floor(min_value*0.9 if min_value> 0 
                                               else min_value*1.1)

        bins, x_fine = select_bins(self.min_value, self.max_value)
        
        if cfg.get('event_definition') == 'rarity':
            for Clim in self.climates:
                if Clim.name == 'factual':
                    FactClim = Clim
                    # the event is redifined using the best guess event 
                    # probability obtained from observations
                    event = calculate_intensity(FactClim.BestGuessGEV, x_fine,
                                                                     event_rp)
            for Clim in self.climates:
                Clim.update_event(event)
        for Clim in self.climates:
            Clim.calculate_intensities(x_fine, event_rp, cfg['CI_quantiles'])
        self.risk_ratios = calculate_risk_ratios(self.climates, cfg['CI_quantiles'])


    def plot_uncert_plot(self, cfg: dict, prov_dic: dict):

        return


    def plot_attribution_plot(self, cfg: dict, prov_dic: dict):

        return
    

    def save_risk_ratios(self, cfg: dict):

        return


def main(cfg):
    """Main function, essentially diagnostic itself."""

    # loading the observational info from ancestor task
    obs_file = io.get_ancestor_file(cfg, 'obs_information_event.yml')
    obs_info = yaml.safe_load(open(obs_file, 'r'))
    logger.info(f"Loaded information on OBS statistics from {obs_file}")

    provenance_dic = {'authors': ['malinina_elizaveta']}
    # add ancestors to provenance

    input_data = cfg['input_data']

    datasets = group_metadata(input_data.values(), 'dataset')

    output_dic = {}

    if cfg['individual_models']:
    # provide statistics for individual models
        for dataset in datasets.keys():
            dataset_info = datasets[dataset]
            logger.info(f"Processing {dataset}")
            DatasetClimates = Climates(dataset_info, dataset, cfg, obs_info)
            for Climate in DatasetClimates.climates:
                output_dic[dataset] = Climate.reform_to_yml()
            output_dic[dataset]['risk_ratios'] = DatasetClimates.save_risk_ratios(
                                                            cfg['CI_quantiles'])

    logger.info(f"Processing Multi-Model-Ensemble")
    MMEClimates = Climates(input_data.values(), 'Multi-Model-Ensemble', cfg, obs_info)    
    MMEClimates.plot_attribution_plot(cfg, provenance_dic)
    MMEClimates.plot_uncert_plot(cfg, provenance_dic)
    for MMEClimate in MMEClimates.climates:
        output_dic['Multi-Model-Ensemble'] = MMEClimate.reform_to_yml()
    output_dic['Multi-Model-Ensemble'] = MMEClimates.save_risk_ratios(cfg['CI_quantiles'])

    model_stats_path = os.path.join(cfg['work_dir'],'models_statistics.yml')
    with open(model_stats_path, 'w') as model_info_yml:
        yaml.dump(output_dic, model_info_yml, sort_keys=False)
    logger.info(f"Saved OBS info into {model_stats_path}")

    logger.info("Successfully completed diagnostic")


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)