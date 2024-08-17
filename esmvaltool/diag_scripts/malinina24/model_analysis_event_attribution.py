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
from esmvaltool.diag_scripts.malinina24.observational_return_periods import StationaryRP

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
    '''

    def __init__(self, input_data: list, group: str, anomaly_calculation: bool,
                                                    cfg: dict, obs_info: dict):
        
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

        self.bootstrap = bootstrap_gev(self.data, obs_info['event'], datasets,
                                                                           cfg)
        if cfg['initial_conditions']: 
            initial_cond = calculate_initial_cond(self.bootstrap)
        else:
            initial_cond = None
        self.BestGuessGEV = StationaryRP(self.data.flatten(), 
                                            obs_info['event'],
                                            weights=self.weights,
                                            initial=initial_cond)


    def reform_to_yml(self):

        out_dic = {}

        return out_dic

class Climates:
    '''
    Class with climates information for a dataset

    Attributes
    ----------
    name : str
        name of the dataset
    climates : list
        list with the Climates classes 
    '''

    def __init__(self, input_data: list, dataset_name: str, cfg: dict, 
                                                            obs_info: dict):

        self.name = dataset_name
        self.climates = []

        groups = list(group_metadata(input_data, 'variable_group').keys())
        if 'anomaly' in groups:
            groups.remove('anomaly'); anomaly_calculation = True

        for group in groups:
            GroupClimate = Climate(input_data, group, anomaly_calculation, cfg, obs_info)
            self.climates.append(GroupClimate)
        self.risk_ratios = {}


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
            output_dic[dataset]['risk_ratios'] = DatasetClimates.save_risk_ratios(cfg)

    logger.info(f"Processing Multi-Model-Ensemble")
    MMEClimates = Climates(input_data.values(), 'Multi-Model-Ensemble', cfg, obs_info)    
    MMEClimates.plot_attribution_plot(cfg, provenance_dic)
    MMEClimates.plot_uncert_plot(cfg, provenance_dic)
    for MMEClimate in MMEClimates.climates:
        output_dic['Multi-Model-Ensemble'] = MMEClimate.reform_to_yml()
    output_dic['Multi-Model-Ensemble'] = MMEClimates.save_risk_ratios(cfg)

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