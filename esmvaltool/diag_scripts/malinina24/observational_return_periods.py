import xarray as xr
# import climextremes as cex
import logging
import numpy as np
import matplotlib.pyplot as plt
import os

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_data
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def obtain_covariate_data(input_data: dict, dataset: str, 
                                            covariate_dataset: str| None):
    '''
    This function obtains covariate data array from the input data

    Parameters
    ----------
    input_data: 
        metadata dictionary from ESMValTool with all available data
    dataset:
        name of the dataset of interest
    covariate_dataset: 
        name of the covariate dataset from recipe, is ignored if
        dataset has a valid entry in obs_covariate group

    Raises
    ------
    ValueError
        if no valid covariate dataset for the dataset was found
    '''
    all_covariates = select_metadata(input_data.values(), 
                                        variable_group='obs_covariate')
    if all_covariates:
        covariate_meta = select_metadata(all_covariates, dataset=dataset)
        if not covariate_meta: 
            covariate_meta = select_metadata(all_covariates, 
                                                dataset=covariate_dataset)
            if not covariate_meta:
                raise ValueError(f"Since the dataset {dataset} does not "
                                    "have a valid 'covariate_dataset' key "
                                    "should be provided in the recipe. "
                                    f"Current entry: {covariate_dataset}")
            path = covariate_meta[0]['filename']
            short_name = input_data[path]['short_name']
            covariate_data = xr.open_dataset(path)[short_name]
        else: 
            covariate_data = None
            logger.info("No covariate data was found")

        return covariate_data

class StationaryRP:
    '''
    Class with statistics and return period from the stationary GEV

    Attributes
    ----------
    rp : float
        return period in years
    shape : float | nan
        shape parameter of the fit
    loc : float | nan
        location parameter of the fit
    scale : float| nan
        scale parameter of the fit
    notion : str
        which GEV shape parameter notion is used. Possible options: 
        scipy or R. This attribute is used simply for checking.
    '''
    def __init__(self, data: xr.DataArray, event: float):
        '''
        Parameters
        ----------
            data : 
                a data array with the timeseries that is used for fitting
            event : 
                the strength of the event for return period calculation
        '''
        stat_gev = cex.fit_gev(data.data, returnValue=event, getParams=True)
        self.rp = np.exp(stat_gev['logReturnPeriod'][0])
        self.shape = -1 *stat_gev['mle'][2]
        self.loc = stat_gev['mle'][0]
        self.scale = stat_gev['mle'][1]
        self.notation= 'scipy'

class NonStationaryRP:
    '''
    Class with statistics and return period from the non-stationary GEV

    The non-stationary GEV fit assumes that one or multiple GEV params
    are dependent on a covariate (eg., GSAT, CO2 etc). While the
    dependency can be non-linear, in climextremes only linear is 
    implemented and only for location parameter. The formula:
    loc = loc_a*covarite + loc_b

    Attributes
    ----------
    rp : float
        return period in years for the idx (e.g., year of interest)
    shape : float | nan
        shape parameter of the fit
    loc_a : float | nan
        slope parameter of the fitted location parameter 
        (loc = loc_a*covariate + loc_b)
    loc_b : float | nan
        intersect parameter of the fitted location parameter 
        (loc = loc_a*covariate + loc_b)
    loc_idx : float | nan
        the location parameter for the idx (e.g., year of interest)
    scale : float | nan
        scale parameter of the fit
    notion : str
        which GEV shape parameter notion is used. Possible options: 
        scipy or R. This attribute is used simply for checking.
    '''
    def __init__(self, data: xr.DataArray, covariate: xr.DataArray, 
                                    event: float, idx: int = -1):
        '''
        Parameters
        ----------
            data : 
                a data array with the timeseries that is used for fitting
            covariate: 
                a data array with a covariate that is used for fitting
                (e.g., GAST, CO2 etc.)
            event : 
                the strength of the event for return period calculation
            idx : 
                index for which the retrun period has to be calculated
                (e.g., year of interest), default is -1 (the last)
        
        Raises
        ------
        ValueError
            If the length of the data and covariate is different.
        '''
        if len(data.data) != len(covariate.data):
            raise ValueError("The lengths of the data and covariate"
                                                        "don't match.")
        stat_gev = cex.fit_gev(data.data, covariate.data, locationFun=1, 
                               returnValue=event, getParams=True)
        self.rp = np.exp(stat_gev['logReturnPeriod'][idx])
        self.shape = -1 *stat_gev['mle'][3]
        self.loc_a = stat_gev['mle'][0]
        self.loc_b = stat_gev['mle'][1]
        self.loc_idx = self.loc_a * covariate.data[idx] + self.loc_b
        self.scale = stat_gev['mle'][2]
        self.notation= 'scipy'


class SingleObsDataset:
    '''
    Class with individual OBS data and appropriate GEV statistics
    
    Attributes
    ----------
    name
    data
    covariate_data
    event
    StationaryRP
    stationary_rp_CI: dict
    NonStationaryRP (optional)
    non_stationary_rp_CI: dict
    '''
    def __init__(self, dataset: str, path: str, cfg: dict):
        '''
        Parameters
        ----------
        dataset : 
            name of the dataset for which the class will be initialized
        path : 
            path to the file which has the data for the dataset
        cfg:
            internal ESMValTool dict with parameters for the diagnostic
        '''
        all_input_data = cfg['input_data']
        short_name = all_input_data[path]['short_name']
        self.data = xr.open_dataset(path)[short_name]

    #     calculate_anomaly
    #         if cfg.get('anomaly') is not None:
    #           get type 
    #           if cfg['anomaly'].get('period') not None:
    #               get sub_array
    #           else:
    #               get group 'anomaly' path for dataset
    #               sub_array = xr.open_datarray(path_anomaly)

        covariate_dataset = cfg.get('covariate_dataset')
        covariate_data = obtain_covariate_data(all_input_data, dataset, 
                                               covariate_dataset)
        self.covariate_data = covariate_data

        # define event, idx

        self.StationaryRP = StationaryRP(self.data, self.event)
        if not self.covariate_data:
            self.NonStationaryRP = NonStationaryRP(self.data, 
                                                   self.covariate_data,
                                                   self.event, idx)

    #   confidence_intervals()


def plot_figure(Dataset: SingleObsDataset, cfg: dict):

    return


def create_yml_data_entry(Dataset: SingleObsDataset): 

    return

def main(cfg):

    obs_data_metadata = select_metadata(cfg['input_data'].values(), 
                                        variable_group='obs_data')
    datasets = group_metadata(obs_data_metadata, 'dataset', sort=True)

    output_dict = {}

    for dataset in datasets.keys(): 
        datapath = datasets[dataset][0]['filename']
        Dataset = SingleObsDataset(dataset, datapath, cfg)
        plot_figure(Dataset, cfg)
        dataset_yml_entry = create_yml_data_entry(Dataset) 
        output_dict[dataset] = dataset_yml_entry
    
    # save yml file with Provenance

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)