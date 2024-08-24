'''
This diagnostic calculates return periods for an extreme event of 
interest from observational datasets. 

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
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata, get_diagnostic_filename, save_data, save_figure
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools
from esmvaltool.diag_scripts.shared import ProvenanceLogger

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))

def select_bins(min_val:  float | int, max_val: float | int):
    '''
    Create the bins and x-values for plotting

    Parameters
    ----------
    min_val :  
        minimum value 
    max_val : 
        maximum value

    Returns
    -------
    bins : np.ndarray
        bins for plotting 
    x_fine : np.ndarray
        x values with fine resolution 
    '''

    if max_val - min_val < 5: 
        bins = np.arange(min_val, max_val+0.1, 0.1).round(1)
        x_fine = np.arange(min_val, max_val+0.01, 0.01).round(2)
    elif 5 <= max_val - min_val < 10: 
        bins = np.arange(min_val, max_val+0.1, 0.5).round(1)
        x_fine = np.arange(min_val, max_val+0.01, 0.05).round(2)
    elif 10 <= max_val - min_val < 30: 
        bins = np.arange(min_val, max_val+0.1, 1).round()
        x_fine = np.arange(min_val, max_val+0.01, 0.1).round(1)
    elif 30 <= max_val - min_val < 50: 
        bins = np.arange(min_val, max_val+0.1, 2).round()
        x_fine = np.arange(min_val, max_val+0.01, 0.2).round(1)
    else: 
        bins = np.arange((min_val - 5)//5, (max_val + 5)//5 + 0.1, 5).round()
        x_fine = np.arange(min_val, max_val+0.01, 0.5).round(1)    

    return bins, x_fine


def obtain_covariate_data(input_data: dict, dataset: str, 
                                            covariate_dataset: str| None, 
                                            smooth_kernel: int | None, 
                                            add_year: bool | None):
    '''
    This function obtains covariate data array from the input data

    Parameters
    ----------
    input_data: 
        metadata dictionary from ESMValTool with all available data
    dataset:
        name of the dataset of extreme of interest
    covariate_dataset: 
        name of the covariate dataset from recipe, is ignored if
        dataset has a valid entry in obs_covariate group
    smooth_kernel:
        the length of the smoothing window
    add_year:
         flag to possibly add a year at the end of the covariate data

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
        if smooth_kernel:
            covariate_data = covariate_data.rolling(time=smooth_kernel, 
                                                    min_periods=1).mean()
        if add_year:
            last_value = covariate_data[-1]
            last_date = datetime.strptime(str(last_value.time.values).split('.'
                                                      )[0],'%Y-%m-%dT%H:%M:%S')
            days = 366 if last_date.year%4 == 0 else 365
            new_time = np.datetime64(last_date + timedelta(days=days))
            last_value.assign_coords(time=new_time)
            covariate_data = xr.concat([covariate_data, 
                                        last_value.assign_coords(
                                        time=new_time)], dim='time')
    else: 
        covariate_data = None
        logger.info("No covariate data was found")

    return covariate_data


def calculate_anomaly(data: xr.DataArray, anomaly_dic: dict, name: str,
                       input_data: dict): 
    '''
    Parameters
    ----------
    data : 
        a data array for which anomalies to be calculated
    anomaly_dict : 
        a dictionary with the anomaly calculation parameters
    name :
        name of the dataset
    input_data : 
        ESMValTool metadata dictionary with all the available data

    Returns
    -------
    anomaly_data:
        xarray DataArray with the anomaly data
    
    Raises
    ------
    ValueError
        if the base period file for the data doesn't exist or
        if the short names don't match
    NotImplementedError
        if the anomly type is not 'anomaly' or 'ratio' 
    '''
    base_period = anomaly_dic.get('period')
    base_year_start = base_period[0] ; base_year_end = base_period[1]
    base_data = data.sel(time=slice(f"{base_year_start}-01-01", 
                                    f"{base_year_end}-12-31"))
    if len(base_data)==0:
        logger.info(f"The base period {base_year_start}-{base_year_end} "
                    f"wasn't available in the data file for dataset {name}. "
                    "Looking for an 'anomaly' group for this dataset.")
        base_meta = select_metadata(input_data.values(), dataset=name,
                                        variable_group='anomaly')
        if not base_meta:
            raise ValueError(f"No data for {name} in 'anomaly' group "
                             "has been found. ")
        else:
            base_path = base_meta[0]['filename']
            base_short_name = input_data[base_path]['short_name']
            if base_data.name != base_short_name:
                raise ValueError(f"The short name ({base_short_name}) of the "
                                 f"base dataset for {name} doesn't match the "
                                 f"original data ({base_data.name}).")
            base_data = xr.open_dataset(base_path)[base_short_name]

    if 'time' in base_data.dims:
        base_data = base_data.mean(dim='time')

    ano_type = anomaly_dic.get('type')

    if ano_type == 'anomaly':
        anomaly_data = data - base_data
    elif ano_type == 'ratio':
        anomaly_data = data / base_data
    else:
        raise NotImplementedError(f"The anomaly type {ano_type} isn't "
                                  f"impleneted. Possible entries: 'anomaly' "
                                  "and 'ratio'.")

    return anomaly_data


def define_event(data: xr.DataArray, ana_year: int | None):
    '''
    Parameters
    ----------
    data:
        an observational dataset data from which event is defined 
    ana_year: 
        a year of interest which is coming from the recipe, optional
    
    Returns
    -------
    event : float 
        event of interest
    year : int
        year of the event of interest 
    idx : int
        index of the event of interest in the data array
    '''
    if ana_year:
        idx = np.argwhere(data.year.data == ana_year).flatten()[0]
    else:
        idx = -1 
    event = data.isel(time=idx).data.item()
    year = data.year.data[idx]

    return float(event), int(year), int(idx)
    

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
    def __init__(self, data: np.ndarray, event: float, 
                 weights: np.ndarray | None = None, 
                 initial: dict | None = None):
        '''
        Parameters
        ----------
            data : 
                an array with the timeseries that is used for fitting
            event : 
                the strength of the event for return period calculation
            weights : 
                an array with the weights, the same shape as data
            initial : 
                an dictionary with the initial conditions for the fit
                keys are {'location', 'scale', 'shape'} 
        
        Raises
        ------
        ValueError 
            if the weights and data are do not have matching shape or
            the keywords in the initial condition dictionary are not 
            correct or initial conditions are not float type
        '''
        stat_gev = cex.fit_gev(data, returnValue=event, getParams=True,
                               weights=weights, initial=initial)
        try:
            self.rp = float(np.exp(stat_gev['logReturnPeriod'][0]))
            self.shape = float(-1 *stat_gev['mle'][2])
            self.loc = float(stat_gev['mle'][0])
            self.scale = float(stat_gev['mle'][1])
        except:
            self.rp = float(np.nan) ; self.shape = float(np.nan)
            self.loc = float(np.nan) ; self.scale = float(np.nan)
        self.notation= 'scipy'
        if weights:
            if data.shape != weights.shape:
                raise ValueError("The shapes of the data and weights for "
                                 "stationary GEV fit have unmatching shapes")
        if initial:
            if list(initial.keys()) != ['location', 'scale', 'shape']:
                raise ValueError("The initial conditions supposed to be "
                                 "['location', 'scale', 'shape'], currently "
                                 f"it is {list(initial.keys())}")
            if not(isinstance(initial['shape'], float)):
                raise ValueError("The type of shape parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['shape'])}")
            if not(isinstance(initial['scale'], float)):
                raise ValueError("The type of scale parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['scale'])}")
            if not(isinstance(initial['location'], float)):
                raise ValueError("The type of location parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['location'])}")

    def obtain_pdf(self, x_gev : np.ndarray):
        '''
        This function calculates pdf for x_gev

        Parameters
        ----------
        x_gev : 
            array with x-es for which pdfs are calculated
        '''
        self.pdf = gev.pdf(x_gev, self.shape, loc=self.loc, scale=self.scale)

        return

    def obtain_rp_curve(self, x_gev : np.ndarray):
        '''
        This function calculates return period curve for x_gev

        Parameters
        ----------
        x_gev : 
            array with x-es for which pdfs are calculated
        '''
        self.rp_curve = 1/gev.sf(x_gev, self.shape, loc=self.loc, scale=self.scale)

        return

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
    cov_value: float
        value for the covariate for the idx (e.g., year of interest)
    notion : str
        which GEV shape parameter notion is used. Possible options: 
        scipy or R. This attribute is used simply for checking.
    '''
    def __init__(self, data: np.ndarray, covariate: np.ndarray, 
                     event: float, idx: int = -1, initial: dict | None = None):
        '''
        Parameters
        ----------
            data : 
                an array with the timeseries that is used for fitting
            covariate: 
                an array with a covariate that is used for fitting
                (e.g., GSAT, CO2 etc.)
            event : 
                the strength of the event for return period calculation
            idx : 
                index for which the retrun period has to be calculated
                (e.g., year of interest), default is -1 (the last)
            initial : 
                an dictionary with the initial conditions for the fit
                keys are {'location', 'scale', 'shape'}
            
        
        Raises
        ------
        ValueError
            If the length of the data and covariate is different or
            the keywords in the initial condition dictionary are not
            correct or initial conditions are not float type
        '''
        if len(data) != len(covariate):
            raise ValueError("The lengths of the data and covariate"
                                                        "don't match.")
        stat_gev = cex.fit_gev(data, covariate, locationFun=1, 
                               returnValue=event, getParams=True)
        try: 
            self.rp = float(np.exp(stat_gev['logReturnPeriod'][idx]))
            self.shape = float(-1 *stat_gev['mle'][3])
            self.loc_a = float(stat_gev['mle'][0])
            self.loc_b = float(stat_gev['mle'][1])
            self.scale = float(stat_gev['mle'][2])
        except: 
            self.rp = float(np.nan)  
            self.loc_a = float(np.nan) ; self.loc_b = float(np.nan)
            self.shape = float(np.nan) ; self.scale = float(np.nan)
        self.loc_idx = float(self.loc_a * covariate[idx] + self.loc_b)
        self.cov_value = float(covariate[idx])
        self.notation= 'scipy'
        if initial:
            if not(isinstance(initial['shape'], float)):
                raise ValueError("The type of shape parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['shape'])}")
            if not(isinstance(initial['scale'], float)):
                raise ValueError("The type of scale parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['scale'])}")
            if not(isinstance(initial['location'], float)):
                raise ValueError("The type of location parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['location'])}")


def obtain_confidence_intervals(data: xr.DataArray, event: float, 
                                conf_int: list[int | float], seed: int=1,  
                                yblock: int = 1, 
                                covariate_data: xr.DataArray | None = None,
                                covariate_value: float | None = None):
    '''
    This function calculates confidence intervals with bootstrap

    Parameters
    ----------
    data: 
        data array with the extremes data
    event: 
        strength of the extreme event in the data units
    conf_int: 
        quantiles for the confidence interval calculation
    seed: 
        optional seed for the start of the bootstrap random sequence
    covariate_data:
        optional data with the covariate for non stationary GEV
    covariate_value:
        optional covariate value for the year of interest

    Returns
    -------
    out_dic:
        dictionary with quantile number as a key and corresponding RP 
        percentile as a value      
    '''
    rps = []

    high = int(np.ceil(len(data)/yblock))
    if len(data)%yblock > 0: 
        logger.warning(f"The data length {len(data)} is not devisible by "
                       f"the averaging block {yblock}, the data for bootstrap "
                       f"will be corrected by choosing {(len(data)//yblock)+1}"
                       f" samples and truncated to {len(data)}.")

    rng = np.random.default_rng(seed)
    for i in range(1000):
        fit_idx = rng.integers(low=0, high=high, size=high)
        fit_idx = np.asarray([np.arange(i, i+yblock) for i in fit_idx]
                                                                    ).flatten()
        if covariate_data is None:
            TmpStat = StationaryRP(data[fit_idx].data, event)
            shape = TmpStat.shape
            scale = TmpStat.scale
            loc = TmpStat.loc
        else: 
            TmpNonStat = NonStationaryRP(data[fit_idx].data, 
                                         covariate_data[fit_idx].data, event)
            shape = TmpNonStat.shape
            scale = TmpNonStat.scale
            loc = TmpNonStat.loc_a*covariate_value + TmpNonStat.loc_b
        try:
            rp = np.around(1/gev.sf(event, shape, loc=loc, scale=scale), 1)
        except:
            rp = np.nan
        rps.append(rp) 
    
    cis = np.nanpercentile(rps, conf_int, 
                           method='closest_observation').round(1)      
    out_dic = {}

    for n, perc in enumerate(conf_int):
        out_dic[perc] = float(cis[n])

    return out_dic
 

class SingleObsDataset:
    '''
    Class with individual OBS data and appropriate GEV statistics
    
    Attributes
    ----------
    name : str
        name of the dataset
    data : xr.DataArray
        data array with the extremes data
    covariate_data: xr.DataArray | None
        data array with the covariate data (or None)
    event : float 
        the strength of the event 
    year : int 
        year of the event of interest  
    StationaryRP : class
        class with info on the stationary GEV fit
    stationary_rp_CI: dict
        dictionary with the confidence interval from stationary GEV
    NonStationaryRP (optional) : class
        class with info on the non-stationary GEV fit
    non_stationary_rp_CI (optional): dict
        dictionary with the confidence interval from non-stationary GEV
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
        self.name = dataset

        all_input_data = cfg['input_data']
        short_name = all_input_data[path]['short_name']
        data = xr.open_dataset(path)[short_name]
        if cfg.get('anomaly'):
            data = calculate_anomaly(data, cfg['anomaly'], dataset)

        self.data = data

        covariate_dataset = cfg.get('covariate_dataset')
        smooth_kernel = cfg.get('smooth_covariate_years')
        add_year = cfg.get('add_covariate_year')
        covariate_data = obtain_covariate_data(all_input_data, dataset, 
                                               covariate_dataset,
                                               smooth_kernel, add_year)
        self.covariate_data = covariate_data

        event, year, idx = define_event(self.data, cfg.get('analysis_year'))
        self.event = event
        self.year = year

        self.StationaryRP = StationaryRP(self.data.data, self.event)
        # TODO add the y_block as for models
        self.stationary_cis = obtain_confidence_intervals(
                                                self.data, self.event,
                                                conf_int=cfg['CI_quantiles'],
                                                seed=cfg.get('bootstrap_seed'),
                                                yblock=cfg.get('yblock'))
        if self.covariate_data is None:
            logger.info(f"For {dataset} only the stationary GEV was calculated")
        else:
            self.NonStationaryRP = NonStationaryRP(self.data.data, 
                                                   self.covariate_data.data,
                                                   self.event, idx)
            self.non_stationary_cis = obtain_confidence_intervals(
                                                self.data, self.event,
                                                conf_int=cfg['CI_quantiles'],
                                                seed=cfg.get('bootstrap_seed'),
                                                covariate_data=covariate_data,
                                                covariate_value=
                                                self.NonStationaryRP.cov_value,
                                                yblock=cfg.get('yblock'))
    
    def reform_to_yml_dic(self): 
        """Reforms class to the output yml friendly dictionary"""
        dataset_dic = {}

        dataset_dic['event'] = self.event; dataset_dic['year'] = self.year
        dataset_dic['units'] = self.data.units
        dataset_dic['stationary_gev'] = { 'shape': self.StationaryRP.shape,
                                        'scale': self.StationaryRP.scale,
                                        'loc': self.StationaryRP.loc,
                                        'RP' : self.StationaryRP.rp,
                                        'CI': self.stationary_cis  }
        if self.covariate_data is None:
            logger.info(f"For {self.name} only stationary entry is created")
        else:
            dataset_dic['non_stationary_gev'] = { 
                'shape': self.NonStationaryRP.shape,
                'scale': self.NonStationaryRP.scale,
                'loc_a': self.NonStationaryRP.loc_a,
                'loc_b': self.NonStationaryRP.loc_b, 
                'loc_idx': self.NonStationaryRP.loc_idx, 
                'cov_value' : self.NonStationaryRP.cov_value,
                'RP': self.NonStationaryRP.rp,
                'CI': self.non_stationary_cis
            }

        return dataset_dic


def plot_figure(Dataset: SingleObsDataset, cfg: dict, prov_dic: dict):
    '''
    Plots figure similar to Figs. 3-4 in Malinina and Gillett (2024)

    Parameters
    ----------
    Dataset:
        class with the data to be plotted
    cfg:
        dictionary with all essentiall parameters
    prov_dic: 
        provenance dictionary
    '''
    prov_dic['caption'] = f"{Dataset.name} GEV statistics for" +\
                                                f" the {Dataset.year} event"
    
    mpl_st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(mpl_st_file)
    color_st = eplot.get_dataset_style(Dataset.name, cfg.get('color_style'))

    label = cfg.get('var_label') if cfg.get('var_label') else Dataset.data.name
    reg = cfg.get('region') if cfg.get('region') else 'region'

    min_val = np.floor(Dataset.data.min()*0.9 if Dataset.data.min()> 0 
                                               else Dataset.data.min()*1.1)
    max_val = np.ceil(Dataset.data.max()*1.1 if Dataset.data.max()> 0 
                                               else Dataset.data.max()*0.9)

    bins, x_fine = select_bins(min_val, max_val)

    if Dataset.covariate_data is None:
        obs_surv_func = gev.sf(x_fine, Dataset.StationaryRP.shape, 
                               loc=Dataset.StationaryRP.loc,
                               scale=Dataset.StationaryRP.scale)
        up_ci_bound = Dataset.non_stationary_cis[
                                        max(Dataset.non_stationary_cis.keys())]
        rp = np.around(Dataset.NonStationaryRP.rp, 1)
    else:
        obs_surv_func = gev.sf(x_fine, Dataset.NonStationaryRP.shape, 
                        loc=Dataset.NonStationaryRP.loc_idx,
                        scale=Dataset.NonStationaryRP.scale)
        up_ci_bound = Dataset.non_stationary_cis[
                                        max(Dataset.non_stationary_cis.keys())]
        rp = np.around(Dataset.NonStationaryRP.rp, 1)

    fig_obs = plt.figure(constrained_layout=False)
    fig_obs.set_size_inches(12., 8.)
    gs = fig_obs.add_gridspec(nrows=2, ncols=5)
    ax_hist = fig_obs.add_subplot(gs[0, :-2])
    ax_surv= fig_obs.add_subplot(gs[0, -2:])
    ax_tseries= fig_obs.add_subplot(gs[1, :])

    ax_hist.hist(Dataset.data.data, bins=bins, edgecolor=color_st['facecolor'],
                          facecolor=color_st['color'], alpha=0.3, density=True)
    ax_hist.scatter(Dataset.event, 0, c=color_st['facecolor'], marker='x',
                                                          s=100, clip_on=False)
    ax_hist.set_xlim(bins[0], bins[-1])
    ax_hist.set_xlabel(f"{label}, {Dataset.data.units}")
    ax_hist.set_ylabel('Number density')
    ax_hist.set_title('(a) Histogram')

    ax_surv.plot(x_fine[obs_surv_func>0.00001], 
                 1/obs_surv_func[obs_surv_func>0.00001], c=color_st['color'])
    ax_surv.scatter(Dataset.event, rp, c=color_st['facecolor'], 
                                            marker='x', s=100, clip_on=False)
    ax_surv.set_title(
                f"(b) {Dataset.name} {label} return period in {Dataset.year}")
    ax_surv.set_xlim(bins[0], bins[-1])
    ax_surv.set_ylim(1,1000 if up_ci_bound < 1000 else 10000)
    ax_surv.set_yscale('log')
    ax_surv.grid(color='silver', axis='both', alpha=0.5)
    ax_surv.set_ylabel('years')
    ax_surv.text(0.15, 0.5, f"{Dataset.name} return period\n{rp} years",
                                transform = ax_surv.transAxes, clip_on=False)

    ax_tseries.plot(Dataset.data.year.data, Dataset.data.data, c=color_st['color'])
    ax_tseries.scatter(Dataset.year, Dataset.event, marker='o', s=100, 
                       facecolors='None', edgecolors=color_st['facecolor'])
    ax_tseries.set_ylim(bins[0], bins[-1])
    ax_tseries.set_xlim(Dataset.data.year.data[0] - 0.5, 
                        Dataset.data.year.data[-1] + 0.5)
    ax_tseries.text(0.7, 0.92, f"{Dataset.name}({Dataset.year}): " + \
                    str(np.around(Dataset.event, 1))+ f" {Dataset.data.units}",
                    transform = ax_tseries.transAxes, clip_on=False)
    ax_tseries.set_title(f"(c) {Dataset.name} {label} timeseries")
    ax_tseries.grid(color='silver', axis='both', alpha=0.5)
    ax_tseries.set_xlabel('time')
    ax_tseries.set_ylabel(f"{label}, {Dataset.data.units}")

    fig_obs.suptitle(f"{Dataset.name} {label} in {reg} and its GEV fit")

    fig_obs.tight_layout()

    fig_path = os.path.join(cfg['plot_dir'], 
                            f"{Dataset.name}_{reg}_{label.lower()}")
    save_figure(fig_path, prov_dic, cfg, fig_obs, close=True)

    return


def main(cfg):
    """Main function, essentially diagnostic itself"""

    provenance_dic = {'authors': ['malinina_elizaveta']}

    obs_data_metadata = select_metadata(cfg['input_data'].values(), 
                                        variable_group='obs_data')
    datasets = group_metadata(obs_data_metadata, 'dataset', sort=True)

    output_dict = {}

    for dataset in datasets.keys(): 
        logger.info(f"Processing {dataset}")
        datapath = datasets[dataset][0]['filename']
        logger.info(f"Using {datapath}")
        Dataset = SingleObsDataset(dataset, datapath, cfg)
        plot_figure(Dataset, cfg, provenance_dic)
        dataset_yml_entry = Dataset.reform_to_yml_dic() 
        output_dict[dataset] = dataset_yml_entry
    
    obs_info_path = os.path.join(cfg['work_dir'],'obs_information_event.yml')
    with open(obs_info_path, 'w') as obs_info_yml:
        yaml.dump(output_dict, obs_info_yml, sort_keys=False)
    logger.info(f"Saved OBS info into {obs_info_path}")

    logger.info('Successfully ran the extremes observational diagnostic')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)