'''
This diagnostic calculates return periods and risk ratios for an
extreme event of interest from models, thus providing actual
attribution.

This diagnostic requires an ancestor task to determine the event and
its rarity from the observations. 

The preprocessed data should be provided in form of spatially averaged
timeseries of the annual maxima/minima of a variable (e.g., tasmax).
While some seasonal data can be used at the preprocessor stage, to the
diagnostic only one value per year shall be provided. Example:
extract JJA season first, calculate annual mean using JJA data only.

This code was developped for the results used in Malinina&Gillett(2024)
and Gillet et al. (2022) (Weather and Climate Extremes). 

The result of this diagnostic is a figure similar to Fig 9 from 
Malinina&Gillett(2024) and a yml-file with the attribution statistics.

Author: Elizaveta Malinina (elizaveta.malinina-rieger@ec.gc.ca)
Initial development: 2021-2022
Last updated: August 2024
'''

import xarray as xr
import climextremes as cex
import logging
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
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
from esmvaltool.diag_scripts.eccc_extremes.observational_return_periods import StationaryRP, select_bins

logger = logging.getLogger(os.path.basename(__file__))


def get_percentiles_dic(data : np.ndarray | list, 
                                percentiles: list[int| float], round : int):
    '''
    This function returns a dictionary with percentiles

    Parameters
    ----------
    data :
        array or list with the data for which percentiles calculated
    percentiles : 
        list with percentiles
    round : 
        the level of rounding of the percentiles
    
    Returns
    -------
    ci_dict : dict
        ditionary with percentiles as keys
    '''

    ci_dict = {}

    for perc in percentiles:
        ci_dict[perc] = float(np.nanpercentile(data, perc, 
                                    method='closest_observation').round(round))

    return ci_dict


def calculate_anomaly(data: np.array, input_data: list, file_meta: dict,
                                            anomaly_type: str = 'anomaly'):
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
        a type of anomaly to calculate, default 'anomaly'
    
    Returns
    -------
    anomaly : 
        numpy array with the anomaly data
    
    Raises
    ------
    NotImplementedError
        if the anomly type is not 'anomaly' or 'ratio' 
    '''

    anomaly_type = 'anomaly' if anomaly_type is None else anomaly_type

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
    # convert to array to use the best counting methods
    datasets = np.asarray(datasets)
    # weights are assigned in a way that each model in a sample gets
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

    # number of unique datasets and realizations in each dataset
    unique_dtsts, reals_dtsts = np.unique(datasets, return_counts=True)

    # number of iterations for bootstrap
    pool_iter_n = np.max([reals_dtsts.max()*100, 1000])

    # the number of the datasets to be included in the pool
    # in case the dataset number is less than 3, using 3 as
    # the minimum ensemble number of to draw conclusions on
    # 3 most likely to be used in the case of only one model
    n_dataset = np.max([len(unique_dtsts), 3])

    n_years = data.shape[1]; n_all_reals = data.shape[0]

    pool_size = int(np.around(n_dataset*n_years/yblock))

    # defining the seed for the random sequence
    rng = np.random.default_rng(seed)

    bootstrap_samples = list()

    for i in range(pool_iter_n):
        sample_data = list() ; sample_datasets = list()
        mod_idxs = rng.integers(0, high=n_all_reals, size=pool_size)
        for mod_idx in mod_idxs:
            y_idx = rng.integers(0, high=n_years-yblock)
            sample_data.append(data[mod_idx, y_idx:y_idx+yblock])
            sample_datasets.append(datasets[mod_idx])
        sample_data = np.asarray(sample_data)
        if cfg.get('model_weighting'):
            sample_weights = determine_weights(sample_datasets, 
                                            sample_data.shape).flatten()
        else:
            sample_weights = None
        SampleGEV = StationaryRP(sample_data.flatten(), event, 
                                                        weights=sample_weights)
        bootstrap_samples.append(SampleGEV)

    return bootstrap_samples


def calculate_intensity(GEV : StationaryRP, event_rp: float,
                                            x_fine : np.ndarray | None = None):
    '''
    Calculate intensity of the event if the event_rp rarity

    Parameters
    ----------
    GEV : 
        class with the stationary GEV info
    event_rp :
        return period of the event for which inetensity is calculated
    x_fine : 
        optional array with the strengths. If None, internal x_gev from
        GEV will be used

    Returns
    -------
    intensity : float 
        the strength of the event of the event_rp probability
    '''

    if x_fine is None:
        rp_func = GEV.rp_curve ; x_fine = GEV.x_gev
    else:
        rp_func = 1/gev.sf(x_fine, GEV.shape, loc=GEV.loc, scale=GEV.scale)

    intensity = float(x_fine[np.argmin(np.abs(rp_func-event_rp))])

    return intensity


class ObsData:
    '''
    Attributes
    ----------
    name : str
        Name of the dataset
    event : float
        Strength of the event of interest
    year: int 
        Year of the event
    event_rp : float
        Return period of the event for the selected type of fit
    rp_ci_max : float
        Upper border of the return period confidence interval
    rp_ci_min : float
        Lower border of the return period confidence interval
    model_event : float (optional)
        Strength of the event of the event_rp rarity in the factual
        climate in the models
    '''
    def __init__(self, obs_info: dict, cfg: dict):
        '''
        Parameters
        ----------
        obs_info : 
            loaded from ancestor file dictionary
        cfg : 
            internal ESMValTool dictionary with input from the recipe
        '''
        self.name = cfg['obs_ref_dataset']
        self.event = obs_info[self.name]['event']
        obs_fit = cfg.get('obs_gev') if cfg.get('obs_gev') else 'stationary'
        self.event_rp = obs_info[self.name][obs_fit+'_gev']['RP']
        self.year = obs_info[self.name]['year']
        max_ci = np.max(list(obs_info[self.name][obs_fit+'_gev']['CI'].keys()))
        min_ci = np.min(list(obs_info[self.name][obs_fit+'_gev']['CI'].keys()))
        self.rp_ci_max = obs_info[self.name][obs_fit+'_gev']['CI'][max_ci]      
        self.rp_ci_min = obs_info[self.name][obs_fit+'_gev']['CI'][min_ci]

    def add_model_event(self, model_event: float):
        '''
        Parameters
        ----------
        model_event : 
            Strength of the event of the event_rp rarity in the factual
            climate in the models
        '''
        self.model_event = model_event  


class Climate:
    '''
    Attributes
    ----------
    name : str
        name of the climate
    start_year : int
        start year of the climate period
    end_year : int
        end year of the climate period (included)
    short_name : str
        name of the variable
    units : str
        units of the data
    data : np.ndarray
        data for the climate, 2D  (N_real, N_years)
    weights : np.ndarray | None
        weights for the data. If array, flattened (N_real*N_years)
    bootstrap : list[StationaryRP]
        list with the boostrap realizations
    BestGuessGEV : StationaryRP
        the best guess of the stationary GEV 
    rp_ci: dict
        dictionary with the confidence intervals
    intensity: float
        strength of the event in of the provided rarity in this climate
    intensity_CI: dict
        dictionary with the confidence interval of the intensities
    '''

    def __init__(self, input_data: list, group: str, anomaly_calculation: bool,
                                                    cfg: dict, event: float):
        '''
        Parameters
        ----------
        input_data : 
            list with the climates metadata including file paths
        group : 
            name of the group/climate
        anomaly_calculation: 
            flag which switches on/off anomaly calculation
        cfg : 
            internal ESMValTool dictionary with input from the recipe
        event : 
            strength of the event of interest
        '''
        
        group_info = select_metadata(input_data, variable_group=group)

        self.name = group
        # since the group is uniform for these keys, saving the first
        self.start_year = list(group_metadata(group_info, 'start_year'
                                                                  ).keys())[0]
        self.end_year = list(group_metadata(group_info, 'end_year').keys())[0]
        self.short_name = list(group_metadata(group_info, 'short_name'
                                                                  ).keys())[0]
        self.units = list(group_metadata(group_info, 'units').keys())[0]

        # obtain data and datasets to calculate weights
        files = group_metadata(group_info, 'filename')
        group_data = list(); datasets = list()
        for file in files.keys():
            file_meta = select_metadata(input_data, filename=file)[0]
            file_dataset = file_meta['dataset']
            data = xr.open_dataset(file)[self.short_name].data
            if anomaly_calculation:
                data = calculate_anomaly(data, input_data, file_meta, 
                                            cfg.get('anomaly_type'))
            group_data.append(data)
            datasets.append(file_dataset)
        group_data = np.asarray(group_data) ; datasets = np.asarray(datasets)

        self.data = group_data
        if cfg.get('model_weighting'):
            self.weights = determine_weights(datasets, group_data.shape).flatten()
        else:
            self.weights = None

        self.bootstrap = bootstrap_gev(self.data, event, datasets, cfg)
        if cfg.get('initial_conditions'): 
            initial_cond = calculate_initial_cond(self.bootstrap)
        else:
            initial_cond = None
        self.BestGuessGEV = StationaryRP(self.data.flatten(), event,
                                                weights=self.weights,
                                                initial=initial_cond)
        self.rp_ci = get_percentiles_dic([b.rp for b in self.bootstrap],
                                          cfg['CI_quantiles'], 1)
    
    def obtain_rps_and_pdfs_curves(self, x_gev):
        '''This function gets RP curves and PDFs for x_gev'''

        self.BestGuessGEV.obtain_pdf(x_gev)
        self.BestGuessGEV.obtain_rp_curve()
        for i in range(len(self.bootstrap)):
            self.bootstrap[i].obtain_pdf(x_gev)
            self.bootstrap[i].obtain_rp_curve()

    def calculate_intensities(self, x_fine : np.ndarray, event_rp : float,
                                                        conf_ints : list):
        '''
        This function calculates strength of the event_rp rarity event

        Parameters
        ----------
        x_fine : 
            array with the x-es used for the strength calculation
        event_rp : 
            return period of the event for which to calculate strength
        conf_ints : 
            list with percentiles for confidence interval calculation
        '''

        self.obtain_rps_and_pdfs_curves(x_fine)

        self.intensity = calculate_intensity(self.BestGuessGEV, event_rp)
        intensities = list()
        for i in range(len(self.bootstrap)):
            intensities.append(calculate_intensity(self.bootstrap[i], 
                                                                event_rp))

        self.intensity_CI = get_percentiles_dic(intensities, conf_ints, 1)

    def update_event(self, event: float):
        '''
        Paramaters
        ----------
        event : 
            strength of the event used for return period calculations
        '''
        self.BestGuessGEV.rp = float(1/gev.sf(event, self.BestGuessGEV.shape, 
                                               loc=self.BestGuessGEV.loc,
                                               scale=self.BestGuessGEV.scale))
        for i in range(len(self.bootstrap)):
            self.bootstrap[i].rp = float(1/gev.sf(event, 
                                                self.bootstrap[i].shape,
                                                loc=self.bootstrap[i].loc,
                                                scale=self.bootstrap[i].scale))
        self.rp_ci = get_percentiles_dic([b.rp for b in self.bootstrap],
                                                    list(self.rp_ci.keys()), 1)


    def reform_to_yml(self):
        '''This function returns the information in yml-dictionary'''

        out_dic = { self.name : {'shape': self.BestGuessGEV.shape,
                                 'scale': self.BestGuessGEV.scale,
                                 'loc': self.BestGuessGEV.loc,
                                 'RP': self.BestGuessGEV.rp,
                                 'RP_CI': self.rp_ci, 
                                 'intensity': self.intensity,
                                 'intensity_CI': self.intensity_CI}}

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
                bootstrap_rrs.append(clim_list[old_clim_idx].bootstrap[i].rp/\
                                    clim_list[clim_idx].bootstrap[i].rp)
            
            risk_ratio_dic[key]['CI'] = get_percentiles_dic(bootstrap_rrs,
                                                                   ci_percs, 2)

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
    risk_ratios : dict
        dictionary with the risk ratios 
    max_value : float
        maximum value in the data for all climates 
    min_value : float
        minimum value in the data for all climates
    '''

    def __init__(self, input_data: list, dataset_name: str, cfg: dict, 
                                                            ObsInfo: ObsData):
        '''
        Parameters
        ----------
        input_data: 
            list with all the metadata and files for the dataset
        dataset_name:
            name of the dataset
        cfg: 
            internal ESMValTool dictionary with the input from recipe
        ObsInfo:
            class with the observational information
        '''

        self.name = dataset_name
        self.climates = []

        groups = list(group_metadata(input_data, 'variable_group').keys())
        if 'anomaly' in groups:
            groups.remove('anomaly'); anomaly_calculation = True

        for group in groups:
            GroupClimate = Climate(input_data, group, anomaly_calculation, cfg,
                                                                 ObsInfo.event)
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
                    # the event is redifined using the best guess event 
                    # probability obtained from observations
                    event = calculate_intensity(Clim.BestGuessGEV, 
                                                            ObsInfo.event_rp, x_fine)
                    ObsInfo.add_model_event(event)
                    break
            for Clim in self.climates:
                Clim.update_event(event)
        for Clim in self.climates:
            Clim.calculate_intensities(x_fine, ObsInfo.event_rp,
                                                        cfg['CI_quantiles'])
        self.risk_ratios = calculate_risk_ratios(self.climates, 
                                                        cfg['CI_quantiles'])


    def plot_attribution_plot(self, cfg: dict, ObsInfo: ObsData, 
                                                        prov_dic: dict):
        '''
        This function plots the Fig 9 in Malinina and Gillett (2024)

        Parameters
        ----------
        cfg : 
            internal ESMValTool dictionary with the input from recipe
        ObsInfo : 
            class with the observational information
        prov_dic : 
            dictionary with provenance parameters
        '''

        prov_dic['caption'] = f"{self.name} (a) probability density function"+\
                              f" ,(b) QQ-plot and (c) return periods."
        mpl_st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
        plt.style.use(mpl_st_file)

        # getting info for the obs plotting
        obs_color = eplot.get_dataset_style(
                                        ObsInfo.name, cfg.get('color_style'))

        reg = cfg['region'] if cfg.get('region') else 'region'
        label = cfg['var_label'] if cfg.get('var_label') else self.climates[0
                                                                    ].short_name
        units = cfg['units'] if cfg.get('units') else self.climates[0
                                                                    ].units
        bins, x_fine = select_bins(self.min_value, self.max_value)

        fig = plt.figure(constrained_layout=False)
        fig.set_size_inches(12., 8.)
        gs = fig.add_gridspec(nrows=2, ncols=6)
        ax_hist = fig.add_subplot(gs[:, :-2])
        ax_qq= fig.add_subplot(gs[0, -2:])
        ax_rps= fig.add_subplot(gs[1, -2:])

        for clim in self.climates:
            clim_leg_label = f"{clim.start_year}-{clim.end_year}"
            clim_color_st = eplot.get_dataset_style(clim.name, 
                                                    cfg.get('color_style'))
            # plotting histogram part
            ax_hist.hist(clim.data.flatten(), weights=clim.weights, bins=bins,
                         edgecolor=clim_color_st['facecolor'],  alpha=0.3,
                         color=clim_color_st['color'], density=True, zorder=2,
                         label=clim_leg_label)
            ax_hist.plot(clim.BestGuessGEV.x_gev, clim.BestGuessGEV.pdf, 
                         color=clim_color_st['color'], zorder=3,
                         label=f"GEV fit {clim_leg_label}")
            # determine the pdf uncert interval
            pdf_max_q = np.nanpercentile([b.pdf for b in clim.bootstrap], 
                                         np.max(cfg['CI_quantiles']), axis=0,
                                         method='closest_observation')
            pdf_min_q = np.nanpercentile([b.pdf for b in clim.bootstrap], 
                                         np.min(cfg['CI_quantiles']), axis=0,
                                         method='closest_observation')
            ax_hist.fill_between(clim.BestGuessGEV.x_gev, pdf_min_q, pdf_max_q,
                                 color=clim_color_st['color'], alpha=0.3, 
                                 linewidth=0, zorder=4)
            # QQ plot part 
            # creating quantile measures for QQ plot
            quantile_measures = np.arange(0, 1.01, 0.01)
            quantile_measures[0] = 0.001 # numpy not accepting 0 for quantile
            # calculating theoretical quantiles of the distribution
            theor_quants = gev(clim.BestGuessGEV.shape, 
                               loc=clim.BestGuessGEV.loc, 
                               scale=clim.BestGuessGEV.scale).ppf(
                                                            quantile_measures)
            # calculating practical quantiles of the data (for QQ plot)
            pract_quants = np.quantile(clim.data.flatten(), quantile_measures)
            ax_qq.scatter(theor_quants, pract_quants, zorder=3, 
                          lw=0.75*mpl.rcParams['lines.linewidth'],
                          label=clim_leg_label, marker=clim_color_st['mark'],
                          edgecolors=clim_color_st['color'], facecolors='None')
            # plotting return periods
            ax_rps.plot(clim.BestGuessGEV.x_gev, clim.BestGuessGEV.rp_curve,
                        color=clim_color_st['color'], zorder=3, 
                        label=clim_leg_label)
            # determine the rps uncert interval
            rps_max_q = np.nanpercentile([b.rp_curve for b in clim.bootstrap], 
                                         np.max(cfg['CI_quantiles']), axis=0,
                                         method='closest_observation')
            rps_min_q = np.nanpercentile([b.rp_curve for b in clim.bootstrap], 
                                         np.min(cfg['CI_quantiles']), axis=0,
                                         method='closest_observation')
            ax_rps.fill_between(clim.BestGuessGEV.x_gev, rps_min_q, rps_max_q,
                                color=clim_color_st['color'], alpha=0.3, 
                                linewidth=0, zorder=4)
            
        # add observational info
        hist_ylims = ax_hist.get_ylim()
        ax_hist.set_ylim(*hist_ylims)
        # get color for the obs dataset 
        ax_hist.vlines(ObsInfo.event, *hist_ylims, color=obs_color['color'],
                       zorder=1, label=f"{ObsInfo.name} ({ObsInfo.year})",
                       lw=mpl.rcParams['lines.linewidth']*1.5)
        ax_rps.vlines(ObsInfo.event, 0.1, ObsInfo.event_rp, 
                      color=obs_color['color'],
                      label=f"{ObsInfo.name} ({ObsInfo.year})", zorder=2)
        ax_rps.hlines(ObsInfo.event_rp, x_fine[0], ObsInfo.event, 
                                        color=obs_color['color'], zorder=2)
        ax_rps.fill_between([x_fine[0], ObsInfo.event], ObsInfo.rp_ci_min, 
                    ObsInfo.rp_ci_max if ObsInfo.rp_ci_max < 10000 else 10000,
                    color=obs_color['color'], zorder=2, alpha=0.15)
        if cfg.get('event_definition') == 'rarity':
            ax_hist.vlines(ObsInfo.model_event, *hist_ylims, 
                           color=obs_color['color'], zorder=1, ls='--',
                           label=f"Factual\n(1 in {np.around(ObsInfo.event_rp, 1)})",
                           lw=mpl.rcParams['lines.linewidth']*1.5)
            ax_rps.vlines(ObsInfo.model_event, 0.1, ObsInfo.event_rp, ls='--',
                          color=obs_color['color'], zorder=2, 
                          label=f"Factual\n(1 in {np.around(ObsInfo.event_rp, 1)})")

        # technical aspects of the plot, limits, captions etc.
        ax_qq.plot(x_fine, x_fine, c='tab:grey', zorder=1)
        # legends
        ax_hist.legend(loc=2, fancybox=False, frameon=False)
        ax_qq.legend(loc=2, fancybox=False, frameon=False, handletextpad=0.01)
        ax_rps.legend(loc=2, fancybox=False, frameon=False)
        # grids
        ax_qq.grid(color='silver', axis='both', alpha=0.5)
        ax_rps.grid(color='silver', axis='both', alpha=0.5)
        # set rp axes ylim
        ax_rps.set_ylim(1, 10000 if ObsInfo.rp_ci_max > 950 else 1000)
        ax_rps.set_yscale('log')
        # set axes x- and ylim
        ax_hist.set_xlim(x_fine[0], x_fine[-1]) 
        ax_qq.set_xlim(x_fine[0], x_fine[-1]) 
        ax_qq.set_ylim(x_fine[0], x_fine[-1]) 
        ax_rps.set_xlim(x_fine[0], x_fine[-1]) 
        # captions and titles
        ax_hist.set_title('(a) Probability density function')
        ax_hist.set_ylabel('Number density')
        ax_hist.set_xlabel(f"{label}, {units}")
        ax_qq.set_title('(b) Quality assessment')
        ax_qq.set_ylabel('Data quantile')
        ax_qq.set_xlabel('GEV quantile')
        ax_rps.set_title('(c) Return period')
        ax_rps.set_ylabel('years')
        ax_rps.set_xlabel(f"{label}, {units}")

        fig.suptitle(f"{label} in {reg}")

        fig.set_dpi(250)

        fig.tight_layout()
        fig_path = os.path.join(cfg['plot_dir'], 
                            f"attribution_{self.name}_{reg}_{label.lower()}")
        save_figure(fig_path, prov_dic, cfg, fig, close=True)

        return


def main(cfg):
    """Main function, essentially diagnostic itself."""

    # loading the observational info from ancestor task
    obs_file = io.get_ancestor_file(cfg, 'obs_information_event.yml')
    obs_info = yaml.safe_load(open(obs_file, 'r'))
    logger.info(f"Loaded information on OBS statistics from {obs_file}")
    ObsInfo = ObsData(obs_info, cfg)

    provenance_dic = {'authors': ['malinina_elizaveta'], 
                      'ancestors': [obs_file] }

    input_data = cfg['input_data']

    datasets = group_metadata(input_data.values(), 'dataset')

    output_dic = {}

    if cfg.get('individual_models'):
    # provide statistics for individual models
        for dataset in datasets.keys():
            dataset_info = datasets[dataset]
            logger.info(f"Processing {dataset}")
            DatasetClimates = Climates(dataset_info, dataset, cfg, ObsInfo)
            output_dic[dataset] = {}
            for Climate in DatasetClimates.climates:
                output_dic[dataset].update(Climate.reform_to_yml())
            output_dic[dataset]['risk_ratios'] = DatasetClimates.risk_ratios

    logger.info(f"Processing Multi-Model-Ensemble")
    MMEClimates = Climates(input_data.values(), 'Multi-Model-Ensemble', cfg, 
                                                                    ObsInfo)
    MMEClimates.plot_attribution_plot(cfg, ObsInfo, provenance_dic)
    output_dic['Multi-Model-Ensemble'] = {}
    for MMEClimate in MMEClimates.climates:
        output_dic['Multi-Model-Ensemble'].update(MMEClimate.reform_to_yml())
    output_dic['Multi-Model-Ensemble']['risk_ratios'] = MMEClimates.risk_ratios

    model_stats_path = os.path.join(cfg['work_dir'],'models_statistics.yml')
    with open(model_stats_path, 'w') as model_info_yml:
        yaml.dump(output_dic, model_info_yml, sort_keys=False)
    logger.info(f"Saved model info into {model_stats_path}")

    logger.info("Successfully completed diagnostic")


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)