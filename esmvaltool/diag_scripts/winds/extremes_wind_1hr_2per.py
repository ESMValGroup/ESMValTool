import csv
import esmvalcore.preprocessor as eprep
import iris
import cf_units
import cftime
import climextremes as cex
import pandas as pd
import logging
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import genextreme as gev
from scipy.stats import kstest, cramervonmises

# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata, group_metadata
import esmvaltool.diag_scripts.shared.plot as eplot
from esmvaltool.diag_scripts.ocean import diagnostic_tools as diagtools

# # This part sends debug statements to stdout
logger = logging.getLogger(os.path.basename(__file__))
# logging.getLogger().addHandler(logging.StreamHandler(sys.stdout))


def obtain_obs_info(ano_obs_cb, abs_obs_cb, groups, cfg):

    obs_gsat_info = groups.pop('obs_gsat')

    ano_obs_arr = ano_obs_cb.data

    gsat_obs_cb = iris.load_cube(obs_gsat_info[0]['filename'])
    gsat_df = pd.DataFrame(gsat_obs_cb.data, columns=['gsat'])

    # smoothing gsat
    kernel_size = cfg['smooth_gsat_years']
    gsat_smooth_df = gsat_df.rolling(kernel_size, min_periods=1).mean()
    gsat_smooth_arr = gsat_smooth_df['gsat'].to_numpy()
    
    ana_year_const = iris.Constraint(time = lambda cell: cell.point.year == cfg['analysis_year'])
    ana_year_value =  float(ano_obs_cb.extract(ana_year_const).data)
    # determine the index of the ana_year_value, if several, assume the latest
    ana_arg = max(np.where(ano_obs_arr == ana_year_value)[0])
    ana_gsat = gsat_smooth_arr[ana_arg]

    # calculate stationary gev 
    obs_stat = cex.fit_gev(ano_obs_arr, returnValue=ana_year_value, getParams=True)
    orig_stat_rp = np.exp(obs_stat['logReturnPeriod'][0]) # it is the only value, we are making it a float

    # calculate non-stationary gev 
    obs_non_stat = cex.fit_gev(ano_obs_arr, gsat_smooth_arr, locationFun=1, returnValue=ana_year_value, getParams=True)
    orig_nonstat_rp = np.exp(obs_non_stat['logReturnPeriod'][ana_arg])

    bootstrap_rps = list()
    bootstrap_st_rps = list()
    rng = np.random.default_rng(501)

    for i in range(1000): 
        for_fit_indices = rng.integers(low=0, high=len(ano_obs_arr), size=len(ano_obs_arr))
        for_fit_indices.sort()
        temp_gev = cex.fit_gev(ano_obs_arr[for_fit_indices], gsat_smooth_arr[for_fit_indices],
                               locationFun=1, initial={'location':float(np.around(obs_stat['mle'][0],2)),
                                            'scale':float(np.around(obs_stat['mle'][1],2)), 
                                            'shape':float(np.around(obs_stat['mle'][2],2))}, getParams=True)
        temp_st_gev = cex.fit_gev(ano_obs_arr[for_fit_indices], initial={'location':float(np.around(obs_stat['mle'][0],2)),
                                            'scale':float(np.around(obs_stat['mle'][1],2)), 
                                            'shape':float(np.around(obs_stat['mle'][2],2))}, getParams=True)
        try:
            temp_loc = temp_gev['mle'][0] + temp_gev['mle'][1]*ana_gsat
            temp_rp = np.around(1/gev.sf(ana_year_value, -1*temp_gev['mle'][3],
                                        loc= temp_loc, scale=temp_gev['mle'][2]), 1)
            bootstrap_rps.append(temp_rp)
        except:
            bootstrap_rps.append(np.nan)
        try:
            temp_st_rp = np.around(np.exp(temp_st_gev['logReturnPeriod'][0]), 1)
            bootstrap_st_rps.append(temp_st_rp)
        except:
            bootstrap_st_rps.append(np.nan)
    
    bootstrap_rps = np.asarray(bootstrap_rps)
    bootstrap_st_rps = np.asarray(bootstrap_st_rps)

    rp_perc = np.nanpercentile(bootstrap_rps, [5,10,50,90,95]).round(1)
    rp_st_perc = np.nanpercentile(bootstrap_st_rps, [5,10,50,90,95]).round(1)

    era_csv = open(os.path.join(cfg['work_dir'], 'gev_era_data_'+cfg['region'].lower()+'_'+cfg['ax_var_label'].lower()+'.csv'), 'w', newline='')
    era_csv_writer = csv.writer(era_csv, delimiter=',')
    era_csv_writer.writerow([str(cfg['analysis_year'])+' ERA5 '+cfg['ax_var_label']+' value '+str(ana_year_value)+ ', ERA5 smoothed GSAT value '+str(ana_gsat)])
    era_csv_writer.writerow(['ERA5 non-stationary GEV params'])
    era_csv_writer.writerow(list(obs_non_stat['mle_names']))
    era_csv_writer.writerow(list(obs_non_stat['mle']))
    era_csv_writer.writerow([str(cfg['analysis_year'])+' ERA5 return period', str(np.around(orig_nonstat_rp,1))])
    era_csv_writer.writerow(['Bootstrapped uncertanties on ERA5 nonstationary return period'])
    era_csv_writer.writerow(['5_perc', '10_perc', '50_perc', '90_perc', '95_perc'])
    era_csv_writer.writerow(rp_perc)
    era_csv_writer.writerow(['ERA5 stationary GEV params'])
    era_csv_writer.writerow(list(obs_stat['mle_names']))
    era_csv_writer.writerow(list(obs_stat['mle']))
    era_csv_writer.writerow(['ERA5 stationary return period', str(np.around(orig_stat_rp,1))])
    era_csv_writer.writerow(['Bootstrapped uncertanties on ERA5 stationary return period'])
    era_csv_writer.writerow(['5_perc', '10_perc', '50_perc', '90_perc', '95_perc'])
    era_csv_writer.writerow(rp_st_perc)
    era_csv.close()

    obs_gev_data={'gev_param_names' : obs_non_stat['mle_names'],
                  'gev_param_values': obs_non_stat['mle'],
                  'abs_obs_cb': abs_obs_cb,
                  'ana_year_value': ana_year_value, 
                  'ana_gsat_value': ana_gsat, 
                  'ana_year_RP': orig_nonstat_rp,
                  'ana_year_RP_CI': rp_perc,
                #   'ana_year_RP': orig_stat_rp,
                #   'ana_year_RP_CI': rp_st_perc,
                  'ano_obs_cb': ano_obs_cb,
                  'smoothed_gsat':gsat_smooth_arr}

    return obs_gev_data


def get_era_wxx(cfg):

    wind_dir = cfg['ERA_wind_loc']
    aux_dir = cfg['auxiliary_data_dir']

    dates_csv = open(os.path.join(cfg['work_dir'], 'dates_max_wind.csv'), 'w', newline='')
    dates_csv_writer = csv.writer(dates_csv, delimiter=',')
    dates_csv_writer.writerow(['date', 'latitude', 'longitude', 'ws'])

    wind_file_list = sorted(os.listdir(wind_dir))

    era_ws_arr = np.zeros(len(wind_file_list)) 
    years = np.zeros(len(wind_file_list))
    for n, wind_f in enumerate(wind_file_list):
        wind_cb = iris.load_cube(os.path.join(wind_dir, wind_f))
        seas_cb = eprep.extract_season(wind_cb, season=cfg['season'])
        reg_cb = eprep.extract_shape(seas_cb, os.path.join(aux_dir,cfg['shape_region']), method='contains', crop=True)
        max_value = reg_cb.data.max()
        max_index = np.unravel_index(reg_cb.data.argmax(), reg_cb.shape)
        max_time = cftime.num2pydate(reg_cb.coord('time').points[max_index[0]], 
                                     reg_cb.coord('time').units.origin, 
                                     reg_cb.coord('time').units.calendar)
        max_lat = reg_cb.coord('latitude').points[max_index[1]]
        max_lon = reg_cb.coord('longitude').points[max_index[2]]
        dates_csv_writer.writerow([str(max_time), max_lat, max_lon, max_value])
        era_ws_arr[n] = max_value ; years[n] = max_time.year

    dates_csv.close()

    tims = [cftime.datetime(y, 8, 15, calendar='gregorian') for y in years]
    tim_dim = iris.coords.DimCoord(cftime.date2num(tims,'days since 1850-01-01', calendar='gregorian'), 
                                   standard_name='time', long_name='time', var_name='time',
                                    units=cf_units.Unit('days since 1850-01-01', calendar='gregorian'))
    tim_dim.guess_bounds()   

    era_big_cube = iris.cube.Cube(era_ws_arr, standard_name='wind_speed', long_name='10m wind speed', var_name='10m_ws', 
                                    units="m s-1", dim_coords_and_dims=[(tim_dim,0)])

    ano_cb = eprep.anomalies(era_big_cube, 'full',
                        reference={'start_year': cfg['reference_period'][0],
                        'start_month': 1, 'start_day':1,
                        'end_year': cfg['reference_period'][1],
                        'end_month': 12, 'end_day':31})

    return ano_cb, era_big_cube


def bootstrap_gev(data_dic): 

    if type(data_dic) == dict:  
        # determining the max length of the model realisation, the number of bootstrap
        # iterations is this value * 100
        max_cblst_len = np.asarray([len(data_dic[model]['data']) for model in data_dic.keys()]).max()
        iter_pool = max_cblst_len *100
        n_real = len(data_dic.keys())
        n_years = np.asarray([data_dic[model]['data'][0].shape for model in data_dic.keys()]).max()
        pool_size = int(np.around(n_real*n_years))
    elif type(data_dic) == iris.cube.CubeList: 
        iter_pool = np.max([len(data_dic) *100, 1000])
        n_years = data_dic[0].shape[0]
        n_real = 3 # this is a number of realisations which shown to be enough to draw conclusions
        if len(data_dic)<n_real: 
            pool_size = len(data_dic) * n_years
        else: 
            pool_size = n_real * n_years

    pool_data = list()
    if type(data_dic) == dict: 
        for model in data_dic.keys():
            for i in range(len(data_dic[model]['data'])):
                pool_data.append(data_dic[model]['data'][i].data)
    elif type(data_dic) == iris.cube.CubeList:
        for cb in data_dic: 
            pool_data.append(cb.data)
    pool_data = np.asarray(pool_data)

    shapes = np.zeros(iter_pool)
    locs = np.zeros(iter_pool)
    scales = np.zeros(iter_pool)

    rng = np.random.default_rng(501)

    for i in range(0, iter_pool):
        sample_data = list()
        mod_idxs = rng.choice(pool_data.shape[0], size = pool_size)
        for mod_idx in mod_idxs: 
            y_idx = rng.choice(n_years)
            sample_data.append(pool_data[mod_idx, y_idx])
        sample_data = np.asarray(sample_data).flatten()
        gev_params = cex.fit_gev(sample_data, getParams=True)
        try:
            shapes[i] = gev_params['mle'][2]; locs[i]=gev_params['mle'][0]
            scales[i] = gev_params['mle'][1]
        except: 
            shapes[i] =np.nan ; locs[i] = np.nan ; scales[i] = np.nan  

    param_dic = {'loc': locs, 'scale': scales, 'shape': shapes}    

    return param_dic


def make_uncert_figures(data_dic, cfg, border):

    colors = {}
    colors['wind_now'] = '#c57900'
    colors['wind_nat'] = '#005000'
    colors['wind_fut'] = '#8036A8'

    exp_list = list(data_dic.keys())  ; exp_list.remove('obs_info') 
    models = list(data_dic[exp_list[0]].keys())

    x_gev = np.arange(border[0],border[1], 0.1)

    tlocs = {'wind_now': 0.04 , 'wind_nat': 0.01,  'wind_fut': 0.08}
    uncert_band = {}
    for exp in exp_list: 
        uncert_band[exp] = {}

    gev_params = ['shape', 'loc', 'scale']

    for model in models: 

        # this is a figure where we will plot single distributions from bootstrap
        fig_single_bootstrap, ax_single_bootstrap = plt.subplots(1)
        fig_single_bootstrap.set_size_inches(12., 8.)

        # this a figure where we plot how the values for GEV params are distributed 
        fig_gev_distr, ax_gev_distr = plt.subplots(len(gev_params))
        fig_gev_distr.set_size_inches(8., 12.)

        for exp in exp_list: 
            gev_dic = data_dic[exp][model].pop('uncert')
            data_dic[exp][model] = data_dic[exp][model].pop('data')
            # here we plot distribution of single GEV params
            for n, gev_param in enumerate(gev_params):
                ax_gev_distr[n].hist(gev_dic[gev_param], bins=50, edgecolor='none',
                        facecolor = colors[exp], alpha=0.3, label = exp, density=True)
                ax_gev_distr[n].set_xlabel(gev_param)
                ax_gev_distr[n].set_ylabel('Number density')
                ax_gev_distr[n].set_title('GEV parameter: ' + gev_param)

            # this is an array where we'll throw all pdfs, to calculate 5/95 perc later 
            all_pdfs = np.zeros((len(x_gev), len(gev_dic['loc'])))
            all_sfs = np.zeros((len(x_gev), len(gev_dic['loc'])))

            # here we plot single distributions
            for i in range(len(gev_dic['loc'])):
                gev_pdf = gev.pdf(x_gev, -1*gev_dic['shape'][i],gev_dic['loc'][i], gev_dic['scale'][i])
                gev_sf = gev.sf(x_gev, -1*gev_dic['shape'][i],gev_dic['loc'][i], gev_dic['scale'][i])
                all_pdfs[:, i] = gev_pdf
                all_sfs[:, i] = gev_sf
                ax_single_bootstrap.plot(x_gev, gev_pdf, color = colors[exp], alpha=0.03)
            
            uncert_band[exp][model] = {'x_gev': x_gev, 'pdf_5th_perc' : np.nanpercentile(all_pdfs, 5, axis = 1), 
                                                'pdf_95th_perc': np.nanpercentile(all_pdfs, 95, axis = 1),
                                                'sf_5th_perc' : np.nanpercentile(all_sfs, 5, axis = 1), 
                                                'sf_95th_perc' : np.nanpercentile(all_sfs, 95, axis = 1),
                                                'loc_min': gev_dic['loc'].min(), 'loc_max': gev_dic['loc'].max(),
                                                'loc_5': np.nanpercentile(gev_dic['loc'], 5, interpolation='nearest'),
                                                'loc_95': np.nanpercentile(gev_dic['loc'], 95, interpolation='nearest'),
                                                'scale_min': gev_dic['scale'].min(), 'scale_max': gev_dic['scale'].max(),
                                                'scale_5': np.nanpercentile(gev_dic['scale'], 5, interpolation='nearest'),
                                                'scale_95': np.nanpercentile(gev_dic['scale'], 95, interpolation='nearest'),
                                                'shape_min': gev_dic['shape'].min(), 'shape_max': gev_dic['shape'].max(),
                                                'shape_5': np.nanpercentile(gev_dic['shape'], 5, interpolation='nearest'),
                                                'shape_95': np.nanpercentile(gev_dic['shape'], 95, interpolation='nearest'),
                                                'return_periods_all':1/all_sfs}
    
            param_str = '                             '+exp
            param_str += '\nshape mean:'+str(np.around(gev_dic['shape'].mean(),3))+', max:'+str(np.around(gev_dic['shape'].max(),3)) + ', min:' + str(np.around(gev_dic['shape'].min(),3)) 
            param_str += '\n loc mean:'+str(np.around(gev_dic['loc'].mean(),3))+', max:'+str(np.around(gev_dic['loc'].max(),3)) + ', min:' + str(np.around(gev_dic['loc'].min(),3)) \
                + '\nscale mean:'+ str(np.around(gev_dic['scale'].mean(),3))+', max:'+str(np.around(gev_dic['scale'].max(),3)) + ', min:' + str(np.around(gev_dic['scale'].min(),3))               
            
            ax_single_bootstrap.text(-0.95*border[0], tlocs[exp], param_str, color = colors[exp])

        ax_gev_distr[0].legend(loc=0, fancybox=False, frameon=False)
        fig_gev_distr.suptitle('Distribution of GEV parameters after bootstrap in ' +model,
                    fontsize = 'x-large')
        fig_gev_distr.set_dpi(250)
        plt.tight_layout()
        fig_gev_distr.savefig(os.path.join(cfg['plot_dir'], 'figure_'+cfg['region'].lower()+'_extreme_distr_param_gev_'+ model + diagtools.get_image_format(cfg)))
        plt.close(fig_gev_distr)

        ax_single_bootstrap.legend(loc=0, fancybox=False, frameon=False)
        ax_single_bootstrap.set_xlim(border[0], border[1])
        ax_single_bootstrap.set_ylim(0, ax_single_bootstrap.get_ylim()[1])
        ax_single_bootstrap.set_xlabel(cfg['ax_var_label']+' anomaly, '+ cfg['var_units'])
        ax_single_bootstrap.set_ylabel('Number density')

        fig_single_bootstrap.suptitle('Estimation of GEV fit uncertainty from '+model+' with Bootstrap method',
                    fontsize = 'x-large')
        fig_single_bootstrap.set_dpi(250)

        fig_single_bootstrap.savefig(os.path.join(cfg['plot_dir'], 'figure_'+cfg['region'].lower()+'_extremes_bootstrap_gev_'+ model + diagtools.get_image_format(cfg)))
        plt.close(fig_single_bootstrap)
                            
    return uncert_band


def make_hist_figure(data_dic, cfg, uncert_band, border):

    obs_info_arr = data_dic['obs_info']
    event = obs_info_arr['ana_year_value']

    era_event_prob = 1/obs_info_arr['ana_year_RP']
    era_event_RP = obs_info_arr['ana_year_RP']

    risk_csv = open(os.path.join(cfg['work_dir'], 'GEV_'+cfg['region']+'_'+cfg['ax_var_label'] +'_risk_data.csv'), 'w', newline='')
    risk_csv_writer = csv.writer(risk_csv, delimiter=',')
    risk_head_row = ['model']

    risk_uncert_csv = open(os.path.join(cfg['work_dir'], 'gev_'+cfg['region']+'_'+cfg['ax_var_label'] +'_uncert_risk_data.csv'), 'w', newline='')
    risk_uncert_csv_writer = csv.writer(risk_uncert_csv, delimiter=',')
    risk_uncert_head_row = ['model','all/nat_r_r_5', 'all/nat_r_r_10', 'all/nat_r_r_50','all/nat_r_r_90', 'all/nat_r_r_95']
    risk_uncert_csv_writer.writerow(risk_uncert_head_row)

    # exp_list = list(data_dic.keys()) ; exp_list.remove('obs_info')       
    # hard code for now to test the approach       
    exp_list = ['wind_now', 'wind_nat']

    quantile_measures = np.arange(0, 1.01, 0.01); quantile_measures[0] = 0.001

    colors = {}
    colors['wind_now'] = '#c57900'
    colors['wind_nat'] = '#005000'
    colors['wind_fut'] = '#8036A8'

    csv_file = open(os.path.join(cfg['work_dir'], 'gev_'+cfg['region']+'_'+cfg['ax_var_label'] +'_parameters.csv'), 'w', newline='')
    gevs_csv_writer = csv.writer(csv_file, delimiter=',')
    head_row = ['model']
    for exp_key in exp_list:
        head_row.extend(['shape_'+exp_key, 'shape_min_'+exp_key, 'shape_5_'+exp_key, 'shape_95_'+exp_key,'shape_max_'+exp_key])
        head_row.extend(['loc_'+exp_key,'loc_min_'+exp_key,'loc_5_'+exp_key,'loc_95_'+exp_key,'loc_max_'+exp_key])
        head_row.extend(['scale_'+exp_key,'scale_min_'+exp_key,'scale_5_'+exp_key,'scale_95_'+exp_key,'scale_max_'+exp_key])
        head_row.extend(['ks_stat_'+exp_key, 'ks_pvalue_'+exp_key, 'cvm_stat_'+exp_key, 'cvm_pvalue_'+exp_key])
        risk_head_row.append(exp_key+'_prob')
        risk_head_row.extend([exp_key+'_return_p', exp_key+'_return_p_5', exp_key+'_return_p_95', exp_key+'_return_p_5',exp_key+'_return_p_95'])
        risk_head_row.extend([exp_key+'_intens', exp_key+'_intens_5',exp_key+'_intens_95'])
    gevs_csv_writer.writerow(head_row)
    risk_csv_writer.writerow(risk_head_row)
    models = data_dic[exp_list[0]].keys()
    for model in models: 
        fig = plt.figure(constrained_layout=False)
        fig.set_size_inches(12., 8.)
        gs = fig.add_gridspec(nrows=2, ncols=6)
        ax_hist = fig.add_subplot(gs[:, :-2])
        ax_qq= fig.add_subplot(gs[0, -2:])
        ax_surv= fig.add_subplot(gs[1, -2:])
        model_row = [model]
        risk_model_row = [model]
        for exp_key in exp_list:
            ens_cubelist = data_dic[exp_key][model]
            distrib_data = []
            weights = []
            for cube in ens_cubelist:
                if (model == 'Multi-Model-Mean')&(cfg['model_weighting']):
                    cube_weight = cube.attributes['ensemble_weight']*cube.attributes['reverse_dtsts_n']
                else:
                    cube_weight = 1
                for point in cube.data:
                    distrib_data.append(np.around(point, 1))
                    weights.append(cube_weight/len(cube.data))
            distrib_data = np.asarray(distrib_data)
            weights = np.asarray(weights)
            if cfg['model_weighting']:
                un_weights = np.unique(weights)
                rev_un_weights = np.asarray(1/un_weights).round(0).astype('int32')
                large_denom = np.gcd.reduce(rev_un_weights)
                dev_weights = rev_un_weights/large_denom
                least_mult = np.lcm.reduce(dev_weights.astype('int32'))
                un_factors = least_mult/dev_weights
                factors = np.zeros(len(weights))
                for n_w, un_wght in enumerate(un_weights):
                    factors[np.where(weights==un_wght)] = un_factors[n_w]
                factors = factors.astype('int32')
                upd_distr_data = list()
                new_weights = list()
                for n_dp, distrib_point in enumerate(distrib_data): 
                    for f in range(factors[n_dp]):
                        upd_distr_data.append(distrib_point)
                        new_weights.append(weights[n_dp]/factors[n_dp])
                distrib_data = np.asarray(upd_distr_data)
                weights = np.asarray(new_weights)
            x_gev = uncert_band[exp_key][model]['x_gev']
            w_distr_par = cex.fit_gev(distrib_data, getParams=True)
            try:
                w_distr_loc = w_distr_par['mle'][0]; w_distr_scale = w_distr_par['mle'][1] 
                w_distr_shape = w_distr_par['mle'][2]
            except:
                w_distr_loc = np.nan ; w_distr_scale = np.nan ; w_distr_shape = np.nan
            w_pdf = gev.pdf(x_gev, -1*w_distr_shape, loc=w_distr_loc, scale=w_distr_scale)
            w_survival = gev.sf(x_gev, -1*w_distr_shape, loc=w_distr_loc, scale=w_distr_scale)
            theor_quants = gev(-1*w_distr_shape, loc=w_distr_loc, scale=w_distr_scale).ppf(quantile_measures)
            ks_res = kstest(distrib_data, gev(-1*w_distr_shape, loc=w_distr_loc, scale=w_distr_scale).cdf)
            cvm_res = cramervonmises(distrib_data, gev(-1*w_distr_shape, loc=w_distr_loc, scale=w_distr_scale).cdf)
            model_row.extend([w_distr_shape, uncert_band[exp_key][model]['shape_min'], uncert_band[exp_key][model]['shape_5']])
            model_row.extend([uncert_band[exp_key][model]['shape_95'], uncert_band[exp_key][model]['shape_max']])
            model_row.extend([w_distr_loc, uncert_band[exp_key][model]['loc_min'], uncert_band[exp_key][model]['loc_5']])
            model_row.extend([uncert_band[exp_key][model]['loc_95'], uncert_band[exp_key][model]['loc_max']])            
            model_row.extend([w_distr_scale, uncert_band[exp_key][model]['scale_min'], uncert_band[exp_key][model]['scale_5']])
            model_row.extend([uncert_band[exp_key][model]['scale_95'], uncert_band[exp_key][model]['scale_max']])            
            model_row.extend([ks_res.statistic, ks_res.pvalue, cvm_res.statistic, cvm_res.pvalue])
            n_bins = np.arange(border[0], border[1]+0.1, 1)
            ax_hist.hist(distrib_data, bins=n_bins, edgecolor=colors[exp_key],
                    facecolor = colors[exp_key], alpha=0.3, label=cfg['name_' + exp_key] , density=True, weights=weights, zorder = 2)
            ax_hist.plot(x_gev, w_pdf, c = colors[exp_key], ls = 'solid', label = 'GEV fit '+cfg['name_' + exp_key], zorder=3)
            pdf_perc_5 = uncert_band[exp_key][model]['pdf_5th_perc']
            pdf_perc_95 = uncert_band[exp_key][model]['pdf_95th_perc']
            sf_perc_5 = uncert_band[exp_key][model]['sf_5th_perc']
            sf_perc_95 = uncert_band[exp_key][model]['sf_95th_perc']
            ax_hist.fill_between(x_gev, pdf_perc_5, pdf_perc_95, color=colors[exp_key], alpha=0.3, linewidth=0, zorder=4)
            ax_surv.fill_between(x_gev, 1/sf_perc_5, 1/sf_perc_95, color=colors[exp_key], alpha=0.3, linewidth=0, zorder=4)
            if exp_key == 'wind_now':
                mock_event_idx = np.argmin(np.abs(w_survival-era_event_prob))
                mock_event_int = x_gev[mock_event_idx]
            intens = x_gev[np.argmin(np.abs(w_survival-era_event_prob))]
            intens_95 = x_gev[np.argmin(np.abs(sf_perc_95-era_event_prob))]
            intens_5 = x_gev[np.argmin(np.abs(sf_perc_5-era_event_prob))]
            mock_event_sf = gev.sf(mock_event_int, -1*w_distr_shape, loc=w_distr_loc, scale=w_distr_scale)
            mock_event_5p = np.nanpercentile(uncert_band[exp_key][model]['return_periods_all'][mock_event_idx,:],5)
            mock_event_95p = np.nanpercentile(uncert_band[exp_key][model]['return_periods_all'][mock_event_idx,:],95)
            try:
                risk_model_row.extend([mock_event_sf, 1/mock_event_sf, mock_event_5p, mock_event_95p, intens, intens_5, intens_95])
            except:
                risk_model_row.extend([np.nan, np.nan, mock_event_5p, mock_event_95p, intens, intens_5, intens_95])
            pract_quants = np.quantile(distrib_data, quantile_measures)
            ax_qq.scatter(theor_quants, pract_quants, edgecolors=colors[exp_key], marker='o', facecolors='None', lw=0.75, label=cfg['name_' + exp_key], zorder=3)
            ax_surv.plot(x_gev, 1/w_survival, color=colors[exp_key], zorder=3)
            # ax_hist.scatter(intens, 0, s=100)
        gevs_csv_writer.writerow(model_row)
        risk_csv_writer.writerow(risk_model_row)

        ylims = ax_hist.get_ylim()
        ax_hist.set_ylim(*ylims)

        # ax_hist.text(border[1]/2.5, ylims[1]*0.65,'  Number of\nrealisations ' +str(len(ens_cubelist)), fontsize='large')
        ax_hist.vlines(event, *ylims, color = 'indianred', linestyle = 'solid', lw=1.5, zorder=1, label = 'ERA5 ('+str(cfg['analysis_year'])+')')
        ax_hist.vlines(mock_event_int, *ylims, color = 'indianred', linestyle = 'dashed', lw=1.5, zorder=1, label = 'mock event')

        ax_surv.vlines(event, 0.1, era_event_RP, linestyle = 'solid', color='indianred', zorder=2,  label = 'ERA5 ('+str(cfg['analysis_year'])+')')
        ax_surv.vlines(mock_event_int, 0.1, era_event_RP, linestyle = 'dashed', color='indianred', zorder=2,  label = 'mock event')
        ax_surv.hlines(era_event_RP, x_gev[0], event,linestyle = 'solid', color='indianred', zorder=2)

        ax_hist.legend(loc=0, fancybox=False, frameon=False)
        ax_qq.legend(loc=2, fancybox=False, frameon=False, handletextpad=0.01)
        ax_qq.plot(border, border, c='tab:grey', zorder=1)
        ax_qq.grid(color='silver', axis='both', alpha=0.5)
        ax_hist.set_xlim(border[0]/1.2, border[1]/1.2)
        ax_qq.set_xlim(border[0]/1.2, border[1]/1.2)
        ax_qq.set_ylim(border[0]/1.2, border[1]/1.2)
        ax_surv.set_ylim(1,10000)
        ax_surv.set_xlim(border[0]/1.2, border[1]/1.2)
        ax_surv.set_yscale('log')
        ax_surv.grid(color='silver', axis='both', alpha=0.5)
        ax_hist.set_title('Probability density function')
        ax_qq.set_title('Quality assessment')
        ax_qq.set_ylabel('Data quantile')
        ax_qq.set_xlabel('GEV quantile')
        ax_surv.set_title('Return period')
        ax_surv.set_ylabel('years')
        ax_surv.set_xlabel(cfg['ax_var_label'] + ' anomaly, ' +cfg['var_units'])
        ax_hist.set_xlabel(cfg['ax_var_label'] + ' anomaly, ' +cfg['var_units'])
        ax_hist.set_ylabel('Number density')

        if model == 'Multi-Model-Mean':
            fig.suptitle(cfg['title_var_label']+' anomalies in '+cfg['region']+' relative to '+ str(cfg['reference_period'][0]) \
            + '-' + str(cfg['reference_period'][1])+ ' calculated from '+str(len(models) - 1) +' CMIP6 models' , fontsize = 'x-large')
        else:
            fig.suptitle(cfg['title_var_label']+' anomalies in '+cfg['region']+' relative to '+ str(cfg['reference_period'][0]) \
            + '-' + str(cfg['reference_period'][1])+ ' calculated from '+model, fontsize = 'x-large')
        fig.set_dpi(250)

        plt.tight_layout()
        
        fig.savefig(os.path.join(cfg['plot_dir'], 'figure_'+cfg['region'].lower()+'_extremes_gev_'+model + diagtools.get_image_format(cfg)))

        all_to_nat = list()
        for i in range(len(uncert_band['wind_now'][model]['return_periods_all'][mock_event_idx,:])):
                # we do RPs, so they are flipped
                all_to_nat.append(uncert_band['wind_nat'][model]['return_periods_all'][mock_event_idx,i]/uncert_band['wind_now'][model]['return_periods_all'][mock_event_idx,i])
        all_to_nat_perc = np.nanpercentile(all_to_nat, [5,10,50,90,95])
        risk_uncert_row = [model]
        risk_uncert_row.extend(all_to_nat_perc)
        risk_uncert_csv_writer.writerow(risk_uncert_row)
    
    risk_uncert_csv.close()

    return


def make_era_dist_figure(obs_info_dic, cfg, border):

    abs_cube = obs_info_dic['abs_obs_cb']
    ana_year_const = iris.Constraint(time = lambda cell: cell.point.year == cfg['analysis_year'])
    ana_year_value = abs_cube.extract(ana_year_const).data; ana_arg = np.max(np.where(abs_cube.data == ana_year_value)[0])

    border[0] = np.floor(abs_cube.data.min()*0.9)
    border[1] = np.ceil(abs_cube.data.max()*1.1)
    x_gev_fine = np.arange(border[0], border[1]+0.1, 0.1)
    n_bins = np.arange(border[0], border[1]+0.1, 1)

    era_param = cex.fit_gev(abs_cube.data, obs_info_dic['smoothed_gsat'], locationFun=1, returnValue=float(ana_year_value), getParams=True)
    era_RP = np.exp(era_param['logReturnPeriod'][ana_arg])
    era_loc = obs_info_dic['ana_gsat_value'] * era_param['mle'][1] + era_param['mle'][0] 
    era_scale = era_param['mle'][2] ; era_shape = era_param['mle'][3] 
    era_sf = gev.sf(x_gev_fine,  -1*era_shape, loc = era_loc, scale = era_scale)

    tims = cf_units.num2pydate(abs_cube.coord('time').points, abs_cube.coord('time').units.origin, abs_cube.coord('time').units.calendar)
    years = np.asarray([t.year for t in tims])

    fig_era = plt.figure(constrained_layout=False)
    fig_era.set_size_inches(12., 8.)
    gs = fig_era.add_gridspec(nrows=2, ncols=5)
    ax_era_hist = fig_era.add_subplot(gs[0, :-2])
    ax_era_surv= fig_era.add_subplot(gs[0, -2:])
    ax_era_tseries= fig_era.add_subplot(gs[1, :])

    ax_era_hist.hist(abs_cube.data, bins=n_bins, edgecolor='indianred', facecolor='indianred', alpha=0.3, density=True)
    ax_era_hist.scatter(ana_year_value, 0, c='indianred', marker='x',lw=2, s=100, label=str(cfg['analysis_year']), clip_on=False, zorder=4)
    ax_era_hist.set_xlim(*border)
    ax_era_hist.set_xlabel(cfg['ax_var_label'] +', ' + cfg['var_units'])
    ax_era_hist.set_ylabel('Number density')
    ax_era_hist.set_title('Probability density function')

    ax_era_surv.plot(x_gev_fine[era_sf>0.00001], 1/era_sf[era_sf>0.00001], c='indianred')
    ax_era_surv.set_title('ERA5 ' +  cfg['ax_var_label'] +' return period in '+ str(cfg['analysis_year']))
    ax_era_surv.set_xlabel(cfg['ax_var_label'] +', ' + cfg['var_units'])
    ax_era_surv.set_ylabel('years')
    ax_era_surv.scatter(ana_year_value, era_RP, c='indianred', marker='x',lw=2, s=100, label= str(cfg['analysis_year']), clip_on=False, zorder=4)
    ax_era_surv.text(border[1]*0.8, 5,'ERA5 return period\n    '+str(np.around(era_RP, 1))+ ' years', color='k')
    ax_era_surv.set_yscale('log')
    ax_era_surv.grid(color='silver', axis='both', alpha=0.5)
    ax_era_surv.set_ylim(1,10000)
    ax_era_surv.set_xlim(*border)

    ax_era_tseries.set_title('ERA5 ' + cfg['ax_var_label'] +' timeseries')
    ax_era_tseries.plot(years, abs_cube.data, c='indianred')
    ax_era_tseries.grid(color='silver', axis='both', alpha=0.5)
    ax_era_tseries.set_xlabel('time')
    ax_era_tseries.set_ylabel(cfg['ax_var_label']+', ' + cfg['var_units'])
    ax_era_tseries.scatter(years[np.argmax(abs_cube.data)], abs_cube.data.max(), edgecolors='indianred', marker='o', facecolors='None', s=100, lw=2)
    ax_era_tseries.set_xlim(years[0]-0.5, years[-1]+0.5)
    ax_era_tseries.set_ylim(*border)
    ax_era_tseries.arrow(years[-1] - len(years)*0.1, abs_cube.data.mean()*0.2 + 0.8*abs_cube.data.max(), len(years)*0.09 ,
                                     0.19*abs_cube.data.max()-abs_cube.data.mean()*0.2, color='k', length_includes_head=True,  
                                                                                        head_width=0.5, head_length=0.5)
    ax_era_tseries.text(years[-1] - len(years)*0.16, abs_cube.data.mean()*0.2 + 0.77*abs_cube.data.max(), 
                        'ERA5 '+ str(cfg['analysis_year'])+': '+str(np.around(abs_cube.data.max(),1))+' '+cfg['var_units'])

    fig_era.suptitle('ERA5 ' + cfg['title_var_label'] + ' in ' + cfg['region'] + ' and its GEV fit', fontsize = 'x-large')

    plt.tight_layout()

    fig_era.savefig(os.path.join(cfg['plot_dir'], 'figure_'+cfg['region'].lower()+'_'+cfg['ax_var_label'] +'_era' + diagtools.get_image_format(cfg)))


    return


def main(cfg):

    input_data = cfg['input_data']

    groups = group_metadata(input_data.values(), 'variable_group', sort=True)

    anomalies = groups.pop('wind_ano')

    wind_groups_l=[]
    for k in groups.keys():
        if 'wind' in k:
            wind_groups_l.append(k)
    
    # while the cmorizing is not implemented
    ano_era_cb, abs_era_cube = get_era_wxx(cfg)

    obs_info = obtain_obs_info(ano_era_cb, abs_era_cube, groups, cfg)

    mins = list(); maxs = list()

    plotting_dic = {}

    fit_param_apr = {}

    for group in wind_groups_l:
        plotting_dic[group] = {}
        group_data = groups[group]
        datasets = group_metadata(group_data, 'dataset')
        ens_cubelist = iris.cube.CubeList()
        for dataset in datasets.keys():
            filepaths = list(group_metadata(datasets[dataset], 'filename').keys())
            n_real = len(filepaths)
            mod_cubelist = iris.cube.CubeList()
            for filepath in filepaths:
                mod_cb = iris.load_cube(filepath)
                file_metadata = select_metadata(datasets[dataset], filename = filepath)
                ens = file_metadata[0]['ensemble']
                anom_cb = iris.load_cube(select_metadata(anomalies, dataset=dataset, ensemble=ens)[0]['filename'])
                mod_cb = mod_cb-anom_cb
                mins.append(mod_cb.collapsed('time', iris.analysis.MIN).data)
                maxs.append(mod_cb.collapsed('time', iris.analysis.MAX).data)
                mod_cb.attributes['ensemble_weight'] = 1 / n_real
                mod_cb.attributes['reverse_dtsts_n'] = 1/ len(datasets)
                ens_cubelist.append(mod_cb)
                mod_cubelist.append(mod_cb)
            plotting_dic[group][dataset] = {'data': mod_cubelist}
            plotting_dic[group][dataset]['uncert'] = bootstrap_gev(mod_cubelist)
        plotting_dic[group]['Multi-Model-Mean'] = {'data' : ens_cubelist}
        plotting_dic[group]['Multi-Model-Mean']['uncert'] = bootstrap_gev(plotting_dic[group]) 
        fit_param_apr[group] = {'loc': np.around(plotting_dic[group]['Multi-Model-Mean']['uncert']['loc'].mean(),3),
                                'scale': np.around(plotting_dic[group]['Multi-Model-Mean']['uncert']['scale'].mean(),3),
                                'shape': np.around(plotting_dic[group]['Multi-Model-Mean']['uncert']['shape'].mean(),3)}
    
    plotting_dic['obs_info'] = obs_info

    min_var = np.asarray(mins).min() ; max_var = np.asarray(maxs).max()  

    border = [np.floor(min_var*1.5), np.ceil(max_var*1.5)]

    st_file = eplot.get_path_to_mpl_style(cfg.get('mpl_style'))
    plt.style.use(st_file)

    uncert_band = make_uncert_figures(plotting_dic, cfg, border)

    make_hist_figure(plotting_dic, cfg, uncert_band, border)

    make_era_dist_figure(obs_info, cfg, border)

    logger.info('Success')


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        main(config)
