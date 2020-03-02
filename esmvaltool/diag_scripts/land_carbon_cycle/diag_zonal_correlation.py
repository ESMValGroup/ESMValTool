"""
Look at this module for guidance how to write your own.

Read the README_PERSONAL_DIAGNOSTIC file associated with this example;

Module for personal diagnostics (example).
Internal imports from exmvaltool work e.g.:

from esmvalcore.preprocessor import regrid
from esmvaltool.diag_scripts.shared.supermeans import get_supermean

Pipe output through logger;

Please consult the documentation for help with esmvaltool's functionalities
and best coding practices.
"""
# place your module imports here:

# operating system manipulations (e.g. path constructions)
import os

# to manipulate iris cubes
import iris
import matplotlib.pyplot as plt

import numpy as np
# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import get_diagnostic_filename, group_metadata, select_metadata, run_diagnostic, extract_variables
from esmvalcore.preprocessor import area_statistics


import extraUtils as xu

# esmvaltool.diag_scripts.shared.extract_variables

def load_variable(metadata, var_name):
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube

### user-defined functions for calculating partial correlations
def percentile_based_outlier(data, threshold=98):
    diff = (100 - threshold) / 2.0
    minval, maxval = np.percentile(data, [diff, 100 - diff])
    return (data < minval) | (data > maxval)


def pCorr(C):
    d1=C[:,0]
    d2=C[:,1]
    d3=C[:,2]
    minPoints=20
    correLation_method = 'pearson'
    if d1.size > minPoints:
        d1out=percentile_based_outlier(d1,threshold=98)
        d1[d1out]=np.nan
        d2out=percentile_based_outlier(d2,threshold=98)
        d2[d2out]=np.nan
        d3out=percentile_based_outlier(d3,threshold=98)
        d3[d3out]=np.nan
        d1,d2,d3=get_common_mask(d1,d2,d3)
        d1=np.ma.masked_invalid(d1).compressed().flatten()
        d2=np.ma.masked_invalid(d2).compressed().flatten()
        d3=np.ma.masked_invalid(d3).compressed().flatten()
        if correLation_method == 'pearson':
            r12,p = xu.calc_pearson_r(d1, d2,outlierPerc=0)
            r13,p = xu.calc_pearson_r(d1, d3,outlierPerc=0)
            r23,p = xu.calc_pearson_r(d2, d3,outlierPerc=0)
        elif correLation_method == 'spearman':
            r12,p = xu.calc_spearman_r(d1, d2,outlierPerc=2)
            r13,p = xu.calc_spearman_r(d1, d3,outlierPerc=2)
            r23,p = xu.calc_spearman_r(d2, d3,outlierPerc=2)
        else:
            print('choose a valid correLation_method [pearson/spearman]')
            exit
        r123=(r12-r13*r23)/np.sqrt((1-r13**2)*(1-r23**2))

    else:
        r123=np.nan
    return(r123)

def get_common_mask(_dat1,_dat2,_dat3):
    '''
    returns an array with 1 where all three data arrays have valid numeric values and zero elsewhere
    '''
    _dat1Mask=np.ma.getmask(np.ma.masked_invalid(_dat1))
    _dat2Mask=np.ma.getmask(np.ma.masked_invalid(_dat2))
    _dat3Mask=np.ma.getmask(np.ma.masked_invalid(_dat3))
    _valMaskA=1-(1-_dat1Mask)*(1-_dat2Mask)*(1-_dat3Mask)
    _valMask=np.ma.nonzero(_valMaskA)
    _dat1[_valMask]=np.nan
    _dat2[_valMask]=np.nan
    _dat3[_valMask]=np.nan
    _dat1=np.ma.masked_invalid(_dat1)
    _dat2=np.ma.masked_invalid(_dat2)
    _dat3=np.ma.masked_invalid(_dat3)
    return(_dat1,_dat2,_dat3)
def calc_zonal_correlation(_dat,_pr,_tas,bandsize=9):
    __dat=np.zeros((np.shape(_dat)[0],2))
    _dat,_pr,_tas=get_common_mask(_dat,_pr,_tas)
    bandsize=9
    for li in range(len(__dat)):
        istart=max(0,li-bandsize)
        iend=min(359,li+bandsize+1)
        _datZone=_dat[istart:iend,:]
        _prZone=_pr[istart:iend,:]
        _tasZone=_tas[istart:iend,:]
        _datZoneC=np.ma.masked_invalid(_datZone).compressed().flatten()
        _prZoneC=np.ma.masked_invalid(_prZone).compressed().flatten()
        _tasZoneC=np.ma.masked_invalid(_tasZone).compressed().flatten()
        pc_vpt=pCorr(np.column_stack((_datZoneC,_prZoneC,_tasZoneC)))
        pc_vtp=pCorr(np.column_stack((_datZoneC,_tasZoneC,_prZoneC)))
        __dat[li,0]=pc_vpt
        __dat[li,1]=pc_vtp
    return(__dat)
def _get_obs_data(cfg):
    """Get all data."""
    if not cfg.get('obs_files'):
        raise ValueError('The observation files needs to be specified in the '
                        'recipe (see recipe description for details)')
    else:
        input_files = [os.path.join(cfg['auxiliary_data_dir'], obs_file) for obs_file in cfg.get('obs_files')]
    all_data = {}
    varList = cfg.get('obs_variables')
    print (input_files,varList)
    for _var in varList:
        variable_constraint = iris.Constraint(cube_func=(lambda c: c.var_name == _var))
        cube = iris.load_cube(input_files,constraint=variable_constraint)
        print(cube)    
        all_data[_var]=cube.data
    return (all_data)
def fisher_z(_rdat):
    _zdat=0.5*(np.log(1+_rdat)-np.log(1-_rdat))
    return(_zdat)

def inverse_fisher_z(_zdat):
    _rdat=(np.exp(2*_zdat)-1)/(np.exp(2*_zdat)+1)
    return(_rdat)
def get_zonal_correlation(cfg):
    """
    A diagnostic function to calculate the total carbon stock from the list of variables
    passed to the diagnostic script.

    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    """
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    # import pdb; pdb.set_trace()

    # iterate over key(dataset) and values(list of vars)
    allModdat={}
    for key, value in my_files_dict.items():
        c_total = load_variable(value, 'ctotal')
        gpp = load_variable(value, 'gpp')
        ctau = (c_total/gpp) / (86400*365)
        pr      = load_variable(value, 'pr')
        tas     = load_variable(value, 'tas')
        zon_corr = calc_zonal_correlation(ctau.data, pr.data,tas.data)
        allModdat[key]=zon_corr
    allObsdat = _get_obs_data(cfg)
    _plot_zonal_correlation(allModdat,allObsdat,cfg)
    return 'I am done with my first ESMValTool diagnostic!'

def cln_ax(x_lab,ax_fs=8,axlw=0.4,rem_list=['top','right']):
    plt.xlim(-1,1)
    plt.ylim(-65,85)
    plt.axhline(y=0,lw=0.48,color='grey')
    plt.axvline(x=0,lw=0.48,color='grey')
    plt.xlabel(x_lab,fontsize=ax_fs)
    ax=plt.gca()
    for loc, spine in ax.spines.items():
        if loc in rem_list:
            spine.set_position(('outward',0)) # outward by 10 points
            spine.set_linewidth(0.)
        else:
            spine.set_linewidth(axlw)
    return
def draw_legend(ax_fs=8):
    leg = plt.legend(loc=(1.00974,.06),fontsize=ax_fs,ncol=1,columnspacing=0.05,fancybox=True,handlelength=0.8)
    leg.get_frame().set_linewidth(0)
    leg.get_frame().set_facecolor('#eeeeee')
    leg.legendPatch.set_alpha(0.45)
    texts = leg.get_texts()
    plt.setp(texts, fontsize=ax_fs*0.9)
    return(leg)

def _plot_zonal_correlation(all_mod_dat,all_obs_dat,cfg):
    """
    makes the line plots of zonal correlations from all models 

    Arguments:
        cfg - nested dictionary of metadata
        all_mod_dat - dictionary of correlations from all models
        all_obs_dat - dictionary of correlations and ranges from observation

    Returns:
        string; makes some time-series plots
    """
    models=list(all_mod_dat.keys())
    nmodels=len(models)
    lats=np.linspace(-89.75,89.75,360,endpoint=True)[::-1]
    ax_fs=8
    fig=plt.figure(figsize=(6,4))
    # tau-tas correlations
    sp1=plt.subplot(1,2,1)
    x_lab='$r_{\\tau-tas,pr}$'
    cln_ax(x_lab)
    # get the observations out of the dictionary
    r_tau_ctotal_tas=all_obs_dat['r_tau_ctotal_tas']
    r_tau_ctotal_tas_5=all_obs_dat['r_tau_ctotal_tas_5']
    r_tau_ctotal_tas_95=all_obs_dat['r_tau_ctotal_tas_95']
    # plot the correlations from observation
    sp1.plot(r_tau_ctotal_tas,lats,color='k',lw=1.1,label='Observation')
    sp1.fill_betweenx(lats, r_tau_ctotal_tas_5,r_tau_ctotal_tas_95, facecolor='grey',alpha=0.40)

    # tau-pr correlations
    sp2=plt.subplot(1,2,2)
    x_lab='$r_{\\tau-pr,tas}$'
    cln_ax(x_lab)

    # get the observations out of the dictionary
    r_tau_ctotal_pr=all_obs_dat['r_tau_ctotal_pr']
    r_tau_ctotal_pr_5=all_obs_dat['r_tau_ctotal_pr_5']
    r_tau_ctotal_pr_95=all_obs_dat['r_tau_ctotal_pr_95']

    # plot the correlations from observation
    sp2.plot(r_tau_ctotal_pr,lats,color='k',lw=1.1,label='Observation')
    sp2.fill_betweenx(lats, r_tau_ctotal_pr_5,r_tau_ctotal_pr_95, facecolor='grey',alpha=0.40)

    ######$$$$$PLOTTING for models
    # define arrays to store the zonal correlation of each model
    rAll_pt=np.zeros((len(lats),nmodels))
    rAll_tp=np.zeros((len(lats),nmodels))
    # loop over models and plot zonal correlations
    for row_m in range(nmodels):
        row_mod=models[row_m]
        mod_dat_row=all_mod_dat[row_mod]
        r_mod_vtp=np.flipud(mod_dat_row[:,1])
        rAll_tp[:,row_m]=r_mod_vtp
        sp1.plot(np.ma.masked_equal(r_mod_vtp,np.nan),lats,lw=0.3,label=row_mod)
        r_mod_vpt=np.flipud(mod_dat_row[:,0])
        rAll_pt[:,row_m]=r_mod_vpt
        sp2.plot(np.ma.masked_equal(r_mod_vpt,np.nan),lats,lw=0.3,label=row_mod)

    # get the normalized mean zonal correlation of all models for tau-tas,pr
    zAll_tp=fisher_z(rAll_tp) # do the fisher's transformation of r
    zAll_tp[np.isinf(zAll_tp)]=np.nan
    # mean of fisher's z
    zmm_ens=np.nanmean(zAll_tp,axis=1)
    zmm_ens_std=np.nanstd(zAll_tp,axis=1)
    r_mmod=inverse_fisher_z(zmm_ens)
    # get the uncertainty ranges
    r_mmod_std_low=inverse_fisher_z(zmm_ens-zmm_ens_std)
    r_mmod_std_hi=inverse_fisher_z(zmm_ens+zmm_ens_std)
    # plot the normalized mean and uncertainty
    sp1.plot(np.ma.masked_equal(r_mmod,np.nan),lats,color='blue',ls='--',lw=1,label='Normalized Mean r')
    sp1.fill_betweenx(lats, np.ma.masked_equal(r_mmod_std_low,np.nan),np.ma.masked_equal(r_mmod_std_hi,np.nan), facecolor='#42d4f4',alpha=0.25)

    # get the normalized mean zonal correlation of all models for tau-pr,tas
    zAll_pt=fisher_z(rAll_pt)
    zAll_pt[np.isinf(zAll_pt)]=np.nan
    zmm_ens=np.nanmean(zAll_pt,axis=1)
    zmm_ens_std=np.nanstd(zAll_pt,axis=1)
    r_mmod=inverse_fisher_z(zmm_ens)
    r_mmod_std_low=inverse_fisher_z(zmm_ens-zmm_ens_std)
    r_mmod_std_hi=inverse_fisher_z(zmm_ens+zmm_ens_std)

    sp2.plot(np.ma.masked_equal(r_mmod,np.nan),lats,color='blue',ls='--',lw=1,label='Normalized Mean r')
    sp2.fill_betweenx(lats, np.ma.masked_equal(r_mmod_std_low,np.nan),np.ma.masked_equal(r_mmod_std_hi,np.nan), facecolor='#42d4f4',alpha=0.25)

    plt.gca().yaxis.set_label_position("right")
    # draw the legend
    leg=draw_legend()
    # generate output path and save
    local_path = cfg['plot_dir']
    t_x=plt.figtext(0.5,0.5,' ',transform=plt.gca().transAxes)
    png_name = 'comparison_zonal_correlation_turnovertime_climate_Carvalhais2014.png'
    plt.savefig(os.path.join(local_path, png_name),bbox_inches='tight',bbox_extra_artists=[t_x,leg],dpi=450)
    plt.close()

    return 'Plotting complete'

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        get_zonal_correlation(config)
