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


# esmvaltool.diag_scripts.shared.extract_variables

def load_variable(metadata, var_name):
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube
            

def _zonal_tau(_datgpp,_datcTotal):
    __dat=np.zeros((np.shape(_datgpp)[0]))
    minPoints=20
    bandsize=9
    for li in range(len(__dat)):
        istart=max(0,li-bandsize)
        iend=min(359,li+bandsize+1)
        _datgppZone=_datgpp[istart:iend,:]#* _areaZone
        _datcTotalZone=_datcTotal[istart:iend,:]#* _areaZone
        # _datgppZone=remove_tail_percentiles(_datgppZone,_outlierPerc=10)
        # _datcTotalZone=remove_tail_percentiles(_datcTotalZone,_outlierPerc=10)
        nValids=np.nansum(1-np.ma.getmask(np.ma.masked_invalid(_datgppZone)))
        if nValids > minPoints:
            __dat[li] = (np.nansum(_datcTotalZone)/np.nansum(_datgppZone)) / (86400*365)
        else:
            __dat[li] = np.nan
            print(nValids,'valid points for tau, setting nan')
    return(__dat)
def _get_zonal_tau(cfg):
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
        # load the cube from data files only
        # using a single variable here so just grab the first (and only)
        # list element

        c_total = load_variable(value, 'ctotal')
        gpp = load_variable(value, 'gpp')
        print(gpp,c_total)

        c_tau = _zonal_tau(gpp.data, c_total.data)
        print(key,value)
        allModdat[key]=c_tau
    allObsdat = _get_obs_data(cfg)
    _plot_zonal_tau(allModdat,allObsdat, cfg)
    return 'I am done with my first ESMValTool diagnostic!'

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

def rem_axLine(rem_list=['top','right'],axlw=0.4):
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
def _plot_zonal_tau(all_mod_dat,all_obs_dat,cfg):
    """
    makes the maps of variables 

    Arguments:
        cfg - nested dictionary of metadata
        cube - the cube to plot
        dataset - name of the dataset to plot

    Returns:
        string; makes some time-series plots

    Note: this function is private; remove the '_'
    so you can make it public.
    """
    models=list(all_mod_dat.keys())
    nmodels=len(models)
    lats=np.linspace(-89.75,89.75,360,endpoint=True)[::-1]
    tau_obs=all_obs_dat['tau_ctotal']
    tau_obs_5=all_obs_dat['tau_ctotal_5']
    tau_obs_95=all_obs_dat['tau_ctotal_95']
    plt.figure(figsize=(3,5))
    sp0=plt.subplot(1,1,1)
    sp0.plot(tau_obs,lats,color='k',lw=1.5,label='Observation')
    sp0.fill_betweenx(lats, tau_obs_5,tau_obs_95, facecolor='grey',alpha=0.40)
    ax_fs=8
    for row_m in range(nmodels):
        row_mod=models[row_m]
        dat_mod_tau=all_mod_dat[row_mod]
        sp0.plot(np.ma.masked_equal(np.flipud(dat_mod_tau),np.nan),lats,lw=0.5,label=row_mod)


    # dat_mod_tau=mm_ens_zonal_tau
    # sp0.plot(np.ma.masked_equal(dat_mod_tau,np.nan),lats,color='blue',lw=lwMainLine,label='Model Ensemble')
    leg=draw_legend()
    valrange_md=(1,1000)
    plt.gca().set_xscale('log')
    plt.xlim(valrange_md[0],valrange_md[1])
    plt.xlim(2,1000)
    plt.ylim(-65,85)
    plt.axhline(y=0,lw=0.48,color='grey')
    x_lab='$\\tau$'

    plt.xlabel(x_lab,fontsize=ax_fs)
    plt.ylabel('Latitude ($^\\circ N$)',fontsize=ax_fs,ma='center')
    rem_axLine(['top','right'])
    local_path = cfg['plot_dir']
    t_x=plt.figtext(0.5,0.5,' ',transform=plt.gca().transAxes)
    png_name = 'comparison_zonal_turnovertime_Carvalhais2014.png'
    plt.savefig(os.path.join(local_path, png_name),bbox_inches='tight',bbox_extra_artists=[t_x,leg],dpi=450)
    plt.close()

    return 'Plotting complete'


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        _get_zonal_tau(config)
