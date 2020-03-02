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
# plotting functions
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.colors

# map library 
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy

import numpy as np
# import internal esmvaltool modules here
from esmvaltool.diag_scripts.shared import get_diagnostic_filename, group_metadata, select_metadata, run_diagnostic, extract_variables
from esmvalcore.preprocessor import area_statistics

import sys, os, os.path
# user-defined function
import extraUtils as xu

# esmvaltool.diag_scripts.shared.extract_variables

def load_variable(metadata, var_name):
    candidates = select_metadata(metadata, short_name=var_name)
    assert len(candidates) == 1
    filename = candidates[0]['filename']
    cube = iris.load_cube(filename)
    return cube
            
def rem_nan(tmp,_fill_val=-9999.):
    whereisNan=np.isnan(tmp)
    tmp[whereisNan]=_fill_val
    whereisNan=np.isinf(tmp)
    tmp[whereisNan]=_fill_val
    return(tmp)


def apply_mask(_dat,_mask_dat,fill_val=np.nan):
    mask_where=np.ma.getmask(np.ma.masked_less(_mask_dat,1.))
    _dat[mask_where]=fill_val
    return(_dat)

def _get_hex_data(_dat1,_dat2,valrange_sc):
    _dat1=np.ma.masked_outside(_dat1,valrange_sc[0],valrange_sc[1]).filled(np.nan)
    _dat2=np.ma.masked_outside(_dat2,valrange_sc[0],valrange_sc[1]).filled(np.nan)
    datMask1=np.ones(np.shape(_dat1))
    datMask2=np.ones(np.shape(_dat1))
    datMask1[_dat1==-9999.]=0
    datMask2[_dat2==-9999.]=0
    datMask=datMask1*datMask2
    _dat1m=apply_mask(_dat1,datMask)
    _dat2m=apply_mask(_dat2,datMask)
    _dat1mc=np.ma.masked_equal(_dat1m,np.nan).compressed()
    _dat2mc=np.ma.masked_equal(_dat2m,np.nan).compressed()
    return(_dat1mc,_dat2mc)

def _get_uncer_mask(_mmdat,_dat_5,_dat_95,_nmodels,fill_val=-9999.,nmodel_reject=2):
    wnan=np.isnan(_mmdat)
    _mmdat[wnan]=-9999
    mmdat_5=_mmdat-_dat_5
    mmdat_5c=np.ma.masked_less(np.ma.masked_greater_equal(mmdat_5,0).filled(1),0).filled(0).sum(0)
    mmdat_95=_mmdat-_dat_95
    mmdat_95c=(-1*np.ma.masked_greater(np.ma.masked_less_equal(mmdat_95,0).filled(-1),0).filled(0)).sum(0)
    mmdat_595=np.ones((2,360,720))
    mmdat_595[0]=mmdat_5c
    mmdat_595[1]=mmdat_95c
    mmInval=mmdat_595.max(0)-mmdat_595.min(0)
    nValid=_nmodels-mmInval-1
    uncer_mask=np.ma.masked_greater(np.ma.masked_less_equal(nValid,nmodel_reject).filled(0),nmodel_reject).filled(1)
    return(uncer_mask)
def get_ratio_colorbarInfo():
    cbInfo_r={}
    border=0.9
    ncolo=128
    # get the colormap
    _bounds_rat=np.concatenate((np.geomspace(0.2,0.25,num=int(ncolo)/4),np.geomspace(0.25,0.33,num=int(ncolo)/4),np.geomspace(0.33,0.5,num=int(ncolo)/4),np.geomspace(0.5,border,num=int(ncolo)/4),np.linspace(border,1/border,num=int(ncolo/4)),np.geomspace(1/border,2,num=int(ncolo)/4),np.geomspace(2,3,num=int(ncolo)/4),np.geomspace(3,4,num=int(ncolo)/4),np.geomspace(4,5,num=int(ncolo)/4)))
    colors1 = plt.cm.Blues(np.linspace(0.15, 0.998, ncolo))[::-1]
    colorsgr=np.tile(np.array([0.8,0.8,0.8,1]),int(ncolo/4)).reshape(int(ncolo/4),-1)
    colors2 = plt.cm.Reds(np.linspace(0.15, 0.998, ncolo))#[::-1]
    # combine them and build a new colormap
    colors1g = np.vstack((colors1, colorsgr))
    colors = np.vstack((colors1g,colors2))
    cm_rat = mpl.colors.LinearSegmentedColormap.from_list('my_colormap', colors)

    cbticks_rat=[0.2,0.25,0.33,0.5,0.9,1.1,2,3,4,5]
    cbticks_lab=['$\\dfrac{1}{5}$','$\\dfrac{1}{4}$','$\\dfrac{1}{3}$','$\\dfrac{1}{2}$','$\\dfrac{1}{1.1}$','$1.1$','$2$','$3$','$4$','$5$']
    cbInfo_r['_bounds_rat']=_bounds_rat
    cbInfo_r['cm_rat']=cm_rat
    cbInfo_r['cbticks_rat']=cbticks_rat
    cbInfo_r['cbticks_lab']=cbticks_lab
    return(cbInfo_r)
def get_dia_colorbarInfo():
    cbInfo_d={}
    cbName='plasma_r'
    _bounds_dia=np.concatenate(([1],np.linspace(8,16,num=10)[:-1],np.linspace(16,32,num=10)[:-1],np.linspace(32,64,num=10)[:-1],np.linspace(64,128,num=10)[:-1],np.linspace(128,256,num=10)[:-1],np.linspace(256,1000,num=2,endpoint=True)))
    cbticks=np.array([1,8,16,32,64,128,256])
    clist_=xu.get_colomap(cbName,_bounds_dia,lowp=0.,hip=1)
    cm_dia = mpl.colors.ListedColormap(clist_)
    cbInfo_d['_bounds_dia']=_bounds_dia
    cbInfo_d['cm_dia']=cm_dia
    cbInfo_d['cbticks']=cbticks
    return(cbInfo_d)


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
        all_data[_var]=np.roll(np.flipud(cube.data),360,axis=1)
    return (all_data)

def get_cTau(cfg):
    """
    A diagnostic function to calculate the ecosystem carbon turnover time from total carbon stock and gpp.

    Arguments:
        cfg - nested dictionary of metadata

    Returns:
        string; runs the user diagnostic

    """
    varNames=list(extract_variables(cfg).keys())
    my_files_dict = group_metadata(cfg['input_data'].values(), 'dataset')
    # import pdb; pdb.set_trace()

    allModdat={}
    for key, value in my_files_dict.items():
        c_total = load_variable(value, 'ctotal')
        gpp = load_variable(value, 'gpp')
        c_tau = (c_total/gpp) / (86400*365)
        c_tau.var_name = 'cTau'
        c_tau.standard_name = None 
        c_tau.long_name = 'Ecosystem Turnover time of carbon'
        c_tau.units = 'years'
        ofilename = get_diagnostic_filename(key+'_cTau', cfg)
        iris.save(c_tau, ofilename)
        allModdat[key]=c_tau.data
        _plot_single_map(c_tau.data,key,cfg)
    allObsdat = _get_obs_data(cfg)
    _plot_single_map(allObsdat['tau_ctotal'],'Carvalhais2014',cfg)
    _plot_multimodel_agreement(allModdat,allObsdat, cfg)
    _plot_matrix_map(allModdat,allObsdat, cfg)
class Map(dict):
    """
    Example:
    m = Map({'first_name': 'Eduardo'}, last_name='Pool', age=24, sports=['Soccer'])
    """
    def __init__(self, *args, **kwargs):
        super(Map, self).__init__(*args, **kwargs)
        for arg in args:
            if isinstance(arg, dict):
                for k, v in arg.iteritems():
                    self[k] = v

        if kwargs:
            for k, v in kwargs.iteritems():
                self[k] = v

    def __getattr__(self, attr):
        return self.get(attr)

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __setitem__(self, key, value):
        super(Map, self).__setitem__(key, value)
        self.__dict__.update({key: value})

    def __delattr__(self, item):
        self.__delitem__(item)

    def __delitem__(self, key):
        super(Map, self).__delitem__(key)
        del self.__dict__[key]
        
def _plot_single_map(_dat,_name,_cfg):
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
    fill_val=np.nan#-9999.
    ax_fs=7.1

    # colorbar sizes
    map_extent_dat = [-180,180, -60, 90]

    # get the colormap for diagonal maps
    cbInfo_d = get_dia_colorbarInfo()
    cbticks=cbInfo_d['cbticks']
    _bounds_dia=cbInfo_d['_bounds_dia']
    cm_dia=cbInfo_d['cm_dia']


    valrange_sc=(2,256)
    varName='$\\tau_{median}$'
    varUnit='years'
    ### move this to a subfunction when the configurations of the figure can be passed as a dict
    fig=plt.figure(figsize=(5,3))
    plot_dat=_dat
    _ax=plt.axes([0.1,0.1,0.9,0.9],projection=ccrs.Robinson(),frameon=False)
    _ax.set_global()
    plt.imshow(np.ma.masked_less(np.roll(np.flipud(plot_dat[60:,:]),360,axis=1),-9999.),extent=map_extent_dat,norm=matplotlib.colors.BoundaryNorm(_bounds_dia, len(_bounds_dia)),cmap=cm_dia,origin='upper',vmin=_bounds_dia[0],vmax=_bounds_dia[-1], transform=ccrs.PlateCarree())
    _ax.coastlines(linewidth=0.4,color='grey')
    plt.gca().outline_patch.set_visible(False)
    plt.title('Ecosystem Turnover Time (years), '+_name,fontsize=0.98*ax_fs)
    _axcol_dia=[0.15,0.06,0.8,0.035]
    cb=xu.mk_colo_tau(_axcol_dia,_bounds_dia,cm_dia,tick_locs=cbticks,cbfs=0.86*ax_fs,cbtitle='',cbrt=90)
    local_path = _cfg['plot_dir']
    png_name = 'global_turnovertime_' + _name + '.png'
    t_x=plt.figtext(0.5,0.5,' ',transform=plt.gca().transAxes)
    plt.savefig(os.path.join(local_path, png_name),bbox_inches='tight',bbox_extra_artists=[t_x],dpi=450)
    plt.close()
    return
def _plot_matrix_map(all_mod_dat,all_obs_dat,cfg):
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
    fcfg = Map()
    fill_val=np.nan#-9999.
    fcfg.fill_val=np.nan#-9999.
    models=list(all_mod_dat.keys())
    models.insert(0,'obs')
    all_mod_dat['obs']=all_obs_dat['tau_ctotal']
    nmodels=len(models)
    #FIGURES SETTINGS AND PARAMETER of the figure(preferably read from configurations file)
    x0=0.02
    y0=1.0
    wp=1./nmodels
    hp=wp
    xsp=0.0
    aspect_data=2160./4320.
    ysp=-0.03
    xsp_sca=wp/3*(aspect_data)
    ysp_sca=hp/3*(aspect_data)
    ax_fs=7.1

    # colorbar sizes
    hcolo=0.085*hp*nmodels/7
    wcolo=0.25
    cb_off_y=0.06158
    map_extent_dat = [-180,180, -60, 90]
    # define the data and information for plotting ratios 
    cbInfo_r = get_ratio_colorbarInfo()
    _bounds_rat = cbInfo_r['_bounds_rat']
    cm_rat = cbInfo_r['cm_rat']
    cbticks_rat = cbInfo_r['cbticks_rat']
    cbticks_lab = cbInfo_r['cbticks_lab']

    # get the colormap for diagonal maps
    cbInfo_d = get_dia_colorbarInfo()
    cbticks=cbInfo_d['cbticks']
    _bounds_dia=cbInfo_d['_bounds_dia']
    cm_dia=cbInfo_d['cm_dia']


    valrange_sc=(2,256)
    varName='$\\tau_{median}$'
    varUnit='years'
    
    fig=plt.figure(figsize=(9,6))
    for row_m in range(nmodels):
        row_mod=models[row_m]
        mod_dat_row=all_mod_dat[row_mod]
        for col_m in range(nmodels):
            col_mod=models[col_m]
            mod_dat_col=all_mod_dat[col_mod]
            print('---'+row_mod+' vs '+col_mod+'---')
            if row_m == col_m:

                _ax=plt.axes([x0+row_m*wp+row_m*xsp,y0-(col_m*hp+col_m*ysp),wp,hp],projection=ccrs.Robinson(),frameon=False)#,sharex=right,sharey=all)
                _ax.set_global()
                plot_dat=mod_dat_row
                gl_mean=np.nanmedian(np.ma.masked_less(mod_dat_row,-999.).filled(np.nan))
                plt.imshow(np.ma.masked_less(np.roll(np.flipud(plot_dat[60:,:]),360,axis=1),-9999.),extent=map_extent_dat,norm=matplotlib.colors.BoundaryNorm(_bounds_dia, len(_bounds_dia)),cmap=cm_dia,origin='upper',vmin=_bounds_dia[0],vmax=_bounds_dia[-1], transform=ccrs.PlateCarree())
                _ax.coastlines(linewidth=0.4,color='grey')
                plt.gca().outline_patch.set_visible(False)
                tit_str=varName+" = "+str(round(gl_mean,2))+' ('+varUnit+')'
            if row_m < col_m:
                _ax=plt.axes([x0+row_m*wp+row_m*xsp+xsp_sca,y0-(col_m*hp+col_m*ysp)+ysp_sca,wp*aspect_data,hp*aspect_data])#,sharex=right,sharey=all)
                xdat,ydat=mod_dat_col,mod_dat_row
                dat1h,dat2h=_get_hex_data(xdat,ydat,valrange_sc)
                _ax.hexbin(dat1h,dat2h,bins='log',mincnt=10, gridsize=40, cmap='viridis_r',linewidths=0)
                plt.ylim(valrange_sc[0],valrange_sc[1]*1.05)
                plt.xlim(valrange_sc[0],valrange_sc[1]*1.05)
                ymin,ymax = plt.ylim()
                xmin,xmax = plt.xlim()
                plt.plot((xmin,xmax),(ymin,ymax),'k',lw=0.1)
                r,p=xu.calc_spearman_r(xdat,ydat)
                tit_str="$R^2$="+str(round(r**2,2))
                plt.title(tit_str,fontsize=ax_fs*0.953,ma='left',y=1.175,va="top")
                print (tit_str)
                if row_m !=0 and col_m != nmodels-1:
                    xu.ax_clr(axfs=ax_fs)
                    xu.rotate_labels(which_ax='x',axfs=ax_fs,rot=90)
                elif row_m == 0 and col_m != nmodels-1:
                    xu.ax_clrX(axfs=ax_fs)
                    xu.rotate_labels(which_ax='x',axfs=ax_fs,rot=90)
                elif col_m == nmodels-1 and row_m !=0 :
                    xu.ax_clrY(axfs=ax_fs)
                    xu.rotate_labels(which_ax='x',axfs=ax_fs,rot=90)
                if row_m == 0 and col_m == nmodels-1:
                    xu.ax_orig(axfs=ax_fs)
                    xu.rotate_labels(which_ax='x',axfs=ax_fs,rot=90)
                    plt.ylabel('Column',fontsize=ax_fs)
                    plt.xlabel('Row',fontsize=ax_fs)
            if row_m > col_m:
                _ax=plt.axes([x0+row_m*wp+row_m*xsp,y0-(col_m*hp+col_m*ysp),wp,hp],projection=ccrs.Robinson(),frameon=False)#,sharex=right,sharey=all)
                _ax.set_global()
                plot_dat=rem_nan(mod_dat_row/mod_dat_col)
                _ax.imshow(np.ma.masked_equal(np.roll(np.flipud(plot_dat[60:,:]),360,axis=1),-9999.),extent=map_extent_dat,norm=matplotlib.colors.BoundaryNorm(_bounds_rat, len(_bounds_rat)),interpolation='none',vmin=_bounds_rat[0],vmax=_bounds_rat[-1],cmap=cm_rat,origin='upper', transform=ccrs.PlateCarree())
                _ax.coastlines(linewidth=0.4,color='grey')
                plt.gca().outline_patch.set_visible(False)
            if row_m == nmodels-1 or row_m == 0:
                if col_mod == 'obs':
                    _title_sp='observation'
                else:
                    _title_sp= col_mod
                plt.ylabel(_title_sp,fontsize=0.809*ax_fs)
                if row_m == nmodels-1:
                    plt.gca().yaxis.set_label_position("right")
            if col_m == 0 or col_m == nmodels-1:
                if row_mod == 'obs':
                    _title_sp='Carvalhais2014'
                else:
                    _title_sp=row_mod
                if col_m == 0:
                    plt.title(_title_sp,fontsize=0.809*ax_fs)
                else:
                    plt.xlabel(_title_sp,fontsize=0.809*ax_fs)

    t_x=plt.figtext(0.5,0.5,' ',transform=plt.gca().transAxes)
    x_colo=0.02 
    y_colo=y0+hp+cb_off_y
    _axcol_dia=[x_colo,y_colo,wcolo,hcolo]
    print (x_colo)
    cb_tit_d='Ecosystem Turnover Time (years)'
    cb=xu.mk_colo_tau(_axcol_dia,_bounds_dia,cm_dia,tick_locs=cbticks,cbfs=0.86*ax_fs,cbtitle=cb_tit_d,cbrt=90)
    x_colo=0.76
    y_colo=y0+hp+cb_off_y
    _axcol_rat=[x_colo,y_colo,wcolo,hcolo]
    cb=xu.mk_colo_cont(_axcol_rat,_bounds_rat,cm_rat,cbfs=0.7*ax_fs,cbrt=90,col_scale='log',cbtitle='Ratio (Column/Row)',tick_locs=cbticks_rat)
    cb.ax.set_xticklabels(cbticks_lab,fontsize=0.86*ax_fs,ha='center',rotation=0)
    local_path = cfg['plot_dir']
    png_name = 'global_comparison_matrix_models_Carvalhais2014.png'
    plt.savefig(os.path.join(local_path, png_name),bbox_inches='tight',bbox_extra_artists=[t_x],dpi=450)
    plt.close()

    return 'Plotting complete'
def _plot_multimodel_agreement(all_mod_dat,all_obs_dat,cfg):
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
    fcfg = Map()
    fill_val=np.nan#-9999.
    fcfg.fill_val=np.nan#-9999.
    models=list(all_mod_dat.keys())
    tau_obs=all_obs_dat['tau_ctotal']
    tau_obs_5=all_obs_dat['tau_ctotal_5']
    tau_obs_95=all_obs_dat['tau_ctotal_95']
    nmodels=len(models)
    #FIGURES SETTINGS AND PARAMETER of the figure(preferably read from configurations file)
    ax_fs=7.1
    cbInfo_r = get_ratio_colorbarInfo()
    _bounds_rat = cbInfo_r['_bounds_rat']
    cm_rat = cbInfo_r['cm_rat']
    cbticks_rat = cbInfo_r['cbticks_rat']
    cbticks_lab = cbInfo_r['cbticks_lab']
    
    # calculate the multimodel bias

    mm_full_tau=np.ones((nmodels,np.shape(tau_obs)[0],np.shape(tau_obs)[1]))*fill_val
    for row_m in range(nmodels):
        row_mod=models[row_m]
        mod_dat_row_tau=all_mod_dat[row_mod]
        mm_full_tau[row_m]=mod_dat_row_tau
    uncer_mask_tau=_get_uncer_mask(mm_full_tau,tau_obs_5,tau_obs_95,nmodels,nmodel_reject=int(nmodels/4),fill_val=np.nan)
    # uncer_mask_tau=np.zeros((360,720))+np.random.random((360,720))
    mm_tau=rem_nan(np.nanmedian(mm_full_tau,axis=0))
    mm_bias_tau=mm_tau/tau_obs
    mm_bias_tau=rem_nan(mm_bias_tau)
    mm_bias_tau[np.isnan(mm_bias_tau)]=0

    lats=np.linspace(-89.75,89.75,360,endpoint=True)[::-1]
    lons=np.linspace(-179.75,179.75,720,endpoint=True)

    latint=lats[1]-lats[0]
    lonint=lons[1]-lons[0]
    map_extent_dat = [-180,180, -90, 90]
    # 
    # x_scat,y_scat=_mp(*np.meshgrid(lons, lats))

    ### move this to a subfunction when the configurations of the figure can be passed as a dict
    fig=plt.figure(figsize=(5,3))
    plot_dat=mm_bias_tau
    _ax=plt.axes([0.1,0.1,0.9,0.9],projection=ccrs.Robinson(),frameon=False)#,sharex=right,sharey=all)
    x,y=np.meshgrid(lons-lonint/2, lats-latint/2)
    _ax.set_global()
    _ax.imshow(np.ma.masked_equal(np.roll(np.flipud(plot_dat[:,:]),360,axis=1),-9999.),extent=map_extent_dat,norm=matplotlib.colors.BoundaryNorm(_bounds_rat, len(_bounds_rat)),interpolation='none',vmin=_bounds_rat[0],vmax=_bounds_rat[-1],cmap=cm_rat,origin='upper', transform=ccrs.PlateCarree())
    # _ax.pcolor(x, y, np.ma.masked_not_equal(uncer_mask_tau,0), alpha=0.,hatch='//////',linewidth=2, transform=ccrs.PlateCarree())
    _ax.pcolormesh(x, y, np.ma.masked_not_equal(np.roll(np.flipud(uncer_mask_tau[:,:]),360,axis=1),0), alpha=0,hatch='..',linewidth=1.9, transform=ccrs.PlateCarree())
    # _ax.contourf(x, y, np.ma.masked_not_equal(np.roll(np.flipud(uncer_mask_tau[:,:]),360,axis=1),0), alpha=0.,hatches='......',linewidth=0.2,extent=map_extent_dat, transform=ccrs.PlateCarree())
    _ax.coastlines(linewidth=0.4,color='grey')
    plt.gca().outline_patch.set_visible(False)
    plt.title('Bias of multimodel median and model agreement',fontsize=0.98*ax_fs)
    _axcol_rat=[0.15,0.06,0.8,0.035]
    cb=xu.mk_colo_cont(_axcol_rat,_bounds_rat,cm_rat,cbfs=0.7*ax_fs,cbrt=90,col_scale='log',cbtitle='',tick_locs=cbticks_rat)
    cb.ax.set_xticklabels(cbticks_lab,fontsize=0.86*ax_fs,ha='center',rotation=0)
    local_path = cfg['plot_dir']
    png_name = 'global_multimodelAgreement_turnovertime_Carvalhais2014.png'
    t_x=plt.figtext(0.5,0.5,' ',transform=plt.gca().transAxes)
    plt.savefig(os.path.join(local_path, png_name),bbox_inches='tight',bbox_extra_artists=[t_x],dpi=450)
    plt.close()

    return 'Plotting complete'


if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        get_cTau(config)
