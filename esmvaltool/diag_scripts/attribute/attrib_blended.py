"""Perform attribution using blended and masked temperatures."""
import logging
import os
import xarray as xr
from pprint import pformat
from scipy.stats import t
import netCDF4

import iris

import numpy,time
import detatt_mk as da
import matplotlib
matplotlib.use('Agg') #Turn off interactive plots.
import matplotlib.pyplot as plt
import ncblendmask_esmval as ncbm

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
   ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


def attrib_warming(beta,betaCI,dec_warming,ci95_dec_warming,ci95_beta_obs=0.):
  attrib_warming=beta*dec_warming
  attrib_range=[attrib_warming-numpy.absolute(beta)*dec_warming*(((beta-betaCI[0])/beta)**2+(ci95_dec_warming/dec_warming)**2+(ci95_beta_obs/beta)**2)**0.5,attrib_warming+numpy.absolute(beta)*dec_warming*(((betaCI[1]-beta)/beta)**2+(ci95_dec_warming/dec_warming)**2+(ci95_beta_obs/beta)**2)**0.5]
  return (attrib_warming,numpy.sort(attrib_range))


def merge_noaa():
#Merge NOAAGlobalTemp with HadCRUT4 over the period 1850-1879, since NOAAGlobalTemp starts in 1880. Note that climatology period for HadCRUT4 is 1961-1990 whereas NOAAGlobalTemp is 1971-2000.
  hadcrut4=xr.open_dataset('/home/rng/data/esmvaltool/HadCRUT.4.6.0.0.median.nc')
  noaa=xr.open_dataset('/home/rng/data/esmvaltool/NOAAGlobalTemp_v5.0.0_gridded_s188001_e202002_c20200308T133325.nc')
  noaa=noaa.rename({'lon': 'longitude','lat': 'latitude','anom': 'temperature_anomaly'})
  noaa=noaa.squeeze(drop=True)
#Reorder longitudes to match HadCRUT4.
  index=numpy.arange(-36,36,1)
  index[0:36]=index[0:36]+72
  noaa.temperature_anomaly.values=noaa.temperature_anomaly.values[:,:,index]
  noaa.longitude.values=hadcrut4.longitude.values
  noaa=noaa.interp_like(hadcrut4) #Done purely to change start time to 1850.
  #Rebase HadCRUT4 relative to 1971-2000.
  had4_base=hadcrut4.temperature_anomaly.values[(1971-1850)*12:(2001-1850)*12,:,:]
  for m in range(12):
    norm = numpy.nanmean(had4_base[m::12,:,:],axis=0)
    hadcrut4.temperature_anomaly.values[m::12,:,:] = hadcrut4.temperature_anomaly.values[m::12,:,:] - norm
  #Concatenate HadCRUT4 1850-1879 and NOAA 1880 onwards.
  noaa.temperature_anomaly.values[0:(1880-1850)*12,:,:]=hadcrut4.temperature_anomaly.values[0:(1880-1850)*12,:,:]
  noaa=noaa.fillna(-1.0e30)
  #Write out to new file with xarray.
  noaa.to_netcdf('/home/rng/data/esmvaltool/NOAA_merged.nc')
  return()


def merge_giss():
#Merge GISSTemp with HadCRUT4 over the period 1850-1879, since GISSTemp starts in 1880.
  hadcrut4=xr.open_dataset('/home/rng/data/esmvaltool/HadCRUT.4.6.0.0.median.nc')
  giss=xr.open_dataset('/home/rng/data/esmvaltool/gistemp1200_GHCNv4_ERSSTv5.nc')
  giss=giss.rename({'lon': 'longitude','lat': 'latitude','tempanomaly': 'temperature_anomaly'})
  giss=giss.interp_like(hadcrut4)
  #Rebase HadCRUT4 relative to 1951-1980.
  had4_base=hadcrut4.temperature_anomaly.values[(1951-1850)*12:(1981-1850)*12,:,:]
  for m in range(12):
    norm = numpy.nanmean(had4_base[m::12,:,:],axis=0)
    hadcrut4.temperature_anomaly.values[m::12,:,:] = hadcrut4.temperature_anomaly.values[m::12,:,:] - norm
  #Concatenate HadCRUT4 1850-1879 and GISSTemp 1880 onwards.
  giss.temperature_anomaly.values[0:(1880-1850)*12,:,:]=hadcrut4.temperature_anomaly.values[0:(1880-1850)*12,:,:]
  giss=giss.fillna(-1.0e30)
  #Write out to new file with xarray.
  giss.to_netcdf('/home/rng/data/esmvaltool/GISS_merged.nc')
  return()
  
def main(cfg):
    plt.ioff() #Turn off interactive plotting.
    font = {'size'   : 5}
    matplotlib.rc('font', **font)
    mm_conv=0.0394
 
    """Carry out an attribution analysis using masked and blended temperature."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    obs=cfg['obs']
    exp_flag=cfg['exp_flag']
    diag_name=cfg['diag_name']
    rcplot=cfg['rcplot']
    simple_uncert=cfg['simple_uncert']
    ensobs=''
#    rcplot=False
#    simple_uncert=False #Flag to add simple uncertainty estimate in attributable warming to plots, considering only uncertainty in regression coeffs, for comparison.
    pool_int_var=True #Flag to pool internal variability estimates.
    if obs=='had5':
        obs_file='/home/rng/data/esmvaltool/HadCRUT.5.0.0.0.analysis.anomalies.ensemble_median.nc'
        ensobs='/home/rng/data/esmvaltool/HadCRUT.5.0.0.0.analysis.anomalies.' #If set apply multi-model analysis using ensemble obs data.
    elif obs=='had4':
        obs_file='/home/rng/data/esmvaltool/HadCRUT.4.6.0.0.median.nc'  #Updated to end of 2019.
        ensobs='/home/rng/data/esmvaltool/HadCRUT.4.6.0.0.anomalies.' #If set apply multi-model analysis using ensemble obs data.
    elif obs=='noaa':
        dummy=merge_noaa()
        obs_file='/home/rng/data/esmvaltool/NOAA_merged.nc'
    elif obs=='giss':
        dummy=merge_giss()
        obs_file='/home/rng/data/esmvaltool/GISS_merged.nc'
    else:
        exit ('Observations not recognised')
    sftlf_file='/home/rng/data/esmvaltool/CNRM-CM6-1-5x5-sftlf.nc'
    
#sftlf_file='/pf/b/b380746/CNRM-CM6-1-5x5-sftlf.nc' #Hard-coded path to sftlf file for CNRM-CM6 on a 5x5 grid. (Can't input through pre-processor at the moment. Update with sftlf for each model through preprocessor later.)
    
    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='ensemble')
    logger.info(
        "Group input data by model and sort by ensemble:"
        "\n%s", pformat(grouped_input_data))
    nmodel=len(grouped_input_data)
# Define colours
    colors=numpy.array([[0,0,0],[196,121,0],[178,178,178],[0,52,102],[0,79,0],[200,0,0],[0,200,0],[0,0,200],[112,160,205],])/256.
    shade_cols=numpy.array([[128,128,128,128],[204,174,113,128],[191,191,191,128],[67,147,195,128],[223,237,195,128],[255,150,150,128],[150,255,150,128],[150,150,255,128],[91,174,178,128]])/256.
    if exp_flag=='GHG':
      experiments=['historical-ssp245','hist-nat','hist-GHG','hist-aer','hist-CO2','hist-stratO3','hist-volc','hist-sol','hist-nat-ssp245-nat','hist-GHG-ssp245-GHG','hist-aer-ssp245-aer']
      label=['OTH','NAT','GHG']
      cols=colors[[3,4,2],:]
#      shade_cols=shade_cols[[3,4,2],:]
    else:
      label=['GHG','NAT','AER']
      cols=colors[[2,4,3],:]        
#      shade_cols=shade_cols[[2,4,3],:]        
      experiments=['historical-ssp245','hist-nat','hist-aer','hist-GHG','hist-CO2','hist-stratO3','hist-volc','hist-sol','hist-nat-ssp245-nat','hist-aer-ssp245-aer','hist-GHG-ssp245-GHG'] 
    nexp=len(experiments)-3 #-3 to account for repetition of hist-nat, hist-nat-ssp245-nat etc.
    print ('Number of experiments', nexp)
    # Loop over variables/datasets in alphabetical order
#   Define variables for D&A analysis
    av_yr=int(diag_name[4:6]) #Last two digits of diag_name are averaging period in yrs.
    years=list(numpy.arange(1850+av_yr/2,2020+av_yr/2,av_yr)) #Used for plotting.
    ldiag=int(170/av_yr) #length of diagnostic, assuming 170 years of data.
    if diag_name[0:4]=='hemi': ldiag=ldiag*2
   
    anom_max=500 #arbitrary max size for number of anomalies.
    nens_max=100
    mean_diag=numpy.zeros((ldiag,nexp,nmodel))
    mean_dec_warming=numpy.zeros((nexp,nmodel))
    ci95_dec_warming=numpy.zeros((nexp,nmodel))
    anom=numpy.zeros((ldiag,anom_max))
    anom_mod=numpy.zeros((anom_max))
    anom_index=0
    ens_sizes=numpy.zeros((nexp,nmodel))
    ensobs_diag=[]
    ensobs_dec_warming=[]
    ci95_beta_obs2=numpy.zeros((2))
    ci95_beta_obs3=numpy.zeros((3))
    
    for mm, dataset in enumerate(grouped_input_data):
#        if mm > 3: continue
        logger.info("*************** Processing model %s", dataset)
        grouped_model_input_data = group_metadata(
            grouped_input_data[dataset], 'exp', sort='ensemble')
        for exp in grouped_model_input_data:
            logger.info("***** Processing experiment %s", exp)
            exp_string = [experiments.index(i) for i in experiments if exp == i]
            experiment = exp_string[0]
            #Label hist-nat-ssp245-nat as hist-nat etc (some models' hist-nat ends in 2014 so is merged with ssp245-nat).
            if experiment > 7: experiment=experiment-7
            print ('*** Experiment',exp,'Index:',experiment)
            grouped_exp_input_data = group_metadata(
              grouped_model_input_data[exp], 'ensemble', sort='variable_group')
            nens=len(grouped_exp_input_data)
            ens_sizes[experiment,mm]=nens
            exp_diags=numpy.zeros((ldiag,nens))
            exp_dec_warming=numpy.zeros(nens)
        
            for ee, ensemble in enumerate(grouped_exp_input_data):
                logger.info("** Processing ensemble %s", ensemble)
                files=[]
                for attributes in grouped_exp_input_data[ensemble]:
                    logger.info("Processing variable %s", attributes['variable_group'])
                    files.append(attributes['filename'])
                logger.info("*************** Files for blend and mask %s", files)
                dec_warming=[]
                obs_dec_warming=[]
                (exp_diags[:,ee],had4_diag)=ncbm.ncblendmask_esmval('max', files[0],files[1],files[2],sftlf_file,obs_file,dec_warming,obs_dec_warming,0,0,diag_name,obs,ensobs,ensobs_diag,ensobs_dec_warming)
                ensobs='' #Set to empty string so that ensemble obs diagnostics are only calculated on the first iteration.
                exp_dec_warming[ee]=dec_warming[0]
#Calculate attributable warming in 2010-2019 vs 1850-1899 GMST.
            mean_diag[:,experiment,mm]=numpy.mean(exp_diags,axis=1)
            mean_dec_warming[experiment,mm]=numpy.mean(exp_dec_warming)
            if nens==1:
                ci95_dec_warming[experiment,mm]=ci95_dec_warming[numpy.nonzero(ci95_dec_warming*(ens_sizes**0.5))].mean() #use mean of CIs already calculated, corrected for ensemble size if ensemble size is 1.
            else:        
                ci95_dec_warming[experiment,mm]=(numpy.std(exp_dec_warming,ddof=1)/((nens)**0.5))*t.ppf(0.95,nens-1)
            if nens>1: anom[:,anom_index:anom_index+nens-1]=(exp_diags[:,0:nens-1]-mean_diag[:,experiment:experiment+1,mm])*((nens/(nens-1))**0.5) #Intra-ensemble anomalies for use as pseudo-control. Only use nens-1 ensemble members to ensure independence, and inflate variance by sqrt (nens / (nens-1)) to account for subtraction of ensemble mean.
            anom_mod[anom_index:anom_index+nens-1]=mm
            anom_index=anom_index+nens-1

    anom=anom[:,0:anom_index]
    anom_mod=anom_mod[0:anom_index]
    print ('anom_mod',anom_mod)

    if len(ensobs_diag) != 0:
    #Carry out multi-model analysis with each ensemble member of ensemble obs dataset if flag is set.
      enssize=len(ensobs_diag)
      beta2=numpy.zeros((2,enssize))
      beta3=numpy.zeros((3,enssize))
      neff=nmodel**2/numpy.sum(1./ens_sizes[0:3,:],axis=1) #Effective ensemble size when using multi-model mean.
      plt.figure(1,figsize=[180*mm_conv,60*mm_conv]) 
      for ee in range(enssize):
        (xr,yr,cn1,cn2)=da.reduce_dim(numpy.mean(mean_diag[:,0:3,:],axis=2),ensobs_diag[ee][:,None],anom[:,list(range(1,anom_index,2))],anom[:,list(range(0,anom_index,2))])
        att_out2=da.tls(xr[:,0:2],yr,cn1,ne=neff[[0,1]],cn2=cn2,flag_2S=1,RCT_flag=False)
        att_out3=da.tls(xr[:,0:3],yr,cn1,ne=neff[[0,1,2]],cn2=cn2,flag_3S=1,RCT_flag=False)
        beta2[:,ee]=numpy.squeeze(att_out2['beta'][:])
        beta3[:,ee]=numpy.squeeze(att_out3['beta'][:])
    #2-way regression coefficients.
        plt.subplot(121)
        plt.plot([ee+0.9,ee+0.9],numpy.transpose(att_out2['betaCI'][0,:]),color=colors[1,:],linewidth=0.5)
        plt.plot([ee+1.1,ee+1.1],numpy.transpose(att_out2['betaCI'][1,:]),color=colors[4,:],linewidth=0.5)
        plt.plot([ee+0.9],att_out2['beta'][0],color=colors[1,:],marker='+')
        plt.plot([ee+1.1],att_out2['beta'][1],color=colors[4,:],marker='+')
#3-way regression coefficients.
        plt.subplot(122)
        plt.plot([ee+0.8,ee+0.8],numpy.transpose(att_out3['betaCI'][2,:]),color=cols[2,:],linewidth=0.5)
        plt.plot([ee+1.0,ee+1.0],numpy.transpose(att_out3['betaCI'][1,:]),color=cols[1,:],linewidth=0.5)    
        plt.plot([ee+1.2,ee+1.2],numpy.transpose(att_out3['betaCI'][0,:]),color=cols[0,:],linewidth=0.5)
        plt.plot([ee+0.8],att_out3['beta'][2],color=cols[2,:],marker='+')
        plt.plot([ee+1.0],att_out3['beta'][1],color=cols[1,:],marker='+')
        plt.plot([ee+1.2],att_out3['beta'][0],color=cols[0,:],marker='+')

      plt.subplot(121)
      plt.plot([0,enssize+1],[1,1],color='black',linewidth=1,ls=':')
      plt.ylabel('Regression coefficients')#,size='x-small')
      plt.plot([0,enssize+1],[0,0],color='black',linewidth=1,ls='--')
      plt.axis([0,enssize+1,-1,3])
      plt.text (-10,3.3,'a',fontsize =7,fontweight='bold', va='center', ha='center')

      plt.subplot(122)
      plt.plot([0,enssize+1],[1,1],color='black',linewidth=1,ls=':')
      plt.plot([0,enssize+1],[0,0],color='black',linewidth=1,ls='--')
      plt.axis([0,enssize+1,-1,3])
      plt.text (-10,3.3,'b',fontsize =7,fontweight='bold', va='center', ha='center')
      plt.savefig('/home/rng/plots/esmvaltool/reg_obsens_'+diag_name+'_'+exp_flag+'_'+obs+'.eps')
      plt.close()

      for experiment in range(2):      
        ci95_beta_obs2[experiment]=numpy.std(beta2[experiment,:],ddof=1)*t.ppf(0.95,enssize-1)
      for experiment in range (3):
        ci95_beta_obs3[experiment]=numpy.std(beta3[experiment,:],ddof=1)*t.ppf(0.95,enssize-1)

      print ('beta2',beta2)
      print ('ci95_beta_obs2', ci95_beta_obs2)
      print ('ci95_beta_obs3', ci95_beta_obs3)

#Calculate uncertainty in obs warming.
#      ens_warming=numpy.zeros(enssize)
#      for ens in range(enssize):
#        ens_warming[ens]=numpy.squeeze(numpy.mean(ensobs_diag[ens][32:34])-numpy.mean(ensobs_diag[ens][0:10])) #Assume 5-yr means, calculate 2010-2019-1850-1899 in GMST.
#      ens_warming=numpy.sort(ens_warming)
#      print ('ens_warming',ens_warming)
#      print ('mean',numpy.mean(ens_warming))
 
      
    att_out={}
    att_out3={}
    model_names=[]
    model_indices=[]
    fig={}
    mm_attrib=0 #Counter over those models used for attribution.
    print ('Number of anomaly segments',anom_index)
    for mm, dataset in enumerate(grouped_input_data):
        if mean_diag[0,1,mm] == 0: #If there is no hist-nat simulation skip over model.
            continue
        model_names.append(dataset)
        model_indices.append(mm)
        print ('Model:',dataset)
        start=time.time()
        if pool_int_var:
          model_anom=anom
          nanom=anom_index
        else:
          model_anom=anom[:,(anom_mod==mm)] #Use only anomalies from one model.
          nanom=numpy.count_nonzero(anom_mod==mm)
        (xr,yr,cn1,cn2)=da.reduce_dim(mean_diag[:,[0,1,2],mm],had4_diag[:,None],model_anom[:,list(range(0,nanom,2))],model_anom[:,list(range(1,nanom,2))])
        time1 =time.time()
        print ('1***',time1-start)
        att_out[dataset]=da.tls(xr[:,0:2],yr,cn1,ne=ens_sizes[[0,1],mm],cn2=cn2,flag_2S=1,RCT_flag=rcplot)
        time2=time.time()
        print ('2***',time2-time1)
        att_out3[dataset]=da.tls(xr[:,0:3],yr,cn1,ne=ens_sizes[[0,1,2],mm],cn2=cn2,flag_3S=1,RCT_flag=rcplot)
# Replace beta confidence interval bounds with large positive and negative numbers if NaNs.
        print ('model, att_out3',dataset,att_out3[dataset])
        print ('model, att_out',dataset,att_out[dataset])
        att_out3[dataset]['betaCI'][numpy.isnan(att_out3[dataset]['betaCI'][:,0]),0]=-1000. 
        att_out3[dataset]['betaCI'][numpy.isnan(att_out3[dataset]['betaCI'][:,1]),1]=1000.
        att_out[dataset]['betaCI'][numpy.isnan(att_out[dataset]['betaCI'][:,0]),0]=-1000. 
        att_out[dataset]['betaCI'][numpy.isnan(att_out[dataset]['betaCI'][:,1]),1]=1000.
# For plotting purposes: If upper end of confidence range < best estimate, replace with large positive number. If lower end of confidence range > best estimate replace with large negative number. 
        for pat in range(3):
            if att_out3[dataset]['betaCI'][pat,0]>att_out3[dataset]['beta'][pat]:
                att_out3[dataset]['betaCI'][pat,0]=-1000.
            if att_out3[dataset]['betaCI'][pat,1]<att_out3[dataset]['beta'][pat]:
                att_out3[dataset]['betaCI'][pat,1]=1000.
        for pat in range(2):
            if att_out[dataset]['betaCI'][pat,0]>att_out[dataset]['beta'][pat]:
                att_out[dataset]['betaCI'][pat,0]=-1000.
            if att_out[dataset]['betaCI'][pat,1]<att_out[dataset]['beta'][pat]:
                att_out[dataset]['betaCI'][pat,1]=1000.
                
                
        time3=time.time()
        print ('3***',time3-time2)
        if mm_attrib == 0:
            ant='ANT'
            nat='NAT'
            ghg=exp_flag #Set to GHG or AER.
            oth='OTH'
        else:
            ant=None
            nat=None
            ghg=None
            oth=None
        if rcplot:
            topleft=321
            topright=322
            bottomleft=325
            bottomright=326
        else:        
            topleft=221
            topright=222
            bottomleft=223
            bottomright=224
#2-way regression coefficients.
        if rcplot:
            plt.figure(0,figsize=[180*mm_conv,180*mm_conv]) 
        else:
            plt.figure(0,figsize=[180*mm_conv,120*mm_conv]) 
        plt.subplot(topleft)
        plt.plot([mm_attrib+0.9,mm_attrib+0.9],numpy.transpose(att_out[dataset]['betaCI'][0,:]),color=colors[1,:],linewidth=2,label=ant)
        plt.plot([mm_attrib+1.1,mm_attrib+1.1],numpy.transpose(att_out[dataset]['betaCI'][1,:]),color=colors[4,:],linewidth=2,label=nat)
        plt.plot([mm_attrib+0.9],att_out[dataset]['beta'][0],color=colors[1,:],marker='+')
        plt.plot([mm_attrib+1.1],att_out[dataset]['beta'][1],color=colors[4,:],marker='+')
        if rcplot:
          plt.subplot(323)
          plt.bar([mm_attrib+1],[att_out[dataset]['rc_pvalue']],color='gray')
#2-way attributable warming.
        plt.subplot(bottomleft)
        [att_warming,att_warming_range]=attrib_warming(att_out[dataset]['beta'][0],att_out[dataset]['betaCI'][0,:],mean_dec_warming[0,mm]-mean_dec_warming[1,mm],(ci95_dec_warming[0,mm]**2+ci95_dec_warming[1,mm]**2)**0.5)
        plt.plot([mm_attrib+0.9,mm_attrib+0.9],att_warming_range,color=colors[1,:],linewidth=2,label=ant)
        plt.plot(mm_attrib+0.9,att_warming,color=colors[1,:],marker='+')
        [att_warming,att_warming_range]=attrib_warming(att_out[dataset]['beta'][1],att_out[dataset]['betaCI'][1,:],mean_dec_warming[1,mm],ci95_dec_warming[1,mm])
        plt.plot  ([mm_attrib+1.1,mm_attrib+1.1],att_warming_range,color=colors[4,:],linewidth=2,label=nat)
        plt.plot(mm_attrib+1.1,att_warming,color=colors[4,:],marker='+')
        if simple_uncert:
          plt.plot([mm_attrib+0.9,mm_attrib+0.9],att_out[dataset]['betaCI'][0,:]*(mean_dec_warming[0,mm]-mean_dec_warming[1,mm]),color=colors[1,:],marker='_',linestyle='None')
          plt.plot  ([mm_attrib+1.1,mm_attrib+1.1],att_out[dataset]['betaCI'][1,:]*mean_dec_warming[1,mm],color=colors[4,:],marker='_',linestyle='None')
#3-way regression coefficients.
        plt.subplot(topright)
        plt.plot([mm_attrib+0.8,mm_attrib+0.8],numpy.transpose(att_out3[dataset]['betaCI'][2,:]),color=cols[2,:],linewidth=2,label=label[2])
        plt.plot([mm_attrib+1.0,mm_attrib+1.0],numpy.transpose(att_out3[dataset]['betaCI'][1,:]),color=cols[1,:],linewidth=2,label=label[1])    
        plt.plot([mm_attrib+1.2,mm_attrib+1.2],numpy.transpose(att_out3[dataset]['betaCI'][0,:]),color=cols[0,:],linewidth=2,label=label[0])
        plt.plot([mm_attrib+0.8],att_out3[dataset]['beta'][2],color=cols[2,:],marker='+')
        plt.plot([mm_attrib+1.0],att_out3[dataset]['beta'][1],color=cols[1,:],marker='+')
        plt.plot([mm_attrib+1.2],att_out3[dataset]['beta'][0],color=cols[0,:],marker='+')

        if rcplot:
            plt.subplot(324)
            plt.bar([mm_attrib+1],[att_out3[dataset]['rc_pvalue']],color='gray')
#3-way attributable warming.
        plt.subplot(bottomright)
        [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][2],att_out3[dataset]['betaCI'][2,:],mean_dec_warming[2,mm],ci95_dec_warming[2,mm])
        plt.plot([mm_attrib+0.8,mm_attrib+0.8],att_warming_range,color=cols[2,:],linewidth=2,label=label[2])
        plt.plot(mm_attrib+0.8,att_warming,color=cols[2,:],marker='+')

        
        [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][1],att_out3[dataset]['betaCI'][1,:],mean_dec_warming[1,mm],ci95_dec_warming[1,mm])
        plt.plot([mm_attrib+1.0,mm_attrib+1.0],att_warming_range,color=cols[1,:],linewidth=2,label=label[1])
        plt.plot(mm_attrib+1.0,att_warming,color=cols[1,:],marker='+')

        
        [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][0],att_out3[dataset]['betaCI'][0,:],mean_dec_warming[0,mm]-mean_dec_warming[1,mm]-mean_dec_warming[2,mm],(ci95_dec_warming[0,mm]**2+ci95_dec_warming[1,mm]**2+ci95_dec_warming[2,mm]**2)**0.5)
        plt.plot([mm_attrib+1.2,mm_attrib+1.2],att_warming_range,color=cols[0,:],linewidth=2,label=label[0])
        plt.plot(mm_attrib+1.2,att_warming,color=cols[0,:],marker='+')
        if simple_uncert:
          plt.plot([mm_attrib+0.8,mm_attrib+0.8],att_out3[dataset]['betaCI'][2,:]*mean_dec_warming[2,mm],color=cols[2,:],marker='_',linestyle='None')
          plt.plot([mm_attrib+1.0,mm_attrib+1.0],att_out3[dataset]['betaCI'][1,:]*mean_dec_warming[1,mm],color=cols[1,:],marker='_',linestyle='None')
          plt.plot([mm_attrib+1.2,mm_attrib+1.2],att_out3[dataset]['betaCI'][0,:]*(mean_dec_warming[0,mm]-mean_dec_warming[1,mm]-mean_dec_warming[2,mm]),color=cols[0,:],marker='_',linestyle='None')
        
        mm_attrib=mm_attrib+1 #Counter over those models used for attribution.
        label=[None,None,None]
#Multi-model analysis.
    dataset='Multi'
    model_names.append(dataset)
#    anom=anom*1.5 #Scale internal variability.
    (xr,yr,cn1,cn2)=da.reduce_dim(numpy.mean(mean_diag[:,0:3,model_indices],axis=2),had4_diag[:,None],anom[:,list(range(1,anom_index,2))],anom[:,list(range(0,anom_index,2))])
    neff=mm_attrib**2/numpy.sum(1./ens_sizes[0:3,model_indices],axis=1) #Effective ensemble size when using multi-model mean.
    print ('*****attrib_blended ens_sizes',ens_sizes[0:3,model_indices])
    print ('*****attrib_blended neff',neff)
    att_out[dataset]=da.tls(xr[:,0:2],yr,cn1,ne=neff[[0,1]],cn2=cn2,flag_2S=1,RCT_flag=rcplot)
    att_out3[dataset]=da.tls(xr[:,0:3],yr,cn1,ne=neff[[0,1,2]],cn2=cn2,flag_3S=1,RCT_flag=rcplot)
    multim_mean_dec_warming=numpy.mean(mean_dec_warming[:,model_indices],axis=1)
#Compute uncertainty in attributable warming based on spread in ratio of GSAT to GMST warming across models.
    if diag_name[0:4]=='hemi':
      nlat=2
      nper=ldiag//2
      mean_diag=numpy.mean(numpy.reshape(mean_diag,(nlat,nper,nexp,nmodel)),axis=0)
    mean_dec_warming_gmst=numpy.squeeze(numpy.mean(mean_diag[32:34,:,:],axis=0)-numpy.mean(mean_diag[0:10,:,:],axis=0)) #Assume 5-yr means, calculate 2010-2019-1850-1899 in GMST.
    multim_mean_dec_warming_gmst=numpy.mean(mean_dec_warming_gmst[:,model_indices],axis=1)
    #Define uncertainty in multi-model mean warming to 2010-2019 based on spread in ratio of GSAT to GMST increase across models.
    print ('mean_dec_warming',mean_dec_warming)
    print ('mean_dec_warming_gmst',mean_dec_warming_gmst)
    print ('ratio',mean_dec_warming/mean_dec_warming_gmst)
    print (numpy.std(mean_dec_warming/mean_dec_warming_gmst,ddof=1,axis=1))
    print (numpy.mean(mean_dec_warming_gmst,axis=1))
    multim_ci95_dec_warming=numpy.std(mean_dec_warming/mean_dec_warming_gmst,ddof=1,axis=1)*t.ppf(0.95,nmodel-1)*numpy.mean(mean_dec_warming_gmst,axis=1)
#    multim_ci95_dec_warming=(numpy.mean(ci95_dec_warming[:,model_indices]**2,axis=1)/len(model_indices))**0.5
    print ('multim_ci95_dec_warming',multim_ci95_dec_warming)
    print ('Ratio of historical GSAT to GMST mean',numpy.mean(mean_dec_warming[0,:])/numpy.mean(mean_dec_warming_gmst[0,:]))

    if pool_int_var: #Only plot multi model results if pooling internal var.

      #2-way regression coefficients.
      plt.subplot(topleft)
      plt.plot([mm_attrib+0.9,mm_attrib+0.9],numpy.transpose(att_out[dataset]['betaCI'][0,:]),color=colors[1,:],linewidth=2,label=ant)
      plt.plot([mm_attrib+1.1,mm_attrib+1.1],numpy.transpose(att_out[dataset]['betaCI'][1,:]),color=colors[4,:],linewidth=2,label=nat)
      plt.plot([mm_attrib+0.9],att_out[dataset]['beta'][0],color=colors[1,:],marker='+')
      plt.plot([mm_attrib+1.1],att_out[dataset]['beta'][1],color=colors[4,:],marker='+')
      if rcplot:
          plt.subplot(323)
          plt.bar([mm_attrib+1],[att_out[dataset]['rc_pvalue']],color='gray')
  #2-way attributable warming.
      plt.subplot(bottomleft)
      print ('Two-way attributable warming')
      [att_warming,att_warming_range]=attrib_warming(att_out[dataset]['beta'][0],att_out[dataset]['betaCI'][0,:],multim_mean_dec_warming[0]-multim_mean_dec_warming[1],(multim_ci95_dec_warming[0]**2+multim_ci95_dec_warming[1]**2)**0.5,ci95_beta_obs2[0])
      plt.plot([mm_attrib+0.9,mm_attrib+0.9],att_warming_range,color=colors[1,:],linewidth=2,label=ant)
      plt.plot(mm_attrib+0.9,att_warming,color=colors[1,:],marker='+')
      print ('ANT:',att_warming_range)
      [att_warming,att_warming_range]=attrib_warming(att_out[dataset]['beta'][1],att_out[dataset]['betaCI'][1,:],multim_mean_dec_warming[1],multim_ci95_dec_warming[1],ci95_beta_obs2[1])
      plt.plot  ([mm_attrib+1.1,mm_attrib+1.1],att_warming_range,color=colors[4,:],linewidth=2,label=nat)
      plt.plot(mm_attrib+1.1,att_warming,color=colors[4,:],marker='+')
      print ('NAT:',att_warming_range)
      if simple_uncert:
        plt.plot([mm_attrib+0.9,mm_attrib+0.9],att_out[dataset]['betaCI'][0,:]*(multim_mean_dec_warming[0]-multim_mean_dec_warming[1]),color=colors[1,:],marker='_',linestyle='None')
        plt.plot  ([mm_attrib+1.1,mm_attrib+1.1],att_out[dataset]['betaCI'][1,:]*multim_mean_dec_warming[1],color=colors[4,:],marker='_',linestyle='None')


  #3-way regression coefficients.
      plt.subplot(topright)
      plt.plot([mm_attrib+0.8,mm_attrib+0.8],numpy.transpose(att_out3[dataset]['betaCI'][2,:]),color=cols[2,:],linewidth=2)
      plt.plot([mm_attrib+1.0,mm_attrib+1.0],numpy.transpose(att_out3[dataset]['betaCI'][1,:]),color=cols[1,:],linewidth=2)    
      plt.plot([mm_attrib+1.2,mm_attrib+1.2],numpy.transpose(att_out3[dataset]['betaCI'][0,:]),color=cols[0,:],linewidth=2)
      plt.plot([mm_attrib+0.8],att_out3[dataset]['beta'][2],color=cols[2,:],marker='+')
      plt.plot([mm_attrib+1.0],att_out3[dataset]['beta'][1],color=cols[1,:],marker='+')
      plt.plot([mm_attrib+1.2],att_out3[dataset]['beta'][0],color=cols[0,:],marker='+')

      if rcplot:
          plt.subplot(324)
          plt.bar([mm_attrib+1],[att_out3[dataset]['rc_pvalue']],color='gray')
  #3-way attributable warming.
      print ('Three-way attributable warming')
      plt.subplot(bottomright)
      [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][2],att_out3[dataset]['betaCI'][2,:],multim_mean_dec_warming[2],multim_ci95_dec_warming[2],ci95_beta_obs3[2])
      plt.plot([mm_attrib+0.8,mm_attrib+0.8],att_warming_range,color=cols[2,:],linewidth=2)
      plt.plot(mm_attrib+0.8,att_warming,color=cols[2,:],marker='+')
      print (exp_flag,att_warming_range)
      [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][1],att_out3[dataset]['betaCI'][1,:],multim_mean_dec_warming[1],multim_ci95_dec_warming[1],ci95_beta_obs3[1])
      plt.plot([mm_attrib+1.0,mm_attrib+1.0],att_warming_range,color=cols[1,:],linewidth=2)    
      plt.plot(mm_attrib+1.0,att_warming,color=cols[1,:],marker='+')
      print ('NAT:',att_warming_range)
      [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][0],att_out3[dataset]['betaCI'][0,:],multim_mean_dec_warming[0]-multim_mean_dec_warming[1]-multim_mean_dec_warming[2],(multim_ci95_dec_warming[0]**2+multim_ci95_dec_warming[1]**2+multim_ci95_dec_warming[2]**2)**0.5,ci95_beta_obs3[0])
      plt.plot([mm_attrib+1.2,mm_attrib+1.2],att_warming_range,color=cols[0,:],linewidth=2)
      plt.plot(mm_attrib+1.2,att_warming,color=cols[0,:],marker='+')
      print ('OTH:',att_warming_range)
      if simple_uncert:
        plt.plot([mm_attrib+0.8,mm_attrib+0.8],att_out3[dataset]['betaCI'][2,:]*multim_mean_dec_warming[2],color=cols[2,:],linestyle='None',marker='_')
        plt.plot([mm_attrib+1.0,mm_attrib+1.0],att_out3[dataset]['betaCI'][1,:]*multim_mean_dec_warming[1],color=cols[1,:],linestyle='None',marker='_')    
        plt.plot([mm_attrib+1.2,mm_attrib+1.2],att_out3[dataset]['betaCI'][0,:]*(multim_mean_dec_warming[0]-multim_mean_dec_warming[1]-multim_mean_dec_warming[2]),color=cols[0,:],linestyle='None',marker='_')

#Calculate attributable warming ranges in GMST.
      print ('Two-way attributable warming - GMST')
      [att_warming,att_warming_range]=attrib_warming(att_out[dataset]['beta'][0],att_out[dataset]['betaCI'][0,:],multim_mean_dec_warming_gmst[0]-multim_mean_dec_warming_gmst[1],0.,ci95_beta_obs2[0])
      print ('ANT:',att_warming_range)
      [att_warming,att_warming_range]=attrib_warming(att_out[dataset]['beta'][1],att_out[dataset]['betaCI'][1,:],multim_mean_dec_warming_gmst[1],0.,ci95_beta_obs2[1])
      print ('NAT:',att_warming_range)
      print ('Three-way attributable warming - GMST')
      [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][2],att_out3[dataset]['betaCI'][2,:],multim_mean_dec_warming_gmst[2],0.,ci95_beta_obs3[2])
      print (exp_flag,att_warming_range)
      [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][1],att_out3[dataset]['betaCI'][1,:],multim_mean_dec_warming_gmst[1],0.,ci95_beta_obs3[1])
      print ('NAT:',att_warming_range)
      [att_warming,att_warming_range]=attrib_warming(att_out3[dataset]['beta'][0],att_out3[dataset]['betaCI'][0,:],multim_mean_dec_warming_gmst[0]-multim_mean_dec_warming_gmst[1]-multim_mean_dec_warming_gmst[2],0.,ci95_beta_obs3[0])
      print ('OTH:',att_warming_range)
        
      
    nmodel_attrib=mm_attrib
      
    obs_warming=obs_dec_warming[0]*numpy.mean(mean_dec_warming[0,:])/numpy.mean(mean_dec_warming_gmst[0,:])
#    if obs=='had5':
#        obs_warming=obs_dec_warming[0]*1.1059 #Equivalent ratio for HadCRUT5
#    elif obs=='had4':
#        obs_warming=obs_dec_warming[0]*1.1629959046011968 #Scale up by ensemble mean GSAT/GMST ratio in CMIP6.
#    obs_warming=obs_dec_warming[0] #Set for plotting GMST warming.
#    elif obs=='noaa': 
#        obs_warming=obs_dec_warming[0]*10  #Fix this.
#    elif obs=='giss': 
#        obs_warming=obs_dec_warming[0]*10  #Fix this.
#    else:
#        exit('Obs not recognised')
    print ('Obs warming',obs_warming)
    print ('Obs warming in GMST',obs_dec_warming[0])
    print ('att_out3',att_out3)
    panel_labels=['a','b','c','d','e','f']
    panel_counter=0
    

    for ff in [topleft,topright]:
        plt.subplot(ff)
        plt.plot([0,nmodel_attrib+2],[1,1],color='black',linewidth=1,ls=':')
        if ff == topleft: plt.ylabel('Regression coefficients')#,size='x-small')
        plt.plot([0,nmodel_attrib+2],[0,0],color='black',linewidth=1,ls='--')
        if pool_int_var:
          plt.axis([0,nmodel_attrib+2,-1,3])
          plt.xticks(list(range(1,nmodel_attrib+2)),[''])
        else:
          plt.axis([0,nmodel_attrib+1,-1,3])
          plt.xticks(list(range(1,nmodel_attrib+1)),[''])
        plt.legend(loc="upper left")
        plt.text (-2.5,3.5,panel_labels[panel_counter],fontsize =7,fontweight='bold', va='center', ha='center')
        panel_counter=panel_counter+1

    if rcplot:
        for ff in [323,324]:
            plt.subplot(ff)
            plt.plot([0,nmodel_attrib+2],[0.05,0.05],color='black',linewidth=1,ls='--')
            plt.plot([0,nmodel_attrib+2],[0.95,0.95],color='black',linewidth=1,ls='--')
            if pool_int_var:
              plt.axis([0,nmodel_attrib+2,0,1])
              plt.xticks(list(range(1,nmodel_attrib+2)),[''])
            else:
              plt.axis([0,nmodel_attrib+1,0,1])
              plt.xticks(list(range(1,nmodel_attrib+1)),[''])
            if ff == 323: plt.ylabel('RCT P-value')#,size='x-small')
            plt.text (-2,1.075,panel_labels[panel_counter],fontsize =7,fontweight='bold', va='center', ha='center')
            panel_counter=panel_counter+1

    if simple_uncert: #Label lower panels a and b.
      panel_counter=0
      
    for ff in [bottomleft,bottomright]:
        plt.subplot(ff)
        plt.plot([0,nmodel_attrib+2],[obs_warming,obs_warming],color='black',linewidth=1,label='Had4 GSAT')
        if ff == bottomleft: plt.ylabel('Attributable change 2010-2019 vs 1850-1900 ($^\circ$C)')#,size='x-small') 
        plt.plot([0,nmodel_attrib+2],[0,0],color='black',linewidth=1,ls='--')
        if pool_int_var:
          plt.axis([0,nmodel_attrib+2,-2,3])
          plt.xticks(list(range(1,nmodel_attrib+2)),model_names,rotation=30.,ha="right")
        else:
          plt.axis([0,nmodel_attrib+1,-2,3])
          plt.xticks(list(range(1,nmodel_attrib+1)),model_names[0:nmodel_attrib],rotation=30.,ha="right")

        plt.text (-2,3.3,panel_labels[panel_counter],fontsize =7,fontweight='bold', va='center', ha='center')
        panel_counter=panel_counter+1


    pool_flag='' if pool_int_var else '_not_pooled'
    uncert_flag='__simple_uncert' if simple_uncert else ''
    plt.savefig('/home/rng/plots/esmvaltool/reg_attrib_'+diag_name+'_'+exp_flag+'_'+obs+pool_flag+uncert_flag+'.eps')
    plt.close()

#Debug plot.
    mod_cols=numpy.array([[0,73,73],[255,255,109],[0,146,146],[255,109,182],[255,182,119],[146,0,0],[0,109,219],[182,109,255],[109,182,255],[182,219,255],[73,0,146],[146,73,0],[219,209,0],[36,255,36]])/256.    
    if diag_name[0:4]=='hemi':
      had4_diag=numpy.mean(numpy.reshape(had4_diag,(nlat,nper)),axis=0)
      anom=numpy.mean(numpy.reshape(anom,(nlat,nper,anom_index)),axis=0)
    plt.plot(years,anom[:,0],color="Gray",label='Pseudo-control')
    for ann in range(1,anom_index):
        plt.plot(years,anom[:,ann],color="Gray")
    for mm, dataset in enumerate(grouped_input_data):
      if mm == 0 or mm == 7:
        plt.plot(years,mean_diag[:,0,mm],color=mod_cols[mm],linewidth=4,label=dataset+'_ALL')
        plt.plot(years,mean_diag[:,1,mm],color=mod_cols[mm],linewidth=4,ls='--')
#Multi-model mean ANT response
    ant=numpy.mean(mean_diag[:,0,:]-mean_diag[:,1,:],axis=1)*att_out['Multi']['beta'][0]
    plt.plot(years,ant,color="Red",linewidth=4,label='Scaled anthropogenic')
    plt.plot(years,had4_diag,color="Black",linewidth=4,label=obs)
    plt.legend(loc="upper left")
    plt.xlabel('Year')
    plt.ylabel('GMST anomaly ($^\circ$C)')
    plt.savefig('/home/rng/plots/esmvaltool/timeseries_ant_nat_'+diag_name+'_'+exp_flag+'_'+obs+'.png')
    plt.close()

    
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
