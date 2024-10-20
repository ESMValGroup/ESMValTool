"""Perform attribution using blended and masked temperatures."""
"""Plots Figure 2, and Extended Data Figures 2-6 and 8 in Gillett et al.."""
import logging
import os
import csv
import xarray as xr
import iris
from pprint import pformat
from scipy.stats import t

import numpy
import esmvaltool.diag_scripts.attribute.detatt_mk as da
import matplotlib
matplotlib.use('Agg') #Turn off interactive plots.
import matplotlib.pyplot as plt
import esmvaltool.diag_scripts.attribute.ncblendmask_esmval as ncbm

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic, select_metadata)


logger = logging.getLogger(os.path.basename(__file__))


def attrib_warming(beta,betaCI,dec_warming,ci90_dec_warming,ci90_beta_obs=0.):
  #Calculate attributable warming and confidence range, considering uncertainties in beta, dec_warming and obs.
  attrib_warming=beta*dec_warming
  attrib_range=[attrib_warming-numpy.absolute(beta)*dec_warming*(((beta-betaCI[0])/beta)**2+(ci90_dec_warming/dec_warming)**2+(ci90_beta_obs/beta)**2)**0.5,
                attrib_warming+numpy.absolute(beta)*dec_warming*(((betaCI[1]-beta)/beta)**2+(ci90_dec_warming/dec_warming)**2+(ci90_beta_obs/beta)**2)**0.5]
  perc_17_83=[attrib_warming-0.58*numpy.absolute(beta)*dec_warming*(((beta-betaCI[0])/beta)**2+(ci90_dec_warming/dec_warming)**2+(ci90_beta_obs/beta)**2)**0.5,
              attrib_warming+0.58*numpy.absolute(beta)*dec_warming*(((betaCI[1]-beta)/beta)**2+(ci90_dec_warming/dec_warming)**2+(ci90_beta_obs/beta)**2)**0.5] #Calculate 17th and 83rd percentiles, assuming a Gaussian.
  attrib_range.sort()
  perc_17_83.sort()
  return (attrib_warming,attrib_range,perc_17_83)
  
def main(cfg):
    plt.ioff() #Turn off interactive plotting.
    font = {'size'   : 5}
    matplotlib.rc('font', **font)
    mm_conv=0.0394
 
    """Carry out an attribution analysis using masked and blended temperature."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    plot_dir=cfg['plot_dir']
    output_file_type=cfg['output_file_type']
    exp_flag=cfg['exp_flag']
    rcplot=cfg['rcplot']
    simple_uncert=cfg['simple_uncert']
    year_block = cfg.get('year_block')
    warming_base = cfg.get('warming_base')
    warming_years = cfg.get('warming_years')
    pool_int_var=True #Flag to pool internal variability estimates.


    obs_file = select_metadata(input_data, variable_group='obs_mean')[0]['filename']
    obs_cb = iris.load_cube(obs_file)
    y_start = obs_cb.coord('time').cell(0).point.year
    y_end = obs_cb.coord('time').cell(-1).point.year

    hadlabel=select_metadata(input_data, variable_group='obs_mean')[0]['dataset']

    ens_obs_files = [m['filename'] for m in select_metadata(input_data, variable_group='obs_ens')]
    ens_obs_cblst = iris.load(ens_obs_files)

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='ensemble')
    grouped_input_data.pop(hadlabel)
    logger.info(
        "Group input data by model and sort by ensemble:"
        "\n%s", pformat(grouped_input_data))
    nmodel=len(grouped_input_data)
    
# Define colours
    colors=numpy.array([[0,0,0],[196,121,0],[178,178,178],[0,52,102],[0,79,0],[200,0,0],[0,200,0],[0,0,200],[112,160,205],])/256.
    if exp_flag=='GHG':
        experiments=['historical-ssp245','hist-nat','hist-GHG','hist-aer','hist-stratO3','hist-nat-ssp245-nat','hist-GHG-ssp245-GHG','hist-aer-ssp245-aer','hist-stratO3-ssp245-stratO3']
        label=['OTH','NAT','GHG']
        cols=colors[[3,4,2],:]
    else:
        label=['GHG','NAT','AER']
        cols=colors[[2,4,3],:]        
        experiments=['historical-ssp245','hist-nat','hist-aer','hist-GHG','hist-stratO3','hist-nat-ssp245-nat','hist-aer-ssp245-aer','hist-GHG-ssp245-GHG','hist-stratO3-ssp245-stratO3'] 
    nexp=len(experiments)-4 #-4 to account for repetition of hist-nat, hist-nat-ssp245-nat etc.
    print ('Number of experiments', nexp)

    #   Define variables for D&A analysis
    av_yr=year_block #Last two digits of diag_name are averaging period in yrs.
    years=list(numpy.arange(y_start+av_yr/2,y_end+1+av_yr/2,av_yr)) #Used for plotting.
    #Added 5 years here.
    nyear=y_end - y_start +1 #Number of years, hard-coded.
    ldiag=int(nyear/av_yr) #length of diagnostic.
   
    anom_max=500 #arbitrary max size for number of anomalies.
    mean_diag=numpy.zeros((ldiag,nexp,nmodel))
    mean_dec_warming=numpy.zeros((nexp,nmodel))
    ci90_dec_warming=numpy.zeros((nexp,nmodel))
    anom=numpy.zeros((ldiag,anom_max))
    anom_mod=numpy.zeros((anom_max))
    anom_index=0
    ens_sizes=numpy.zeros((nexp,nmodel))
    ensobs_diag=[]
    ensobs_dec_warming=[]
    ci90_beta_obs2=numpy.zeros((2))
    ci90_beta_obs3=numpy.zeros((3))
    mean_gmst_comp_warming=numpy.zeros((nyear,nexp,nmodel))
    ci90_gmst_comp_warming=numpy.zeros((nyear,nexp,nmodel))

    #Loop over models, then datasets, then ensemble members.
    for mm, dataset in enumerate(grouped_input_data):
        logger.info("*************** Processing model %s", dataset)
        grouped_model_input_data = group_metadata(
            grouped_input_data[dataset], 'exp', sort='ensemble')
        for exp in grouped_model_input_data:
            logger.info("***** Processing experiment %s", exp)
            exp_string = [experiments.index(i) for i in experiments if exp == i]
            experiment = exp_string[0]
            #Label hist-nat-ssp245-nat as hist-nat etc (some models' hist-nat ends in 2014 so is merged with ssp245-nat).
            if experiment > 4:
               experiment=experiment-4
            print ('*** Experiment',exp,'Index:',experiment)
            grouped_exp_input_data = group_metadata(
              grouped_model_input_data[exp], 'ensemble', sort='variable_group')
            nens=len(grouped_exp_input_data)
            ens_sizes[experiment,mm]=nens
            exp_diags=numpy.zeros((ldiag,nens))
            exp_dec_warming=numpy.zeros(nens)
            exp_gmst_comp_warming=numpy.zeros((nyear,nens))
       
            for ee, ensemble in enumerate(grouped_exp_input_data):
                logger.info("** Processing ensemble %s", ensemble)
                files=[]
                for attributes in grouped_exp_input_data[ensemble]:
                    logger.info("Processing variable %s", attributes['variable_group'])
                    files.append(attributes['filename'])
                logger.info("*************** Files for blend and mask %s", files)
                dec_warming=[]
                obs_dec_warming=[]
                gmst_comp_warming=[]
                
                #Calculate the diagnostic used for attribution from an individual simulation.
                (exp_diags[:,ee],had4_diag)=ncbm.ncblendmask_esmval(files[0],obs_cb,dec_warming,obs_dec_warming,0,gmst_comp_warming,year_block,ens_obs_cblst,ensobs_diag,ensobs_dec_warming, warming_years=warming_years, warming_base=warming_base)
                ens_obs_cblst='' #Set to empty string so that ensemble obs diagnostics are only calculated on the first iteration.
                exp_dec_warming[ee]=dec_warming[0]
                exp_gmst_comp_warming[:,ee]=gmst_comp_warming[0]
            #Average diagnostic and warming in 2010-2019 vs 1850-1899 GSAT over ensemble members.
            if mm==1 and experiment ==1: #Select GISS nat simulations.
                plt.figure(2,figsize=[180*mm_conv,60*mm_conv]) 
                for ee in range(nens):
                  plt.plot(years,exp_diags[:,ee],color='black')
                plt.savefig(plot_dir+'/giss_nat.pdf')
                plt.close()

            mean_diag[:,experiment,mm]=numpy.mean(exp_diags,axis=1)
            mean_dec_warming[experiment,mm]=numpy.mean(exp_dec_warming)
            mean_gmst_comp_warming[:,experiment,mm]=numpy.mean(exp_gmst_comp_warming,axis=1)
            if nens==1:
                ci90_dec_warming[experiment,mm]=ci90_dec_warming[numpy.nonzero(ci90_dec_warming*(ens_sizes**0.5))].mean() #use mean of CIs already calculated, corrected for ensemble size if ensemble size is 1.
            else:        
                ci90_dec_warming[experiment,mm]=(numpy.std(exp_dec_warming,ddof=1)/((nens)**0.5))*t.ppf(0.95,nens-1)
                ci90_gmst_comp_warming[:,experiment,mm]=(numpy.std(exp_gmst_comp_warming,axis=1,ddof=1)/((nens)**0.5))*t.ppf(0.95,nens-1)
            if nens>1: 
                anom[:,anom_index:anom_index+nens-1]=(exp_diags[:,0:nens-1]-mean_diag[:,experiment:experiment+1,mm])*((nens/(nens-1))**0.5) #Intra-ensemble anomalies for use as pseudo-control. Only use nens-1 ensemble members to ensure independence, and inflate variance by sqrt (nens / (nens-1)) to account for subtraction of ensemble mean.
            anom_mod[anom_index:anom_index+nens-1]=mm
            anom_index=anom_index+nens-1

    anom=anom[:,0:anom_index]
    anom_mod=anom_mod[0:anom_index]

    if len(ensobs_diag) != 0:
    #Carry out multi-model analysis with each ensemble member of ensemble obs dataset if flag is set.
      enssize=len(ensobs_diag)
      beta2=numpy.zeros((2,enssize))
      beta3=numpy.zeros((3,enssize))
      neff=nmodel**2/numpy.sum(1./ens_sizes[0:3,:],axis=1) #Effective ensemble size when using multi-model mean.
      #Plot Extended Data Fig 6.
      plt.figure(1,figsize=[180*mm_conv,60*mm_conv]) 
      for ee in range(enssize):
        #Apply attribution analysis.
        (xr,yr,cn1,cn2)=da.reduce_dim(numpy.mean(mean_diag[:,0:3,:],axis=2),ensobs_diag[ee][:,None],
                                      anom[:,list(range(1,anom_index,2))],anom[:,list(range(0,anom_index,2))])
        yr = yr.data
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
      plt.savefig(plot_dir+'/reg_obsens_'+exp_flag+'.'+output_file_type)
      plt.close()

      #Calculate 5-95% range for betas based on obs uncertainties to use in multi-model analysis.
      for experiment in range(2):      
        ci90_beta_obs2[experiment]=numpy.std(beta2[experiment,:],ddof=1)*t.ppf(0.95,enssize-1)
      for experiment in range (3):
        ci90_beta_obs3[experiment]=numpy.std(beta3[experiment,:],ddof=1)*t.ppf(0.95,enssize-1)
        
    #Set up main figure.
    if rcplot:
      plt.figure(0,figsize=[180*mm_conv,180*mm_conv]) 
    else:
      plt.figure(0,figsize=[180*mm_conv,120*mm_conv]) 

    #Main attribution analysis.      
    att_out={}
    att_out3={}
    model_names=[]
    model_indices=[]
    mm_attrib=0 #Counter over those models used for attribution.
    for mm, dataset in enumerate(grouped_input_data):
        if mean_diag[0,1,mm] == 0: #If there is no hist-nat simulation skip over model.
            continue
        model_names.append(dataset)
        model_indices.append(mm)
        print ('Model:',dataset)
        if pool_int_var:
            model_anom=anom
            nanom=anom_index
        else:
            model_anom=anom[:,(anom_mod==mm)] #Use only anomalies from one model.
            nanom=numpy.count_nonzero(anom_mod==mm)
        #Apply attribution analysis.
        (xr,yr,cn1,cn2)=da.reduce_dim(mean_diag[:,[0,1,2],mm],had4_diag[:,None],
                                      model_anom[:,list(range(0,nanom,2))],model_anom[:,list(range(1,nanom,2))])
        yr = yr.data
        att_out[dataset]=da.tls(xr[:,0:2],yr,cn1,ne=ens_sizes[[0,1],mm],cn2=cn2,flag_2S=1,RCT_flag=rcplot)
        att_out3[dataset]=da.tls(xr[:,0:3],yr,cn1,ne=ens_sizes[[0,1,2],mm],cn2=cn2,flag_3S=1,RCT_flag=rcplot)
# Replace beta confidence interval bounds with large positive and negative numbers if NaNs.
        print ('model, att_out3',dataset,att_out3[dataset])
        print ('model, att_out',dataset,att_out[dataset])
        att_out3[dataset]['betaCI'][numpy.isnan(att_out3[dataset]['betaCI'][:,0]),0]=-1000. 
        att_out3[dataset]['betaCI'][numpy.isnan(att_out3[dataset]['betaCI'][:,1]),1]=1000.
        att_out[dataset]['betaCI'][numpy.isnan(att_out[dataset]['betaCI'][:,0]),0]=-1000. 
        att_out[dataset]['betaCI'][numpy.isnan(att_out[dataset]['betaCI'][:,1]),1]=1000.
# For plotting purposes: If upper end of confidence range < best estimate, replace with large positive number. 
# If lower end of confidence range > best estimate replace with large negative number. 
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
        [att_warming,att_warming_range,dummy]=attrib_warming(att_out[dataset]['beta'][0],att_out[dataset]['betaCI'][0,:],
                                                             mean_dec_warming[0,mm]-mean_dec_warming[1,mm],
                                                             (ci90_dec_warming[0,mm]**2+ci90_dec_warming[1,mm]**2)**0.5)
        plt.plot([mm_attrib+0.9,mm_attrib+0.9],att_warming_range,color=colors[1,:],linewidth=2,label=ant)
        plt.plot(mm_attrib+0.9,att_warming,color=colors[1,:],marker='+')
        [att_warming,att_warming_range,dummy]=attrib_warming(att_out[dataset]['beta'][1],att_out[dataset]['betaCI'][1,:],
                                                             mean_dec_warming[1,mm],ci90_dec_warming[1,mm])
        plt.plot  ([mm_attrib+1.1,mm_attrib+1.1],att_warming_range,color=colors[4,:],linewidth=2,label=nat)
        plt.plot(mm_attrib+1.1,att_warming,color=colors[4,:],marker='+')
        if simple_uncert:
          plt.plot([mm_attrib+0.9,mm_attrib+0.9],att_out[dataset]['betaCI'][0,:]*(mean_dec_warming[0,mm]-mean_dec_warming[1,mm]),
                   color=colors[1,:],marker='_',linestyle='None')
          plt.plot  ([mm_attrib+1.1,mm_attrib+1.1],att_out[dataset]['betaCI'][1,:]*mean_dec_warming[1,mm],
                     color=colors[4,:],marker='_',linestyle='None')
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
        [att_warming,att_warming_range,dummy]=attrib_warming(att_out3[dataset]['beta'][2],
                                                             att_out3[dataset]['betaCI'][2,:],
                                                             mean_dec_warming[2,mm],ci90_dec_warming[2,mm])
        plt.plot([mm_attrib+0.8,mm_attrib+0.8],att_warming_range,color=cols[2,:],linewidth=2,label=label[2])
        plt.plot(mm_attrib+0.8,att_warming,color=cols[2,:],marker='+')

        
        [att_warming,att_warming_range,dummy]=attrib_warming(att_out3[dataset]['beta'][1],
                                                             att_out3[dataset]['betaCI'][1,:],
                                                             mean_dec_warming[1,mm],ci90_dec_warming[1,mm])
        plt.plot([mm_attrib+1.0,mm_attrib+1.0],att_warming_range,color=cols[1,:],linewidth=2,label=label[1])
        plt.plot(mm_attrib+1.0,att_warming,color=cols[1,:],marker='+')
        
        [att_warming,att_warming_range,dummy]=attrib_warming(att_out3[dataset]['beta'][0],
                                                             att_out3[dataset]['betaCI'][0,:],
                                                             mean_dec_warming[0,mm]-mean_dec_warming[1,mm]-mean_dec_warming[2,mm],
                                                             (ci90_dec_warming[0,mm]**2+ci90_dec_warming[1,mm]**2+ci90_dec_warming[2,mm]**2)**0.5)
        plt.plot([mm_attrib+1.2,mm_attrib+1.2],att_warming_range,color=cols[0,:],linewidth=2,label=label[0])
        plt.plot(mm_attrib+1.2,att_warming,color=cols[0,:],marker='+')
        if simple_uncert:
          plt.plot([mm_attrib+0.8,mm_attrib+0.8],att_out3[dataset]['betaCI'][2,:]*mean_dec_warming[2,mm],
                   color=cols[2,:],marker='_',linestyle='None')
          plt.plot([mm_attrib+1.0,mm_attrib+1.0],att_out3[dataset]['betaCI'][1,:]*mean_dec_warming[1,mm],
                   color=cols[1,:],marker='_',linestyle='None')
          plt.plot([mm_attrib+1.2,mm_attrib+1.2],att_out3[dataset]['betaCI'][0,:]*(mean_dec_warming[0,mm]-mean_dec_warming[1,mm]-mean_dec_warming[2,mm]),
                   color=cols[0,:],marker='_',linestyle='None')
        
        mm_attrib=mm_attrib+1 #Counter over those models used for attribution.
        label=[None,None,None]
#Multi-model analysis.
    dataset='Multi'
    model_names.append(dataset)
    (xr,yr,cn1,cn2)=da.reduce_dim(numpy.mean(mean_diag[:,0:3,model_indices],axis=2),
                                  had4_diag[:,None],anom[:,list(range(1,anom_index,2))],
                                  anom[:,list(range(0,anom_index,2))])
    yr =yr.data
    neff=mm_attrib**2/numpy.sum(1./ens_sizes[0:3,model_indices],axis=1) #Effective ensemble size when using multi-model mean.
    att_out[dataset]=da.tls(xr[:,0:2],yr,cn1,ne=neff[[0,1]],cn2=cn2,flag_2S=1,RCT_flag=rcplot)
    att_out3[dataset]=da.tls(xr[:,0:3],yr,cn1,ne=neff[[0,1,2]],cn2=cn2,flag_3S=1,RCT_flag=rcplot)
    multim_mean_dec_warming=numpy.mean(mean_dec_warming[:,model_indices],axis=1)
    multim_mean_gmst_comp_warming=numpy.mean(mean_gmst_comp_warming[:,:,model_indices],axis=2)
#Compute uncertainty in attributable warming based on spread in ratio of GSAT to GMST warming across models.
    mean_dec_warming_gmst=numpy.squeeze(numpy.mean(mean_diag[int(warming_years[0]/year_block):int(warming_years[1]/year_block)+1,:,:],
                                                   axis=0)-numpy.mean(mean_diag[int(warming_base[0]/year_block):int(warming_base[1]/year_block)+1,:,:],axis=0)) #Assume 5-yr means, calculate 2010-2019-1850-1899 in GMST.
    multim_ci90_dec_warming=(numpy.mean(ci90_dec_warming**2,axis=1))**0.5 #RMS dec warming uncertainty across models.
    multim_ci90_gmst_comp_warming=numpy.zeros((nyear,nexp))
    for yy in range(nyear):
      multim_ci90_gmst_comp_warming[yy,:]= (numpy.mean(ci90_gmst_comp_warming[yy,:,:]**2,axis=1))**0.5 #RMS dec warming uncertainty across models.   

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
      [att_warming,att_warming_range,dummy]=attrib_warming(att_out[dataset]['beta'][0],
                                                           att_out[dataset]['betaCI'][0,:],
                                                           multim_mean_dec_warming[0]-multim_mean_dec_warming[1],
                                                           (multim_ci90_dec_warming[0]**2+multim_ci90_dec_warming[1]**2)**0.5,ci90_beta_obs2[0])
      print ('main',att_out[dataset]['beta'][0],att_out[dataset]['betaCI'][0,:],multim_mean_dec_warming[0]-multim_mean_dec_warming[1],
             (multim_ci90_dec_warming[0]**2+multim_ci90_dec_warming[1]**2)**0.5,ci90_beta_obs2[0])
      fitted_ant=att_out[dataset]['beta'][0]*numpy.mean((mean_diag[:,0,:]-mean_diag[:,1,:]),axis=1)
      fitted_nat=att_out[dataset]['beta'][1]*numpy.mean((mean_diag[:,1,:]),axis=1)
      plt.plot([mm_attrib+0.9,mm_attrib+0.9],att_warming_range,color=colors[1,:],linewidth=2,label=ant)
      plt.plot(mm_attrib+0.9,att_warming,color=colors[1,:],marker='+')
      print ('ANT:',att_warming,att_warming_range)
      [att_warming,att_warming_range,dummy]=attrib_warming(att_out[dataset]['beta'][1],
                                                           att_out[dataset]['betaCI'][1,:],multim_mean_dec_warming[1],
                                                           multim_ci90_dec_warming[1],ci90_beta_obs2[1])
      plt.plot  ([mm_attrib+1.1,mm_attrib+1.1],att_warming_range,color=colors[4,:],linewidth=2,label=nat)
      plt.plot(mm_attrib+1.1,att_warming,color=colors[4,:],marker='+')
      print ('NAT:',att_warming,att_warming_range)
      if simple_uncert:
        plt.plot([mm_attrib+0.9,mm_attrib+0.9],
                 att_out[dataset]['betaCI'][0,:]*(multim_mean_dec_warming[0]-multim_mean_dec_warming[1]),
                 color=colors[1,:],marker='_',linestyle='None')
        plt.plot  ([mm_attrib+1.1,mm_attrib+1.1],
                   att_out[dataset]['betaCI'][1,:]*multim_mean_dec_warming[1],
                   color=colors[4,:],marker='_',linestyle='None')


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
      [att_warming,att_warming_range,dummy]=attrib_warming(att_out3[dataset]['beta'][2],att_out3[dataset]['betaCI'][2,:],
                                                           multim_mean_dec_warming[2],multim_ci90_dec_warming[2],ci90_beta_obs3[2])
      fitted_ghg=att_out3[dataset]['beta'][2]*numpy.mean((mean_diag[:,2,:]),axis=1)
      fitted_oth=att_out3[dataset]['beta'][0]*numpy.mean((mean_diag[:,0,:]-mean_diag[:,1,:]-mean_diag[:,2,:]),axis=1)
      plt.plot([mm_attrib+0.8,mm_attrib+0.8],att_warming_range,color=cols[2,:],linewidth=2)
      plt.plot(mm_attrib+0.8,att_warming,color=cols[2,:],marker='+')
      print (exp_flag,att_warming,att_warming_range)
      [att_warming,att_warming_range,dummy]=attrib_warming(att_out3[dataset]['beta'][1],
                                                           att_out3[dataset]['betaCI'][1,:],
                                                           multim_mean_dec_warming[1],
                                                           multim_ci90_dec_warming[1],ci90_beta_obs3[1])
      plt.plot([mm_attrib+1.0,mm_attrib+1.0],att_warming_range,color=cols[1,:],linewidth=2)    
      plt.plot(mm_attrib+1.0,att_warming,color=cols[1,:],marker='+')
      print ('NAT:',att_warming,att_warming_range)
      [att_warming,att_warming_range,dummy]=attrib_warming(att_out3[dataset]['beta'][0],
                                                           att_out3[dataset]['betaCI'][0,:],
                                                           multim_mean_dec_warming[0]-multim_mean_dec_warming[1]-multim_mean_dec_warming[2],
                                                           (multim_ci90_dec_warming[0]**2+multim_ci90_dec_warming[1]**2+multim_ci90_dec_warming[2]**2)**0.5,
                                                           ci90_beta_obs3[0])
      plt.plot([mm_attrib+1.2,mm_attrib+1.2],att_warming_range,color=cols[0,:],linewidth=2)
      plt.plot(mm_attrib+1.2,att_warming,color=cols[0,:],marker='+')
      print ('OTH:',att_warming,att_warming_range)
      if simple_uncert:
        plt.plot([mm_attrib+0.8,mm_attrib+0.8],att_out3[dataset]['betaCI'][2,:]*multim_mean_dec_warming[2],
                 color=cols[2,:],linestyle='None',marker='_')
        plt.plot([mm_attrib+1.0,mm_attrib+1.0],att_out3[dataset]['betaCI'][1,:]*multim_mean_dec_warming[1],
                 color=cols[1,:],linestyle='None',marker='_')    
        plt.plot([mm_attrib+1.2,mm_attrib+1.2],
                 att_out3[dataset]['betaCI'][0,:]*(multim_mean_dec_warming[0]-multim_mean_dec_warming[1]-multim_mean_dec_warming[2]),
                 color=cols[0,:],linestyle='None',marker='_')

      
    nmodel_attrib=mm_attrib
      
    obs_warming=obs_dec_warming[0]*numpy.mean(mean_dec_warming[0,:])/numpy.mean(mean_dec_warming_gmst[0,:])

    print ('Obs warming',obs_warming)
    print ('Obs warming in GMST',obs_dec_warming[0])
    #Finish off plots.
    panel_labels=['a','b','c','d','e','f']
    panel_counter=0
    

    for ff in [topleft,topright]:
        plt.subplot(ff)
        plt.plot([0,nmodel_attrib+2],[1,1],color='black',linewidth=1,ls=':')
        if ff == topleft: plt.ylabel('Regression coefficients')#,size='x-small')
        plt.plot([0,nmodel_attrib+2],[0,0],color='black',linewidth=1,ls='--')
        if pool_int_var:
            plt.axis([0,nmodel_attrib+2,-1,3])
#          plt.xticks(list(range(1,nmodel_attrib+2)),[''])
        else:
            plt.axis([0,nmodel_attrib+1,-1,3])
#          plt.xticks(list(range(1,nmodel_attrib+1)),[''])
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
#              plt.xticks(list(range(1,nmodel_attrib+2)),[''])
            else:
              plt.axis([0,nmodel_attrib+1,0,1])
#              plt.xticks(list(range(1,nmodel_attrib+1)),[''])
            if ff == 323: plt.ylabel('RCT P-value')#,size='x-small')
            plt.text (-2,1.075,panel_labels[panel_counter],fontsize =7,fontweight='bold', va='center', ha='center')
            panel_counter=panel_counter+1

    if simple_uncert: #Label lower panels a and b.
        panel_counter=0
      
    for ff in [bottomleft,bottomright]:
        plt.subplot(ff)
        plt.plot([0,nmodel_attrib+2],[obs_warming,obs_warming],color='black',linewidth=1,label='Had4 GSAT')
        if ff == bottomleft: 
           plt.ylabel(f'Attributable change {warming_years[0]}-{warming_years[1]} vs {warming_base[0]}-{warming_base[1]} ($^\circ$C)')#,size='x-small') 
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
    plt.savefig(plot_dir+'/reg_attrib_'+exp_flag+'_'+pool_flag+uncert_flag+'.'+output_file_type)
    plt.close()
    plt.figure(2,figsize=[180*mm_conv,60*mm_conv]) 
    for aa in range(anom_index):
        plt.plot(years,anom[:,aa])
    plt.savefig(plot_dir+'/anom.pdf')
    plt.close()
    plt.figure(2,figsize=[180*mm_conv,60*mm_conv]) 
#    for mm in range(3):
    plt.plot(years,mean_diag[:,1,0],color='black')
    plt.plot(years,mean_diag[:,1,1],color='green')
    plt.plot(years,mean_diag[:,1,2],color='gray')
    plt.savefig(plot_dir+'/mean_diag.pdf')
    plt.close()

    plt.figure(2,figsize=[180*mm_conv,60*mm_conv]) 
#    for mm in range(3):
    plt.plot(years,fitted_ant-numpy.mean(fitted_ant[0:10]),color='red')
    plt.plot(years,fitted_nat-numpy.mean(fitted_nat[0:10]),color='green')
    plt.plot(years,fitted_ghg-numpy.mean(fitted_ghg[0:10]),color='gray')
    plt.plot(years,fitted_oth-numpy.mean(fitted_oth[0:10]),color='blue')
    plt.plot(years,had4_diag[:,None]-numpy.mean(had4_diag[0:10,None]),color='black')
    plt.savefig(plot_dir+'/fitted_model.pdf')
    plt.close()

#Calculate annual mean timeseries for gmst_comp attributable warming, ANT, NAT, GHG, OTH.
    
    with open(plot_dir+'/cmip6_fitted_comp_gmst.csv', mode='w') as file:
        data_writer=csv.writer(file,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        data_writer.writerow(['Year', 'ANT 5', '17', '83', '95', '50', 'NAT 5', '17', '83', '95', '50','GHG 5', '17', '83', '95', '50', 'OTH 5', '17', '83', '95', '50'])
        data_writer.writerow(['CMIP6 DAMIP models spatially-complete GMST attributable components'])
        for yy in range(nyear):
          [ant_warming,ant_warming_range,ant_17_83]=attrib_warming(att_out[dataset]['beta'][0],
                                                                    att_out[dataset]['betaCI'][0,:],
                                                                    multim_mean_gmst_comp_warming[yy,0]-multim_mean_gmst_comp_warming[yy,1],
                                                                    (multim_ci90_gmst_comp_warming[yy,0]**2+multim_ci90_gmst_comp_warming[yy,1]**2)**0.5,ci90_beta_obs2[0])
          [nat_warming,nat_warming_range,nat_17_83]=attrib_warming(att_out[dataset]['beta'][1],att_out[dataset]['betaCI'][1,:],
                                                                    multim_mean_gmst_comp_warming[yy,1],
                                                                    multim_ci90_gmst_comp_warming[yy,1],ci90_beta_obs2[1])
          [ghg_warming,ghg_warming_range,ghg_17_83]=attrib_warming(att_out3[dataset]['beta'][2],
                                                                    att_out3[dataset]['betaCI'][2,:],
                                                                    multim_mean_gmst_comp_warming[yy,2],
                                                                    multim_ci90_gmst_comp_warming[yy,2],ci90_beta_obs3[2])
          [oth_warming,oth_warming_range,oth_17_83]=attrib_warming(att_out3[dataset]['beta'][0],
                                                                    att_out3[dataset]['betaCI'][0,:],
                                                                    multim_mean_gmst_comp_warming[yy,0]-multim_mean_gmst_comp_warming[yy,1]-multim_mean_gmst_comp_warming[yy,2],
                                                                    (multim_ci90_gmst_comp_warming[yy,0]**2+multim_ci90_gmst_comp_warming[yy,1]**2+multim_ci90_gmst_comp_warming[yy,2]**2)**0.5,
                                                                    ci90_beta_obs3[0])
          data_writer.writerow([y_start+yy,float(ant_warming_range[0]),float(ant_17_83[0]),float(ant_17_83[1]),float(ant_warming_range[1]),float(ant_warming),\
                                        float(nat_warming_range[0]),float(nat_17_83[0]),float(nat_17_83[1]),float(nat_warming_range[1]),float(nat_warming),\
                                        float(ghg_warming_range[0]),float(ghg_17_83[0]),float(ghg_17_83[1]),float(ghg_warming_range[1]),float(ghg_warming),\
                                        float(oth_warming_range[0]),float(oth_17_83[0]),float(oth_17_83[1]),float(oth_warming_range[1]),float(oth_warming)])
                                

      

    
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
