"""Carry out cross-validation of attribution results using blended and masked temperatures."""
import logging
import os
from pprint import pformat
from scipy.stats import t

import iris

import numpy
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

def attrib_warming(beta,betaCI,dec_warming,ci95_dec_warming):
  print ('attrib_warming inputs:',beta,betaCI,dec_warming,ci95_dec_warming)
  attrib_warming=beta*dec_warming
  attrib_range=[attrib_warming-numpy.absolute(beta)*dec_warming*(((beta-betaCI[0])/beta)**2+(ci95_dec_warming/dec_warming)**2)**0.5,attrib_warming+numpy.absolute(beta)*dec_warming*(((betaCI[1]-beta)/beta)**2+(ci95_dec_warming/dec_warming)**2)**0.5]
  print ('attrib_warming, attrib_range',attrib_warming,numpy.sort(attrib_range))
  print ('betaCI*dec_warming',betaCI*dec_warming)
  return (attrib_warming,numpy.sort(attrib_range))



def main(cfg):
    plt.ioff() #Turn off interactive plotting.
    font = {'size'   : 5}
    matplotlib.rc('font', **font)
    mm_conv=0.0394
    plt.figure(figsize=[88*mm_conv,88*mm_conv])

    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    had4_file='/home/rng/data/esmvaltool/HadCRUT.4.6.0.0.median.nc'
    sftlf_file='/home/rng/data/esmvaltool/CNRM-CM6-1-5x5-sftlf.nc'
#sftlf_file='/pf/b/b380746/CNRM-CM6-1-5x5-sftlf.nc' #Hard-coded path to sftlf file for CNRM-CM6 on a 5x5 grid. (Can't input through pre-processor at the moment. Update with sftlf for each model through preprocessor later.)
    
    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='ensemble')
    logger.info(
        "Group input data by model and sort by ensemble:"
        "\n%s", pformat(grouped_input_data))
    print (type (grouped_input_data))
    print (len(grouped_input_data))
    nmodel=len(grouped_input_data)
# Colours
    colors=numpy.array([[0,0,0],[196,121,0],[178,178,178],[0,52,102],[0,79,0],[200,0,0],[0,200,0],[0,0,200],[112,160,205],])/256.
    shade_cols=numpy.array([[128,128,128,128],[204,174,113,128],[191,191,191,128],[67,147,195,128],[223,237,195,128],[255,150,150,128],[150,255,150,128],[150,150,255,128],[91,174,178,128]])/256.
    exp_flag='AER' #Use hist-aer for 3-way regression.
#    exp_flag='GHG' #Use hist-GHG
    if exp_flag=='GHG':
        experiments=['historical-ssp245','hist-nat','hist-GHG','hist-aer','hist-CO2','hist-stratO3','hist-volc','hist-sol','hist-nat-ssp245-nat','hist-GHG-ssp245-GHG','hist-aer-ssp245-aer']
        cols=[3,4,2]
    else:
        experiments=['historical-ssp245','hist-nat','hist-aer','hist-GHG','hist-CO2','hist-stratO3','hist-volc','hist-sol','hist-nat-ssp245-nat','hist-aer-ssp245-aer','hist-GHG-ssp245-GHG'] 
        cols=[2,4,3]        
    nexp=len(experiments)-3
    print ('Number of experiments', nexp)
    # Loop over variables/datasets in alphabetical order
#   Define variables for D&A analysis
    diag_name='gmst05'
    av_yr=int(diag_name[4:6]) #Last two digits of diag_name are averaging period in yrs.
    years=list(numpy.arange(1850+av_yr/2,2020+av_yr/2,av_yr)) #Used for plotting.
    ldiag=int(170/av_yr) #length of diagnostic, assuming 170 years of data.
    if diag_name[0:4]=='hemi': ldiag=ldiag*2
#    diag_name='gmst10'
#    years=list(numpy.arange(1855,2025,10)) #Used for plotting.
#    ldiag=17 #length of diagnostic,hard-coded for the moment.

    anom_max=500 #arbitrary max size for number of anomalies.
    nens_max=100
    mean_diag=numpy.zeros((ldiag,nexp,nmodel))
    mean_dec_warming=numpy.zeros((nexp,nmodel))
    ci95_dec_warming=numpy.zeros((nexp,nmodel))
    pseudo_obs=numpy.zeros((ldiag,nens_max,nmodel))
    anom=numpy.zeros((ldiag,anom_max))
    anom_index=0
    ens_sizes=numpy.zeros((nexp,nmodel),dtype=int)

    for mm, dataset in enumerate(grouped_input_data):
#        if mm > 5: continue
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
                had4_dec_warming=[]
                
                (exp_diags[:,ee],had4_diag)=ncbm.ncblendmask_esmval('max', files[0],files[1],files[2],sftlf_file,had4_file,dec_warming,had4_dec_warming,0,0,diag_name)
                if exp=='historical-ssp245': pseudo_obs[:,ee,mm]=exp_diags[:,ee]
                exp_dec_warming[ee]=dec_warming[0]
            mean_diag[:,experiment,mm]=numpy.mean(exp_diags,axis=1)
            mean_dec_warming[experiment,mm]=numpy.mean(exp_dec_warming)
            if nens==1:
                ci95_dec_warming[experiment,mm]=ci95_dec_warming[numpy.nonzero(ci95_dec_warming*(ens_sizes**0.5))].mean() #use mean of CIs already calculated, corrected for ensemble size.
            else:        
                ci95_dec_warming[experiment,mm]=(numpy.std(exp_dec_warming,ddof=1)/((nens)**0.5))*t.ppf(0.95,nens-1)
                print ('t.ppf',t.ppf(0.95,nens-1))
            if nens>1: anom[:,anom_index:anom_index+nens-1]=(exp_diags[:,0:nens-1]-mean_diag[:,experiment:experiment+1,mm])*((nens/(nens-1))**0.5) #Intra-ensemble anomalies for use as pseudo-control. Only use nens-1 ensemble members to ensure independence, and inflate variance by sqrt (nens / (nens-1)) to account for subtraction of ensemble mean.
            anom_index=anom_index+nens-1
    anom=anom[:,0:anom_index]
    mean_dec_warming[0,:]=mean_dec_warming[0,:]-mean_dec_warming[1,:]-mean_dec_warming[2,:] #Replace historical with OTH.
    ci95_dec_warming[0,:]=(ci95_dec_warming[0,:]**2+ci95_dec_warming[1,:]**2+ci95_dec_warming[2,:]**2)**0.5
    att_out={}
    model_names=[]
    model_indices=[]
    fig={}
    mm_attrib=0 #Counter over those models used for attribution.
    for mm, dataset in enumerate(grouped_input_data):
        if mean_diag[0,1,mm] == 0: #If there is no hist-nat simulation skip over model.
            continue
        model_names.append(dataset)
        model_indices.append(mm)
        mm_attrib=mm_attrib+1
    nin=[0,0,0]
    nout=0
    if exp_flag=='GHG':
      label=['OTH','NAT','GHG']
#      cols=['orange','blue','green']
    else:
      label=['GHG','NAT','AER']
#      cols=['green','blue','orange']        
    print ('******* mean_dec_warming',mean_dec_warming[0:3,model_indices])
    print ('******* ci95_dec_warming',ci95_dec_warming[0:3,model_indices])
    print ('******* ens_sizes',ens_sizes[0:3,model_indices])
#    anom=anom*1.5 #Inflate standard deviation of internal variability by 50%.
    for mm in range(mm_attrib):
#        if mm > 5: continue
        other_model_indices=model_indices.copy()
        other_model_indices.pop(mm) #All models but the one used for pseudo-obs.
        mms=model_indices[mm] #The model used for pseudo-obs.
        print ('Model indices',model_indices)
        print ('Other model indices',other_model_indices)
        other_model_mean=numpy.mean(mean_diag[:,0:3,other_model_indices],axis=2)
        neff=(mm_attrib-1)**2/numpy.sum(1./ens_sizes[0:3,other_model_indices],axis=1) #Effective ensemble size when using multi-model mean.
        other_mean_dec_warming=numpy.mean(mean_dec_warming[:,other_model_indices],axis=1)
#        other_mean_ci95=numpy.mean(ci95_dec_warming[:,other_model_indices],axis=1)/len(other_model_indices)
#Compute uncertainty in attributable warming based on spread in ratio of GSAT to GMST warming across models.
        if diag_name[0:4]=='hemi':
          nlat=2
          nper=ldiag//2
          mean_diag=numpy.mean(numpy.reshape(mean_diag,(nlat,nper,nexp,nmodel)),axis=0)
        mean_dec_warming_gmst=numpy.squeeze(numpy.mean(mean_diag[32:34,:,:],axis=0)-numpy.mean(mean_diag[0:10,:,:],axis=0)) #Assume 5-yr means, calculate 2010-2019-1850-1899 in GMST.
        mean_dec_warming_gmst[0,:]=mean_dec_warming_gmst[0,:]-mean_dec_warming_gmst[1,:]-mean_dec_warming_gmst[2,:] #Replace historical with OTH.
        #Define uncertainty in multi-model mean warming to 2010-2019 based on spread in ratio of GSAT to GMST increase across models.
        other_mean_ci95=numpy.std(mean_dec_warming/mean_dec_warming_gmst,ddof=1,axis=1)*t.ppf(0.95,nmodel-1)*numpy.mean(mean_dec_warming_gmst,axis=1)
        print ('mean_dec_warming/mean_dec_warming_gmst',mean_dec_warming/mean_dec_warming_gmst)
        print ('numpy.std(mean_dec_warming/mean_dec_warming_gmst,ddof=1,axis=1)*t.ppf(0.95,nmodel-1)',numpy.std(mean_dec_warming/mean_dec_warming_gmst,ddof=1,axis=1)*t.ppf(0.95,nmodel-1))
        print ('numpy.mean(mean_dec_warming_gmst,axis=1)',numpy.mean(mean_dec_warming_gmst,axis=1))
        print ('numpy.mean(mean_dec_warming_gmst,axis=1)',numpy.mean(mean_dec_warming_gmst,axis=1))
        
        for ee in range(ens_sizes[0,model_indices[mm]]):
            offset =0.005*ee-0.0025*(ens_sizes[0,model_indices[mm]]-1)
            pobs=pseudo_obs[:,ee,mms]
            print ('mms,ee',mms,ee)
#            print ('pobs',pobs)
#            print ('other_model_mean',other_model_mean)
            (xr,yr,cn1,cn2)=da.reduce_dim(other_model_mean,pobs[:,None],anom[:,list(range(1,anom_index,2))],anom[:,list(range(0,anom_index,2))])
#            print ('shapes',numpy.shape(other_model_mean),numpy.shape(pseudo_obs[:,ee,mm]),numpy.shape(xr),numpy.shape(yr))
            att_out=da.tls(xr,yr,cn1,cn2=cn2,ne=neff,flag_3S=1,RCT_flag=0) 
#            print ('betaCI',att_out['betaCI'])
#            print ('mean_dec_warming',mean_dec_warming[0:3,mms])
#            print ('other_mean_dec_warming[0:3]',other_mean_dec_warming[0:3])
            print ('Simulated AER warming',mean_dec_warming[2,mms],'AER-attrib warming',numpy.transpose(att_out['betaCI'][2,:]*other_mean_dec_warming[2]))
            print ('Simulated NAT warming',mean_dec_warming[1,mms],'NAT-attrib warming',numpy.transpose(att_out['betaCI'][1,:]*other_mean_dec_warming[1]))
            print ('Simulated GHG warming',mean_dec_warming[0,mms],'GHG-attrib warming',numpy.transpose(att_out['betaCI'][0,:]*other_mean_dec_warming[0]))
            for experiment in range(0,3):
                [att_warming,att_range]=attrib_warming(att_out['beta'][experiment],att_out['betaCI'][experiment,:],other_mean_dec_warming[experiment],other_mean_ci95[experiment])
                plt.plot([mean_dec_warming[experiment,mms]+offset,mean_dec_warming[experiment,mms]+offset],att_range,color=shade_cols[cols[experiment],:],zorder=1,linewidth=0.4)
                plt.plot([mean_dec_warming[experiment,mms]],[att_warming],color=colors[cols[experiment],:],marker='x',markersize=2,label=label[experiment],zorder=2)
                if att_warming-((att_warming-att_range[0])**2+ci95_dec_warming[experiment,mms]**2)**0.5 < mean_dec_warming[experiment,mms] < att_warming+((att_range[1]-att_warming)**2+ci95_dec_warming[experiment,mms]**2)**0.5: nin[experiment]=nin[experiment]+1
#                print (att_warming,mean_dec_warming[experiment,mms],ci_att_warming,ci95_dec_warming[experiment,mms])
            label=['','','']
    print ('Fraction within CI')
    print (label)
    print (nin/(numpy.sum(ens_sizes[0,model_indices])))
    print (nin,numpy.sum(ens_sizes[0,model_indices]))


    plt.axis([-1.5,3,-1.5,3])
    plt.xlabel('Simulated change in GSAT ($^\circ$C)',fontsize=7)
    plt.ylabel('Reconstructed change in GSAT ($^\circ$C)',fontsize=7)
    plt.plot([-2,3],[-2,3],color='black',linewidth=1,ls=':')             
    plt.legend(loc="upper left")
    plt.savefig('/home/rng/plots/esmvaltool/xval_'+diag_name+'_'+exp_flag+'.pdf')
    plt.close()
     
    

if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
