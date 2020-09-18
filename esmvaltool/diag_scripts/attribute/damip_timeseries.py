"""Plots timeseries of masked and blended GMST from CMIP6 models, for comparison with obs."""
import logging
import os
from pprint import pformat

import iris

import numpy
import csv
import matplotlib
from scipy import stats
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import ncblendmask_esmval as ncbm

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic,
                                            select_metadata, sorted_metadata)
from esmvaltool.diag_scripts.shared._base import (
    ProvenanceLogger, get_diagnostic_filename, get_plot_filename)
from esmvaltool.diag_scripts.shared.plot import quickplot

logger = logging.getLogger(os.path.basename(__file__))


def get_provenance_record(attributes, ancestor_files):
    """Create a provenance record describing the diagnostic data and plot."""
    caption = ("Average {long_name} between {start_year} and {end_year} "
               "according to {dataset}.".format(**attributes))

    record = {
        'caption': caption,
        'statistics': ['mean'],
        'domains': ['global'],
        'plot_type': 'zonal',
        'authors': [
            'ande_bo',
            'righ_ma',
        ],
        'references': [
            'acknow_project',
        ],
        'ancestors': ancestor_files,
    }
    return record


def main(cfg):
    matplotlib.use('Agg')
    plt.ioff() #Turn off interactive plotting.
    """Compute the time average for each input dataset."""
    # Get a description of the preprocessed data that we will use as input.
    input_data = cfg['input_data'].values()
    obs=cfg['obs']
#    obs='had4'
    if obs=='had5':
        obs_file='/home/rng/data/esmvaltool/HadCRUT.5.0.0.0.analysis.anomalies.ensemble_median.nc'
        hadlabel='HadCRUT5'
    elif obs=='had4':
        obs_file='/home/rng/data/esmvaltool/HadCRUT.4.6.0.0.median.nc'  #Updated to end of 2019.
        hadlabel='HadCRUT4'
    else:
        exit('Obs not recognised')
    sftlf_file='/home/rng/data/esmvaltool/CNRM-CM6-1-5x5-sftlf.nc' #Hard-coded path to sftlf file for CNRM-CM6 on a 5x5 grid. (Can't input through pre-processor at the moment. Update with sftlf for each model through preprocessor later.)


    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='ensemble')
    logger.info(
        "Group input data by model and sort by ensemble:"
        "\n%s", pformat(grouped_input_data))
    type (grouped_input_data)
    print (len(grouped_input_data))
    nmodel=len(grouped_input_data)
    experiments=['historical-ssp245','hist-GHG','hist-aer','hist-nat','hist-volc','hist-sol','hist-stratO3','hist-CO2','hist-GHG-ssp245-GHG','hist-aer-ssp245-aer','hist-nat-ssp245-nat']
    labels=['Anthropogenic and natural forcings','Greenhouse gases','Aerosols','Natural forcings']
    mod_cols=['olive','teal','navy','orange','yellow','lime','green','cyan','blue','purple','magenta','pink','coral']
    cols=numpy.array([[0,0,0],[196,121,0],[178,178,178],[0,52,102],[0,79,0],[200,0,0],[0,200,0],[0,0,200],[112,160,205]])/256.
    shade_cols=numpy.array([[128,128,128,128],[204,174,113,128],[191,191,191,128],[67,147,195,128],[223,237,195,128],[255,150,150,128],[150,255,150,128],[150,150,255,128],[91,174,178,128]])/256.
    nexp=len(experiments)-3 #Subtract three to account for repetition of hist-GHG, hist-nat etc.
    print ('Number of experiments', nexp)
    # Loop over variables/datasets in alphabetical order
    diag_name='gmst01' #Annual mean GMST.
    ldiag=170 #length of diagnostic,hard-coded for the moment.
    years=list(range(1850,2020,1)) #Used for plotting.
    anom_max=500 #arbitrary max size for number of anomalies.
    mean_diag=numpy.zeros((ldiag,nexp,nmodel))
    mean_gmst_comp_warming=numpy.zeros((ldiag,nexp,nmodel))
    mean_ann_warming=numpy.zeros((ldiag,nexp,nmodel))
    mm_ann_warming=numpy.zeros((ldiag,nexp))
    range_ann_warming=numpy.zeros((ldiag,nexp,2)) # 5-95% range.
    nensmax=50
    msval=1e20
    all_ann_warming=numpy.full((ldiag,nexp,nmodel,nensmax),msval)
    all_ann_warming_comp=numpy.full((ldiag,nexp,nmodel,nensmax),msval)
    all_ann_warming_gsat=numpy.full((ldiag,nexp,nmodel,nensmax),msval)
    ens_sizes=numpy.zeros((nexp,nmodel))
    model_names=[]

    plt.figure(figsize=[7,9])
    plt.subplot(211)
    
    
    for mm, dataset in enumerate(grouped_input_data):
        logger.info("*************** Processing model %s", dataset)
        model_names.append(dataset)
        lbl=dataset
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
            exp_ann_warming=numpy.zeros((ldiag,nens))
            exp_gmst_comp_warming=numpy.zeros((ldiag,nens))
            
        
            for ee, ensemble in enumerate(grouped_exp_input_data):
                logger.info("** Processing ensemble %s", ensemble)
                files=[]
                for attributes in grouped_exp_input_data[ensemble]:
                    logger.info("Processing variable %s", attributes['variable_group'])
                    files.append(attributes['filename'])
                logger.info("*************** Files for blend and mask %s", files)
                dec_warming=[]
                had4_dec_warming=[]
                ann_warming=[]
                gmst_comp_warming=[]
                (exp_diags[:,ee],had4_diag)=ncbm.ncblendmask_esmval('max', files[0],files[1],files[2],sftlf_file,obs_file,dec_warming,had4_dec_warming,ann_warming,gmst_comp_warming,diag_name,obs)
                exp_diags[:,ee]=exp_diags[:,ee]-numpy.mean(exp_diags[0:(1901-1850),ee]) #Take anomalies relative to 1850-1900.
                had4_diag=had4_diag-numpy.mean(had4_diag[0:(1901-1850)])
                exp_ann_warming[:,ee]=ann_warming[0]
                exp_gmst_comp_warming[:,ee]=gmst_comp_warming[0]
                if exp=="historical-ssp245":
                    alpha_ens=1. if ee==0 else 0.2
                    plt.plot(years,exp_diags[:,ee],color=mod_cols[mm],linewidth=1,label=lbl,zorder=1-ee,alpha=alpha_ens)
                    lbl=""
            mean_diag[:,experiment,mm]=numpy.mean(exp_diags,axis=1)
            mean_ann_warming[:,experiment,mm]=numpy.mean(exp_ann_warming,axis=1)
            mean_gmst_comp_warming[:,experiment,mm]=numpy.mean(exp_gmst_comp_warming,axis=1)
            all_ann_warming[:,experiment,mm,0:nens]=exp_diags #Use GMST.
            all_ann_warming_comp[:,experiment,mm,0:nens]=exp_gmst_comp_warming #Use spatially complete GMST.
            all_ann_warming_gsat[:,experiment,mm,0:nens]=exp_ann_warming #Use GSAT.

    with open('/home/rng/plots/esmvaltool/cmip6_gmst.csv', mode='w') as file:
        data_writer=csv.writer(file,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        data_writer.writerow(['CMIP6 DAMIP models HadCRUT4-masked blended GMST (Cowtan et al., 2015) and globally-complete GSAT'])
        for experiment in range(nexp):
            data_writer.writerow(['Experiment:',experiments[experiment]])
            for mm, dataset in enumerate(grouped_input_data):
                data_writer.writerow([dataset])
                for ee in range(int(ens_sizes[experiment,mm])):
                    data_writer.writerow(['Ensemble member',ee])
                    data_writer.writerow(['Year, GMST_complete, GMST_HadCRUT4_masked, GSAT'])
                    for yy in range(ldiag):
                        data_writer.writerow([years[yy],all_ann_warming[yy,experiment,mm,ee],all_ann_warming_comp[yy,experiment,mm,ee],all_ann_warming_gsat[yy,experiment,mm,ee]])

            
    slope_gsat,intercept,r_value,p_value,std_err=stats.linregress(years[0:45],numpy.mean(mean_ann_warming[1970-1850:2015-1850,0,:],axis=1))
    slope_gmst,intercept,r_value,p_value,std_err=stats.linregress(years[0:45],numpy.mean(mean_diag[1970-1850:2015-1850,0,:],axis=1))
    
    print ('Ratio of 1970-2014 warming trends',slope_gsat/slope_gmst)
    print ('slope_gsat',slope_gsat,'slope_gmst',slope_gmst)
    print ('Ratio of GSAT to GMST warming',numpy.mean(mean_ann_warming[2010-1850:2020-1850,0,:])/numpy.mean(mean_diag[2010-1850:2020-1850,0,:]))
    print ('Ratio of GSAT to GMST_comp warming',numpy.mean(mean_ann_warming[2010-1850:2020-1850,0,:])/numpy.mean(mean_gmst_comp_warming[2010-1850:2020-1850,0,:]))
    print ('Ratio of GSAT to GMST_comp warming by model',numpy.mean(mean_ann_warming[2010-1850:2020-1850,0,:],axis=0)/numpy.mean(mean_gmst_comp_warming[2010-1850:2020-1850,0,:],axis=0))
    print ('Ratio of GSAT to GMST warming by model',numpy.mean(mean_ann_warming[2010-1850:2020-1850,0,:],axis=0)/numpy.mean(mean_diag[2010-1850:2020-1850,0,:],axis=0))
    print ('Ratio of GSAT to GMST warming by simulation',numpy.mean(all_ann_warming_gsat[2010-1850:2020-1850,0,:,:],axis=0)/numpy.mean(all_ann_warming[2010-1850:2020-1850,0,:,:],axis=0))
    print ('Mean warming by experiment and model',numpy.mean(mean_ann_warming[2010-1850:2020-1850,:,:],axis=0))
           
    denom=numpy.mean(all_ann_warming[2010-1850:2020-1850,0,:,:],axis=0)
#    denom[denom==msval]=0. #Assign missing values to zero.
    print ('denom',denom)
    ratio_by_model=numpy.mean(all_ann_warming_gsat[2010-1850:2020-1850,0,:,:],axis=0)/denom
    copy_ratio_by_model=numpy.reshape(ratio_by_model[:,:],nmodel*nensmax)
    ratio_by_model[denom==msval]=numpy.nan
    print ('Ratio across models **',ratio_by_model)
    print ('Standard deviation across simulations',numpy.nanstd(ratio_by_model))
    print ('Mean across simulations',numpy.nanmean(ratio_by_model))
#Equivalent calculation for spatially-complete GMST.    
    denom=numpy.mean(all_ann_warming_comp[2010-1850:2020-1850,0,:,:],axis=0)
#    denom[denom==msval]=0. #Assign missing values to zero.
    print ('comp_denom',denom)
    ratio_by_model_comp=numpy.mean(all_ann_warming_gsat[2010-1850:2020-1850,0,:,:],axis=0)/denom
    copy_ratio_by_model_comp=numpy.reshape(ratio_by_model_comp[:,:],nmodel*nensmax)
    ratio_by_model_comp[denom==msval]=numpy.nan

    
    plt.plot(years,had4_diag,color='black',linewidth=2,label=hadlabel)
    plt.plot(years,numpy.mean(mean_diag[:,0,:],axis=1),color=cols[1],linewidth=2,label='Model mean GMST')
    plt.plot(years,numpy.mean(mean_ann_warming[:,0,:],axis=1),color='red',linewidth=2,label='Model mean GSAT')
    plt.plot([1850,2020],[0,0],color='black',linewidth=1,ls='--',zorder=0)
    plt.axis([1850,2020,-1,2])
    plt.xlabel('Year')
    plt.ylabel('Global mean temperature anomaly ($^\circ$C)')
    plt.legend(loc="upper left",ncol=2)
    print ('Ens sizes')
    print (ens_sizes)
    for experiment in range(nexp):
        wts=numpy.zeros((nmodel,nensmax))
        for mm in range(nmodel):
            wts[mm,0:int(ens_sizes[experiment,mm])]=1./ens_sizes[experiment,mm]
        wts=numpy.reshape(wts,nmodel*nensmax)/numpy.sum(wts)
#        print ('wts',wts)
        if experiment==0:
           sort_ratio=numpy.sort(copy_ratio_by_model)
           print ('sort_ratio',sort_ratio)
           sort_index=numpy.argsort(copy_ratio_by_model)
           cdf=numpy.cumsum(wts[sort_index])
           range_ratio=[sort_ratio[cdf>=0.05][0],sort_ratio[cdf>=0.95][0]]
           print ('5-95% range of ratio for GSAT/HadCRUT4-masked GMST',range_ratio)
#Equivalent calculation for spatially-complete GMST
           sort_ratio=numpy.sort(copy_ratio_by_model_comp)
           print ('sort_ratio',sort_ratio)
           sort_index=numpy.argsort(copy_ratio_by_model_comp)
           cdf=numpy.cumsum(wts[sort_index])
           range_ratio=[sort_ratio[cdf>=0.05][0],sort_ratio[cdf>=0.95][0]]
           print ('5-95% range of ratio for GSAT/spatially-complete GMST',range_ratio)
        for yy in range(ldiag):
                year_warming=numpy.reshape(all_ann_warming[yy,experiment,:,:],nmodel*nensmax)
                sort_warming=numpy.sort(year_warming)
                sort_index=numpy.argsort(year_warming)
                cdf=numpy.cumsum(wts[sort_index])
                range_ann_warming[yy,experiment,:]=[sort_warming[cdf>=0.05][0],sort_warming[cdf>=0.95][0]]
                mm_ann_warming[yy,experiment]=numpy.sum(year_warming*wts)
    plt.subplot(212)
    zzs=[3,1,0,2]
    for experiment in range(4):
#        offset=experiment*-2.
        offset=0
        plt.fill_between(years,range_ann_warming[:,experiment,0]+offset,range_ann_warming[:,experiment,1]+offset,color=shade_cols[experiment+1,:],zorder=zzs[experiment])
#        plt.plot([1850,2025],[offset,offset],color='black',linewidth=1)
        plt.plot(years,mm_ann_warming[:,experiment]+offset,color=cols[experiment+1,:],linewidth=2,label=labels[experiment],zorder=zzs[experiment]+4)

#    plt.plot(years,had4_diag*1.1719758280521986,color='black',linewidth=1,label='Observations',zorder=8)
    plt.plot(years,had4_diag,color='black',linewidth=2,label=hadlabel,zorder=8)
    plt.axis([1850,2020,-1,2])
    plt.plot([1850,2020],[0,0],color='black',linewidth=1,ls='--',zorder=0)
    plt.xlabel('Year')
    plt.ylabel('Global mean surface temperature anomaly ($^\circ$C)')
    plt.legend(loc="upper left")

    plt.savefig('/home/rng/plots/esmvaltool/Fig1_'+obs+'.png')

    plt.close()
    print ('ANT timeseries')
    print (mm_ann_warming[:,0]-mm_ann_warming[:,3])
#    print ('Years with negative ANT')
#    print (years[numpy.reshape(mm_ann_warming[:,0]-mm_ann_warming[:,3],ldiag)<0.])

    fig=    plt.figure(figsize=[7,14])
    ax1=fig.add_subplot(111)
    for experiment in range(nexp):
        offset=experiment*-1.5
        plt.fill_between(years,range_ann_warming[:,experiment,0]+offset,range_ann_warming[:,experiment,1]+offset,color=shade_cols[experiment+1,:])
        plt.plot([1850,2025],[offset,offset],color='black',linewidth=1)
        plt.plot(years,mm_ann_warming[:,experiment]+offset,color=cols[experiment+1,:],linewidth=1,label=experiments[experiment])
        plt.text(1860,offset+0.4,experiments[experiment])
    plt.plot(years,had4_diag,color='black',linewidth=1,label=hadlabel,zorder=8)
#    ax1.legend(loc="upper left",ncol=2)
    ax1.set_xlim(1850,2020)
    ax1.set_xlabel('Year')
    ax1.set_ylim(-11,2)
    plt.yticks(numpy.arange(26)*0.5-11,['','','-1.0','-0.5','0.0','0.5','1.0','','-1.0','-0.5','0.0','0.5','1.0','','-1.0','-0.5','0.0','0.5','1.0','','-1.0','-0.5','0.0','0.5','1.0','',''])
    ax1.set_ylabel('Global mean surface temperature change ($^\circ$C)')
    ax2=ax1.twinx()
    ax2.set_ylim(-11,2)
    plt.yticks(numpy.arange(26)*0.5-11,['-0.5',' 0.0',' 0.5',' 1.0','','-1.0','-0.5',' 0.0',' 0.5',' 1.0','','-1.0','-0.5',' 0.0',' 0.5',' 1.0','','-1.0','-0.5',' 0.0',' 0.5',' 1.0','','','','',''])
    plt.savefig('/home/rng/plots/esmvaltool/supplement_timeseries'+obs+'.png')

    plt.close()

#Plots of ratios of GSAT to GMST. 
    plt.figure(figsize=[7,9])
    plt.subplot(211)
    plt.ylabel('Ratio of GSAT to HadCRUT4-masked GMST',size='small') 
    plt.axis([0,nmodel+1,1.0,1.3])
    plt.xticks(list(range(1,nmodel+1)),model_names, size='xx-small',rotation=30.,ha="right")
    for mm, dataset in enumerate(grouped_input_data):
        print ('mm',mm)
        print (ens_sizes[0,mm])
        ens_size=ens_sizes[0,mm]
        print (type(ens_size))
        print (type(ens_size.astype(int)))
        for ee in range(ens_size.astype(int)):
          plt.plot(mm+1,numpy.mean(all_ann_warming_gsat[2010-1850:2020-1850,0,mm,ee],axis=0)/numpy.mean(all_ann_warming[2010-1850:2020-1850,0,mm,ee],axis=0),color=mod_cols[mm],marker='+')
    plt.subplot(212)
    plt.ylabel('Ratio of GSAT to globally-complete GMST',size='small') 
    plt.axis([0,nmodel+1,1.0,1.3])
    plt.xticks(list(range(1,nmodel+1)),model_names, size='x-small',rotation=30.,ha="right")
    for mm, dataset in enumerate(grouped_input_data):
        ens_size=ens_sizes[0,mm]
        for ee in range(ens_size.astype(int)):
          plt.plot(mm+1,numpy.mean(all_ann_warming_gsat[2010-1850:2020-1850,0,mm,ee],axis=0)/numpy.mean(all_ann_warming_comp[2010-1850:2020-1850,0,mm,ee],axis=0),color=mod_cols[mm],marker='+')
    plt.savefig('/home/rng/plots/esmvaltool/gmst_to_gsat_ratios'+obs+'.png')
    plt.close()

#Plots of GSAT/GMST ratios versus warming. 
    plt.figure(figsize=[5,7])
    plt.subplot(211)
    plt.ylabel('Ratio of GSAT to HadCRUT4-masked GMST',size='small') 
    plt.xlabel('GSAT warming in 2010-2019 vs 1850-1900',size='small') 
    plt.axis([0,2,1.0,1.3])
    for mm, dataset in enumerate(grouped_input_data):
        print ('mm, x,y',mm,numpy.mean(all_ann_warming_gsat[2010-1850:2020-1850,0,mm,:]),numpy.mean(all_ann_warming_gsat[2010-1850:2020-1850,0,mm,:])/numpy.mean(all_ann_warming[2010-1850:2020-1850,0,mm,:]))
        plt.plot(numpy.mean(mean_ann_warming[2010-1850:2020-1850,0,mm]),numpy.mean(mean_ann_warming[2010-1850:2020-1850,0,mm])/numpy.mean(mean_diag[2010-1850:2020-1850,0,mm]),color=mod_cols[mm],marker='+')
    plt.subplot(212)
    plt.ylabel('Ratio of GSAT to globally-complete GMST',size='small') 
    plt.xlabel('GSAT warming in 2010-2019 vs 1850-1900',size='small') 
    plt.axis([0,2,1.0,1.1])
    for mm, dataset in enumerate(grouped_input_data):
        plt.plot(numpy.mean(mean_ann_warming[2010-1850:2020-1850,0,mm]),numpy.mean(mean_ann_warming[2010-1850:2020-1850,0,mm])/numpy.mean(mean_gmst_comp_warming[2010-1850:2020-1850,0,mm]),color=mod_cols[mm],marker='+')
    plt.savefig('/home/rng/plots/esmvaltool/gmst_to_gsat_ratios_vs_warming'+obs+'.png')

    plt.close()
  
    
   



 
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
