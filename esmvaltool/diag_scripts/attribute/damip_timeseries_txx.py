"""Plots Figures 1, EDF Figures 1 and 7 in Gillett et al."""
import logging
import os
from pprint import pformat
import numpy
import csv
import iris
import matplotlib
from scipy import stats
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import esmvaltool.diag_scripts.attribute.ncblendmask_esmval_txx as ncbm
import random
import math

from esmvaltool.diag_scripts.shared import (group_metadata, run_diagnostic, select_metadata)

logger = logging.getLogger(os.path.basename(__file__))


def main(cfg):
    matplotlib.use('Agg')
    plt.ioff() #Turn off interactive plotting.
    # Get a description of the preprocessed data that we will use as input, as well as other flags.
    input_data = cfg['input_data'].values()
    plot_dir=cfg['plot_dir']
    output_file_type=cfg['output_file_type']
    year_block = cfg.get('year_block')
    warming_base = cfg.get('warming_base')
    warming_years = cfg.get('warming_years')

    sh_name = select_metadata(input_data, variable_group='models')[0]['short_name']
    obs_file = os.path.join(cfg['auxiliary_data_dir'], f'{sh_name}_gridded.nc')
    # obs_file = select_metadata(input_data, variable_group='obs_mean')[0]['filename']
    # obs_cb = iris.load_cube(obs_file)

    obs_cb = iris.load_cube(obs_file)
    y_start = obs_cb.coord('time').cell(0).point.year
    y_end = obs_cb.coord('time').cell(-1).point.year

    hadlabel=f'{sh_name}_gridded' # change 
    # hadlabel=select_metadata(input_data, variable_group='obs_mean')[0]['dataset'] 

    ens_obs_cblst = ''

    grouped_input_data = group_metadata(
        input_data, 'dataset', sort='ensemble')
    logger.info(
        "Group input data by model and sort by ensemble:"
        "\n%s", pformat(grouped_input_data))
    # grouped_input_data.pop(hadlabel) # to delete
    nmodel=len(grouped_input_data)

#Initialise variables.
    experiments=['historical-ssp245','hist-GHG','hist-aer','hist-nat','hist-stratO3','hist-GHG-ssp245-GHG','hist-aer-ssp245-aer','hist-nat-ssp245-nat','hist-stratO3-ssp245-stratO3']
    labels=['Anthropogenic and natural forcings','Greenhouse gases','Aerosols','Natural forcings']
    nexp=len(experiments)-4 #Subtract three to account for repetition of hist-GHG, hist-nat, hist-aer.
    nyear=y_end - y_start +1 #Number of years, hard-coded.
    ldiag=int(nyear/year_block) #length of diagnostic.
    years=numpy.arange(y_start,y_end+1,year_block) #Used for plotting.
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
    ensobs_diag=[]
    ensobs_dec_warming=[]


#Set up figure including colours.
    font = {'size'   : 5}
    matplotlib.rc('font', **font)
    mm_conv=0.0394
    mod_cols=numpy.array([[0,73,73],[255,255,109],[0,146,146],[255,109,182],[255,182,119],[146,0,0],[0,109,219],[182,109,255],[109,182,255],[182,219,255],[73,0,146],[146,73,0],[219,209,0],[36,255,36]])/256.    
    cols=numpy.array([[0,0,0],[196,121,0],[178,178,178],[0,52,102],[0,79,0],[200,0,0],[0,200,0],[0,0,200],[112,160,205]])/256.
    shade_cols=numpy.array([[128,128,128,128],[204,174,113,128],[191,191,191,128],[67,147,195,128],[223,237,195,128],[255,150,150,128],[150,255,150,128],[150,150,255,128],[91,174,178,128]])/256.

    plt.figure(figsize=[88*mm_conv,113*mm_conv], dpi=200)
    plt.subplot(211)
    
#Loop over models, then experiments, then ensemble members.    
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
            #Label hist-nat-ssp245-nat as hist-nat, hist-ghg-ssp245-ghg as hist-ghg etc
            #(some models' hist-nat ends in 2014 so is merged with ssp245-nat).
            if experiment > 4: 
                experiment=experiment-4
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
                obs_dec_warming=[]
                ann_warming=[]
                gmst_comp_warming=[]
                #Calculate masked and blended GMST for individual simulation.
                (exp_diags[:,ee],obs_diag)=ncbm.ncblendmask_esmval(files[0],obs_cb,dec_warming,obs_dec_warming,ann_warming,gmst_comp_warming,year_block,ens_obs_cblst,ensobs_diag,ensobs_dec_warming, warming_years=warming_years, warming_base=warming_base)
                ens_obs_cblst='' #Set to empty string so that ensemble obs diagnostics are only calculated on the first iteration.
                #Take anomalies relative to 1850-1900.
                print ('exp_diags[:,ee]',exp_diags[:,ee])
                exp_diags[:,ee]=exp_diags[:,ee]-numpy.mean(exp_diags[0:(warming_base[1]+1-y_start),ee]) 
                obs_diag=obs_diag-numpy.mean(obs_diag[0:(warming_base[1]+1-y_start)])
                exp_ann_warming[:,ee]=ann_warming[0]
                exp_gmst_comp_warming[:,ee]=gmst_comp_warming[0]
                #Plot first ensemble member of historical.
                if exp=="historical-ssp245" and ee==0:
                    alpha_ens=1. if ee==0 else 0.2
                    ls='dashed' if mm > 7 else 'solid'
                    plt.plot(years,exp_diags[:,ee],color=mod_cols[mm],linewidth=0.5,label=lbl,zorder=1-ee,alpha=alpha_ens,linestyle=ls)
                    lbl=""
            mean_diag[:,experiment,mm]=numpy.mean(exp_diags,axis=1)
            mean_ann_warming[:,experiment,mm]=numpy.mean(exp_ann_warming,axis=1)
            mean_gmst_comp_warming[:,experiment,mm]=numpy.mean(exp_gmst_comp_warming,axis=1)
            all_ann_warming[:,experiment,mm,0:nens]=exp_diags #Use GMST.
            all_ann_warming_comp[:,experiment,mm,0:nens]=exp_gmst_comp_warming #Use spatially complete GMST.
            all_ann_warming_gsat[:,experiment,mm,0:nens]=exp_ann_warming #Use GSAT.
            
    #Write GMST and GSAT timeseries to CSV file for other applications (not needed for Gillett et al. plots).
    with open(plot_dir+'/cmip6_gmst.csv', mode='w') as file:
        data_writer=csv.writer(file,delimiter=',',quotechar='"', quoting=csv.QUOTE_MINIMAL)
        data_writer.writerow(['CMIP6 DAMIP models HadCRUT5-masked blended GMST (Cowtan et al., 2015) and globally-complete GSAT'])
        for experiment in range(nexp):
            data_writer.writerow(['Experiment:',experiments[experiment]])
            for mm, dataset in enumerate(grouped_input_data):
                data_writer.writerow([dataset])
                for ee in range(int(ens_sizes[experiment,mm])):
                    data_writer.writerow(['Ensemble member',ee])
                    data_writer.writerow(['Year, GMST_complete, GMST_HadCRUT5_masked, GSAT'])
                    for yy in range(ldiag):
                        data_writer.writerow([years[yy],all_ann_warming[yy,experiment,mm,ee],all_ann_warming_comp[yy,experiment,mm,ee],all_ann_warming_gsat[yy,experiment,mm,ee]])


#Calculate ratio of GSAT to GMST warming for individual simulations.
    denom=numpy.mean(all_ann_warming[(warming_years[0]-y_start):(warming_years[1] +1 -y_start),0,:,:],axis=0)
    ratio_by_model=numpy.mean(all_ann_warming_gsat[(warming_years[0]-y_start):(warming_years[1] +1 -y_start),0,:,:],axis=0)/denom
    copy_ratio_by_model=numpy.reshape(ratio_by_model[:,:],nmodel*nensmax)
#Equivalent calculation for spatially-complete GMST.    
    denom=numpy.mean(all_ann_warming_comp[(warming_years[0]-y_start):(warming_years[1] +1 -y_start),0,:,:],axis=0)
    ratio_by_model_comp=numpy.mean(all_ann_warming_gsat[(warming_years[0]-y_start):(warming_years[1] +1 -y_start),0,:,:],axis=0)/denom
    copy_ratio_by_model_comp=numpy.reshape(ratio_by_model_comp[:,:],nmodel*nensmax)

    
    plt.plot(years,obs_diag,color='black',linewidth=1,label=hadlabel)
    plt.plot(years,numpy.mean(mean_diag[:,0,:],axis=1),color=cols[1],linewidth=1,label='Model mean')
    print ('Mean GMST',numpy.mean(mean_diag[:,0,:],axis=1))
    plt.plot([y_start,2025],[0,0],color='black',linewidth=0.5,ls='--',zorder=0)
    plt.axis([y_start,2025,-4, 4])
    plt.xlabel('Year')
    plt.ylabel(f'Canada {sh_name} anomaly ($^\circ$C)')
    plt.legend(loc="upper left",ncol=2, fancybox=False, frameon=False)
    for experiment in range(nexp):
        wts=numpy.zeros((nmodel,nensmax))
        for mm in range(nmodel):
            wts[mm,0:int(ens_sizes[experiment,mm])]=1./ens_sizes[experiment,mm]
        wts=numpy.reshape(wts,nmodel*nensmax)/numpy.sum(wts)
        if experiment==0:
#Calculate ratio of GSAT to obs-masked GMST warming 5-95% range.
           sort_ratio=numpy.sort(copy_ratio_by_model)
           sort_index=numpy.argsort(copy_ratio_by_model)
           cdf=numpy.cumsum(wts[sort_index])
           range_ratio=[sort_ratio[cdf>=0.05][0],sort_ratio[cdf>=0.95][0]]
           print ('5-95% range of ratio for GSAT/obs-masked GMST',range_ratio)
#Equivalent calculation for spatially-complete GMST
           sort_ratio=numpy.sort(copy_ratio_by_model_comp)
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
    plt.text (1825,2.25,'a',fontsize =7,fontweight='bold', va='center', ha='center')        
    plt.subplot(212)
    zzs=[3,1,0,2]
    for experiment in range(4):
        offset=0
        plt.fill_between(years,range_ann_warming[:,experiment,0]+offset,range_ann_warming[:,experiment,1]+offset,
                         color=shade_cols[experiment+1,:],zorder=zzs[experiment])
        plt.plot(years,mm_ann_warming[:,experiment]+offset,color=cols[experiment+1,:],linewidth=1,
                 label=labels[experiment],zorder=zzs[experiment]+4)

    plt.plot(years,obs_diag,color='black',linewidth=1,label=hadlabel,zorder=8)
    plt.axis([y_start,2025,-4, 4])
    plt.plot([y_start,2025],[0,0],color='black',linewidth=0.5,ls='--',zorder=0)
    plt.xlabel('Year')
    plt.ylabel(f'Canada {sh_name} anomaly ($^\circ$C)')
    plt.legend(loc="upper left", fancybox=False, frameon=False,ncol=2)
    plt.text (1825,2.25,'b',fontsize =7,fontweight='bold', va='center', ha='center')        
    plt.savefig(plot_dir+'/timeseries.'+output_file_type)
    plt.close()

#Plot Extended Data Fig 1 showing all DAMIP GMST timeseries.
    fig=    plt.figure(figsize=[88*mm_conv,176*mm_conv], dpi=200)
    ax1=fig.add_subplot(111)
    for experiment in range(nexp):
        offset=experiment*-1.5
        plt.fill_between(years,range_ann_warming[:,experiment,0]+offset,range_ann_warming[:,experiment,1]+offset,color=shade_cols[experiment+1,:])
        plt.plot([y_start,2025],[offset,offset],color='black',linewidth=0.5)
        plt.plot(years,mm_ann_warming[:,experiment]+offset,color=cols[experiment+1,:],linewidth=0.5,label=experiments[experiment])
        plt.text(1860,offset+0.4,experiments[experiment])
    plt.plot(years,obs_diag,color='black',linewidth=0.5,label=hadlabel,zorder=8)
    ax1.set_xlim(1850,2025)
    ax1.set_xlabel('Year')
    ax1.set_ylim(-11,2)
    plt.yticks(numpy.arange(27)*0.5-11,['','','-1.0','-0.5','0.0','0.5','1.0','','-1.0','-0.5','0.0','0.5',
                                        '1.0','','-1.0','-0.5','0.0','0.5','1.0','','-1.0','-0.5',
                                        '0.0','0.5','1.0','',''])
    ax1.set_ylabel('Global mean surface temperature change ($^\circ$C)')
    ax2=ax1.twinx()
    ax2.set_ylim(-11,2)
    plt.yticks(numpy.arange(27)*0.5-11,['-0.5',' 0.0',' 0.5',' 1.0','','-1.0','-0.5',' 0.0',' 0.5',' 1.0',
                                        '','-1.0','-0.5',' 0.0',' 0.5',' 1.0','','-1.0','-0.5',' 0.0',
                                        ' 0.5',' 1.0','','','','',''])
    plt.savefig(plot_dir+'/supplement_timeseries.'+output_file_type)
    plt.close()

#Plot of ratios of GSAT to GMST (Extended Data Fig 7). 
    plt.figure(figsize=[88*mm_conv,113*mm_conv])
    plt.subplot(211)
    plt.ylabel('Ratio of GSAT to HadCRUT5-masked GMST')
    plt.axis([0,nmodel+1,1.0,1.3])
    plt.xticks(list(range(1,nmodel+1)),model_names,rotation=30.,ha="right")
    for mm, dataset in enumerate(grouped_input_data):
        ens_size=ens_sizes[0,mm]
        for ee in range(ens_size.astype(int)):
         plt.plot(mm+1,numpy.mean(all_ann_warming_gsat[(warming_years[0]-y_start):(warming_years[1] +1 -y_start),0,mm,ee],
                                  axis=0)/numpy.mean(all_ann_warming[(warming_years[0]-y_start):(warming_years[1]+1-y_start),0,mm,ee],
                                                     axis=0),color=mod_cols[mm],marker='+')
    plt.text (-2,1.31,'a',fontsize =7,fontweight='bold', va='center', ha='center')
    plt.tight_layout()

    plt.subplot(212)
    plt.ylabel('Ratio of GSAT to globally-complete GMST') 
    plt.axis([0,nmodel+1,1.0,1.3])
    plt.xticks(list(range(1,nmodel+1)),model_names,rotation=30.,ha="right")
    for mm, dataset in enumerate(grouped_input_data):
        ens_size=ens_sizes[0,mm]
        for ee in range(ens_size.astype(int)):
          plt.plot(mm+1,numpy.mean(all_ann_warming_gsat[(warming_years[0]-y_start):(warming_years[1] +1 -y_start),0,mm,ee],axis=0)/numpy.mean(all_ann_warming_comp[(warming_years[0]-y_start):(warming_years[1] +1 -y_start),0,mm,ee],axis=0),color=mod_cols[mm],marker='+')
    plt.text (-2,1.31,'b',fontsize =7,fontweight='bold', va='center', ha='center')
    plt.tight_layout()
    plt.savefig(plot_dir+'/gmst_to_gsat_ratios.'+output_file_type)
    plt.close()

#Calculate uncertainty in GMST and GSAT warming in obs (as reported in Gillett et al.).
 

    plt.figure(2,figsize=[180*mm_conv,60*mm_conv], dpi=200) 
    plt.plot(years,all_ann_warming[:,3,1,0],color='green')
    plt.plot(years,all_ann_warming[:,3,1,1],color='green')
    plt.plot(years,all_ann_warming[:,3,1,2],color='green')
    plt.plot(years,all_ann_warming[:,0,1,0],color='black')
    plt.plot(years,all_ann_warming[:,0,1,1],color='black')
    plt.plot(years,all_ann_warming[:,0,1,2],color='black')
    plt.savefig(plot_dir+'/giss_timeseries.'+output_file_type)
    plt.close()

#print ('plot_dir',plot_dir)    
 
if __name__ == '__main__':

    with run_diagnostic() as config:
        main(config)
