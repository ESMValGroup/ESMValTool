import numpy as np
import glob
import yaml
import netCDF4 as nc4
import xarray as xr
import nctoolkit as nc
from datetime import datetime as dt
from dateutil.relativedelta import relativedelta

# this deals with one single monthly value in 
# CNRM-ESM2-1 ssp370 r4i1p1f2 that for some reason
# is == 1e+20 in the medium and tight layouts (but not in loose)
def interpolate_invalid(data):
    print('INTERPOLATING INVALID')
    data1=data.copy()
    invalid_array=data>100.0
    invalid_idx=np.where(data>100.0)
    for ii in range(len(invalid_idx)):
        data1[invalid_idx[ii]]=((data[invalid_idx[ii]+1]-data[invalid_idx[ii]-1])/2)+data[invalid_idx[ii]-1]
    return data1

#1. Read input file with the list of processed models
FLIST='/home/users/gig/CRACAB/FIND_FILES/output_match_hist_and_ssp.yml' #ALL MODELS
#FLIST='/home/users/gig/CRACAB/FIND_FILES/output_best_models.yml' #BEST MODELS

# reads into a list, each element of which is a dict
with open(FLIST,'rb') as f:
    Stuff=yaml.safe_load(f)

# type of preproc
PREPROC='medium_None'

INDIR='/home/users/gig/CRACAB/ENSEMBLE_RUN_monthly_all/esmvaltool_output/*/plots/diag_aeu/AMOC_timeseries/'

# reconstruct indexes of hist and ssp time-periods for monthly data
year_start_hist=1950
year_end_hist=2015
year_start_ssp=2015
year_end_ssp=2100

time_hist=np.linspace(1950.,2015.,num=12*(2015-1950)+1)
time_ssp= np.linspace(2015.,2100.,num=12*(2100-2015)+1)

i_start_hist=np.where(time_hist==year_start_hist)[0][0]
i_stop_hist=np.where(time_hist==year_end_hist)[0][0]
i_start_ssp=np.where(time_ssp==year_start_ssp)[0][0]
i_stop_ssp=np.where(time_ssp==year_end_ssp)[0][0]

#2. Load data to a dictionary, average data from the same model in the process
#   (1 Model 1 Vote)

#2.1 Preallocate Data dictionary
dataset_list=[x['dataset'] for x in Stuff]
Data=dict.fromkeys(dataset_list)
for dd in Data.keys():
    Data[dd]={}

for ss in range(len(Stuff)):
    curr_dset=Stuff[ss]['dataset']
    curr_exp=Stuff[ss]['exp']
    Data[curr_dset][curr_exp]=np.array([])

#2.2 Read AEUs and write them to the dictionary
#    average of files with same dataset and exp
#    i.e. average all members
count_hist_ok=0
count_hist_no=0
count_ssp_ok=0
count_ssp_no=0
for dset in Data.keys():
    for exp in Data[dset].keys():
        if exp=='historical':
            istart=i_start_hist
            istop=i_stop_hist
        else:
            istart=i_start_ssp
            istop=i_stop_ssp
        curr_flist=glob.glob(INDIR+'Timeseries_aeu_'+dset+'_'+exp+'_*'+PREPROC+'.nc')
        curr_TS=0.0
        if len(curr_flist)!=0: #if there are files at all!
            for ff in range(len(curr_flist)):
                D=nc4.Dataset(curr_flist[ff],'r')
                #compute monthly means
                curr_TS_x=np.array(D.variables['uo'])[istart:istop]
                if any(curr_TS_x>1000.0):
                    curr_TS_x=interpolate_invalid(curr_TS_x) #SEE def interpolate_invalid(data)
                #curr_TS_y=np.mean(curr_TS_x.reshape(-1,12),axis=0)
                curr_TS_y=curr_TS_x.copy()
                curr_TS=curr_TS+curr_TS_y
                D.close()
            curr_TS=curr_TS/len(curr_flist) #mean of ensemble members for 1 model, 1 scenario
            print(dset+' '+exp+' '+str(len(curr_flist))+' members')
            Data[dset][exp]=curr_TS
            if exp=='historical':
                count_hist_ok=count_hist_ok+1
            else:
                count_ssp_ok=count_ssp_ok+1
        else:
            print(dset+' '+exp+' FAILED TO RUN')
            if exp=='historical':
                count_hist_no=count_hist_no+1
            else:
                count_ssp_no=count_ssp_no+1

print('TOT HISTORICAL: '+str(count_hist_ok))
print('TOT SSP: '+str(count_ssp_ok))
print('FAILED HISTORICAL: '+str(count_hist_no))
print('FAILED SSP: '+str(count_ssp_no))

##2.3 Make another dictionary with the average of the scenarios
#EXPs=['historical','ssp126','ssp245','ssp370','ssp585']
#Data_av=dict.fromkeys(EXPs,0.0)
#exp_count=dict.fromkeys(EXPs,0)
#for dset in Data.keys():
#    for exp in Data[dset].keys():
#        if exp!='ssp119' and len(Data[dset][exp])!=0:
#            Data_av[exp]=Data_av[exp]+Data[dset][exp]
#            exp_count[exp]=exp_count[exp]+1
#
#for ee in EXPs:
#    if exp_count[ee]!=0:
#        Data_av[ee]=Data_av[ee]/exp_count[ee]
    


AR_hist=np.array(np.zeros(65*12))
for i in Data.keys():
    if 'historical' in Data[i]:
        if len(Data[i]['historical'])!=0:
            AR_hist=np.vstack([AR_hist,Data[i]['historical']])

AR_hist=AR_hist[1:,:]

AR_ssp126=np.array(np.zeros(85*12))
for i in Data.keys():
    if 'ssp126' in Data[i]:
        if len(Data[i]['ssp126'])!=0:
            AR_ssp126=np.vstack([AR_ssp126,Data[i]['ssp126']])

AR_ssp126=AR_ssp126[1:,:]

AR_ssp245=np.array(np.zeros(85*12))
for i in Data.keys():
    if 'ssp245' in Data[i]:
        if len(Data[i]['ssp245'])!=0:
            AR_ssp245=np.vstack([AR_ssp245,Data[i]['ssp245']])

AR_ssp245=AR_ssp245[1:,:]

AR_ssp370=np.array(np.zeros(85*12))
for i in Data.keys():
    if 'ssp370' in Data[i]:
        if len(Data[i]['ssp370'])!=0:
            AR_ssp370=np.vstack([AR_ssp370,Data[i]['ssp370']])

AR_ssp370=AR_ssp370[1:,:]

AR_ssp585=np.array(np.zeros(85*12))
for i in Data.keys():
    if 'ssp585' in Data[i]:
        if len(Data[i]['ssp585'])!=0:
            AR_ssp585=np.vstack([AR_ssp585,Data[i]['ssp585']])

AR_ssp585=AR_ssp585[1:,:]



def datearray(startdate,enddate,sinceepoch):
    dates=[]
    while startdate <= enddate:
        dates.append(startdate)
        startdate=startdate+relativedelta(months=1)
    dates=np.array([nc4.date2num(i,sinceepoch) for i in dates])
    return dates

sinceepoch='seconds since 1950-01-01 00:00:00'
years_hist=datearray(dt.strptime("1950-01-16", "%Y-%m-%d"),dt.strptime("2014-12-16", "%Y-%m-%d"),sinceepoch)
years_scen=datearray(dt.strptime("2015-01-16", "%Y-%m-%d"),dt.strptime("2099-12-16", "%Y-%m-%d"),sinceepoch)


XR_hist=xr.Dataset(data_vars=dict(flow_mean_historical=(['time'],np.mean(AR_hist,axis=0)),
                                  flow_std_historical= (['time'],np.std(AR_hist,axis=0)),
                                  flow_min_historical= (['time'],np.min(AR_hist,axis=0)),
                                  flow_max_historical= (['time'],np.max(AR_hist,axis=0))),
                   coords=dict(time=years_hist),
                   attrs=dict(units='Sv',
                              description='monthly averaged Atlantic Equatorial Undercurrent (AEU) statistics, from averaged CMIP6 ensemble, historical simulations (1950-2015), values calculated as positive (east-west) velocity times cell area at 23W, between 3S and -3N and  between 50 and 200m. values are derived by first averaging all ensemble members for each CMIP6 model and then taking the average of all models.',
                              time_units='seconds since 1950-01-01 00:00:00',
                              author='Giovanni Galli, Plymouth Marine Laboratory, 2022'))

XR_ssp=xr.Dataset(data_vars=dict(flow_mean_ssp126=(['time'],np.mean(AR_ssp126,axis=0)),
                                 flow_std_ssp126= (['time'],np.std(AR_ssp126,axis=0)),
                                 flow_min_ssp126= (['time'],np.min(AR_ssp126,axis=0)),
                                 flow_max_ssp126= (['time'],np.max(AR_ssp126,axis=0)),
                                 flow_mean_ssp245=(['time'],np.mean(AR_ssp245,axis=0)),
                                 flow_std_ssp245= (['time'],np.std(AR_ssp245,axis=0)),
                                 flow_min_ssp245= (['time'],np.min(AR_ssp245,axis=0)),
                                 flow_max_ssp245= (['time'],np.max(AR_ssp245,axis=0)),
                                 flow_mean_ssp370=(['time'],np.mean(AR_ssp370,axis=0)),
                                 flow_std_ssp370= (['time'],np.std(AR_ssp370,axis=0)),
                                 flow_min_ssp370= (['time'],np.min(AR_ssp370,axis=0)),
                                 flow_max_ssp370= (['time'],np.max(AR_ssp370,axis=0)),
                                 flow_mean_ssp585=(['time'],np.mean(AR_ssp585,axis=0)),
                                 flow_std_ssp585= (['time'],np.std(AR_ssp585,axis=0)),
                                 flow_min_ssp585= (['time'],np.min(AR_ssp585,axis=0)),
                                 flow_max_ssp585= (['time'],np.max(AR_ssp585,axis=0))),
                   coords=dict(time=years_scen),
                   attrs=dict(units='Sv',
                              description='monthly averaged Atlantic Equatorial Undercurrent (AEU) statistics, from averaged CMIP6 ensemble, scenario simulations (SSP1 2.6, SSP2 4.5, SSP3 7.0, SSP5 8.5, 2015-2100), values calculated as positive (east-west) velocity times cell area at 23W, between 3S and -3N and  between 50 and 200m. values are derived by first averaging all ensemble members for each CMIP6 model and then taking the average of all models.',
                              time_units='seconds since 1950-01-01 00:00:00',
                              author='Giovanni Galli, Plymouth Marine Laboratory, 2022'))


OUTDIR='/home/users/gig/CRACAB/FINAL_PLOT/'
OUTFILE_HIST=OUTDIR+'AEU_monthly_flow_historical_1950-2014.nc'
OUTFILE_SSP=OUTDIR+'AEU_monthly_flow_scenarios_2015-2100.nc'

XR_hist.to_netcdf(OUTFILE_HIST)
XR_ssp.to_netcdf(OUTFILE_SSP)

