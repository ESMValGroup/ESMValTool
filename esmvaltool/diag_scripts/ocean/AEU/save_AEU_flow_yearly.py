import numpy as np
import glob
import yaml
import netCDF4 as nc4
import xarray as xr
import nctoolkit as nc

# NB, a bunch od staff will be passed with argparse:
#     FLIST, PREPROC, 

#1. Read input file with the list of processed models
FLIST='/home/users/gig/CRACAB/FIND_FILES/output_match_hist_and_ssp.yml' #ALL MODELS
#FLIST='/home/users/gig/CRACAB/FIND_FILES/output_best_models.yml' #BEST MODELS


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

# reads into a list, each element of which is a dict
with open(FLIST,'rb') as f:
    Stuff=yaml.safe_load(f)

# type of preproc
#PREPROC='tight_None'
#PREPROC='loose_10_years'
#PREPROC='loose_None'
PREPROC='medium_None'
#PREPROC='medium_1_years'
#PREPROC='medium_10_years'
#PREPROC='tight_10_years'
#PREPROC='tight_1_years'
#PREPROC='tight_None'


INDIR='/home/users/gig/CRACAB/ENSEMBLE_RUN_monthly_all/esmvaltool_output/*/plots/diag_aeu/AMOC_timeseries/'
INDIR_PRE='/home/users/gig/CRACAB/ENSEMBLE_RUN_monthly_all/esmvaltool_output/*/preproc/diag_aeu/aeu/'

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

days_month=np.array([31,28,31,30,31,30,31,31,30,31,30,31],dtype=np.float32)[np.newaxis,:]
days_year=365.0 #sorry Feb 29!


#2.2 Read AEUs and write them to the dictionary
#    average of files with same dataset and exp

# Count ensemble members:
# 1. unique models
count_hist_ok=0
count_hist_no=0
count_ssp_ok=0
count_ssp_no=0
count_models=dict.fromkeys(['historical','ssp119','ssp126','ssp245','ssp370','ssp585'],0)
# 2. total members
count_members=dict.fromkeys(['historical','ssp119','ssp126','ssp245','ssp370','ssp585'],0)
# 3. total models
count_unique_models=dict.fromkeys(Data.keys(),0)

for dset in Data.keys():
    for exp in Data[dset].keys():
        curr_flist=sorted(glob.glob(INDIR+'Timeseries_aeu_'+dset+'_'+exp+'_*'+PREPROC+'.nc'))
        curr_TS=0.0
        if len(curr_flist)!=0: #if there are files at all!
            for ff in range(len(curr_flist)):
                D=nc4.Dataset(curr_flist[ff],'r')
                curr_uo=np.array(D.variables['uo'])
                if any(curr_uo>1000.): # interpolate invalid
                    curr_uo=interpolate_invalid(curr_uo)
                if any(curr_uo>1000.):
                    print('STILL INVALID VALUES ON '+curr_flist[ff])
                # make yearly average of monthly values
                curr_uo=np.sum(days_month*(curr_uo.reshape(-1,12)),axis=1)/days_year
                curr_TS=curr_TS+curr_uo
                D.close()
            curr_TS=curr_TS/len(curr_flist)
            print(dset+' '+exp+' '+str(len(curr_flist))+' members')
            Data[dset][exp]=curr_TS
            count_members[exp]=count_members[exp]+len(curr_flist)
            count_models[exp]=count_models[exp]+1
            count_unique_models[dset]=1
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

# OUTPUTS:
#FILE1:
# flow_mean_historical
# flow_std_historical
# flow_min_historical
# flow_max_historical

#FILE2:
# flow_mean_ssp126
# flow_std_ssp126
# flow_min_ssp126
# flow_max_ssp126
# etc.

AR_hist=np.array(np.zeros(65))
for i in Data.keys():
    if 'historical' in Data[i]:
        if len(Data[i]['historical'])!=0:
            AR_hist=np.vstack([AR_hist,Data[i]['historical']])

AR_hist=AR_hist[1:,:]

AR_ssp126=np.array(np.zeros(86))
for i in Data.keys():
    if 'ssp126' in Data[i]:
        if len(Data[i]['ssp126'])!=0:
            AR_ssp126=np.vstack([AR_ssp126,Data[i]['ssp126']])

AR_ssp126=AR_ssp126[1:,:]

AR_ssp245=np.array(np.zeros(86))
for i in Data.keys():
    if 'ssp245' in Data[i]:
        if len(Data[i]['ssp245'])!=0:
            AR_ssp245=np.vstack([AR_ssp245,Data[i]['ssp245']])

AR_ssp245=AR_ssp245[1:,:]

AR_ssp370=np.array(np.zeros(86))
for i in Data.keys():
    if 'ssp370' in Data[i]:
        if len(Data[i]['ssp370'])!=0:
            AR_ssp370=np.vstack([AR_ssp370,Data[i]['ssp370']])

AR_ssp370=AR_ssp370[1:,:]

AR_ssp585=np.array(np.zeros(86))
for i in Data.keys():
    if 'ssp585' in Data[i]:
        if len(Data[i]['ssp585'])!=0:
            AR_ssp585=np.vstack([AR_ssp585,Data[i]['ssp585']])

AR_ssp585=AR_ssp585[1:,:]


years_hist=np.arange(1950,2015,dtype=np.float32)
years_scen=np.arange(2015,2101,dtype=np.float32)


XR_hist=xr.Dataset(data_vars=dict(flow_mean_historical=(['time'],np.mean(AR_hist,axis=0)),
                                  flow_std_historical= (['time'],np.std(AR_hist,axis=0)),
                                  flow_min_historical= (['time'],np.min(AR_hist,axis=0)),
                                  flow_max_historical= (['time'],np.max(AR_hist,axis=0))),
                   coords=dict(time=years_hist),
                   attrs=dict(units='Sv',
                              description='annual averaged Atlantic Equatorial Undercurrent (AEU) statistics, from averaged CMIP6 ensemble, historical simulations (1950-2015), values calculated as positive (east-west) velocity times cell area at 23W, between 3S and -3N and  between 50 and 200m. values are derived by first averaging all ensemble members for each CMIP6 model and then taking the average of all models.',
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
                              description='annual averaged Atlantic Equatorial Undercurrent (AEU) statistics, from averaged CMIP6 ensemble, scenario simulations (SSP1 2.6, SSP2 4.5, SSP3 7.0, SSP5 8.5, 2015-2100), values calculated as positive (east-west) velocity times cell area at 23W, between 3S and -3N and  between 50 and 200m. values are derived by first averaging all ensemble members for each CMIP6 model and then taking the average of all models.',
                              author='Giovanni Galli, Plymouth Marine Laboratory, 2022'))


OUTDIR='/home/users/gig/CRACAB/FINAL_PLOT/'
OUTFILE_HIST=OUTDIR+'AEU_flow_historical_1950-2014.nc'
OUTFILE_SSP=OUTDIR+'AEU_flow_scenarios_2015-2100.nc'

XR_hist.to_netcdf(OUTFILE_HIST)
XR_ssp.to_netcdf(OUTFILE_SSP)

#ds_hist=nc.from_xarray(XR_hist)
#ds_hist.assign_coords(time=('time',years_hist))
#ds_hist.flow_mean_historical.assign_attrs({'name':'flow_mean_historical',
#                                           'coordinates':'time',
#                                           'units':'Sv'})
#ds_hist.flow_std_historical.assign_attrs({'name':'flow_std_historical',
#                                           'coordinates':'time',
#                                           'units':'Sv'})
#ds_hist.flow_max_historical.assign_attrs({'name':'flow_max_historical',
#                                           'coordinates':'time',
#                                           'units':'Sv'})
#ds_hist.flow_min_historical.assign_attrs({'name':'flow_min_historical',
#                                           'coordinates':'time',
#                                           'units':'Sv'})
