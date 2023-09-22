import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import glob
import yaml
import netCDF4 as nc4
import pickle
import scipy.interpolate
import xarray as xr
import nctoolkit as nc

# NB, a bunch od staff will be passed with argparse:
#     FLIST, PREPROC, 

#1. Read input file with the list of processed models
FLIST='/home/users/gig/CRACAB/FIND_FILES/output_match_hist_and_ssp_succesful_diagnostic.yml' #ALL MODELS
#FLIST='/home/users/gig/CRACAB/FIND_FILES/output_best_models.yml'

#OUTFILE='/home/users/gig/CRACAB/AEU_Ensemble_Averages_uo_transect.pkl' #ALL MODELS
#OUTFILE='/home/users/gig/CRACAB/AEU_Ensemble_Averages_uo_transect_best_models.pkl'

# reads into a list, each element of which is a dict
with open(FLIST,'rb') as f:
    Stuff=yaml.safe_load(f)


#INDIR='/home/users/gig/CRACAB/ENSEMBLE_RUN_3/esmvaltool_output/*/preproc/diag_aeu/aeu/'
INDIR='/home/users/gig/CRACAB/ENSEMBLE_RUN_monthly_all/esmvaltool_output/*/preproc/diag_aeu/aeu/'

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


Coords=dict.fromkeys(dataset_list)
for dd in Coords.keys():
    Coords[dd]={}

for ss in range(len(Stuff)):
    curr_dset=Stuff[ss]['dataset']
    curr_exp=Stuff[ss]['exp']
    Coords[curr_dset][curr_exp]={'lat':np.array([0]),'depth':np.array([0])}

#2.2 Read AEUs and write them to the dictionary
#    average of files with same dataset and exp
count_hist_ok=0
count_hist_no=0
count_ssp_ok=0
count_ssp_no=0
for dset in Data.keys():
    for exp in Data[dset].keys():
        curr_flist=glob.glob(INDIR+'CMIP6_'+dset+'_Omon_'+exp+'*.nc')
        curr_uo=0.0
        if len(curr_flist)!=0: #if there are files at all!
            for ff in range(len(curr_flist)):
                D=nc4.Dataset(curr_flist[ff],'r')
                curr_uo=curr_uo+np.array(D.variables['uo'])
                if ff==0:
                    Coords[dset][exp]['depth']=np.array(D.variables['lev'])
                    Coords[dset][exp]['lat']=np.array(D.variables['lat'])
                D.close()
            curr_uo=curr_uo/len(curr_flist)
            print(dset+' '+exp+' '+str(len(curr_flist))+' members')
            Data[dset][exp]=curr_uo
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


#2.3 Make another dictionary with the average of the scenarios
#    NB, models have a different vertical grid, so need to regrid them

# Load transect mask
D=nc4.Dataset('/home/users/gig/CRACAB/transect_mask.nc','r')
Lat=np.array(D.variables['lat'])
Depth=np.array(D.variables['depth'])
D.close()

Y,Z=np.meshgrid(Lat,Depth)

# REGRID
EXPs=['historical','ssp119','ssp126','ssp245','ssp370','ssp585']
#Data_av=dict.fromkeys(EXPs,0.0)
#exp_count=dict.fromkeys(EXPs,0)
for dset in Data.keys():
    for exp in Data[dset].keys():
        if len(Data[dset][exp])!=0:
            print(dset+' '+exp)
            curr_uo=Data[dset][exp]
            new_uo=np.zeros([curr_uo.shape[0],Y.shape[0],Y.shape[1]])
            Yi,Zi=np.meshgrid(Coords[dset][exp]['lat'],Coords[dset][exp]['depth'])
            for tt in range(curr_uo.shape[0]):
                new_uo[tt,:,:]=scipy.interpolate.griddata((Yi.flatten(),Zi.flatten()),curr_uo[tt,:,:].flatten(),(Y,Z),method='nearest') 
            Data[dset][exp]=new_uo
            #exp_count[exp]=exp_count[exp]+1

#for ee in EXPs:
#    if exp_count[ee]!=0:
#        Data_av[ee]=Data_av[ee]/exp_count[ee]

# YEARLY AVERAGE
days_month=np.array([31,28,31,30,31,30,31,31,30,31,30,31],dtype=np.float32)[np.newaxis,:,np.newaxis,np.newaxis]
days_year=365.0

Data_Y=dict.fromkeys(dataset_list)
for dd in Data_Y.keys():
    Data_Y[dd]={}

for dset in Data_Y.keys():
    for exp in Data[dset].keys():
        if len(Data[dset][exp])!=0:
            C=Data[dset][exp]
            Data_Y[dset][exp]=np.sum(days_month*(C.reshape(-1,12,C.shape[1],C.shape[2])),axis=1)/days_year

## VERTICAL PROFILE
#y_equator=(9,10) #the two cells equally closer to the equator, north and south
#for dset in Data_Y.keys():
#    for exp in Data_Y[dset].keys():
#        if len(Data_Y[dset][exp])!=0:
#            Data_Y[dset][exp]=np.mean(Data_Y[dset][exp][:,:,y_equator[0]:y_equator[1]+1],axis=2)


# to xarrays
AR_hist=np.stack([Data_Y[i]['historical'] for i in Data_Y.keys() if 'historical' in Data_Y[i].keys()])
AR_ssp126=np.stack([Data_Y[i]['ssp126'] for i in Data_Y.keys() if 'ssp126' in Data_Y[i].keys()])
AR_ssp245=np.stack([Data_Y[i]['ssp245'] for i in Data_Y.keys() if 'ssp245' in Data_Y[i].keys()])
AR_ssp370=np.stack([Data_Y[i]['ssp370'] for i in Data_Y.keys() if 'ssp370' in Data_Y[i].keys()])
AR_ssp585=np.stack([Data_Y[i]['ssp585'] for i in Data_Y.keys() if 'ssp585' in Data_Y[i].keys()])

years_hist=np.arange(1950,2015,dtype=np.float32)
years_scen=np.arange(2015,2101,dtype=np.float32)
depth=Z[:,0]

XR_hist=xr.Dataset(data_vars=dict(velocity_mean_historical=(['time','depth','lat'],np.mean(AR_hist,axis=0)),
                                  velocity_std_historical= (['time','depth','lat'],np.std(AR_hist,axis=0)),
                                  velocity_min_historical= (['time','depth','lat'],np.min(AR_hist,axis=0)),
                                  velocity_max_historical= (['time','depth','lat'],np.max(AR_hist,axis=0))),
                   coords=dict(time=years_hist,
                               depth=Depth,
                               lat=Lat),
                   attrs=dict(units='m/s',
                              description='annual averaged Atlantic Equatorial Undercurrent (AEU) velocity and related statistics (mean value, std, min, max), from averaged CMIP6 ensemble, historical simulations (1950-2015), values calculated between 3 deg N and 3 deg S and between 0 and 500m. values are derived by first averaging all ensemble members for each CMIP6 model and then taking the average of all models. Positive velocities are east-to-west',
                              author='Giovanni Galli, Plymouth Marine Laboratory, 2022'))

XR_ssp=xr.Dataset(data_vars=dict(velocity_mean_ssp126=(['time','depth','lat'],np.mean(AR_ssp126,axis=0)),
                                 velocity_std_ssp126= (['time','depth','lat'],np.std(AR_ssp126,axis=0)),
                                 velocity_min_ssp126= (['time','depth','lat'],np.min(AR_ssp126,axis=0)),
                                 velocity_max_ssp126= (['time','depth','lat'],np.max(AR_ssp126,axis=0)),
                                 velocity_mean_ssp245=(['time','depth','lat'],np.mean(AR_ssp245,axis=0)),
                                 velocity_std_ssp245= (['time','depth','lat'],np.std(AR_ssp245,axis=0)),
                                 velocity_min_ssp245= (['time','depth','lat'],np.min(AR_ssp245,axis=0)),
                                 velocity_max_ssp245= (['time','depth','lat'],np.max(AR_ssp245,axis=0)),
                                 velocity_mean_ssp370=(['time','depth','lat'],np.mean(AR_ssp370,axis=0)),
                                 velocity_std_ssp370= (['time','depth','lat'],np.std(AR_ssp370,axis=0)),
                                 velocity_min_ssp370= (['time','depth','lat'],np.min(AR_ssp370,axis=0)),
                                 velocity_max_ssp370= (['time','depth','lat'],np.max(AR_ssp370,axis=0)),
                                 velocity_mean_ssp585=(['time','depth','lat'],np.mean(AR_ssp585,axis=0)),
                                 velocity_std_ssp585= (['time','depth','lat'],np.std(AR_ssp585,axis=0)),
                                 velocity_min_ssp585= (['time','depth','lat'],np.min(AR_ssp585,axis=0)),
                                 velocity_max_ssp585= (['time','depth','lat'],np.max(AR_ssp585,axis=0))),
                   coords=dict(time=years_scen,
                               depth=Depth,
                               lat=Lat),
                   attrs=dict(units='m/s',
                              description='annual averaged Atlantic Equatorial Undercurrent (AEU) velocity and related statistics (mean value, std, min, max), from averaged CMIP6 ensemble, scenario simulations (SSP1 2.6, SSP2 4.5, SSP3 7.0, SSP5 8.5, 2015-2100), values calculated between 3 deg N and 3 deg S and between 0 and 500m. values are derived by first averaging all ensemble members for each CMIP6 model and then taking the average of all models. Positive velocities are east-to-west',
                              author='Giovanni Galli, Plymouth Marine Laboratory, 2022'))


OUTDIR='/home/users/gig/CRACAB/FINAL_PLOT/'
OUTFILE_HIST=OUTDIR+'AEU_velocity_transect_historical_1950-2014.nc'
OUTFILE_SSP=OUTDIR+'AEU_velocity_transect_scenarios_1950-2014.nc'

XR_hist.to_netcdf(OUTFILE_HIST)
XR_ssp.to_netcdf(OUTFILE_SSP)
   



 
