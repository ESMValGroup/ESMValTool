import numpy as np
from scipy.io import loadmat
import h5py
from matplotlib.dates import num2date, set_epoch
from datetime import datetime as dt
import time

# loads data from Brandt et al. 2021, Fig 1a (annual flow timeseries),
# and Fig 1d (vertical profile),
# computes monthly average/max/min seasonality and raw data
# 

#filepath='/users/modellers/gig/Documents/CRACAB/Brandt_Scripts/Brandt_et_al_materials_v2/figure01/subplot_a/plotting_vars_subplot_a/plotting_vars_fig_1_a.mat'
#filepath='/home/users/gig/CRACAB/Brandt_Data/plotting_vars_fig_1_a.mat'

def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch
    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)
    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration
    return date.year + fraction

def year_average(monthly_data,timeline):
    days_month=np.array([31,28,31,30,31,30,31,31,30,31,30,31],dtype=np.float32)[np.newaxis,:]
    month0=timeline[0].month
    monthend=timeline[-1].month
    ext_data=np.append(np.nan*np.ones(month0-1),monthly_data)
    ext_data=np.append(ext_data,np.nan*np.ones(12-monthend))
    nyears=int(len(ext_data)/12)
    time_yr=np.arange(timeline[0].year,1+timeline[-1].year)
    mat_data=ext_data.reshape([nyears,12])
    days_year=np.sum(~np.isnan(mat_data)*days_month,axis=1)[:,np.newaxis] #sorry Feb 29!
    yr_av=np.nansum(mat_data*days_month)/days_year
    yr_av=np.array([np.nansum((mat_data*days_month)[x,:])/days_year[x,0] for x in range(nyears)])
    return yr_av,time_yr

def load_Brandt_Fig1a(filepath):
    set_epoch('-1-12-31T00:00')
    Data=loadmat(filepath)
    D1_mean=Data['data'][0,0]['euc_transport_mm_mean'][0,0]
    D1_anom=Data['data'][0,0]['euc_transport_mm_anom'][:,0] #raw monthly data '01-May-2005', '01-Oct-2019'
    time_aeu=Data['data'][0,0]['time_euc_mm'][:,0]
    timedate=num2date(time_aeu)
    timemonth=np.array([x.month for x in timedate])
    timeyear=np.array([toYearFraction(x) for x in timedate])
    D1_aeu=D1_mean+D1_anom
    D1_aeu_year_av,time_year_av=year_average(D1_aeu,timedate)
    D1_seas_av=np.array([np.mean(x) for x in [D1_aeu[y] for y in [np.where(timemonth==z) for z in np.arange(1,13)]]])
    D1_seas_max=np.array([np.max(x) for x in [D1_aeu[y] for y in [np.where(timemonth==z) for z in np.arange(1,13)]]])
    D1_seas_min=np.array([np.min(x) for x in [D1_aeu[y] for y in [np.where(timemonth==z) for z in np.arange(1,13)]]])
    return timeyear,D1_aeu,time_year_av,D1_aeu_year_av,D1_seas_av,D1_seas_max,D1_seas_min

def load_Brandt_Fig1d(filepath):
    Data={}
    with h5py.File(filepath,'r') as f:
        for k,v in f.items():
            Data[k]=np.array(v)
    Z=Data['P'][0:int(Data['depth_range'][0,0]),0]
    uo=np.nanmean(Data['U_monthly_median'],axis=0)
    return uo,Z
