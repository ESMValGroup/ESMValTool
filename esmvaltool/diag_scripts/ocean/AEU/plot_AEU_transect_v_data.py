import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle
import scipy.io

OUTFILE='/home/users/gig/CRACAB/FINAL_PLOT/AEU_uo_data_v_cmip6.png'

# CMIP6 average data
with open('/home/users/gig/CRACAB/AEU_Ensemble_Averages_uo_transect.pkl','rb') as f: #ALL_MODELS
#with open('/home/users/gig/CRACAB/AEU_Ensemble_Averages_uo_transect_best_models.pkl','rb') as f: #BEST_MODELS
    CMIP6_data=pickle.load(f)

Yc=CMIP6_data['Lat']
Zc=CMIP6_data['Depth']

uo_cmip6=CMIP6_data['historical']

cmip6_time=np.arange(1950,2015)
Brandt_time=(2008,2018)

uo_cmip6_av=np.mean(uo_cmip6[np.where(cmip6_time==Brandt_time[0])[0][0]:,:,:],axis=0)


# AEU reconstructed from ship measurements (Brandt et al. 2020), 
# retrieved from https://zenodo.org/record/4478285#.YWRAIOrML0M
Brandt_Data=scipy.io.loadmat('/home/users/gig/CRACAB/BRANDT_AEU/uvmean_ta_central.mat')

u_mean=np.array(Brandt_Data['u_mean'])/100.0 #cm/s -> m/s
pos_grid=np.array(Brandt_Data['pos_grid'])
depth_grid=np.array(Brandt_Data['depth_grid'])

Yb,Zb=np.meshgrid(pos_grid[0,:],depth_grid[0,:])

# plot parameters
cmax=np.max([np.max(u_mean),np.max(uo_cmip6_av)])
cmin=-cmax

latmin=-3.0
latmax=3.0
depthmin=0.0
depthmax=300.0


# Plot the two things
figwidth=4.92
figheight=3.98

lbwh=np.loadtxt('/home/users/gig/CRACAB/FINAL_PLOT/lbwh_trsct_v_dta.txt')

lbwh[:,0]=lbwh[:,0]/figwidth
lbwh[:,1]=lbwh[:,1]/figheight
lbwh[:,2]=lbwh[:,2]/figwidth
lbwh[:,3]=lbwh[:,3]/figheight

ffs=8
matplotlib.rcParams.update({'font.size': ffs})

Fig=plt.figure(figsize=(figwidth,figheight))

ax0=plt.axes(lbwh[1,:]) #NB: I swapped the axes!
im0=ax0.pcolormesh(Yb,Zb,u_mean,vmin=cmin,vmax=cmax,cmap='seismic')

ax0.set_title('observations (2008-2018)')
ax0.set_xlim(latmin,latmax)
ax0.set_ylim(depthmin,depthmax)
ax0.invert_yaxis()
ax0.set_ylabel('m')
ax0.set_xlabel('lat')
ax0.text(-2.8,25,'f)')

ax1=plt.axes(lbwh[0,:]) #NB: I swapped the axes!
im1=ax1.pcolormesh(Yc,Zc,uo_cmip6_av,vmin=cmin,vmax=cmax,cmap='seismic')

ax1.set_title('historical (2008-2015)')
ax1.set_xlim(latmin,latmax)
ax1.set_ylim(depthmin,depthmax)
ax1.invert_yaxis()
ax1.set_ylabel('m')
#ax1.set_xlabel('lat')
ax1.text(-2.8,25,'e)')


ax_cbr=plt.axes(lbwh[2,:])
cbar=Fig.colorbar(im1,cax=ax_cbr,orientation='vertical')
ax_cbr.set_ylabel(u'$AEU\ u\ velocity\ [m s^{-1}]$')

#Fig.suptitle(u'$AEU\ u\ velocity\ [m s^{-1}]$',fontweight='bold')

print('OUTFILE: '+OUTFILE)
plt.savefig(OUTFILE)
