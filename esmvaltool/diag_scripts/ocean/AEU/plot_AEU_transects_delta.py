import numpy as np
from numpy.ma.core import masked_where as MW
import matplotlib
import matplotlib.pyplot as plt
import pickle
import scipy.io

OUTFILE='/home/users/gig/CRACAB/FINAL_PLOT/AEU_transects_hist_v_ssp_delta.png'

# CMIP6 average data
with open('/home/users/gig/CRACAB/AEU_Ensemble_Averages_uo_transect.pkl','rb') as f:
    CMIP6_data=pickle.load(f)

Yc=CMIP6_data['Lat']
Zc=CMIP6_data['Depth']

uo_hist_all=CMIP6_data['historical']
#uo_119_all=CMIP6_data['ssp119']
uo_126_all=CMIP6_data['ssp126']
uo_245_all=CMIP6_data['ssp245']
uo_370_all=CMIP6_data['ssp370']
uo_585_all=CMIP6_data['ssp585']

hist_time=np.arange(1950,2015)
ssp_time=np.arange(2015,2101)

years_hist=(2000,2010)
years_ssp=(2040,2050)

i_hist=[np.where(hist_time==years_hist[x])[0][0] for x in range(2)]
i_ssp=[np.where(ssp_time==years_ssp[x])[0][0] for x in range(2)]

uo_hist=np.mean(uo_hist_all[i_hist[0]:i_hist[1],:,:],axis=0)
#uo_119=np.mean(uo_119_all[i_ssp[0]:i_ssp[1],:,:],axis=0)
uo_126=np.mean(uo_126_all[i_ssp[0]:i_ssp[1],:,:],axis=0)
uo_245=np.mean(uo_245_all[i_ssp[0]:i_ssp[1],:,:],axis=0)
uo_370=np.mean(uo_370_all[i_ssp[0]:i_ssp[1],:,:],axis=0)
uo_585=np.mean(uo_585_all[i_ssp[0]:i_ssp[1],:,:],axis=0)

#delta_119=MW(uo_119<0.0,uo_119)-uo_hist
delta_126=MW(uo_126<0.0,uo_126)-uo_hist
delta_245=MW(uo_245<0.0,uo_245)-uo_hist
delta_370=MW(uo_370<0.0,uo_370)-uo_hist
delta_585=MW(uo_585<0.0,uo_585)-uo_hist

#Brandt_time=(2008,2018)

#uo_cmip6_av=np.mean(uo_cmip6[np.where(cmip6_time==Brandt_time[0])[0][0]:,:,:],axis=0)


## AEU reconstructed from ship measurements (Brandt et al. 2020), 
## retrieved from https://zenodo.org/record/4478285#.YWRAIOrML0M
#Brandt_Data=scipy.io.loadmat('/home/users/gig/CRACAB/BRANDT_AEU/uvmean_ta_central.mat')
#
#u_mean=np.array(Brandt_Data['u_mean'])/100.0 #cm/s -> m/s
#pos_grid=np.array(Brandt_Data['pos_grid'])
#depth_grid=np.array(Brandt_Data['depth_grid'])
#
#Yb,Zb=np.meshgrid(pos_grid[0,:],depth_grid[0,:])

# plot parameters
#cmax_hist=np.max(uo_hist)
#cmin_hist=-cmax_hist

cmax_delta=np.max([np.max(delta_126),np.max(delta_245),np.max(delta_370),np.max(delta_585)])
cmin_delta=np.min([np.min(delta_126),np.min(delta_245),np.min(delta_370),np.min(delta_585)])
cmax_delta=np.max([cmax_delta, np.abs(cmin_delta)])
cmin_delta=-cmax_delta

latmin=-3.0
latmax=3.0
depthmin=0.0
depthmax=300.0

colors={'historical': (35, 87, 145), # new darkish blue
        'ssp126': (0, 173, 207), # light sky
        'ssp245': (247, 148, 32), # orange
        'ssp370': (231, 29, 37), # red
        'ssp585': (149, 27, 30), } # dark maroon red
colors = {k: (v1/256., v2/256.,v3/256.) for k, (v1,v2,v3) in colors.items()}

# Plot the two things
figwidth=6.73
figheight=3.98

lbwh=np.loadtxt('/home/users/gig/CRACAB/FINAL_PLOT/lbwh_transects_delta.txt')

lbwh[:,0]=lbwh[:,0]/figwidth
lbwh[:,1]=lbwh[:,1]/figheight
lbwh[:,2]=lbwh[:,2]/figwidth
lbwh[:,3]=lbwh[:,3]/figheight

ffs=8
matplotlib.rcParams.update({'font.size': ffs})

Fig=plt.figure(figsize=(figwidth,figheight))

#ax0=plt.axes(lbwh[0,:])
#im0=ax0.pcolormesh(Yc,Zc,uo_hist,vmin=cmin_hist,vmax=cmax_hist,cmap=BrBG)
#
#ax0.set_title('historical (2000-2010)')
#ax0.set_xlim(latmin,latmax)
#ax0.set_ylim(depthmin,depthmax)
#ax0.invert_yaxis()
#ax0.set_ylabel('depth [m]')
#ax0.set_xlabel('lat')


#ax1=plt.axes(lbwh[1,:])
#im1=ax1.pcolormesh(Yc,Zc,delta_119,vmin=cmin_delta,vmax=cmax_delta,cmap=BrBG)
#
#ax1.set_title('ssp119 (2040-2050)')
#ax1.set_xlim(latmin,latmax)
#ax1.set_ylim(depthmin,depthmax)
#ax1.invert_yaxis()
#ax1.set_xlabel('lat')


ax2=plt.axes(lbwh[0,:])
im2=ax2.pcolormesh(Yc,Zc,delta_126,vmin=cmin_delta,vmax=cmax_delta,cmap='BrBG')

ax2.set_title('SSP1-2.6',color=colors['ssp126'],fontweight='bold')
ax2.set_xlim(latmin,latmax)
ax2.set_ylim(depthmin,depthmax)
ax2.invert_yaxis()
ax2.set_ylabel('m')
ax2.set_xticklabels([])
ax2.text(-2.5,25.0,'g)')


ax3=plt.axes(lbwh[1,:])
im3=ax3.pcolormesh(Yc,Zc,delta_245,vmin=cmin_delta,vmax=cmax_delta,cmap='BrBG')

ax3.set_title('SSP2-4.5',color=colors['ssp245'],fontweight='bold')
ax3.set_xlim(latmin,latmax)
ax3.set_ylim(depthmin,depthmax)
ax3.invert_yaxis()
ax3.set_xticklabels([])
ax3.set_yticklabels([])
ax3.text(-2.5,25.0,'h)')


ax4=plt.axes(lbwh[2,:])
im4=ax4.pcolormesh(Yc,Zc,delta_370,vmin=cmin_delta,vmax=cmax_delta,cmap='BrBG')

ax4.set_title('SSP3-7.0',color=colors['ssp370'],fontweight='bold')
ax4.set_xlim(latmin,latmax)
ax4.set_ylim(depthmin,depthmax)
ax4.invert_yaxis()
ax4.set_xlabel('lat')
ax4.set_ylabel('m')
ax4.text(-2.5,25.0,'i)')


ax5=plt.axes(lbwh[3,:])
im5=ax5.pcolormesh(Yc,Zc,delta_585,vmin=cmin_delta,vmax=cmax_delta,cmap='BrBG')

ax5.set_title('SSP5-8.5',color=colors['ssp585'],fontweight='bold')
ax5.set_xlim(latmin,latmax)
ax5.set_ylim(depthmin,depthmax)
ax5.invert_yaxis()
ax5.set_xlabel('lat')
ax5.set_yticklabels([])
ax5.text(-2.5,25.0,'j)')


ax_cbr1=plt.axes(lbwh[4,:])
cbar1=Fig.colorbar(im2,cax=ax_cbr1,orientation='vertical')
ax_cbr1.set_ylabel(u'$AEU\ \Delta u\ velocity\ [m s^{-1}]$')

#ax_cbr2=plt.axes(lbwh[7,:])
#cbar2=Fig.colorbar(im1,cax=ax_cbr2,orientation='horizontal')

#Fig.suptitle(u'$AEU\ \Delta u\ velocity\ [m s^{-1}]$',fontweight='bold')

print('OUTFILE: '+OUTFILE)
plt.savefig(OUTFILE)
