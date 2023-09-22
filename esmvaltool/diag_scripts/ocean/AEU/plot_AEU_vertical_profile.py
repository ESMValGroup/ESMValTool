import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pickle
import scipy.io
import Brandt_flow_data

OUTFILE='FINAL_PLOT/AEU_vertical_profile_hist_v_ssp.png'

Brandt_file='Brandt_Data/plotting_vars_fig_1_d.mat'
Brandt_uo,Brandt_z=Brandt_flow_data.load_Brandt_Fig1d(Brandt_file)

# CMIP6 average data
with open('AEU_Ensemble_Averages_uo_transect.pkl','rb') as f: #ALL MODELS
#with open('AEU_Ensemble_Averages_uo_transect_best_models.pkl','rb') as f: #BEST MODELS
    CMIP6_data=pickle.load(f)

Yc=CMIP6_data['Lat']
Zc=CMIP6_data['Depth']
Zc_1d=Zc[:,0]

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

y_equator=(9,10) #the two cells equally closer to the equator, north and south

uo_hist=np.mean(uo_hist_all[i_hist[0]:i_hist[1],:,:],axis=0)
#uo_119=np.mean(uo_119_all[i_ssp[0]:i_ssp[1],:,:],axis=0)
uo_126=np.mean(uo_126_all[i_ssp[0]:i_ssp[1],:,:],axis=0)
uo_245=np.mean(uo_245_all[i_ssp[0]:i_ssp[1],:,:],axis=0)
uo_370=np.mean(uo_370_all[i_ssp[0]:i_ssp[1],:,:],axis=0)
uo_585=np.mean(uo_585_all[i_ssp[0]:i_ssp[1],:,:],axis=0)

uo_hist=np.mean(uo_hist[:,y_equator[0]:y_equator[1]],axis=1)
#uo_119=np.mean(uo_119[:,y_equator[0]:y_equator[1]],axis=1)
uo_126=np.mean(uo_126[:,y_equator[0]:y_equator[1]],axis=1)
uo_245=np.mean(uo_245[:,y_equator[0]:y_equator[1]],axis=1)
uo_370=np.mean(uo_370[:,y_equator[0]:y_equator[1]],axis=1)
uo_585=np.mean(uo_585[:,y_equator[0]:y_equator[1]],axis=1)

delta_uo_126=uo_126-uo_hist
delta_uo_245=uo_245-uo_hist
delta_uo_370=uo_370-uo_hist
delta_uo_585=uo_585-uo_hist

#Brandt_time=(2008,2018)

#uo_cmip6_av=np.mean(uo_cmip6[np.where(cmip6_time==Brandt_time[0])[0][0]:,:,:],axis=0)


## AEU reconstructed from ship measurements (Brandt et al. 2020), 
## retrieved from https://zenodo.org/record/4478285#.YWRAIOrML0M
#Brandt_Data=scipy.io.loadmat('BRANDT_AEU/uvmean_ta_central.mat')
#
#u_mean=np.array(Brandt_Data['u_mean'])/100.0 #cm/s -> m/s
#pos_grid=np.array(Brandt_Data['pos_grid'])
#depth_grid=np.array(Brandt_Data['depth_grid'])
#
#Yb,Zb=np.meshgrid(pos_grid[0,:],depth_grid[0,:])

# plot parameters
cmax=np.max([np.max(uo_hist),np.max(uo_126),np.max(uo_245),np.max(uo_370),np.max(uo_585)])
cmin=-cmax

depthmin=0.0
depthmax=300.0

colors={'historical': (35, 87, 145), # new darkish blue
        'ssp126': (0, 173, 207), # light sky
        'ssp245': (247, 148, 32), # orange
        'ssp370': (231, 29, 37), # red
        'ssp585': (149, 27, 30), } # dark maroon red
colors = {k: (v1/256., v2/256.,v3/256.) for k, (v1,v2,v3) in colors.items()}
color_factor = 0.75
colors_dark = {k: (color_factor*v1, color_factor*v2, color_factor* v3) for k, (v1,v2,v3) in colors.items()}


# Plot the two things
figwidth=3.35
figheight=1.85

lbwh=np.loadtxt('FINAL_PLOT/lbwh_profile.txt')

lbwh[:,0]=lbwh[:,0]/figwidth
lbwh[:,1]=lbwh[:,1]/figheight
lbwh[:,2]=lbwh[:,2]/figwidth
lbwh[:,3]=lbwh[:,3]/figheight


ffs=8
matplotlib.rcParams.update({'font.size': ffs})

Fig=plt.figure(figsize=(figwidth,figheight))

ax0=plt.axes(lbwh[0,:])
#im0=ax0.pcolormesh(Yc,Zc,uo_hist,vmin=cmin,vmax=cmax,cmap='seismic')
im0=ax0.plot(uo_hist,Zc_1d,color=colors['historical'])
#im1=ax0.plot(uo_119,Zc_1d,color='c')
im2=ax0.plot(uo_126,Zc_1d,color=colors['ssp126'])
im3=ax0.plot(uo_245,Zc_1d,color=colors['ssp245'])
im4=ax0.plot(uo_370,Zc_1d,color=colors['ssp370'])
im5=ax0.plot(uo_585,Zc_1d,color=colors['ssp585'])

im6=ax0.plot(np.zeros(len(Zc_1d)),Zc_1d,'--',color=[0.5,0.5,0.5])

im7=ax0.plot(Brandt_uo,Brandt_z,color='k')

ax0.set_title('Profile')
#ax0.set_xlim(latmin,latmax)
ax0.set_ylim(depthmin,depthmax)
ax0.invert_yaxis()
ax0.set_ylabel(u'$m$')
ax0.set_xlabel(u'$uo\ [ms^{-1}]$')
ax0.text((0.1*(ax0.get_xlim()[1]-ax0.get_xlim()[0])+ax0.get_xlim()[0]),40,'c)')

ax1=plt.axes(lbwh[1,:])
im2d=ax1.plot(delta_uo_126,Zc_1d,color=colors['ssp126'])
im3d=ax1.plot(delta_uo_245,Zc_1d,color=colors['ssp245'])
im4d=ax1.plot(delta_uo_370,Zc_1d,color=colors['ssp370'])
im5d=ax1.plot(delta_uo_585,Zc_1d,color=colors['ssp585'])
im6d=ax1.plot(np.zeros(len(Zc_1d)),Zc_1d,'--',color=[0.5,0.5,0.5])

ax1.set_title('Difference')
ax1.set_ylim(depthmin,depthmax)
ax1.invert_yaxis()
ax1.set_yticklabels([])
ax1.set_xlabel(u'$\Delta uo\ [ms^{-1}]$')
ax1.text((0.1*(ax1.get_xlim()[1]-ax1.get_xlim()[0])+ax1.get_xlim()[0]),40,'d)')


print('OUTFILE: '+OUTFILE)
plt.savefig(OUTFILE)
