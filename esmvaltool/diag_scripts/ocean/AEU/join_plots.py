import matplotlib.pyplot as plt
from matplotlib.cbook import get_sample_data
import numpy as np

INDIR='/home/users/gig/CRACAB/FINAL_PLOT/PLOTS_ALL_MODELS/'

OUTFILE=INDIR+'AEU_alltogether.png'

plot0=INDIR+'AEU_Ensemble_medium_None.png'
plot1=INDIR+'AEU_Ensemble_Seasonality_medium_None.png'
plot2=INDIR+'AEU_vertical_profile_hist_v_ssp.png'
plot3=INDIR+'legend.png'
plot4=INDIR+'AEU_uo_data_v_cmip6.png'
plot5=INDIR+'AEU_transects_hist_v_ssp_delta.png'

figwidth=11.65
figheight=8.27

f=open('/home/users/gig/CRACAB/FINAL_PLOT/lbwh_alltogether.txt')
lbwh=[]
for ll in f:
    dd=ll.split()
    lbwh.append([np.float32(i) for i in dd])

f.close()
lbwh=np.array(lbwh)

lbwh[:,0]=lbwh[:,0]/figwidth
lbwh[:,1]=lbwh[:,1]/figheight
lbwh[:,2]=lbwh[:,2]/figwidth
lbwh[:,3]=lbwh[:,3]/figheight


fig=plt.figure(figsize=(figwidth,figheight))


im0=plt.imread(get_sample_data(plot0))
ax0=plt.axes(lbwh[0])
ax0.imshow(im0)
ax0.axis('off')

im1=plt.imread(get_sample_data(plot1))
ax1=plt.axes(lbwh[1])
ax1.imshow(im1)
ax1.axis('off')

im2=plt.imread(get_sample_data(plot2))
ax2=plt.axes(lbwh[2])
ax2.imshow(im2)
ax2.axis('off')

im3=plt.imread(get_sample_data(plot3))
ax3=plt.axes(lbwh[3])
ax3.imshow(im3)
ax3.axis('off')

im4=plt.imread(get_sample_data(plot4))
ax4=plt.axes(lbwh[4])
ax4.imshow(im4)
ax4.axis('off')

im5=plt.imread(get_sample_data(plot5))
ax5=plt.axes(lbwh[5])
ax5.imshow(im5)
ax5.axis('off')

fig.suptitle('AEU flow\n Historical (2000-2010) vs SSP (2040-2050)')

print('saving: '+OUTFILE)
plt.savefig(OUTFILE,dpi=300)

