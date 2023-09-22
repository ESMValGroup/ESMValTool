import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

OUTFILE='/home/users/gig/CRACAB/FINAL_PLOT/legend.png'

colors={'historical': (35, 87, 145), # new darkish blue
        'ssp126': (0, 173, 207), # light sky
        'ssp245': (247, 148, 32), # orange
        'ssp370': (231, 29, 37), # red
        'ssp585': (149, 27, 30), } # dark maroon red
colors = {k: (v1/256., v2/256.,v3/256.) for k, (v1,v2,v3) in colors.items()}

custom_legend=[Line2D([0],[0],color=colors['historical'],lw=2),
               Line2D([0],[0],color=colors['ssp126'],lw=2),
               Line2D([0],[0],color=colors['ssp245'],lw=2),
               Line2D([0],[0],color=colors['ssp370'],lw=2),
               Line2D([0],[0],color=colors['ssp585'],lw=2),
               Line2D([0],[0],color='black',lw=2),
               Line2D([0],[0],color='black',lw=4,alpha=0.5)]

legend_entries=['historical','SSP1-2.6','SSP2-4.5','SSP3-7.0','SSP5-8.5','observations','obs. range']

# read lbwh
figwidth=1.85
figheight=3.19

f=open('/home/users/gig/CRACAB/FINAL_PLOT/lbwh_legend.txt')
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

ffs=8
matplotlib.rcParams.update({'font.size': ffs})

Fig=plt.figure(figsize=(figwidth, figheight))
ax=plt.axes(lbwh[0,:])

plt.axis('off')
plt.legend(custom_legend,legend_entries,frameon=False,bbox_to_anchor=lbwh[0,:],loc='center')
#plt.legend(custom_legend,legend_entries,frameon=False,loc='center')
#plt.legend(custom_legend,legend_entries,frameon=False)

print('OUTFILE:  '+OUTFILE)

plt.savefig(OUTFILE,dpi=300)
