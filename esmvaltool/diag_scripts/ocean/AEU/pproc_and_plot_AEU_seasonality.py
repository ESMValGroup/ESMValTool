import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import glob
import yaml
import netCDF4 as nc4
import Brandt_flow_data

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
Brandt_file='/home/users/gig/CRACAB/Brandt_Data/plotting_vars_fig_1_a.mat'

B_time,B_aeu,B_time_yr_av,B_aeu_yr_av,B_seas_av,B_seas_max,B_seas_min=Brandt_flow_data.load_Brandt_Fig1a(Brandt_file)

# reads into a list, each element of which is a dict
with open(FLIST,'rb') as f:
    Stuff=yaml.safe_load(f)

# type of preproc
PREPROC='medium_None'

INDIR='/home/users/gig/CRACAB/ENSEMBLE_RUN_monthly_all/esmvaltool_output/*/plots/diag_aeu/AMOC_timeseries/'

# reconstruct indexes of hist and ssp time-periods for monthly data
year_start_hist=2000
year_end_hist=2010
year_start_ssp=2040
year_end_ssp=2050

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
                curr_TS_y=np.mean(curr_TS_x.reshape(-1,12),axis=0)
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

#2.3 Make another dictionary with the average of the scenarios
EXPs=['historical','ssp126','ssp245','ssp370','ssp585']
Data_av=dict.fromkeys(EXPs,0.0)
exp_count=dict.fromkeys(EXPs,0)
for dset in Data.keys():
    for exp in Data[dset].keys():
        if exp!='ssp119' and len(Data[dset][exp])!=0:
            Data_av[exp]=Data_av[exp]+Data[dset][exp]
            exp_count[exp]=exp_count[exp]+1

for ee in EXPs:
    if exp_count[ee]!=0:
        Data_av[ee]=Data_av[ee]/exp_count[ee]
    


#3 Plot Results

#colors=['blue','green','gold','orange','red'] #colors of each scenario, in order (Lee's colors)
colors={'historical': (35, 87, 145), # new darkish blue
        'ssp126': (0, 173, 207), # light sky
        'ssp245': (247, 148, 32), # orange
        'ssp370': (231, 29, 37), # red
        'ssp585': (149, 27, 30), } # dark maroon red
colors = {k: (v1/256., v2/256.,v3/256.) for k, (v1,v2,v3) in colors.items()}
color_factor = 0.75
colors_dark = {k: (color_factor*v1, color_factor*v2, color_factor* v3) for k, (v1,v2,v3) in colors.items()}

years_hist=np.arange(2000,2010)
years_scen=np.arange(2040,2050)

# read lbwh
figwidth=3.35
figheight=1.57

f=open('/home/users/gig/CRACAB/FINAL_PLOT/lbwh_season.txt')
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


for exp_n in range(len(EXPs)):
    exp=EXPs[exp_n]
    ax.plot(np.arange(12)+0.5,Data_av[exp],color=colors[exp])

ax.fill_between(np.arange(12)+0.5,B_seas_min,B_seas_max,color='k',alpha=0.5,linewidth=0.0)

ax.plot(np.arange(12)+0.5,B_seas_av,color='k')

#custom_legend=[Line2D([0],[0],color='blue',lw=1),
#               Line2D([0],[0],color='green',lw=1),
#               Line2D([0],[0],color='gold',lw=1),
#               Line2D([0],[0],color='orange',lw=1),
#               Line2D([0],[0],color='red',lw=1),
#               Line2D([0],[0],color='k',lw=1)]
#
#plt.legend(custom_legend,EXPs.append('observations'))

ax.set_ylabel('AEU [Sv]')
ax.set_title(PREPROC)
ax.set_title('Climatology')

ax.set_xlim(0,12)

months=('J','F','M','A','M','J','J','A','S','O','N','D')
NoLogo=np.array(['','','','','','','','','','','',''])

ax.set_xticks(np.arange(12),minor=False)
ax.set_xticks(np.arange(12)+0.5,minor=True)
ax.set_xticklabels(NoLogo,minor=False)
ax.set_xticklabels(months,minor=True)
ax.tick_params(which='minor',length=0)
ax.text(0.5,25.0,'b)')

plt.grid(axis='y',linestyle='--',linewidth=0.5)

OUTFILE='/home/users/gig/CRACAB/FINAL_PLOT/AEU_Ensemble_Seasonality_'+PREPROC+'.png'
print('OUTFILE:  '+OUTFILE)

plt.savefig(OUTFILE,dpi=300)

