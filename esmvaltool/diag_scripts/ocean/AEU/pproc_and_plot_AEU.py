import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import glob
import yaml
import netCDF4 as nc4
import Brandt_flow_data

# NB, a bunch od staff will be passed with argparse:
#     FLIST, PREPROC, 

#1. Read input file with the list of processed models
FLIST='/home/users/gig/CRACAB/FIND_FILES/output_match_hist_and_ssp.yml' #ALL MODELS
#FLIST='/home/users/gig/CRACAB/FIND_FILES/output_best_models.yml' #BEST MODELS

Brandt_file='/home/users/gig/CRACAB/Brandt_Data/plotting_vars_fig_1_a.mat'
B_time,B_aeu,B_time_yr_av,B_aeu_yr_av,B_seas_av,B_seas_max,B_seas_min=Brandt_flow_data.load_Brandt_Fig1a(Brandt_file)


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

#INDIR='/home/users/gig/CRACAB/ENSEMBLE_RUN_3/esmvaltool_output/*/plots/diag_aeu/AMOC_timeseries/'
#INDIR_PRE='/home/users/gig/CRACAB/ENSEMBLE_RUN_3/esmvaltool_output/*/preproc/diag_aeu/aeu/'

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
        # Uncomment to plot just one ensemble member (here r1i1p1f1)
        #curr_flist_preproc=sorted(glob.glob(INDIR_PRE+'CMIP6_'+dset+'_Omon_'+exp+'_*_uo_*.nc'))
        #curr_flist=[curr_flist[x] for x in range(len(curr_flist)) if 'r1i1p1f1' in curr_flist_preproc[x]]
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

#2.3 Make another dictionary with the average of the scenarios
#    Exclude ssp119 because little data
#EXPs=['historical','ssp119','ssp126','ssp245','ssp370','ssp585']
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
    
#2.5 compute mean and std of scenarios
ALL_DATA=dict.fromkeys(EXPs)

for ee,exp in enumerate(EXPs):
    ALL_DATA[exp]=np.zeros([exp_count[exp],len(Data['UKESM1-0-LL'][exp])])
    count=0
    for dset in Data.keys():
        if exp in Data[dset]:
            if len(Data[dset][exp])!=0:
                ALL_DATA[exp][count,:]=Data[dset][exp]
                count=count+1

print('mean,std historical: '+str(np.mean(ALL_DATA['historical']))+', '+str(np.std(ALL_DATA['historical'])))
print('mean,std ssp126: '+str(np.mean(ALL_DATA['ssp126']))+', '+str(np.std(ALL_DATA['ssp126'])))
print('mean,std ssp245: '+str(np.mean(ALL_DATA['ssp245']))+', '+str(np.std(ALL_DATA['ssp245'])))
print('mean,std ssp370: '+str(np.mean(ALL_DATA['ssp370']))+', '+str(np.std(ALL_DATA['ssp370'])))
print('mean,std ssp585: '+str(np.mean(ALL_DATA['ssp585']))+', '+str(np.std(ALL_DATA['ssp585'])))


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

years_hist=np.arange(1950,2015)
years_scen=np.arange(2015,2101)

# read lbwh
figwidth=6.50
figheight=3.39

f=open('/home/users/gig/CRACAB/FINAL_PLOT/lbwh_timeseries.txt')
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

for dset in Data.keys():
    for exp_n in range(len(EXPs)):
        exp=EXPs[exp_n]
        try:
            curr_data=Data[dset][exp]
            if exp=='historical':
                ax.plot(years_hist,curr_data,color=colors_dark[exp],alpha=0.2)
            else:
                ax.plot(years_scen,curr_data,color=colors_dark[exp],alpha=0.2)
        except:
            pass

for exp_n in range(len(EXPs)):
    exp=EXPs[exp_n]
    if exp=='historical':
        ax.plot(years_hist,Data_av[exp],color='k',lw=2.5)
        ax.plot(years_hist,Data_av[exp],color=colors[exp],lw=2)
    else:
        ax.plot(years_scen,Data_av[exp],color='k',lw=2.5)
        ax.plot(years_scen,Data_av[exp],color=colors[exp],lw=2)

# plot Brandt Data
ax.plot(B_time_yr_av,B_aeu_yr_av,color='k')

#custom_legend=[Line2D([0],[0],color='b',lw=1),
#               Line2D([0],[0],color='c',lw=1),
#               Line2D([0],[0],color='g',lw=1),
#               Line2D([0],[0],color='olive',lw=1),
#               Line2D([0],[0],color='r',lw=1),
#               Line2D([0],[0],color='purple',lw=1)]
#
#plt.legend(custom_legend,EXPs)

#plot shaded historical and ssp periods
ymin=ax.get_ylim()[0]
ymax=ax.get_ylim()[1]
ax.fill_between([2000,2010],ymin,ymax,color='k',alpha=0.2,linewidth=0.0)
ax.fill_between([2040,2050],ymin,ymax,color='purple',alpha=0.2,linewidth=0.0)

custom_legend=[Line2D([0],[0],color='k',lw=4,alpha=0.2),
               Line2D([0],[0],color='purple',lw=4,alpha=0.2)]

plt.legend(custom_legend,['Historical period','SSP period'],frameon=False)

ax.set_ylabel('AEU [Sv]')
ax.set_title('Time Series')

ax.set_xlim(1950,2100)
ax.set_ylim(ymin,ymax)

ax.text(1960,22,'a)')

plt.grid(axis='y',linestyle='--',linewidth=0.5)

OUTFILE='/home/users/gig/CRACAB/FINAL_PLOT/AEU_Ensemble_'+PREPROC+'.png'
print('OUTFILE:  '+OUTFILE)

plt.savefig(OUTFILE,dpi=300)

