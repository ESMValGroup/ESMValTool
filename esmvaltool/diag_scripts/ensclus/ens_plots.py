#*********************************
#           ens_plots            *
#*********************************

# Standard packages
import os
import sys
from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import math

def ens_plots(dir_OUTPUT,dir_PLOT,name_outputs,numclus,field_to_plot):
    '''
    \nGOAL:
    Plot the chosen field for each ensemble
    NOTE:
    '''
    
    # User-defined libraries
    from read_netcdf import read_N_2Dfields

    tit=field_to_plot
    print('Number of clusters: {0}'.format(numclus))

    print(name_outputs)
    varname=name_outputs.split("_")[0]
    kind=name_outputs.split("_")[-1]

    #____________Reading the netCDF file of N 2Dfields of anomalies, saved by ens_anom.py
    ifile=os.path.join(dir_OUTPUT,'ens_anomalies_{0}.nc'.format(name_outputs))
    vartoplot, varunits, lat, lon = read_N_2Dfields(ifile)
    print('vartoplot dim: (numens x lat x lon)={0}'.format(vartoplot.shape))
    numens=vartoplot.shape[0]

    #____________Load labels
    namef=os.path.join(dir_OUTPUT,'labels_{0}.txt'.format(name_outputs))
    labels=np.loadtxt(namef,dtype=int)
    print(labels)

    mi=vartoplot.min()
    ma=vartoplot.max()

    if field_to_plot=='anomalies':
        # compute range colorbar for anomalies
        delta=0.05
        if abs(math.floor(mi*100)/100)<math.ceil(ma*100)/100:
            rangecbarmin=-math.ceil(ma*100)/100
            rangecbarmax=math.ceil(ma*100)/100+delta
        else:
            rangecbarmin=math.floor(mi*100)/100
            rangecbarmax=abs(math.floor(mi*100)/100)+delta
    else:
        # compute range colorbar for climatologies
        delta=0.2
        rangecbarmin=math.floor(mi)
        rangecbarmax=math.ceil(ma)+delta
    
    clevels=np.arange(rangecbarmin,rangecbarmax,delta)
    #clevels=np.arange(2,44,delta)
    #clevels=np.arange(-0.7,0.75,delta)
    #clevels=np.arange(0,6.2,delta)
    
    colors=['b','g','r','c','m','y','DarkOrange','grey']
    m = Basemap(projection='cyl',llcrnrlat=min(lat),urcrnrlat=max(lat),llcrnrlon=min(lon),urcrnrlon=max(lon),resolution='c')
    
    x=int(np.ceil(np.sqrt(numens*1.6)))
    y=int(np.ceil(numens/x))
    print(x,y)

    fig = plt.figure(figsize=(24,14))
    for nens in range(numens):
        #print('//////////ENSEMBLE MEMBER {0}'.format(nens))
        ax = plt.subplot(x, y, nens+1)
    
        m.drawcoastlines()
        # use meshgrid to create 2D arrays
        x_i,y_i=np.meshgrid(lon,lat)
        xi,yi=m(x_i,y_i)
        
        # Plot Data
        if field_to_plot=='anomalies':
            map_plot=m.contourf(xi,yi,vartoplot[nens],clevels,cmap=plt.cm.RdBu_r)
        else:
            map_plot=m.contourf(xi,yi,vartoplot[nens],clevels)
        #print('min={0}'.format(vartoplot[nens].min()))
        #print('max={0}\n'.format(vartoplot[nens].max()))
    
        # Add Title
        subtit = nens
        title_obj=plt.title(subtit, fontsize=32, fontweight='bold')
        for nclus in range(numclus):
            if nens in np.where(labels==nclus)[0]:
                title_obj.set_backgroundcolor(colors[nclus])
    
    cax = plt.axes([0.1, 0.03, 0.8, 0.03]) #horizontal
    cb=plt.colorbar(map_plot,cax=cax, orientation='horizontal')#, labelsize=18)
    cb.ax.tick_params(labelsize=18)
    
    plt.suptitle(kind+' '+varname+' '+tit+' ('+varunits+')', fontsize=45, fontweight='bold')
    #plt.suptitle(kind+' (2038-2068) JJA 2m temperature '+tit+' (degC)', fontsize=45, fontweight='bold')
    #plt.suptitle(kind+' (1979-2008) JJA 2m temperature '+tit+' (degC)', fontsize=45, fontweight='bold')
    #plt.suptitle(kind+' (1979-2008) JJA precipitation rate '+tit+' (mm/day)', fontsize=45, fontweight='bold')
    
    plt.subplots_adjust(top=0.85)
    top    = 0.89  # the top of the subplots of the figure
    bottom = 0.12    # the bottom of the subplots of the figure
    left   = 0.02    # the left side of the subplots of the figure
    right  = 0.98  # the right side of the subplots of the figure
    hspace = 0.36   # the amount of height reserved for white space between subplots
    wspace = 0.14    # the amount of width reserved for blank space between subplots
    plt.subplots_adjust(left=left, bottom=bottom, right=right, top=top, wspace=wspace, hspace=hspace)
    
    # plot the selected fields
    namef=os.path.join(dir_PLOT,'{0}_{1}.eps'.format(field_to_plot,name_outputs))
    fig.savefig(namef)#bbox_inches='tight')
    print('An eps figure for the selected fields is saved in {0}'.format(dir_OUTPUT))
    print('____________________________________________________________________________________________________________________')

    return


#========================================================

if __name__ == '__main__':
    print('This program is being run by itself')
    
    print('**************************************************************')
    print('Running {0}'.format(sys.argv[0]))
    print('**************************************************************')
    dir_OUTPUT    = sys.argv[1]          # OUTPUT DIRECTORY
    dir_PLOT      = sys.argv[2]          # OUTPUT PLOT DIRECTORY
    name_outputs  = sys.argv[3]          # name of the outputs
    numclus       = int(sys.argv[4])     # number of clusters
    field_to_plot = sys.argv[5]          #field to plot ('climatologies', 'anomalies', '75th_percentile', 'mean', 'maximum', 'std', 'trend')

    ens_plots(dir_OUTPUT,dir_PLOT,name_outputs,numclus,field_to_plot)

else:
    print('ens_plots is being imported from another module')

