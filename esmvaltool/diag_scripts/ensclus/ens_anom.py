#*********************************
#            ens_anom            *
#*********************************

# Standard packages
import numpy as np
import sys
import os
from scipy import stats

def ens_anom(filenames,dir_OUTPUT,name_outputs,varname,numens,season,area,extreme):
    '''
    \nGOAL: Computation of the ensemble anomalies based on the desired value from the input variable
    (it can be the percentile, mean, maximum, standard deviation or trend)
    OUTPUT: NetCDF files of ensemble mean of climatology, selected value and anomaly maps.
    '''
    
    # User-defined packages
    from read_netcdf import read3Dncfield, save_N_2Dfields
    from sel_season_area import sel_season, sel_area

    print('***********************************OUTPUT***********************************')
    print('The name of the output files will be <variable>_{0}.txt'.format(name_outputs))
    print('Number of ensemble members: {0}'.format(numens))

    #____________Reading the netCDF file of 3Dfield, for all the ensemble members
    var_ens=[]
    for ens in range(numens):
        ifile=filenames[ens]
        #print('ENSEMBLE MEMBER %s' %ens)
        var, varunits, lat, lon, dates, time_units = read3Dncfield(ifile)
    
        #____________Convertion from kg m-2 s-1 to mm/day
        if varunits=='kg m-2 s-1':
            var=var*86400  #there are 86400 seconds in a day
            varunitsnew='mm/day'
        else:
            varunitsnew=varunits
    
        #____________Selecting a season (DJF,DJFM,NDJFM,JJA)
        var_season,dates_season=sel_season(var,dates,season)
        
        #____________Selecting only [latS-latN, lonW-lonE] box region
        var_area,lat_area,lon_area=sel_area(lat,lon,var_season,area)
        
        var_ens.append(var_area)
    
    if varunitsnew=='mm/day':
        print('\nPrecipitation rate units are converted from kg m-2 s-1 to mm/day')
    
    print('The variable is {0} ({1})'.format(varname,varunitsnew))
    print('Original var shape: (time x lat x lon)={0}'.format(var.shape))
    print('var shape after selecting season {0}: (time x lat x lon)={1}'.format(season,var_season.shape))
    print('var shape after selecting season {0} and area {1}: (time x lat x lon)={2}'.format(season,area,var_area.shape))
    print('Check the number of ensemble members: {0}'.format(len(var_ens)))

    if extreme=='mean':
        #____________Compute the time mean over the entire period, for each ensemble member
        # MEAN
        varextreme_ens=[np.mean(var_ens[i],axis=0) for i in range(numens)]

    elif len(extreme.split("_"))==2:
        #____________Compute the chosen percentile over the period, for each ensemble member
        # PERCENTILE
        q=int(extreme.partition("th")[0])
        varextreme_ens=[np.percentile(var_ens[i],q,axis=0) for i in range(numens)]
    
    elif extreme=='maximum':
        #____________Compute the maximum value over the period, for each ensemble member
        # MAXIMUM
        varextreme_ens=[np.max(var_ens[i],axis=0) for i in range(numens)]
    
    elif extreme=='std':
        #____________Compute the standard deviation over the period, for each ensemble member
        # STANDARD DEVIATION
        varextreme_ens=[np.std(var_ens[i],axis=0) for i in range(numens)]
    
    elif extreme=='trend':
        #____________Compute the linear trend over the period, for each ensemble member
        # TREND
        # Reshape grid to 2D (time, lat*lon)  -> Y
        #Y=[var_ens[i].reshape(var_ens[0].shape[0],var_ens[0].shape[1]*var_ens[0].shape[2])for i in range(numens)]
        #print('Reshaped (time, lat*lon) variable: ',Y[0].shape)
        trendmap=np.empty((var_ens[0].shape[1],var_ens[0].shape[2]))
        trendmap_ens=[]
        for i in range(numens):
            for la in range(var_ens[0].shape[1]):
                for lo in range(var_ens[0].shape[2]):
                    slope, intercept, r_value, p_value, std_err = stats.linregress(range(var_ens[0].shape[0]),var_ens[i][:,la,lo])
                    trendmap[la,lo]=slope
            trendmap_ens.append(trendmap)
        varextreme_ens = trendmap_ens

    print(len(varextreme_ens),varextreme_ens[0].shape)
    varextreme_ens_np=np.array(varextreme_ens)
    print(varextreme_ens_np.shape)
    print('\n------------------------------------------------------------')    
    print('Anomalies are computed with respect to the {0}'.format(extreme))
    print('------------------------------------------------------------\n')    

    #____________Compute and save the anomalies with respect to the ensemble
    ens_anomalies=varextreme_ens_np-np.mean(varextreme_ens_np,axis=0)
    varsave='ens_anomalies'
    ofile=os.path.join(dir_OUTPUT,'ens_anomalies_{0}.nc'.format(name_outputs))
    #print(ofile)
    print('Save the anomalies with respect to the ensemble:')
    print('ens_anomalies shape: (numens x lat x lon)={0}'.format(ens_anomalies.shape))
    save_N_2Dfields(lat_area,lon_area,ens_anomalies,varsave,varunitsnew,ofile)

    #____________Compute and save the climatology
    vartimemean_ens=[np.mean(var_ens[i],axis=0) for i in range(numens)]
    ens_climatologies=np.array(vartimemean_ens)
    varsave='ens_climatologies'
    ofile=os.path.join(dir_OUTPUT,'ens_climatologies_{0}.nc'.format(name_outputs))
    #print(ofile)
    print('Save the climatology:')
    save_N_2Dfields(lat_area,lon_area,ens_climatologies,varsave,varunitsnew,ofile)
    
    #____________Save the extreme
    ens_extreme=varextreme_ens_np
    varsave='ens_extreme'
    ofile=os.path.join(dir_OUTPUT,'ens_extreme_{0}.nc'.format(name_outputs))
    #print(ofile)
    print('Save the extreme:')
    save_N_2Dfields(lat_area,lon_area,ens_extreme,varsave,varunitsnew,ofile)
    

    return

#========================================================

if __name__ == '__main__':
    print('This program is being run by itself')
    
    print('**************************************************************')
    print('Running {0}'.format(sys.argv[0]))
    print('**************************************************************')
    filenames     = sys.argv[1].split()  # input file names
    dir_OUTPUT    = sys.argv[2]          # OUTPUT DIRECTORY
    name_outputs  = sys.argv[3]          # name of the outputs
    varname       = sys.argv[4]          # variable name
    numens        = int(sys.argv[5])     # number of ensemble members
    season        = sys.argv[6]          # seasonal average
    area          = sys.argv[7]          # regional average
    extreme       = sys.argv[8]          # chosen extreme to investigate

    ens_anom(filenames,dir_OUTPUT,name_outputs,varname,numens,season,area,extreme)

else:
    print('ens_anom is being imported from another module')

