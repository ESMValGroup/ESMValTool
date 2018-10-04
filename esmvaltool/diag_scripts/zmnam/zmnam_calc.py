import numpy as np
import netCDF4 as nc4
import sys
from scipy import signal

def butter_filter(data, fs, lowcut=None, order= 2):
    
    if lowcut != None  : filttype = 'lowpass'   

    # Sampling determines Nyquist frequency
    nyq = 0.5 * fs

    if filttype == 'lowpass' : low  = lowcut  / nyq ; freqs = low
    
    b, a = signal.butter(order, freqs, btype=filttype)
    w, h = signal.freqs(b, a)
    y = signal.filtfilt(b, a, data,axis=0)

    return y    

def zmnam_calc(indir,outdir,src_props):

    # TODO use src_props to name files
    #print('calc: ',src_props)

    deg_to_r = np.pi/180. 
    lat_weighting = True

    # Note: daily/monthly means have been
    # already subtracted from daily/monthly files

    #=================
    # Open daily data 
    #=================
    da_fname = 'tmp_gh_da_an_zm_hem.nc'
    in_file = nc4.Dataset(indir+da_fname,"r")
    #
    time_dim  = in_file.variables['time'][:]
    time_nam  = in_file.variables['time'].long_name
    time_uni  = in_file.variables['time'].units
    time_cal  = in_file.variables['time'].calendar
    time      = np.array(time_dim[:],dtype='d')
    startdate = nc4.num2date(time[0],time_uni,time_cal)
    #print startdate
    date = nc4.num2date(time,in_file.variables['time'].units,\
    in_file.variables['time'].calendar)
    #
    lev      = np.array(in_file.variables['plev'][:], dtype='d') 
    lev_nam  = in_file.variables['plev'].long_name
    lev_uni  = in_file.variables['plev'].units
    lev_pos  = in_file.variables['plev'].positive
    lev_axi  = in_file.variables['plev'].axis
    #print lev
    #
    lat      = np.array(in_file.variables['lat'][:], dtype='d')
    lat_nam  = in_file.variables['lat'].long_name
    lat_uni  = in_file.variables['lat'].units
    lat_axi  = in_file.variables['lat'].axis
    #
    #print lat
    lon      = np.array(in_file.variables['lon'][:], dtype='d') 
    lon_nam  = in_file.variables['lon'].long_name
    lon_uni  = in_file.variables['lon'].units
    lon_axi  = in_file.variables['lon'].axis
    #print lon
    #
    print('squeezing')
    zg_da    = np.squeeze(np.array(in_file.variables['zg'][:], dtype='d'))
    #zg_da = np.array(in_file.variables['zg'][:], dtype='d')
    #    zg_da = np.array(in_file.variables['zg'][:])
    print('closing file')
    #    zg_da = np.squeeze(zg_da)
    in_file.close()

    n_tim = len(time_dim)
    print('end infile close')



    #==================
    # Index calculation
    #==================

    # Lowpass filter
    zg_da_lp = butter_filter(zg_da,1,lowcut=1./90,order=2)
    print('end filtering')

    # Outputs: stored by level
    # EOFs, eigenvalues, daily and monthly PCs
    eofs   = np.zeros((len(lev),len(lat)),dtype='d')
    eigs   = np.zeros(len(lev),dtype='d')
    pcs_da = np.zeros((n_tim,len(lev)),dtype='d')

    # Calendar-independent monthly mean
    sta_mon = [] # first day of the month
    mid_mon = [] # 15th of the month
    end_mon = [] # last day of the month; add +1 when slicing

    mo = 999
    idate = 0

    while idate < len(date):
        #

        # Save first day of the month
        if date[idate].month != mo : mo = date[idate].month ; sta_mon.append(idate)

        # Save month mid-day
        if date[idate].day == 15 : mid_mon.append(idate)
 
        # Save last day of the month
        if ((idate == len(date)-1) or \
        (date[idate].month == mo and date[idate+1].month != mo)):
            #
            end_mon.append(idate)
        #
        idate += 1
    #


    pcs_mo = np.zeros((len(date[mid_mon]),len(lev)),dtype='d')

    # Perform analysis by level
    for i_lev in np.arange(len(lev)):
        print(i_lev)
        # Latitude weighting
        if lat_weighting == True:
            for j_lat in np.arange(len(lat)):
                zg_da_lp[:,i_lev,j_lat] *= np.sqrt(abs(np.cos(lat[j_lat]*deg_to_r)))
           
        zg_da_lp_an = zg_da_lp[:,i_lev,:]-np.mean(zg_da_lp[:,i_lev,:],axis=0)
        cov = np.dot(zg_da_lp_an.T,zg_da_lp_an)/(n_tim-1)

        # Print covariance matrix shape and check symmetry
        print(np.shape(cov))
        print("Symmetric?", (cov.transpose() == cov).all())

        # Compute eigenvectors and eigenvalues
        eigenval, eigenvec = np.linalg.eig(cov)

        sum_eigenval = np.sum(eigenval)

        eigenval_norm = eigenval[:]/sum_eigenval

        # Largest eigenvalue
        max_eigenval = max(range(len(eigenval_norm)), key = lambda x: eigenval_norm[x])

        print("Maximum eigenvalue", eigenval_norm[max_eigenval])

        print("Eigen-analysis, vec - var")
        print(np.shape(eigenvec), "-", np.shape(eigenval))

        # PC calculation
        pc = np.dot(zg_da_lp_an[:,:],eigenvec)

        print("PCs"); print(np.shape(pc))


        # Latitude de-weighting
        if lat_weighting == True:
            for i_lat in np.arange(len(lat)):
                eigenvec[i_lat,:] /= np.sqrt(abs(np.cos(lat[i_lat]*deg_to_r)))

        # Retain leading standardized PC & EOF
        # Avoid standardization here? 
        lead_pc_mean = np.mean(pc[:,max_eigenval],axis=0)
        lead_pc_std  = np.std(pc[:,max_eigenval],ddof=1,axis=0)
        lead_pc  = (pc[:,max_eigenval]-lead_pc_mean)/lead_pc_std
        lead_eof = eigenvec[:,max_eigenval]

        max_lat = max(range(len(lat)), key = lambda x: lat[x])
        min_lat = min(range(len(lat)), key = lambda x: lat[x])

        if lead_eof[max_lat] > lead_eof[min_lat]:
            lead_pc  *= -1
            lead_eof *= -1

        lead_pc_mo = np.zeros(len(date[mid_mon]),dtype='d')
        time_mo = np.zeros(len(date[mid_mon]),dtype='d')

        for k_mo in range(len(date[mid_mon])):
            lead_pc_mo[k_mo] = np.mean(lead_pc[sta_mon[k_mo]:end_mon[k_mo]+1])
            time_mo[k_mo] = time[mid_mon[k_mo]]

        # Store PC/EOF for this level (no time dependent)
        eigs[i_lev] = eigenval_norm[max_eigenval]
        eofs[i_lev,:] = lead_eof
        pcs_da[:,i_lev] = lead_pc
        pcs_mo[:,i_lev] = lead_pc_mo  


    # Save output files

    # (1) daily PCs
    #file_out = nc4.Dataset(outdir+'daily_pc.nc', mode='w',format = 'NETCDF3_CLASSIC')
    file_out = nc4.Dataset(outdir+'_'.join(src_props)+'_pc_da.nc', \
    mode='w',format = 'NETCDF3_CLASSIC')
    file_out.title = 'Zonal mean annular mode (1)'
    file_out.contact = 'F. Serva (federico.serva@artov.isac.cnr.it); \
    C. Cagnazzo (c.cagnazzo@isac.cnr.it)'
    #
    file_out.createDimension('time', None) #np.size(time_dim))
    file_out.createDimension('plev',  np.size(lev))
    file_out.createDimension('lat',  np.size(lat))
    file_out.createDimension('lon',  np.size(lon))
    #
    time_var = file_out.createVariable('time','d',('time',))
    #time_var.setncattr('dimensions',time_dim) #dimensions = time_dim
    time_var.setncattr('long_name',time_nam) #long_name  = time_nam
    time_var.setncattr('units',time_uni) #units      = time_uni
    time_var.setncattr('calendar',time_cal) #calendar   = time_cal
    time_var[:] = time_dim[:]
    #
    lev_var = file_out.createVariable('plev','d',('plev',))
    lev_var.setncattr('long_name',lev_nam) #long_name = time_nam
    lev_var.setncattr('units',lev_uni) #units     = time_uni
    lev_var.setncattr('positive',lev_pos) #positive  = time_pos
    lev_var.setncattr('axis',lev_axi) #axis      = time_axi
    lev_var[:] = lev[:]
    #
    pcs_da_var = file_out.createVariable('PC_da','d',('time','plev',))
    pcs_da_var.setncattr('long_name','Daily zonal mean annular mode PC') #long_name = 'Daily zonal mean NAM PC'
    pcs_da_var.setncattr('comment','Reference: Baldwin and Thompson (2009), doi:10.1002/qj.479')
    pcs_da_var[:] = pcs_da[:,:]

    file_out.close()

    # (2) monthly PCs
    #file_out = nc4.Dataset(outdir+'monthly_pc.nc', mode='w',format = 'NETCDF3_CLASSIC')
    file_out = nc4.Dataset(outdir+'_'.join(src_props)+'_pc_mo.nc', \
    mode='w',format = 'NETCDF3_CLASSIC') 
    file_out.title = 'Zonal mean annular mode (2)'
    file_out.contact = 'F. Serva (federico.serva@artov.isac.cnr.it); \
    C. Cagnazzo (c.cagnazzo@isac.cnr.it)'
    #
    file_out.createDimension('time', None) #np.size(time_dim))
    file_out.createDimension('plev',  np.size(lev))
    #
    time_var = file_out.createVariable('time','d',('time',))
    #time_var.setncattr('dimensions',time_dim) #dimensions = time_dim
    time_var.setncattr('long_name',time_nam) #long_name  = time_nam
    #time_var.setncattr('units','days since '+str(startdate)) #units      = time_uni
    #time_var.setncattr('calendar','365_day') #calendar   = time_cal
    time_var.setncattr('units',time_uni)
    time_var.setncattr('calendar',time_cal)
    time_var[:] = time_mo #np.arange(12*n_tim/365,dtype='i')
    #
    lev_var = file_out.createVariable('plev','d',('plev',))
    lev_var.setncattr('long_name',lev_nam) #long_name = time_nam
    lev_var.setncattr('units',lev_uni) #units     = time_uni
    lev_var.setncattr('positive',lev_pos) #positive  = time_pos
    lev_var.setncattr('axis',lev_axi) #axis      = time_axi
    lev_var[:] = lev[:]
    #
    pcs_mo_var = file_out.createVariable('PC_mo','d',('time','plev',))
    pcs_mo_var.setncattr('long_name','Monthly zonal mean annular mode PC') #long_name = 'Daily zonal mean NAM PC'
    pcs_mo_var.setncattr('comment','Reference: Baldwin and Thompson (2009), doi:10.1002/qj.479')
    pcs_mo_var[:] = pcs_mo[:,:]

    file_out.close()

    # (3) EOFs and explained variance
    #file_out = nc4.Dataset(outdir+'eofs.nc', mode='w',format = 'NETCDF3_CLASSIC')
    file_out = nc4.Dataset(outdir+'_'.join(src_props)+'_eofs.nc', \
    mode='w',format = 'NETCDF3_CLASSIC')
 
    file_out.title = 'Zonal mean annular mode (3)'
    file_out.contact = 'F. Serva (federico.serva@artov.isac.cnr.it); \
    C. Cagnazzo (c.cagnazzo@isac.cnr.it)'
    #
    file_out.createDimension('time', None) #np.size(time_dim))
    file_out.createDimension('plev',  np.size(lev))
    file_out.createDimension('lat',  np.size(lat))
    file_out.createDimension('lon',  np.size(lon))
    #
    time_var = file_out.createVariable('time','d',('time',))
    #time_var.setncattr('dimensions',time_dim) #dimensions = time_dim
    time_var.setncattr('long_name',time_nam) #long_name  = time_nam
    time_var.setncattr('units',time_uni) #units      = time_uni
    time_var.setncattr('calendar',time_cal) #calendar   = time_cal
    time_var[:] = 0 #time_dim[:]
    #
    lev_var = file_out.createVariable('plev','d',('plev',))
    lev_var.setncattr('long_name',lev_nam) #long_name = time_nam
    lev_var.setncattr('units',lev_uni) #units     = time_uni
    lev_var.setncattr('positive',lev_pos) #positive  = time_pos
    lev_var.setncattr('axis',lev_axi) #axis      = time_axi
    lev_var[:] = lev[:]
    #
    lat_var = file_out.createVariable('lat','d',('lat',))
    lat_var.setncattr('long_name',lat_nam) #long_name = time_nam
    lat_var.setncattr('units',lat_uni) #units     = time_uni
    lev_var.setncattr('axis',lat_axi) #axis      = time_axi
    lat_var[:] = lat[:]
    #
    lon_var = file_out.createVariable('lon','d',('lon',))
    lon_var.setncattr('long_name',lon_nam) #long_name = time_nam
    lon_var.setncattr('units',lon_uni) #units     = time_uni
    lon_var.setncattr('axis',lon_axi) #axis      = time_axi
    lon_var[:] = lon[:]
    #
    eofs_var = file_out.createVariable('EOF','d',('plev','lat'))
    eofs_var.setncattr('long_name','Zonal mean annular mode EOF') 
    eofs_var.setncattr('comment','Reference: Baldwin and Thompson (2009), doi:10.1002/qj.479')
    eofs_var[:] = eofs[:,:]
    #
    eigs_var = file_out.createVariable('eigenvalues','d',('plev'))
    eigs_var.setncattr('long_name','Zonal mean annular mode EOF explained variance') 
    eigs_var.setncattr('comment','Reference: Baldwin and Thompson (2009), doi:10.1002/qj.479')
    eigs_var[:] = eigs[:]
    #
    file_out.close()

    return

