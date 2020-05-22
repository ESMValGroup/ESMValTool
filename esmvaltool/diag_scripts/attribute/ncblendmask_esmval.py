import sys, numpy, scipy.stats, math
import netCDF4
# Following two functions for blending and masking modified from Cowtan 2015.
# Calculate blended temperatures using general methods
# Usage:
#  python ncblendmask.py <mode> tas.nc tos.nc sic.nc sftof.nc [Had4.nc] > blend.temp
#  <mode> is one of xxx, mxx, xax, max, xxf, mxf, xaf, maf
#  see README for more details
# Nathan Gillett - Adapted from ncblendmask-nc4.py from Cowtan 2015
# http://www-users.york.ac.uk//~kdc3/papers/robust2015/methods.html

# cell areas, used for calculating area weighted averages
def areas( grid ):
  area = grid*[0.0]
  for i in range(grid):
    area[i] = ( ( math.sin(math.radians(180.0*(i+1)/grid-90.0)) -
                  math.sin(math.radians(180.0*(i  )/grid-90.0)) ) /
                math.sin(math.radians(180.0/grid)) )
  return area


def ncblendmask_esmval(options,sic_file,tas_file,tos_file,sftlf_file,obs_file,dec_warming,obs_dec_warming,ann_warming,gmst_comp_warming,diag_name,obs='had4',ensobs='',ensobs_diag=[],ensobs_dec_warming=[]):
# MAIN PROGRAM

# m = mask
# a = blend anomalies
# f = fix ice
# (use x for none)

  # read tas.nc
  nc = netCDF4.Dataset(tas_file, "r")
  lats1 = nc.variables["lat"][:]
  lons1 = nc.variables["lon"][:]
#  year=nc.variables["year"][:]
#  y0=year[0]#NPG - Added since existing y0 definition below did not work on ESMValTool preprocessed files.
  y0=1850 #Assume 1850 start time.
  tas = numpy.ma.filled(nc.variables["tas"][:,:,:],-1.0e30)
  nc.close()

  # read tos.nc
  nc = netCDF4.Dataset(tos_file, "r")
  lats2 = nc.variables["lat"][:]
  lons2 = nc.variables["lon"][:]
  tos = numpy.ma.filled(nc.variables["tos"][:,:,:],-1.0e30)
#  y0 = int(nc.variables["time"][:][0]/10000)
  nc.close()

  # read sic.nc
  nc = netCDF4.Dataset(sic_file, "r")
  lats3 = nc.variables["lat"][:]
  lons3 = nc.variables["lon"][:]
  #Use siconca if it exists, otherwise use siconc.
  if 'siconca' in nc.variables:
    sic = numpy.ma.filled(nc.variables["siconca"][:,:,:],-1.0e30)
  else:
    sic = numpy.ma.filled(nc.variables["siconc"][:,:,:],-1.0e30)  
    nc.close()

  # read sftlf.nc (NPG - Changed from sftof, because of better data availability for sftlf).
  nc = netCDF4.Dataset(sftlf_file, "r")
  lats4 = nc.variables["lat"][:]
  lons4 = nc.variables["lon"][:]
  sftof = 1-numpy.ma.filled(nc.variables["sftlf"][:,:],-1.0e30) #NPG - added '1-' to use lf.
  nc.close()



  if 'm' in options:
    # read HadCRUT4 data as mask
    nc = netCDF4.Dataset(obs_file, "r")
    lats5 = nc.variables["latitude"][:]
    lons5 = nc.variables["longitude"][:]
    if obs=='had5':
        enssize=200
        obs_tas = nc.variables["tas_median"][:,:,:]
    #Make it work with HadCRUT5 - repeat last year in obs_tas
        obs_tas=numpy.concatenate((obs_tas,obs_tas[2016:2028,:,:]))
    else:
        enssize=100
        obs_tas = nc.variables["temperature_anomaly"][:,:,:]
    #Pad with missing values to match length of tas from model.
    if tas.shape[0]>obs_tas.shape[0]:
      obs_tas = numpy.concatenate((obs_tas,numpy.full((tas.shape[0]-obs_tas.shape[0],obs_tas.shape[1],obs_tas.shape[2]),fill_value=-1e30)))
    elif tas.shape[0]<obs_tas.shape[0]:
    #Or clip obs to match length of simulations if obs are longer than simulations.
      obs_tas=obs_tas[0:tas.shape[0],:,:]
    cvgmsk = numpy.ma.filled(obs_tas,-1.0e30)
    nc.close()
    #Simple regridding to agree with ESMValTool output, HadCRUT4 longitudes start from -177.5.
    regrid_index=list(range(int(lons5.shape[0]*0.5),lons5.shape[0]))+list(range(int(lons5.shape[0]*0.5)))
    lons5=lons5[regrid_index]
    obs_tas=obs_tas[:,:,regrid_index]
    cvgmsk=cvgmsk[:,:,regrid_index]
    #Print missing value fraction
    print('Missing value frction',1.0*numpy.count_nonzero(cvgmsk<-100.)/numpy.size(cvgmsk),file=sys.stderr)

  sic = sic[0:tas.shape[0],:,:]

  # dates
  dates = (numpy.arange(tas.shape[0])+0.5)/12.0 + y0
  print (dates,file=sys.stderr)

  # force missing cells to be open water/land and scale if stored as percentage
  sic[sic<  0.0] = 0.0
  sic[sic>100.0] = 0.0
  if numpy.max(sic)>90.0: sic = 0.01*sic

  sftof[sftof<  0.0] = 0.0
  sftof[sftof>100.0] = 0.0
  if numpy.max(sftof)>90.0: sftof = 0.01*sftof

  print ("tos ", numpy.min(tos), numpy.max(tos), numpy.mean(tos),file=sys.stderr)

  print ("sic ", numpy.min(sic), numpy.max(sic), numpy.mean(sic),file=sys.stderr)
  print ("sftof ", numpy.min(sftof), numpy.max(sftof), numpy.mean(sftof),file=sys.stderr)

  # optional fixed ice mode
  if 'f' in options:
    # mask all cells with any ice post 1961
    for m0 in range(0,len(dates),12):
      if dates[m0] > 1961: break
      print (m0, dates[m0],file=sys.stderr)
    for i in range(sic.shape[1]):
      for j in range(sic.shape[2]):
        for m in range(12):
          cmax = sic[m0+m::12,i,j].max()
          if cmax > 0.01:
            sic[m::12,i,j] = 1.0

  # combine land/ice masks
  for m in range(sic.shape[0]):
    sic[m,:,:] = (1.0-sic[m,:,:])*sftof

  printmask=0
  if printmask==1:
    # print mask
    s = ""
    sicmax = numpy.max(sic)
    for i in range(sic.shape[1]-1,0,-sic.shape[1]//25):
      for j in range(0,sic.shape[2],sic.shape[2]//50):
        s += ".123456789#"[int(10*sic[-1,i,j]/sicmax)]
      s += "\n"
    print (s, "\n",file=sys.stderr)
    # print tos mask
    s = ""
    for i in range(tos.shape[1]-1,0,-tos.shape[1]//25):
      for j in range(0,tos.shape[2],tos.shape[2]//50):
        s += "#" if -500 < tos[-1,i,j] < 500 else "."
      s += "\n"
    print (s, "\n",file=sys.stderr)
    # print cvg mask
    if 'm' in options:
      s = ""
      for i in range(cvgmsk.shape[1]-1,0,-cvgmsk.shape[1]//25):
        for j in range(0,cvgmsk.shape[2],cvgmsk.shape[2]//50):
#          s += "#" if -100 < cvgmsk[-1,i,j] < 500 else "."
          s += "#" if -100 < cvgmsk[400,i,j] < 500 else "."
        s += "\n"
      print (s, "\n",file=sys.stderr)

  # deal with missing tos through sic
  for m in range(sic.shape[0]):
    sic[m,tos[m,:,:]<-500.0] = 0.0
    sic[m,tos[m,:,:]> 500.0] = 0.0

  # baseline and blend in the desired order
  if 'a' in options:

    # prepare missing
    for m in range(sic.shape[0]):
#      tos[m,tos[m,:,:]<-500.0] = numpy.nan
      tos[m,abs(tos[m,:,:])> 500.0] = numpy.nan 

    # baseline
    mask = numpy.logical_and( dates > 1961, dates < 1991 )
    base = tas[mask,:,:]
    for m in range(12):
      norm = numpy.mean(base[m::12,:,:],axis=0)
      tas[m::12,:,:] = tas[m::12,:,:] - norm
    base = tos[mask,:,:]
    for m in range(12):
      norm = numpy.nanmean(base[m::12,:,:],axis=0)
      tos[m::12,:,:] = tos[m::12,:,:] - norm
    # blend
    for m in range(sic.shape[0]):
      tos[m,:,:] = tas[m,:,:]*(1.0-sic[m,:,:])+tos[m,:,:]*(sic[m,:,:])

  else:

    # blend
    for m in range(sic.shape[0]):
      tos[m,:,:] = tas[m,:,:]*(1.0-sic[m,:,:])+tos[m,:,:]*(sic[m,:,:])
    # baseline
    mask = numpy.logical_and( dates > 1961, dates < 1991 )
    base = tas[mask,:,:]
    for m in range(12):
      norm = numpy.mean(base[m::12,:,:],axis=0)
      tas[m::12,:,:] = tas[m::12,:,:] - norm
    base = tos[mask,:,:]
    for m in range(12):
      norm = numpy.mean(base[m::12,:,:],axis=0)
      tos[m::12,:,:] = tos[m::12,:,:] - norm


  # deal with any remaining nans
  for m in range(sic.shape[0]):
    msk = numpy.isnan(tos[m,:,:])
    tos[m,msk] = tas[m,msk]
  # calculate area weights
  w = numpy.zeros_like(tas)
  wm = numpy.zeros_like(tas)
  a = areas(sftof.shape[0])
  for m in range(w.shape[0]):
#    for i in range(w.shape[1]):
      for j in range(w.shape[2]):
        w[m,:,j] = a[:]

  wm=w.copy()
  if 'm' in options: wm[ cvgmsk[0:wm.shape[0],:,:] < -100 ] = 0.0
  print (w[0,:,:],file=sys.stderr)
  print (wm[0,:,:],file=sys.stderr)
  # calculate diagnostic
  diag=calc_diag(tos,wm,diag_name) #Diagnostic for attribution analysis.
  dec_warming.append(calc_dec_warming(tas,w)) #Diagnose SAT warming with global coverage for attributable trends.
  obs_dec_warming.append(calc_dec_warming(obs_tas,wm))
  
  if ann_warming!=0:
    ann_warming.append(calc_ann_warming(tas,w)) #Calculate ann warming.
  if gmst_comp_warming!=0:
    gmst_comp_warming.append(calc_ann_warming(tos,w)) #Calculate warming in globally-complete blended data.
  obs_diag=calc_diag(obs_tas[0:tos.shape[0],:,:],wm,diag_name)

  #Repeat obs diagnostics for each member of ensemble observational dataset if ensobs is set. Assume missing data mask is the same as for main obs dataset.
  if ensobs != '':
    #Assume 100 member ensemble observations dataset.
    for ens in range(1,enssize+1):
      nc = netCDF4.Dataset(ensobs+str(ens)+'.nc', "r")
      if obs=='had5':
        obs_tas = nc.variables["tas"][:,:,:]
    #Make it work with HadCRUT5 - repeat last year in obs_tas
        obs_tas=numpy.concatenate((obs_tas,obs_tas[2016:2028,:,:]))
      else:
        obs_tas = nc.variables["temperature_anomaly"][:,:,:]
      nc.close()
      obs_tas=obs_tas[:,:,regrid_index] #Regrid to match esmvaltool output.
      ensobs_dec_warming.append(calc_dec_warming(obs_tas[0:tas.shape[0],:,:],wm))
      ensobs_diag.append(calc_diag(obs_tas[0:tas.shape[0],:,:],wm,diag_name))
  return (diag,obs_diag)

def calc_diag(tos,wm,diag_name):
  av_per=int(diag_name[4:6])*12 #Last two digits of diag_name are averaging period in yrs.
  #compute diagnostic based on masked/blended temperatures.
  if diag_name[0:4]=='gmst':
    nlat=1
  elif diag_name[0:4]=='hemi':
    nlat=2
  else:
    print ('Diagnostic ',diag_name,' not supported')
    exit ()
  nper=math.ceil(tos.shape[0]/av_per) #Round up number of averaging periods.
  diag=numpy.zeros((nlat,nper))
  # calculate temperatures
  for m in range(nper):
    for l in range(nlat):
      diag[l,m]=numpy.sum( wm[m*av_per:(m+1)*av_per,l*tos.shape[1]//nlat:(l+1)*tos.shape[1]//nlat,:] * tos[m*av_per:(m+1)*av_per,l*tos.shape[1]//nlat:(l+1)*tos.shape[1]//nlat,:] ) / numpy.sum( wm[m*av_per:(m+1)*av_per,l*tos.shape[1]//nlat:(l+1)*tos.shape[1]//nlat,:] )
  diag=diag-numpy.mean(diag,axis=1,keepdims=True) #Take anomalies over whole period.
  diag=numpy.reshape(diag,nper*nlat)
  return diag

def calc_dec_warming(tas,w):
  gmt_mon=numpy.zeros(tas.shape[0])
  # calculate 2010-2019 mean relative to 1850-1900, assuming data starts in 1850.
  # If last decade is incomplete, just compute mean from available data.
  for m in range(tas.shape[0]):
      s = numpy.sum( w[m,:,:] )
      gmt_mon[m] = numpy.sum( w[m,:,:] * tas[m,:,:] ) / s
  return (numpy.nanmean(gmt_mon[(2010-1850)*12:(2020-1850)*12])-numpy.mean(gmt_mon[0:(1901-1850)*12]))

def calc_ann_warming(tas,w):
  nyr=math.ceil(tas.shape[0]/12) #Round up number of years. 
  diag=numpy.zeros(nyr)
  gsat_mon=numpy.zeros(tas.shape[0])
  # calculate temperatures
  for m in range(tas.shape[0]):
    s = numpy.sum( w[m,:,:] )
    gsat_mon[m] = numpy.sum( w[m,:,:] * tas[m,:,:] ) / s
#    print (gmst_mon)
  for m in range(nyr):
    diag[m]=numpy.mean(gsat_mon[m*12:(m+1)*12]) #Note - will calculate average over incomplete final year.
  diag=diag-numpy.mean(diag[0:(1901-1850)]) #Take anomalies relative to 1850-1901.
  return (diag)
