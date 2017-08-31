;;#############################################################################
;; REFORMAT SCRIPT FOR THE OBSERVATION "MODIS-CFMIP" (satellite IWP
;; climatology, originally created for CFMIP)
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    References:
;;       King, M. D., et al. (2003), Cloud and aerosol properties,
;;          precipitable water, and profiles of temperature and water vapor
;;          from MODIS, IEEE Trans. Geosci. Remote Sens., 41, 442-458,
;;          doi:10.1109/TGRS.2002.808226.
;;       Pincus, R., S. Platnick, S. A. Ackerman, R. S. Hemler,
;;          Robert J. Patrick Hofmann (2012), Reconciling simulated and
;;          observed views of clouds: MODIS, ISCCP, and the limits of
;;          instrument simulators. J. Climate, 25, 4699-4720, 
;;          doi: 10.1175/JCLI-D-11-00267.1.
;;
;; Last access
;;    20160502
;;
;; Download and processing instructions
;;    Download: ftp://laadsweb.nascom.nasa.gov/NetCDF/L3_Monthly/V02/
;;    Processing: this script (reformat_obs_MODIS.ncl)
;;
;; Caveats
;;    These MODIS IWP data are in-cloud values whereas the clivi data from the
;;    CMIP5 models are grid-box averages --> e.g. multiplication with cloud
;;    fraction needed!
;;
;; Modification history
;;    20150417-A_laue_ax: written.
;;
;; ############################################################################

begin

    ;; Source name
    OBSNAME = "MODIS-CFMIP"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T2Ms"

    ;; Period
    YEAR1 = 2003
    YEAR2 = 2014

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_MODIS-CFMIP.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    initialize = 1

    ;; Loop over years
    do yr = YEAR1, YEAR2
       do mo = 1, 12
          month = sprinti("%0.2i", mo)
          fname = INDIR + "MCD08_M3_NC." + yr + "." + month + ".C051.V02.nc"
          info_output("*** "+ fname + " ***", verbosity, 1)

          f = addfile(fname, "r")

          if (initialize.eq.1) then ; do this only once...
             lat = (/f->lat/)
             lon = (/f->lon/)

             nx = dimsizes(lon)
             ny = dimsizes(lat)
             nt = 12

             time = new((/nt/), double)
             iwp = new((/nt,ny,nx/), "float", "No_FillValue")

             lat@standard_name = "latitude"
             lat@units = "degrees_north"
             lat@axis = "Y"
             lat@long_name = "latitude"
             lat!0 = "lat"

             lon@standard_name = "longitude"
             lon@units = "degrees_east"
             lon@axis = "X"
             lon@long_name = "longitude"
             lon!0 = "lon"

             time@standard_name = "time"
             time@units = "days since 1950-01-01 00:00:00"
             time@axis = "T"
             time@long_name = "time"
             time@calendar = "standard"
             if (isatt(time, "_FillValue")) then
                delete(time@_FillValue)
             end if

;             ; convert longitudes from -180...180 to 0...360
;
;             do i = 0, nx - 1
;                if (lon(i).lt.0.0) then
;                   lon(i) = lon(i) + 360.0
;                end if
;             end do

             iwp!0 = "time"
             iwp!1 = "lat"
             iwp!2 = "lon"
             iwp&lat = lat
             iwp&lon = lon
             iwp@units = "kg m-2"
             iwp@long_name = "ice water path"

             initialize = 0
          end if

          iwp(mo-1,:,:) = (/f->Ice_Path_Mean/)
          ; convert units from [g/m2] to [kg/m2] (accounts for missing values)
          iwp(mo-1,:,:) = iwp(mo-1,:,:) / 1000.0
          time(mo-1) = cd_inv_calendar(yr, mo, 1, 0, 0, 0, \
                       "days since 1950-01-01 00:00:00", 0)

       end do ; month loop

       iwp&time = time

       gAtt = True
       gAtt@history     = "Created on " + systemfunc("date")
       gAtt@host        = systemfunc("echo $HOST")
       gAtt@user        = systemfunc("echo $USER")
       gAtt@period      = yr
       gAtt@tier        = TIER
       gAtt@field       = FIELD
       gAtt@reference1  = "Pincus et al. (2012), J. Clim."
       gAtt@reference2  = "King et al. (2003), IEEE Trans. Geosci. " \
                          + "Remote Sens."
       gAtt@title       = "CFMIP MODIS ice water path (IWP) climatology"
       gAtt@conventions = "CF-1.1/CMOR"

       ;; Outfile
       fout = OUTDIR + "OBS_" + OBSNAME + "_sat_v2_" + FIELD + "_clivi_" + \
           yr + "01-" + yr + "12.nc"
       ;; Write variable
       write_nc(fout, "clivi", iwp, gAtt)
       delete(gAtt)
    end do ; year loop

    delete(time)
    delete(iwp)

    leave_msg(diag_script, "", 1)

end

