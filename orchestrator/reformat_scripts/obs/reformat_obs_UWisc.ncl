;;#############################################################################
;; REFORMAT SCRIPT FOR UWisc OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier: 3 (restricted data set)
;;
;; Source
;;    Reference: O'Dell, C. W., F. J. Wentz, and R. Bennartz (2008),
;;               Cloud liquid water path from satellite-based passive microwave
;;               observations: A new climatology over the global oceans,
;;               J. Clim., 21, 1721-1739, doi:10.1175/2007JCLI1958.1.
;;
;; Last access
;;    2008
;;
;; Download and processing instructions
;;    Download: contact Ralf Bennartz, University of Wisconsin-Madison,
;;              Wisconsin, USA (email: Bennartz 'at' aos.wisc.edu)
;;    Processing: this script (reformat_obs_UWisc.ncl)
;;
;; Caveats
;;    None
;;
;; Modification history
;;    20150415-A_laue_ax: written.
;;
;; ############################################################################

begin

    ;; Source name
    OBSNAME = "UWisc"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T2Ms"

    ;; Period
    YEAR1 = 1988
    YEAR2 = 2007

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_UWisc.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over years
    do yr = YEAR1, YEAR2
        fname = INDIR + "UWisc_LWPMEAN_" + yr + "_v2.nc"
        info_output("*** "+ fname + " ***", verbosity, 1)

        f = addfile(fname, "r")
        lwp = (/f->LWP/)
        lwperr = (/f->LWP_ERROR/)
        lat = (/f->lat/)
        lon = (/f->lon/)
        month = (/f->month/)

;; does not work --> set attributes of dimension variables manually
;;       time = time_attrib(time)
;;       lat = lat_attrib(lat)
;;       lon = lon_attrib(lon)

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

        n = dimsizes(month)
        time = new((/n/), double)

        time@standard_name = "time"
        time@units = "days since 1950-01-01 00:00:00"
        time@axis = "T"
        time@long_name = "time"
        time@calendar = "standard"
        if (isatt(time, "_FillValue")) then
            delete(time@_FillValue)
        end if

        do i = 0, n-1
            time(i) = cd_inv_calendar(yr, tointeger(month(i)), 1, 0, 0, 0, \
                "days since 1950-01-01 00:00:00", 0)
        end do

        lwp!0 = "time"
        lwp!1 = "lat"
        lwp!2 = "lon"
        lwp&time = time
        lwp&lat = lat
        lwp&lon = lon

        lwp@units = "kg m-2"
        lwp@long_name = "liquid water path"

        lwperr!0 = "time"
        lwperr!1 = "lat"
        lwperr!2 = "lon"
        lwperr&time = time
        lwperr&lat = lat
        lwperr&lon = lon

        lwperr@units = "kg m-2"
        lwperr@long_name = "liquid water path error"

        ; convert units from [g/m2] to [kg/m2] (accounts for missing values)

        lwp = lwp / 1000.0
        lwperr = lwperr / 1000.0

        gAtt = True
        gAtt@history     = "Created on " + systemfunc("date")
        gAtt@host        = systemfunc("echo $HOST")
        gAtt@user        = systemfunc("echo $USER")
        gAtt@period      = yr
        gAtt@tier        = TIER
        gAtt@field       = FIELD
        gAtt@source      = "Data provided by Ralf Bennartz (Vanderbilt)"
        gAtt@reference   = "O'Dell, C. W. et al., J. Clim., " + \
            "21, 1721-1739, doi:10.1175/2007JCLI1958.1, 2008."
        gAtt@title       = OBSNAME + \
            " satellite data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile (lwp)
        fout = OUTDIR + "OBS_" + OBSNAME + "_sat_v2_" + FIELD + "_clwvi_" + \
            yr + "01-" + yr + "12.nc"
       ;; Write variable
       write_nc(fout, "clwvi", lwp, gAtt)
       f = addfile(fout, "w")
;       f->clwvi_err = lwperr
       delete(f)

       ;; Outfile (lwperr)
       fout = OUTDIR + "OBS_" + OBSNAME + "_sat_v2_" + FIELD + "_lwpStderr_" + \
           yr + "01-" + yr + "12.nc"

       ;; Write variable
       write_nc(fout, "lwpStderr", lwperr, gAtt)
       f = addfile(fout, "w")
       delete(f)

       gAtt = True
       gAtt@tier        = 2
       gAtt@title       = "dummy field"
       gAtt@reference   = ""

       ;; Outfile (dummy "clivi")
       fout = OUTDIR + "OBS_" + OBSNAME + "_sat_v2_" + FIELD + "_clivi_" + \
              yr + "01-" + yr + "12.nc"
       ;; Write variable
       lwp = -999.0
       lwp@long_name = "dummy field"
       write_nc(fout, "clivi", lwp, gAtt)

       delete(gAtt)
    end do

    delete(time)

    leave_msg(diag_script, "", 1)

end

