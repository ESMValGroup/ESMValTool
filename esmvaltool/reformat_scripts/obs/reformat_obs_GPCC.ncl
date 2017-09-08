;;#############################################################################
;; REFORMAT SCRIPT FOR GPCC REANALYSIS DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    ftp://ftp-anon.dwd.de/pub/data/gpcc/html/fulldata_v6_doi_download.html
;;
;; Last access
;;    20150401
;;
;; Download and processing instructions
;;    Select NetCDF format and the 1901 - 2010 precipitation per gridcell.
;;
;; Caveats
;;
;; Modification history
;;    20150401-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "GPCC"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Selected variable (standard name)
    VAR = "pr"

    ;; Name in the raw data
    NAME = "p"

    ;; Period
    YEAR1 = 1901
    YEAR2 = 2010

    ;; Initialize global variable
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_GPCC.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    f = addfile(INDIR + "full_data_v6_precip_25.nc", "r")
    output = f->$NAME$

    ;; Calculate days per month
    date = cd_calendar(output&time, 0)
    dpm = days_in_month(toint(date(:, 0)), toint(date(:, 1)))
    dpmc = conform(output, dpm, 0)

    ;; Check time range
    if (dimsizes(date(:, 0)).ne.12 * (YEAR2 - YEAR1 + 1)) then
        error_msg("f", diag_script, "", "incorrect number of timesteps")
    end if

    ;; Convert units [mm/month] --> [kg/m2/s]
    output = output / (24 * 3600 * dpmc)

    ;; Set variable attributes
    tmp = var_attrib(output, VAR)
    delete(output)
    output = tmp
    delete(tmp)

    ;; Format time coordinate
    ctime = time_attrib(output&time)
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ;; Format latitude coordinate
    output!1 = "lat"
    if (isMonotonic(output&lat) .eq. 0) then
        error_msg("f", diag_script, "", \
                  "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(output&lat) .eq. -1) then
        output = output(:, ::-1, :)
    end if
    clat = lat_attrib(output&lat)
    delete(output&lat)
    output&lat = clat
    delete(clat)

    ;; Format longitude coordinate
    output!2 = "lon"
    if (isMonotonic(output&lon) .eq. 0) then
        error_msg("f", diag_script, "", \
                  "non-monotonic longitude coordinate")
    end if
    if (any(output&lon.lt.0.)) then
        output = lonFlip(output)
    end if
    clon = lon_attrib(output&lon)
    delete(output&lon)
    output&lon = clon
    delete(clon)

    info_output("  Range: " + min(output) + "; " + max(output), verbosity, 1)

    ;; Set global attributes
    gAtt = True
    gAtt@history       = "Created on " + systemfunc("date")
    gAtt@host          = systemfunc("echo $HOST")
    gAtt@user          = systemfunc("echo $USER")
    gAtt@period        = YEAR1 + "-" + YEAR2
    gAtt@field         = FIELD
    gAtt@tier          = TIER
    gAtt@source        = "ftp://ftp-anon.dwd.de/pub/data/gpcc/html/" + \
        "fulldata_v6_doi_download.html"
    gAtt@reference = "Schneider, U. et al., GPCC Full Data Reanalysis " + \
        "Version 6.0 at 2.5°: Monthly Land-Surface Precipitation from " + \
        "Rain-Gauges built on GTS-based and Historic Data, " + \
        "doi: 10.5676/DWD_GPCC/FD_M_V6_250, 2011."
    gAtt@title = OBSNAME + \
        " reanalysis data reformatted for the ESMValTool"
    gAtt@conventions = "CF/CMOR"

    ;; Outfile
    fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_v6_" + FIELD + \
        "_" + VAR + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ;; Write variable
    write_nc(fout, VAR, output, gAtt)
    delete(gAtt)
    delete(output)

    leave_msg(diag_script, "", 1)

end
