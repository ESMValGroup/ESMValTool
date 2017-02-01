;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-AEROSOL SATELLITE DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.icare.univ-lille1.fr/drupal/archive/?dir=CCI-Aerosols/
;;
;; Last access
;;    20151124
;;
;; Download and processing instructions
;;    Login: cci Password: cci
;;    Download the data from ATSR2_SU_v4.2/L3_MONTHLY/ (1998-2002)
;;    and from AATSR_SU_v4.2/L3_MONTHLY/ (2003-2011).
;;    Other years are not considered since they are not complete.
;;
;; Caveats
;;
;; Modification history
;;    20151124-A_righ_ma: switched to monthly raw data (now available).
;;    20150424-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ESACCI-AEROSOL"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = "/export/pa_data01/ESMVal/obs/RAW/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = "/export/pa_data02/ESMVal/obs/Tier" + TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1998
    YEAR2 = 2011

    ;; Selected variable (standard name)
    VARS = "od550aer"

    ;; Name in the raw data
    NAME = "AOD550_MEAN"

    ;; Field
    FIELD = "T2Ms"

end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = 1
    diag_script = "reformat_obs_ESACCI-AEROSOL.ncl"

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    info_output("", verbosity, 1)
    info_output("Processing " + VARS + " (" + FIELD + ") ", verbosity, 1)

    time = create_timec(YEAR1, YEAR2)
    date = cd_calendar(time, 1)

    ;; Create timeseries
    do yy = YEAR1, YEAR2
        do mm = 1, 12

            ldate = yy + sprinti("%0.2i", mm)
            info_output("Processing " + ldate, verbosity, 1)

            ;; Read file
            str = "find " + INDIR + " -type f -name '" + ldate + "*MONTH*.nc'"
            fname = systemfunc(str)

            ;; No files found
            if (all(ismissing(fname))) then
                continue
            end if

            ;; Extract data
            f = addfile(fname, "r")
            xx = f->$NAME$

            ;; Assign to global array
            if (.not.isdefined("OUTPUT")) then
                dims = array_append_record(dimsizes(time), dimsizes(xx), 0)
                OUTPUT = new(dims, float)
                OUTPUT!0 = "time"
                OUTPUT&time = time
                OUTPUT!1 = "lat"
                OUTPUT&lat = f->latitude
                OUTPUT!2 = "lon"
                OUTPUT&lon = f->longitude
            end if
            OUTPUT(ind(toint(ldate).eq.date), :, :) = (/xx/)

            delete(fname)

        end do
    end do

    ;; Set fill value
    OUTPUT = where(OUTPUT.eq.-999, OUTPUT@_FillValue, OUTPUT)

    ;; Set variable attributes
    tmp = var_attrib(OUTPUT, VARS)
    delete(OUTPUT)
    OUTPUT = tmp
    delete(tmp)

    ;; Format time coordinate
    ctime = time_attrib(OUTPUT&time)
    delete(OUTPUT&time)
    OUTPUT&time = ctime
    delete(ctime)

    ;; Format latitude coordinate
    OUTPUT!1 = "lat"
    if (isMonotonic(OUTPUT&lat) .eq. 0) then
        error_msg("f", diag_script, "", "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(OUTPUT&lat) .eq. -1) then
        OUTPUT = OUTPUT(:, ::-1, :)
    end if
    clat = lat_attrib(OUTPUT&lat)
    delete(OUTPUT&lat)
    OUTPUT&lat = clat
    delete(clat)

    ;; Format longitude coordinate
    OUTPUT!2 = "lon"
    if (isMonotonic(OUTPUT&lon) .eq. 0) then
        error_msg("f", diag_script, "", "non-monotonic longitude coordinate")
    end if
    if (any(OUTPUT&lon.lt.0.)) then
        OUTPUT = lonFlip(OUTPUT)
    end if
    clon = lon_attrib(OUTPUT&lon)
    delete(OUTPUT&lon)
    OUTPUT&lon = clon
    delete(clon)

    info_output("  Range: " + min(OUTPUT) + "; " + max(OUTPUT), verbosity, 1)

    ;; Set global attributes
    gAtt = True
    gAtt@comment       = "Combined dataset ERS2-ATSR2 (1998-2002) and " + \
        "ENVISAT-AATSR (2003-2011), based on the University of Swansea " + \
        "algorithm (monthly mean L3 data)"
    gAtt@history       = "Created on " + systemfunc("date")
    gAtt@host          = systemfunc("echo $HOST")
    gAtt@user          = systemfunc("echo $USER")
    gAtt@period        = YEAR1 + "-" + YEAR2
    gAtt@field         = FIELD
    gAtt@tier          = TIER
    gAtt@source       = "http://www.icare.univ-lille1.fr/drupal/archive/" + \
        "?dir=CCI-Aerosols/"
    gAtt@reference     = "Not available"
    gAtt@title = OBSNAME + " satellite data reformatted for the ESMValTool"
    gAtt@conventions = "CF/CMOR"

    ;; Outfile
    fout = OUTDIR + "OBS_" + OBSNAME + "_sat_SU-v4.2_" + FIELD + \
        "_" + VARS + ".nc"

    ;; Write variable
    write_nc(fout, VARS, OUTPUT, gAtt)
    delete(gAtt)
    delete(OUTPUT)

end
