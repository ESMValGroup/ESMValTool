;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-CLOUD SATELLITE DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/
;;
;; Last access
;;    20160519
;;
;; Download and processing instructions
;;    Download the data from:
;;    cloud/data/L3C/merisaatsr_envisat/phase1_v2.0/yyyy/mm/
;;
;; Caveats
;;
;; Modification history
;;    20160519-A_laue_ax: written (reformat_obs_ESACCI-AEROSOL.ncl).
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ESACCI-CLOUD"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 2007
    YEAR2 = 2009

    ;; Selected variable (standard name)
    VARS = (/"clt", "clwvi", "clivi", \
             "lwpStderr", "iwpStderr", "cltStderr"/)

    ;; Name in the raw data
    NAME = (/"cc_total", "cwp", "iwp", \
             "lwp_uncertainty", "iwp_uncertainty", "cc_total_uncertainty"/)

    ;; Field
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_ESACCI-CLOUD.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    do vID = 0, dimsizes(VARS) - 1

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ", \
                    verbosity, 1)

        time = create_timec(YEAR1, YEAR2)
        date = cd_calendar(time, 1)

        ;; Create timeseries
        do yy = YEAR1, YEAR2
            syear = sprinti("%i", yy)
            do mm = 1, 12

                smonth = sprinti("%0.2i", mm)
                info_output("Processing " + syear + "-" + smonth, verbosity, 1)

                ;; Read file
                fname = INDIR + syear + "/" + smonth + "/" + syear + smonth + \
                    "-ESACCI-L3C_CLOUD-CLD_" + \
                    "PRODUCTS-MERISAATSR_ENVISAT-fv2.0.nc"

                ;; No files found
                if (all(ismissing(fname))) then
                    continue
                end if

                ;; Extract data
                f = addfile(fname, "r")
                xx = f->$NAME(vID)$

                ;; if needed, correct units
                if (any(NAME(vID).eq. \
                    (/"cc_total", "cc_total_uncertainty"/))) then
                    xx = xx * 100.0  ; convert to %
                end if
                if (any(NAME(vID).eq.\
                    (/"cwp", "iwp", "lwp_uncertainty", \
                      "iwp_uncertainty"/))) then
                    xx = xx * 0.001  ; convert to kg/m2
                end if

                ;; Assign to global array
                if (.not.isdefined("OUTPUT")) then
                    dims = dimsizes(xx)
                    dims(0) = dimsizes(time)
                    OUTPUT = new(dims, float)
                    OUTPUT!0 = "time"
                    OUTPUT&time = time
                    OUTPUT!1 = "lat"
                    OUTPUT&lat = f->lat
                    OUTPUT!2 = "lon"
                    OUTPUT&lon = f->lon
                end if

                OUTPUT(ind(toint(yy * 100 + mm).eq.date), :, :) = (/xx/)

                delete(fname)
                delete(f)

            end do
        end do

        ;; Set fill value
        OUTPUT = where(OUTPUT.eq.-999, OUTPUT@_FillValue, OUTPUT)

        ;; Set variable attributes
        tmp = var_attrib(OUTPUT, VARS(vID))
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
            error_msg("f", diag_script, \
                      "", "non-monotonic latitude coordinate")
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
            error_msg("f", diag_script, \
                      "", "non-monotonic longitude coordinate")
        end if
        if (any(OUTPUT&lon.lt.0.)) then
            OUTPUT = lonFlip(OUTPUT)
        end if
        clon = lon_attrib(OUTPUT&lon)
        delete(OUTPUT&lon)
        OUTPUT&lon = clon
        delete(clon)

        info_output("  Range: " + min(OUTPUT) + "; " + \
                    max(OUTPUT), verbosity, 1)

        ;; Set global attributes
        gAtt = True
        gAtt@comment = "ESA Cloud CCI Retrieval Products, " + \
            "MERISAATSR ENVISAT (monthly mean L3 data v2.0)"
        gAtt@history = "Created on " + systemfunc("date")
        gAtt@host    = systemfunc("echo $HOST")
        gAtt@user    = systemfunc("echo $USER")
        gAtt@period  = YEAR1 + "-" + YEAR2
        gAtt@field   = FIELD
        gAtt@tier    = TIER
        gAtt@source  = "http://www.esa-cloud-cci.org/"
        gAtt@reference  = "Hollmann, R. et al. (2015): ESA Cloud Climate " + \
            "Change Initiative (ESA Cloud_cci) data: ESACCI-L3C_CLOUD-CLD_" + \
            "PRODUCTS-MERISAATSR_ENVISAT-fv2.0 via Centre for " + \
            "Environmental Data Analysis, 2016"
        gAtt@title = OBSNAME + " satellite data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_sat_L3C-AATSR-ENVISAT_" + \
            FIELD + "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), OUTPUT, gAtt)
        delete(gAtt)
        delete(OUTPUT)

    end do

    leave_msg(diag_script, "", 1)

end
