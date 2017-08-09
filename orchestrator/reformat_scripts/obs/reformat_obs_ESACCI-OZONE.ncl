;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-OZONE SATELLITE DATA (total column ozone)
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.esa-ozone-cci.org/?q=node/160
;;
;; Last access
;;    20160224
;;
;; Download and processing instructions
;;    Login: not needed Password: press enter
;;    Download the data from TC_L3_MRG Level 3.
;;
;; Caveats
;;
;; Modification history
;;    20160224-A_wenz_sa: written based on reformat_obs_ESACCI-AEROSOL.ncl.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ESACCI-OZONE"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1997
    YEAR2 = 2010

    ;; Selected variable (standard name)
    VARS = (/"toz", "tozStderr"/)

    ;; Name in the raw data
    NAME = (/"atmosphere_mole_content_of_ozone", \
        "atmosphere_mole_content_of_ozone_standard_error"/)

    ;; Field
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = 1
    diag_script = "reformat_obs_ESACCI-OZONE.ncl"

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
            do mm = 1, 12

                ldate = yy + sprinti("%0.2i", mm)
                info_output("Processing " + ldate, verbosity, 1)

                ;; File name
                fname = systemfunc("ls " + INDIR + "ESACCI-OZONE-L3S-TC-" + \
                    "MERGED-DLR_1M-" + ldate + "??-fv0100.nc")

                ;; Check
                if (all(ismissing(fname))) then
                    error_msg("f", diag_script, \
                      "", "no file found for date " + ldate)
                end if

                ;; Extract data
                f = addfile(fname(0), "r")
                xx = f->$NAME(vID)$
                xx@_FillValue = 1e+20
                xx@missing_value = xx@_FillValue
                xx = where(xx.lt.0., xx@_FillValue, xx)

                ;; Assign to global array
                if (.not.isdefined("OUTPUT")) then
                    dims = array_append_record(dimsizes(time), dimsizes(xx), 0)
                    OUTPUT = new(dims, typeof(xx))
                    OUTPUT!0 = "time"
                    OUTPUT&time = time
                    OUTPUT!1 = "lat"
                    OUTPUT&lat = f->latitude
                    OUTPUT!2 = "lon"
                    OUTPUT&lon = f->longitude
                end if
                OUTPUT(ind(toint(ldate).eq.date), :, :) = (/xx/)

                delete(fname)
                delete(xx)

            end do
        end do

        ;; Set fill value
        ;OUTPUT = where(OUTPUT.eq.-9.9, OUTPUT@_FillValue, OUTPUT)

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
        gAtt@comment = "L3 total ozone columns merged using the" + \
                       " algorithm GTO-ECV as part of the ESA Ozone CCI"
        gAtt@history = "Created on " + systemfunc("date")
        gAtt@host    = systemfunc("echo $HOST")
        gAtt@user    = systemfunc("echo $USER")
        gAtt@period  = YEAR1 + "-" + YEAR2
        gAtt@field   = FIELD
        gAtt@tier    = TIER
        gAtt@source  = "http://www.esa-ozone-cci.org/?q=node/160"
        gAtt@reference     = "GTO-ECV reference: Loyola et al.," + \
                       " Global long-term monitoring of the ozone layer" + \
                       " - a prerequisite for predictions, International" + \
                       " Journal of Remote Sensing, vol. 30, no. 15, pp. 4295-4318, 2009."
        gAtt@title = OBSNAME + " satellite data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_sat_L3_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), OUTPUT, gAtt)
        delete(gAtt)
        delete(OUTPUT)

    end do

end
