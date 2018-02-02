;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-OZONE SATELLITE DATA (vertical ozone profiles)
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.esa-ozone-cci.org
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
    YEAR1 = 2007
    YEAR2 = 2008

    ;; Selected variable (standard name)
    VARS = (/"tro3prof", "tro3profStderr"/)

    ;; Name in the raw data
    NAME = (/"merged_ozone_vmr", "uncertainty_of_merged_ozone"/)

    ;; Field
    FIELD = "T2Mz"

end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/logging.ncl"
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

                ;; Read file
                str = "find " + INDIR + " -type f -name 'ESACCI-OZONE-*LP*" \
                      + ldate + "*fv0002.nc'"
                fname = systemfunc(str)

                ;; No files found
                if (all(ismissing(fname))) then
                    continue
                end if

                ;; Extract data
                f = addfile(fname(0), "r")
                xx = f->$NAME(vID)$
                ;xx@_FillValue = 1e+20
                ;xx@missing_value = xx@_FillValue
                ;xx = where(xx.lt.0., xx@_FillValue, xx)

                ;; Assign to global array
                if (.not.isdefined("OUTPUT")) then
                    dims = array_append_record(dimsizes(time), dimsizes(xx), 0)
                    OUTPUT = new(dims, typeof(xx))
                    OUTPUT!0 = "time"
                    OUTPUT&time = time
                    OUTPUT!1 = "plev"
                    OUTPUT&plev = f->air_pressure
                    OUTPUT!2 = "lat"
                    OUTPUT&lat = f->latitude_centers
                end if
                OUTPUT(ind(toint(ldate).eq.date), :, :) = (/xx/)

                delete(fname)

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

        ;; Format longitudepressure level coordinate
        OUTPUT!1 = "plev"
        if (isMonotonic(OUTPUT&plev) .eq. 0) then
            error_msg("f", diag_script, \
                      "", "non-monotonic pressure levels coordinate")
        end if
        cplev = plev_attrib(OUTPUT&plev)
        delete(OUTPUT&plev)
        OUTPUT&plev = cplev
        delete(cplev)

        ;; Format latitude coordinate
        OUTPUT!2 = "lat"
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
        gAtt@source  = "http://www.esa-ozone-cci.org"
        gAtt@reference     = "NaN"
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
