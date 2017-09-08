;;#############################################################################
;; REFORMAT SCRIPT FOR NCEP REANALYISIS DATA ON PRESSURE LEVELS
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
;;
;; Last access
;;    20150914
;;
;; Download and processing instructions
;;    Specify the time range of the data as YEAR1-YEAR2 below, considering only
;;    complete years (Jan to Dec).
;;
;; Caveats
;;    Field and time range is expected to be the same for all variables.
;;
;; Modification history
;;    20150914-A_laue_ax: modified for daily fields.
;;    20140128-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "NCEP"

    ;; Tier
    TIER = 2

    ;; Period
    YEAR1 = 1998
    YEAR2 = 2008

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Selected variable (standard name)
    VARS = (/"ua", "va", "pr", "rlut"/)

    ;; Name in the raw data
    NAME = (/"uwnd", "vwnd", "prate", "ulwrf"/)

    ;; suffix for raw data filenames
    SUFFIX = (/"", "", ".sfc.gauss", ".ntat.gauss" /)

    ;; Expected units (according to CMOR standard)
    EXPUNITS = (/"m/s", "m/s", "Kg/m^2/s", "W/m^2"/)

    ;; Fields
    FLDARR = (/"T3D", "T3D", "T2D", "T2D"/)

    ;; Initialize global variable
    FIELD = ""

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_NCEP-daily.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        ;; Set field
        FIELD = FLDARR(vID)

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ", \
                    verbosity, 1)

        do yr = YEAR1, YEAR2

            fname = INDIR + NAME(vID) + SUFFIX(vID) + "." + yr + ".nc"
            f = addfile(fname, "r")
            tmp = f->$NAME(vID)$

            fill_start = num(ismissing(tmp))

            ;; Extract time range (consider only complete years)
            date = cd_calendar(tmp&time, 0)
            idx1 = 0  ;ind(date(:, 1).eq.1 .and. date(:, 0).eq.YEAR1)
            idx2 = dimsizes(tmp&time) - 1  ;ind(date(:, 1).eq.12 .and. date(:, 0).eq.YEAR2)
            if (isStrSubset(FIELD, "T3")) then
                output = tmp(idx1:idx2, :, :, :)
            else
                output = tmp(idx1:idx2, :, :)
            end if
            delete(tmp)
            delete(date)

            ;; Check time range consistency
            date = cd_calendar(output&time, 0)
            if (min(date(:, 0)).lt.YEAR1 .or. max(date(:, 0)).gt.YEAR2) then
                error_msg("f", diag_script, "", "given time range (" + \
                          YEAR1 + "-" + YEAR2 + ")" + \
                          " not consistent with data (" + \
                          min(date(:, 0)) + "-" + max(date(:, 0)) + ")")
            end if
            delete(date)

            ; Unpack variable according to metadata information
            if (isatt(output, "scale_factor") .or. \
                isatt(output, "add_offset")) then
                tmp = output * output@scale_factor + output@add_offset
                copy_VarMeta(output, tmp)
                delete(output)
                output = tmp
                delete(tmp)
            end if

            ;; Convert units
            if (output@units.ne.EXPUNITS(vID)) then
                error_msg("f", diag_script, "", \
                          "Possibly wrong units for " + VARS(vID))
            end if

            info_output("  Range: " + min(output) + "; " + max(output) + \
                        " [" + output@units + "]", verbosity, 1)

            ;; Set variable attributes
            tmp = var_attrib(output, VARS(vID))
            delete(output)
            output = tmp
            delete(tmp)

            ;; Format time coordinate
;            ctime = time_attrib(output&time)
;            delete(output&time)
;            output&time = ctime
;            delete(ctime)

            output&time@standard_name = "time"
            output&time@units = "days since 1950-01-01 00:00:00"
            output&time@axis = "T"
            output&time@long_name = "time"
            output&time@calendar = "standard"
            if (isatt(output&time, "_FillValue")) then
                delete(output&time@_FillValue)
            end if

            output&time(idx1) = cd_inv_calendar(yr, 1, 1, 0, 0, 0, \
                                "days since 1950-01-01 00:00:00", 0)

            do i = idx1 + 1, idx2
                output&time(i) = output&time(i-1) + 1
            end do

            if (isStrSubset(FIELD, "T3")) then
                ;; Format vertical coordinate
                output!1 = "plev"
                if (isMonotonic(output&plev) .eq. 0) then
                    error_msg("f", diag_script, "", \
                              "Non-monotonic vertical coordinate")
                end if
                if (isMonotonic(output&plev).eq.1) then
                    output = output(:, ::-1, :, :)
                end if
                if (output&plev@units.eq."millibar") then
                    output&plev = output&plev * 100
                else
                    error_msg("f", diag_script, "", \
                              "Unexpected units for plev " + output&plev@units)
                end if
                cplev = plev_attrib(output&plev)
                delete(output&plev)
                output&plev = cplev
                delete(cplev)
            end if

            ;; Format latitude coordinate
            if (isStrSubset(FIELD, "T3")) then
                output!2 = "lat"
                if (isMonotonic(output&lat) .eq. 0) then
                    error_msg("f", diag_script, "", \
                              "Non-monotonic latitude coordinate")
                end if
                if (isMonotonic(output&lat) .eq. -1) then
                    output = output(:, :, ::-1, :)
                end if
            else
                output!1 = "lat"
                if (isMonotonic(output&lat) .eq. 0) then
                    error_msg("f", diag_script, "", \
                              "non-monotonic latitude coordinate")
                end if
                if (isMonotonic(output&lat) .eq. -1) then
                    output = output(:, ::-1, :)
                end if
            end if
            clat = lat_attrib(output&lat)
            delete(output&lat)
            output&lat = clat
            delete(clat)

            ;; Format longitude coordinate
            if (isStrSubset(FIELD, "T3")) then
                output!3 = "lon"
            else
                output!2 = "lon"
            end if
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

            info_output("  Range: " + min(output) + "; " + max(output), \
                        verbosity, 1)

            ;; Set global attributes
            gAtt = True
            gAtt@history = "Created on " + systemfunc("date")
            gAtt@host    = systemfunc("echo $HOST")
            gAtt@user    = systemfunc("echo $USER")
            gAtt@period  = yr
            gAtt@field   = FIELD
            gAtt@tier    = TIER
            gAtt@source  = "http://www.esrl.noaa.gov/psd/data/gridded/" + \
                "data.ncep.reanalysis.html"
            gAtt@reference = "Kalnay, E. et al., B. Am. Meteorol. Soc., " + \
                "77, 437-471, " + \
                "doi:10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2, 1996."
            gAtt@title = OBSNAME + \
                " reanalysis data reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Check fill values
            fill_end = num(ismissing(output))
            if (fill_start.ne.fill_end) then
                error_msg("f", diag_script, "", \
                          "Missing values lost during conversion")
            end if
            delete(fill_start)
            delete(fill_end)

            ;; Outfile
            fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
                "_" + VARS(vID) + "_" + yr + "01-" + yr + "12.nc"

            ;; Write variable
            write_nc(fout, VARS(vID), output, gAtt)
            delete(gAtt)
            delete(output)
        end do  ; yr-loop
    end do ; varriable-loop

    leave_msg(diag_script, "", 1)

end
