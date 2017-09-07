;;#############################################################################
;; REFORMAT SCRIPT FOR OAFlux OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://oaflux.whoi.edu/data.html
;;
;;    Please acknowledge the global ocean heat flux and evaporation products
;;    were provided by the WHOI OAFlux project (http://oaflux.whoi.edu) funded
;;    by the NOAA Climate Observations and Monitoring (COM) program.
;;
;; Last access
;;    20150722
;;
;; Download and processing instructions
;;    Connect to ftp.whoi.edu
;;    Browse to /pub/science/oaflux/data_v3/monthly/turbulence
;;    Download and gunzip all the lh_oaflux_YYYY.nc.gz and sh_oaflux_YYYY.nc.gz
;;    Set the variable to be processed as shortname below, choosing among:
;;      hl - Surface latent heat flux (hfls)
;;      sl - Surface sensible heat flux (hfss)
;;      qa - Surface sensible heat flux (huss)
;;      ta - Near-Surface Air Temperature (tas)
;;      ts - Surface Temperature (ts)
;;      ws - Near-Surface Wind Speed (sfcWind)
;;
;; Caveats
;;    The data type in the WHOI data set (short) is inconsistent with the
;;    missing_value type (float) - this produces a warning from NCL
;;
;; Modification history
;;    20160502-A_righ_ma: repalcing input argument with hard setting to allow
;;                        use withing namelist_reformat_obs.xml.
;;    20150722-A_righ_ma: added @reference attribute and pep8.
;;    20150421-A_eval_ma: written.
;;
;;#############################################################################

function get_settings(res[1]:logical,
                      req_variable[1]:string)
;; Settings for the various variables availabe in the WHOI data set
begin

    ;; Default settings
    res@year1 = 1958
    res@year2 = 2013
    res@orig_missing_value = 32766.
    res@variable_offset = 0.

    res@errname = "err"
    res@attrCopy_err2err = (/"comments"/)
    res@attrCopy_var2err = (/"units"/)

    res@raw_data_folder = getenv("ESMValTool_RAWOBSPATH") + \
        "/Tier2/WHOI-OAFlux/"
    res@processed_data_folder = getenv("ESMValTool_OBSPATH") + \
        "/Tier2/WHOI-OAFlux/"

    ;; Variable specific settings
    ;; Monthly mean surface latent heat flux, positive upward"
    if (req_variable .eq. "hl") then
        res@filename_prefix = "lh_oaflux"

        res@orig_vars = "lhtfl"
        res@vars = "hfls"
        res@scaling_factor = 0.1

    ;; Monthly mean surface sensible heat flux, positive upward"
    else if (req_variable .eq. "sl") then
        res@filename_prefix = "sh_oaflux"

        res@orig_vars = "shtfl"
        res@vars = "hfss"
        res@scaling_factor = 0.1

    ;; Monthly mean specific humidity at 2m
    else if (req_variable .eq. "qa") then
        res@filename_prefix = "qa_oaflux"

        res@orig_vars = "hum2m"
        res@vars = "huss"
        res@scaling_factor = 0.01

    ;; Monthly mean air temperature at 2m
    else if (req_variable .eq. "ta") then
        res@filename_prefix = "ta_oaflux"

        res@orig_vars = "tmp2m"
        res@vars = "tas"
        res@scaling_factor = 0.01
        res@variable_offset = 273.15

    ;; Monthly mean sea surface temperature
    else if (req_variable .eq. "ts") then
        res@filename_prefix = "ts_oaflux"

        res@orig_vars = "tmpsf"
        res@vars = "ts"
        res@scaling_factor = 0.01
        res@variable_offset = 273.15

    ;; Monthly mean neutral wind speed at 10m
    else if (req_variable .eq. "ws") then
        res@filename_prefix = "ws_oaflux"

        res@orig_vars = "wnd10"
        res@vars = "sfcWind"
        res@scaling_factor = 0.01

    end if
    end if
    end if
    end if
    end if
    end if

    return(res)
end

function append_history(infile[1]:file,
                        to_add[1]:string)
;; Creates (or appends) a history string
local newline, f, datestamp, history
begin

    newline = str_get_nl()
    datestamp = systemfunc("date -u")

    ;; Check existing history attribute
    if (isatt(infile, "history")) then
        history = infile@history
    else
        history = ""
    end if

    ;; Create or append history
    if (strlen(history) .gt. 0) then
        history = datestamp + ": " + to_add + newline + history
    else
        history = datestamp + ": " + to_add
    end if
    return(history)

end

function get_output_filename(res[1]:logical,
                             year[1]:integer)
;; Define the CMIP5 DRS output filename
local outfile_string
begin

    outfile_string = res@processed_data_folder + "/" + \
        "OBS_WHOI-OAFlux_ocean_1_" + res@vars + "_" + \
        tostring(year) + "01-" + tostring(year) + "12.nc"

    if (isfilepresent(outfile_string)) then
        system("rm -f " + outfile_string)
    end if

    return(outfile_string)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity =  stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_OAFlux.ncl"
    enter_msg(diag_script, "", 1)

    ;; Variable
    ;;   hl - Surface latent heat flux (hfls)
    ;;   sl - Surface sensible heat flux (hfss)
    ;;   qa - Surface sensible heat flux (huss)
    ;;   ta - Near-Surface Air Temperature (tas)
    ;;   ts - Surface Temperature (ts)
    ;;   ws - Near-Surface Wind Speed (sfcWind)
    shortname = "hl"

    res = True
    res = get_settings(res, shortname)

    ;; Create output dir
    if (.not. isfilepresent_esmval(res@processed_data_folder)) then
        system("mkdir -p " + res@processed_data_folder)
    end if

    ;; Loop over all years (WHOI data is provided in yearly files)
    do yy = res@year1, res@year2

        ;; Read input variable and fix _FillValue
        ;; (see caveat abover regarding the warning)
        infile = addfile(res@raw_data_folder + "/"\
                         + res@filename_prefix + "_"\
                         + tostring(yy) + ".nc", "r")
        invar = infile->$res@orig_vars$
        invar_coords = invar
        invar := tofloat(invar)
        invar := where(invar .eq. 32766., default_fillvalue("float"), invar)
        invar@_FillValue = default_fillvalue("float")
        invar = invar * res@scaling_factor
        invar = invar + res@variable_offset
        copy_VarCoords(invar_coords, invar)

        ;; Reformat the time coordinate
        time = new(12, double)
        cal_res = 0
        cal_res@calendar = "proleptic_gregorian"
        do mm = 1, 12
            time(mm - 1) = \
                cd_inv_calendar(yy, mm, 15, 0, 0, 0, TUNITS, cal_res)
        end do

        FIELD = "T2Ms"
        time@YEAR1 = yy
        time@YEAR2 = yy
        time@FIELD = FIELD
        time := time_attrib(time)

        ;; Create output variable
        outvar = new((/12, dimsizes(invar&lat), dimsizes(invar&lon)/), float)

        ;; Set variable attributes
        outvar := var_attrib(invar, res@vars)

        delete(outvar&time)
        outvar!0 = "time"
        outvar&time = time

        ;; If error estimate is defined, copy it along as well
        if (isatt(res, "errname")) then
            if (isfilevar(infile, res@errname)) then
                errvar = infile->$res@errname$
                errvar_coords = errvar
                errvar := tofloat(errvar)
                errvar := where(errvar .eq. 32766., outvar@_FillValue, errvar)
                errvar@_FillValue = outvar@_FillValue
                errvar = errvar * res@scaling_factor
                errvar = errvar + res@variable_offset
                copy_VarCoords(errvar_coords, errvar)

                ;; Copy along manually defined attributes from
                ;; input error estimate
                if (isatt(res, "attrCopy_err2err")) then
                    do idx_attr = 0, dimsizes(res@attrCopy_err2err) - 1
                        attribute = res@attrCopy_err2err(idx_attr)
                        errvar@$attribute$  = errvar_coords@$attribute$
                    end do
                end if

                ;; Copy along manually defined attributes from input variable
                if (isatt(res, "attrCopy_var2err")) then
                    do idx_attr = 0, dimsizes(res@attrCopy_var2err) - 1
                        attribute = res@attrCopy_var2err(idx_attr)
                        errvar@$attribute$  = outvar@$attribute$
                    end do
                end if
            end if
        end if

        if (isatt(res, "errname")) then
            delete(errvar&time)
            errvar!0 = "time"
            errvar&time = time
        end if

        ;; Reformat latitude coordinate
        outvar!1 = "lat"
        outvar&lat = invar&lat
        if (isMonotonic(outvar&lat) .eq. 0) then
            error_msg("f", diag_script, "",\
                      "non-monotonic latitude coordinate")
        end if
        if (isMonotonic(outvar&lat) .eq. -1) then
            outvar = outvar(:, ::-1, :)
        end if
        clat := lat_attrib(outvar&lat)
        outvar&lat := clat

        if (isatt(res, "errname")) then
            delete(errvar&lat)
            errvar!1 = "lat"
            errvar&lat = clat
        end if

        ;; Reformat longitude coordinate
        outvar!2 = "lon"
        outvar&lon = invar&lon
        if (isMonotonic(outvar&lon) .eq. 0) then
            error_msg("f", diag_script, "",\
                      "non-monotonic longitude coordinate")
        end if
        if (any(outvar&lon .lt. 0.)) then
            outvar = lonFlip(outvar)
        end if
        clon := lon_attrib(outvar&lon)
        outvar&lon := clon

        if (isatt(res, "errname")) then
            delete(errvar&lon)
            errvar!2 = "lon"
            errvar&lon = clon
        end if

        ;; Write reformatted variable/error estimate to file
        outfile_string = get_output_filename(res, yy)
        outfile = addfile(outfile_string, "c")
        filedimdef(outfile, "time", -1, True)

        outfile->$res@vars$ = outvar
        if (isatt(res, "errname")) then
            outfile->$res@errname$ = errvar
        end if

        delete(time)
        delete(errvar)
        delete(errvar_coords)
        delete(invar)
        delete(invar_coords)
        delete(outvar)

        ;; Copy along global attributes
        global_attnames = getvaratts(infile)
        do i = 0, dimsizes(global_attnames) - 1
            old_attr_prefix = "original-attr_"
            attr = old_attr_prefix + global_attnames(i)
            ;; 'history' attribute is handled below
            if (.not. (global_attnames(i) .eq. "history")) then
                outfile@$attr$ =  infile@$global_attnames(i)$
            end if
        end do

        ;; Set global attributes
        outfile@history = \
            append_history(infile, "data reformatted to be used in ESMValTool")
        outfile@host    = systemfunc("hostname")
        outfile@user    = systemfunc("whoami")
        outfile@field   = FIELD
        outfile@tier    = 2
        outfile@source  = "ftp://ftp.whoi.edu/pub/science/oaflux/data_v3/" + \
            "monthly/turbulence/"
        outfile@reference = "The global ocean heat flux and evaporation " + \
            "products were provided by the WHOI OAFlux project " + \
            "(http://oaflux.whoi.edu) funded by the NOAA Climate " + \
            "Observations and Monitoring (COM) program."
        outfile@title = "WHOI OAFlux reanalysis data reformatted for the " + \
            "ESMValTool"
        outfile@conventions = "CF/CMOR"

    end do

    leave_msg(diag_script, "", 1)

end
