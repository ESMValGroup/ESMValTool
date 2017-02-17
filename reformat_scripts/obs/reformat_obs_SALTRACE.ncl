;;#############################################################################
;; REFORMAT SCRIPT FOR SALTRACE OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Not available
;;
;; Last access
;;    Data provided by B. Weinzierl.
;;
;; Download and processing instructions
;;    Data has been untabified and converted with dos2unix.
;;
;; Caveats
;;    The data for the 25 flights are collected and vertically distributed
;;    in bins of 1-km using the provided (measured) pressure values.
;;    Statistics is then calculated for each bin.
;;    The variable massconc_rbc_ngkg_corr is used, which includes a correction
;;    to account for the mass outside the SP2 size range (80-480 nm), assuming
;;    a lognormal distribution.
;;
;; Modification history
;;    20150316-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "SALTRACE"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T0M"

    ;; Standard names
    NAME = "mmrbc"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_SALTRACE.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create vertical bins
    alt_i = 1000. * ispan(0, 13, 1)
    alt_m = 1000. * ispan(0, 12, 1) + 500.
    tmp = stdatmus_z2tdp(alt_i)
    press_i = tmp(2, :)
    delete(tmp)
    tmp = stdatmus_z2tdp(alt_m)
    press_m = tmp(2, :)
    delete(tmp)
    delete(alt_i)
    delete(alt_m)

    ;; Create data array
    NLEVS = dimsizes(press_m)
    NCOLS = 12
    OUTPUT = new((/NLEVS, NCOLS/), float)
    OUTPUT!0 = "plev"
    OUTPUT!1 = "column"
    OUTPUT&plev = press_m * 100.
    OUTPUT&column = (/"N", "min", "max", "mean", "stddev", "5", "10", "25", \
                      "median", "75", "90", "95"/)
    OUTPUT&plev@long_name     = "pressure"
    OUTPUT&plev@positive      = "down"
    OUTPUT&plev@axis          = "Z"
    OUTPUT&plev@units         = "Pa"
    OUTPUT&plev@standard_name = "air_pressure"

    ;; Loop over all data files
    infile = systemfunc("ls " + INDIR + OBSNAME + "*R3.nas")
    do ff = 0, dimsizes(infile) - 1

        ;; Initialize
        DATA = OUTPUT

        ;; Read data
        head = readAsciiHead(infile(ff), "UTSECONDS  ")
        skip = dimsizes(head)
        data = readAsciiTable(infile(ff), 1, "string", skip)
        delete(head)
        delete(skip)

        ;; Get campaign/flight name
        CAMPAIGN = str_sub_str(str_get_field( \
            bname(infile(ff)), 3, "_"), "R3.nas", "")

        info_output("Processing " + CAMPAIGN, verbosity, 1)

        ;; Set time period
        date = toint(str_get_field(data(:, 0), 2, " "))
        date1 = tostring(min(date))
        date2 = tostring(max(date))
        delete(date)

        PERIOD = str_get_cols(date1, 0, 3) + "/" + \
            str_get_cols(date1, 4, 5) + "/" + \
            str_get_cols(date1, 6, 7) + "-" + \
            str_get_cols(date2, 0, 3) + "/" + \
            str_get_cols(date2, 4, 5) + "/" + \
            str_get_cols(date2, 6, 7)

        delete(date1)
        delete(date2)

        ;; Set area
        lat = tofloat(str_get_field(data(:, 0), 5, " "))
        lat@_FillValue = -999.
        lat = where(lat.lt.-1000., lat@_FillValue, lat)
        MINLAT  = min(lat)
        MAXLAT = max(lat)
        delete(lat)

        lon = tofloat(str_get_field(data(:, 0), 6, " "))
        lon@_FillValue = -999.
        lon = where(lon.lt.-1000., lon@_FillValue, lon)
        MINLON = min(lon)
        MAXLON = max(lon)
        MINLON = where(MINLON.lt.0., MINLON + 360., MINLON)
        MAXLON = where(MAXLON.lt.0., MAXLON + 360., MAXLON)
        delete(lon)

        ;; Read measured pressure values
        press = tofloat(str_get_field(data(:, 0), 9, " "))
        press@_FillValue = -999.
        press = where(press.le.0., press@_FillValue, press)

        ;; Read measured data (MASSCONC_RBC_STP_CORR) and fix missing values
        var = tofloat(str_get_field(data(:, 0), 14, " "))
        var@_FillValue = -999.
        var = where(var.le.0., var@_FillValue, var)

        ;; Convert units
        var = var * 1.e-12  ; [ng/kg] --> [kg/kg]
        
        ;; Exlcude points without a valid pressure value
        valid = ind(.not.ismissing(press))
        tmp = var(valid)
        delete(var)
        var = tmp
        delete(tmp)
        tmp = press(valid)
        delete(press)
        press = tmp
        delete(tmp)
        delete(valid)

        ntot = dimsizes(var)
        ntot_check = 0

        ;; Loop over verical bins, assign data and calculate stats
        do pp = 0, NLEVS - 1

            if (pp.eq.0) then
                pmax = 2000.  ; account for p > 1013.25 at the ground
            else
                pmax = press_i(pp)
            end if
            pmin = press_i(pp + 1)

            idx = ind(press.le.pmax .and. press.gt.pmin)
            if (all(ismissing(idx))) then
                delete(idx)
                continue
            end if

            subdata = var(idx)

            ;; Add to total number of data for later check
            ntot_check = ntot_check + dimsizes(idx)

            ;; Exclude levels with less than 10 points
            if (num(.not.ismissing(subdata)).lt.10) then
                delete(subdata)                
                delete(idx)
                continue
            end if

            ;; Calculate all relevant statistics
            do st = 0, dimsizes(DATA&column) - 1
                DATA(pp, st) = distrib_stats(subdata, DATA&column(st))
            end do

            delete(idx)
            delete(subdata)

        end do

        ;; Check for  data loss
        if (ntot.ne.ntot_check) then
            error_msg("f", diag_script, "", "data loss")
        end if
        delete(data)
        delete(press)
        delete(var)

        ;; Set variable attributes (CMOR)
        tmp = var_attrib(DATA, NAME)
        delete(DATA)
        DATA = tmp
        delete(tmp)

        ;; Append global attributes
        gAtt = True
        gAtt@comment2    = "Data binned in 1-km altitude slices, bins " + \
            "with less than 10 valid points excluded"
        gAtt@comment1    = "Data are corrected for the undetected mass " + \
            "outside the SP2 size range (80-480 nm), assuming a " + \
            "lognormal distribution"
        gAtt@history     = "Created on " + systemfunc("date")
        gAtt@host        = systemfunc("echo $HOST")
        gAtt@user        = systemfunc("echo $USER")
        gAtt@lonmax      = where(MAXLON.lt.0., MAXLON + 360., MAXLON)
        gAtt@lonmin      = where(MINLON.lt.0., MINLON + 360., MINLON)
        gAtt@latmax      = MAXLAT
        gAtt@latmin      = MINLAT
        gAtt@info        = CAMPAIGN
        gAtt@period      = PERIOD
        gAtt@field       = FIELD
        gAtt@tier        = TIER
        gAtt@source      = "Data provided by Bernadett Weinzierl (DLR)"
        gAtt@reference   = "Not available"
        gAtt@title       = OBSNAME + \
            " vertical profiles reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        d1 = str_get_field(str_get_field(PERIOD, 1, "-"), 1, "/") + \
            str_get_field(str_get_field(PERIOD, 1, "-"), 2, "/")
        d2 = str_get_field(str_get_field(PERIOD, 2, "-"), 1, "/") + \
            str_get_field(str_get_field(PERIOD, 2, "-"), 2, "/")
        fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_" + \
            CAMPAIGN + "_" + FIELD + "_" + NAME + "_" + d1 + "-" + d2 + ".nc"

        ;; Write variable
        write_nc_profile(fout, NAME, DATA, gAtt)
        delete(DATA)

    end do

    leave_msg(diag_script, "", 1)

end
