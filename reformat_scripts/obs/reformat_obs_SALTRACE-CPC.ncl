;;#############################################################################
;; REFORMAT SCRIPT FOR SALTRACE-CPC OBSERVATIONAL DATA
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
;;
;;    The variable ambient_N5 is used, which includes a correction
;;    to account for pressure-dependent loss of (flow) efficiency relative to 
;;    CPC2 (N14) using a broken power law parametrization with a break point
;;    around 400hPa, based on lab measurements performed during and after 
;;    SALTRACE.
;;
;; Modification history
;;    20170823-A_beer_ci: written, based on SALTRACE-SP2.
;;
;;##############################################################################

begin

    ;; Source name
    OBSNAME = "SALTRACE"
    SUBNAME = "CPC"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/" + SUBNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T0M"

    ;; Standard names
    VARS = (/"conccnd5", "conccnd14", "conccnd10"/)

    ;; Names in raw data
    NAMES = (/"AMBIENT_N5", "AMBIENT_N14", "AMBIENT_N10TD"/)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_SALTRACE-CPC.ncl"
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

    ;; Loop over variables
    do varID = 0, dimsizes(VARS) - 1

        ;; Initialize
        OUTPUT = OUTPUT@_FillValue

        ;; Loop over all data files
        infile = systemfunc("ls " + INDIR + "*v01.ames")
        do ff = 0, dimsizes(infile) - 1

            ;; Initialize
            DATA = OUTPUT

            ;; Read data
            head = readAsciiHead(infile(ff), "UTSECONDS  ")
            skip = dimsizes(head)
            data = readAsciiTable(infile(ff), 1, "string", skip)
            delete(skip)
            head := head(dimsizes(head)-1)

            ;; Find column number for this variable
            ncols = str_fields_count(head, " ")
            cols = new(ncols, string)
            do cc = 1, ncols
                cols(cc - 1) = str_get_field(head, cc, " ")
            end do
            idcol = ind(NAMES(varID).eq.cols) + 1
            if (dimsizes(idcol).gt.1) then
                error_msg("f", diag_script, "", "Multiple entries for " + \
                          "variable " + NAMES(varID) + " in the input data")
            end if
            if (ismissing(idcol)) then
                error_msg("f", diag_script, "", "Variable " + NAMES(varID) + \
                          " not found in the input data")
            end if
            delete(head)
            delete(cols)

            ;; Get campaign/flight name
            CAMPAIGN = str_sub_str(str_get_field( \
                bname(infile(ff)), 1, "_"), "cpc_v01.ames", "")

            info_output("Processing " + CAMPAIGN, verbosity, 1)

            ;; Set time period
            date = "20" + \
                str_get_cols(str_get_field(bname(infile(ff)), 1, "_"), 0, 5)
            PERIOD = str_get_cols(date, 0, 3) + "/" + \
                str_get_cols(date, 4, 5) + "/" + \
                str_get_cols(date, 6, 7) + "-" + \
                str_get_cols(date, 0, 3) + "/" + \
                str_get_cols(date, 4, 5) + "/" + \
                str_get_cols(date, 6, 7)
            delete(date)

            ;; Set area
            lat = tofloat(str_get_field(data(:, 0), 2, " "))
            lat@_FillValue = -999.
            lat = where(lat.lt.-1000., lat@_FillValue, lat)
            MINLAT  = min(lat)
            MAXLAT = max(lat)
            delete(lat)

            lon = tofloat(str_get_field(data(:, 0), 3, " "))
            lon@_FillValue = -999.
            lon = where(lon.lt.-1000., lon@_FillValue, lon)
            MINLON = min(lon)
            MAXLON = max(lon)
            MINLON = where(MINLON.lt.0., MINLON + 360., MINLON)
            MAXLON = where(MAXLON.lt.0., MAXLON + 360., MAXLON)
            delete(lon)

            ;; Read measured pressure values
            press = tofloat(str_get_field(data(:, 0), 6, " "))
            press@_FillValue = -999.
            press = where(press.le.0., press@_FillValue, press)

            ;; Read measured data (AMBIENT_N5) and fix missing values
            var = tofloat(str_get_field(data(:, 0), 9, " "))
            var@_FillValue = -999.
            var = where(var.le.0., var@_FillValue, var)

            ;; Convert units
            var = var * 1.e6    ; [1/cm3] --> [1/m3]

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

            ;; Loop over vertical bins, assign data and calculate stats
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
            tmp = var_attrib(DATA, VARS(varID))
            delete(DATA)
            DATA = tmp
            delete(tmp)

            ;; Append global attributes
            gAtt = True
            gAtt@comment2    = "Data binned in 1-km altitude slices, bins " + \
                "with less than 10 valid points excluded"
            gAtt@comment1    = "Data are corrected for pressure-dependent " + \
                "loss of (flow) efficiency relative to CPC2 (N14) using a " + \
                "broken power law parametrization"
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
            gAtt@title       = OBSNAME + "-" + SUBNAME + \
                " vertical profiles reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Outfile
            d1 = str_get_field(str_get_field(PERIOD, 1, "-"), 1, "/") + \
                str_get_field(str_get_field(PERIOD, 1, "-"), 2, "/")
            d2 = str_get_field(str_get_field(PERIOD, 2, "-"), 1, "/") + \
                str_get_field(str_get_field(PERIOD, 2, "-"), 2, "/")
            fout = OUTDIR + "OBS_" + OBSNAME + "-" + SUBNAME + "_campaign_" + \
                CAMPAIGN + "_" + FIELD + "_" + VARS(varID) + "_" + \
                d1 + "-" + d2 + ".nc"

            ;; Write variable
            write_nc_profile(fout, VARS(varID), DATA, gAtt)
           delete(DATA)

        end do
    end do

    leave_msg(diag_script, "", 1)

end
