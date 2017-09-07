;;#############################################################################
;; REFORMAT SCRIPT FOR ACCESS OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Not available.
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
;; Modification history
;;    20160225-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ACCESS-2"

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
    NAME = (/"conccnd10", "conccnd5"/)

    ;; Corresponding columns
    COLID = (/13, 19/)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_ACCESS-2.ncl"
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

    ;; Loop over variables
    do varid = 0, dimsizes(NAME) - 1

        ;; Create data array
        NLEVS = dimsizes(press_m)
        NCOLS = 12
        OUTPUT = new((/NLEVS, NCOLS/), float)
        OUTPUT!0 = "plev"
        OUTPUT!1 = "column"
        OUTPUT&plev = press_m * 100.
        OUTPUT&column = (/"N", "min", "max", "mean", "stddev", \
                          "5", "10", "25", "median", "75", "90", "95"/)
        OUTPUT&plev@long_name     = "pressure"
        OUTPUT&plev@positive      = "down"
        OUTPUT&plev@axis          = "Z"
        OUTPUT&plev@units         = "Pa"
        OUTPUT&plev@standard_name = "air_pressure"

        ;; Loop over all data files
        infile = systemfunc("ls " + INDIR + "*v02.ames")
        do ff = 0, dimsizes(infile) - 1
            
            ;; Initialize
            DATA = OUTPUT

            ;; Read data
            head = readAsciiHead(infile(ff), "UTSECONDS  ")
            skip = dimsizes(head)
            data = readAsciiTable(infile(ff), 1, "string", skip)

            delete(skip)

            ;; Get campaign/flight name
            CAMPAIGN = str_sub_str(str_get_field( \
                bname(infile(ff)), 1, "_"), ".ames", "")

            info_output("Processing " + CAMPAIGN, verbosity, 1)

            ;; Set time period
            date = toint(str_get_field(data(:, 0), 1, " "))
            dunit = "seconds after " + str_get_field(head(6), 1, " ") + "-" + \
                str_get_field(head(6), 2, " ") + "-" + \
                str_get_field(head(6), 3, " ") + " 00:00:00"
            delete(head)

            date1 = min(date)
            date1@units = dunit
            date1 := tostring(cd_calendar(date1,-2))
            date2 = max(date)
            date2@units = dunit
            date2 := tostring(cd_calendar(date2,-2))
            delete(dunit)
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

            ;; Read measured data and fix missing values
            var = tofloat(str_get_field(data(:, 0), COLID(varid), " "))
            var@_FillValue = -999.
            var = where(var.le.0., var@_FillValue, var)

            ;; Convert units
            var = var * 1.e6  ; [cm-3] --> [m-3]
        
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
            tmp = var_attrib(DATA, NAME(varid))
            delete(DATA)
            DATA = tmp
            delete(tmp)

            ;; Append global attributes
            gAtt = True
            gAtt@comment1    = "Data binned in 1-km altitude slices, bins " + \
                "with less than 10 valid points excluded"
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
                CAMPAIGN + "_" + FIELD + "_" + NAME(varid) + \
                "_" + d1 + "-" + d2 + ".nc"

            ;; Write variable
            write_nc_profile(fout, NAME(varid), DATA, gAtt)
            delete(DATA)

        end do

    end do

    leave_msg(diag_script, "", 1)

end
