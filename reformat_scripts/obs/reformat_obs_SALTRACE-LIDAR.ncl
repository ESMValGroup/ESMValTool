;;#############################################################################
;; REFORMAT SCRIPT FOR SALTRACE-LIDAR OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Gross et al., Atmos. Chem. Phys., 15, 11067-11080, 
;;    doi:10.5194/acp-15-11067-2015, 2015.
;;    Gross et al., Atmos. Chem. Phys., 16, 11535-11546,
;;    doi:10.5194/acp-16-11535-2016, 2016.
;;
;; Last access
;;    Data provided by S. Gross.
;;
;; Download and processing instructions
;;    Data has been untabified and converted with dos2unix.
;;    Commas have been replaced with dots.
;;    Underscores have been added to the header (line 1) of each file, to
;;    correctly identify column names (since these are different in each file).
;;
;; Caveats
;;    The data for the each flight are collected and vertically distributed
;;    in bins of 500 m using the provided (measured) pressure values.
;;    Statistics is then calculated for each bin.
;;    Data below 200 m are not used, since POLIS Lidar has full overlap only
;;    above 200m (Gross et al., ACP, 2016)
;;
;; Modification history
;;    20170823-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "SALTRACE"
    SUBNAME = "LIDAR"

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
    VARS = (/"concdust", "ec532dust"/)

    ;; Names in raw data
    NAMES = (/"dust_conc", "alpha_dust"/)    

    ;; Conversion factor to CMOR units
    CONV = (/1.e-3, 1.e-3/)

    ;; Location (see Gross et al., ACP, 2015, Section 2.1)
    MINLAT = 13.14
    MAXLAT = MINLAT
    MINLON = 360. - 59.62
    MAXLON = MINLON

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_SALTRACE-LIDAR.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create vertical bins
    alt_i = ispan(0, 22000, 500)
    alt_m = ispan(0, 21500, 500) + 250.

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
        infile = systemfunc("ls " + INDIR + "*.dat")

        do ff = 0, dimsizes(infile) - 1

            ;; Initialize
            DATA = OUTPUT

            ;; Read data
            head = readAsciiHead(infile(ff), "height")
            skip = dimsizes(head)
            data = readAsciiTable(infile(ff), 1, "string", skip)
            delete(skip)

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

            ;; Mark missing values with negative
            data = str_sub_str(data, "--", "-999")

            ;; Get campaign/flight name
            CAMPAIGN = str_sub_str(str_get_field( \
                  bname(infile(ff)), 1, "_"), ".nas", "")

            info_output("Processing " + CAMPAIGN, verbosity, 1)

            ;; Set time period
            date1 = str_get_field(CAMPAIGN, 1, "-")
            date2 = date1

            PERIOD = "20" + str_get_cols(date1, 0, 1) + "/" + \
                str_get_cols(date1, 2, 3) + "/" + \
                str_get_cols(date1, 4, 5) + "-" + \
                "20" + str_get_cols(date2, 0, 1) + "/" + \
                str_get_cols(date2, 2, 3) + "/" + \
                str_get_cols(date2, 4, 5)
            delete(date1)
            delete(date2)

            ;; Read measured altitude and convert to pressure
            height = tofloat(str_get_field(data(:, 0), 1, " "))
            tmp = stdatmus_z2tdp(height)
            press = tmp(2,:)
            press@_FillValue = -999
            press = where(press.le.0., press@_FillValue, press)
            delete(tmp)

            ;; Indexes for data below 200 m to be excluded
            idx200 = ind(height.ge.200.)
            delete(height)

            ;; Read data and fix missing values
            var = tofloat(str_get_field(data(:, 0), idcol, " "))
            var@_FillValue = -999
            var = where(var.le.0., var@_FillValue, var)

            ;; Convert units
            var = var * CONV(varID)
        
            ;; Exlcude points below 200 m
            press := press(idx200)
            var := var(idx200)

            ntot = dimsizes(var)
            ntot_check = 0

            ;; Loop over verical bins, assign data and calculate stats
            do pp = 0, NLEVS - 1

                pmax = press_i(pp)
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
            gAtt@comment2    = "Data binned in 500-m altitude slices, bins " + \
                "with less than 10 valid points excluded"
            gAtt@comment1    = "Data below 200 m are not used, since " + \
                "POLIS Lidar has full overlap only above 200m"
            gAtt@history     = "Created on " + systemfunc("date")
            gAtt@host        = systemfunc("echo $HOST")
            gAtt@user        = systemfunc("echo $USER")
            gAtt@lonmax      = MAXLON
            gAtt@lonmin      = MINLON
            gAtt@latmax      = MAXLAT
            gAtt@latmin      = MINLAT
            gAtt@info        = CAMPAIGN
            gAtt@period      = PERIOD
            gAtt@field       = FIELD
            gAtt@tier        = TIER
            gAtt@source      = "Data provided by Silke Gross (DLR)"
            gAtt@reference2  = "Gross et al., Atmos. Chem. Phys., 16, " + \
                "11535-11546, doi:10.5194/acp-16-11535-2016, 2016."
            gAtt@reference1  = "Gross et al., Atmos. Chem. Phys., 15, " + \
                "11067-11080, doi:10.5194/acp-15-11067-2015, 2015."
            gAtt@title       = OBSNAME + "-" + SUBNAME + \
                " vertical profiles reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Outfile
            d1 = str_get_field(str_get_field(PERIOD, 1, "-"), 1, "/") + \
                str_get_field(str_get_field(PERIOD, 1, "-"), 2, "/")
            d2 = str_get_field(str_get_field(PERIOD, 2, "-"), 1, "/") + \
                str_get_field(str_get_field(PERIOD, 2, "-"), 2, "/")
            fout = OUTDIR + "OBS_" + OBSNAME + "-" + SUBNAME + "_campaign_" + \
                CAMPAIGN + "_" + FIELD + "_" + VARS(varID) + \
                "_" + d1 + "-" + d2 + ".nc"

            ;; Write variable
            write_nc_profile(fout, VARS(varID), DATA, gAtt)
            delete(DATA)

        end do
    end do

    leave_msg(diag_script, "", 1)

end
