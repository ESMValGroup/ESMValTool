;;#############################################################################
;; REFORMAT SCRIPT FOR Emmons OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    Not available
;;
;; Last access
;;
;; Download and processing instructions
;;
;; Caveats
;;    Only certain units strings are recognized and converted to mol/mol
;;    (CMOR). If the units are not recognized the corresponding input file is
;;    skipped.
;;
;; Modification history
;;    20140127-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "Emmons"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/profiles/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Field
    FIELD = "T0M"

    ;; Selected variables
    VARIABLES = (/"hcooh", "acet", "ch2o", "ch3ooh", \
                  "acetone", "methanol", "c2h4", "c3h6", "hno3", \
                  "no2", "nox", "so2", "c2h2", "c2h6", "c3h8", \
                  "pan", "ch4", "co2", "co", "h2o", "no", \
                  "o3"/)

    ;; Corresponding standard names
    NAMES = (/"vmrhcooh", "vmrch3coch3", "vmrch2o", "vmrch3ooh", \
              "vmrch3coch3", "vmrch3oh", "vmrc2h4", "vmrc3h6", "vmrhno3", \
              "vmrno2", "vmrnox", "vmrso2", "vmrc2h2", "vmrc2h6", "vmrc3h8", \
              "vmrpan", "vmrch4", "vmrco2", "vmrco", "vmrh2o", "vmrno", \
              "tro3"/)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_Emmons.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Generate data array
    NLEVS = 16
    NCOLS = 10
    OUTPUT = new((/NLEVS, NCOLS/), float)
    OUTPUT!0 = "plev"
    OUTPUT!1 = "column"
    tmp = stdatmus_z2tdp(1000. * ispan(0, 15, 1))
    OUTPUT&plev = todouble(100 * tmp(2, :))
    OUTPUT&column = (/"N", "min", "max", "mean", \
                      "stddev", "5%", "25%", "median", "75%", "95%"/)
    OUTPUT&plev@long_name     = "pressure"
    OUTPUT&plev@positive      = "down"
    OUTPUT&plev@axis          = "Z"
    OUTPUT&plev@units         = "Pa"
    OUTPUT&plev@standard_name = "air_pressure"
    delete(tmp)

    ;; Loop over variables
    do vID = 0, dimsizes(VARIABLES) - 1

        info_output("Processing " + VARIABLES(vID) + " (" + NAMES(vID) + ")", \
                    verbosity, 1)

        ;; Loop over files
        flist = systemfunc("find " + INDIR + "*_" + VARIABLES(vID) + "_*.stat")

        do ff = 0, dimsizes(flist) - 1

            info_output("  " + flist(ff), verbosity, 1)

            fname = systemfunc("basename " + flist(ff))

            location = str_get_field(str_sub_str(fname, ".stat", ""), 4, "_")

            ;; Copy metadata
            output = OUTPUT

            ;; Get header info
            head = readAsciiHead(flist(ff), "# ialt")
            campaign = str_get_field(head(0), 2, " ") + "-" + \
                str_get_field(head(0), 3, " ")
            units = str_get_field(str_get_field(head(1), 2, "["), 1, "]")
            timerange = str_get_field(head(2), 2, " ")
            region = str_squeeze(str_sub_str(str_sub_str(\
                str_sub_str(head(3), "LAT=", ""), "LON=", ","), "#", ""))
            minlat = tofloat(str_get_field(region, 1, ","))
            maxlat = tofloat(str_get_field(region, 2, ","))
            minlon = tofloat(str_get_field(region, 3, ","))
            maxlon = tofloat(str_get_field(region, 4, ","))

            ;; Get columns info
            column = str_sub_str(head(6), "#", "")
            ncols = str_fields_count(column, " ")
            colnames = new(ncols, string)
            do ii = 0, ncols - 1
                colnames(ii) = str_squeeze(str_get_field(column, ii + 1, ","))
            end do

            ;; Check if file has data
            if (toint(systemfunc("wc -l <" + flist(ff))) .eq. \
                dimsizes(head)) then
                delete(head)
                continue
            end if

            ;; Read data
            data = readAsciiTable(flist(ff), 1, "string", dimsizes(head))
            delete(head)
            nrows = dimsizes(data(:, 0))

            ;; Read altitude and convert to integer to use as index
            alt = tofloat(str_get_field(data(:, 0), 1, " "))
            if (any(toint(alt).ne.alt)) then
                error_msg("f", "reformat_obs_Emmons.ncl", "", \
                          "non integer altitude value in " + flist(ff))
            else
                tmp = toint(alt)
                delete(alt)
                alt = tmp
                delete(tmp)
            end if

            ;; Loop over other columns and dump into array
            do jj = 1, ncols - 1
                if (colnames(jj).ne.output&column(jj - 1)) then
                    error_msg("f", "reformat_obs_Emmons.ncl", "", \
                              "unexpected column name in" + flist(ff))
                end if
                output(alt, jj - 1) = \
                    (/tofloat(str_get_field(data(:, 0), jj + 1, " "))/)
            end do

            ;; Set variable attributes (CMOR)
            tmp = var_attrib(output, NAMES(vID))
            delete(output)
            output = tmp
            delete(tmp)

            ;; Convert units (skip first column = N measurements)
            if (all(output@units.ne.(/"mole mole-1", "1e-9"/))) then
                error_msg("f", "reformat_obs_Emmons.ncl", "", \
                          "cannot convert to " + output@units)
            end if

            if (any(units.eq.(/"ppmv", "ppm", "PPMV", "PPM"/))) then
                if (output@units.eq."1e-9") then
                    output(:, 1:) = output(:, 1:) * 1.e3
                else
                    output(:, 1:) = output(:, 1:) / 1.e6
                end if
            else if (any(units.eq.(/"ppbv", "ppb", "PPBV", "PPB"/))) then
                if (output@units.ne."1e-9") then
                    output(:, 1:) = output(:, 1:) / 1.e9
                end if
            else if (any(units.eq.(/"pptv", "ppt", "PPTV", "PPT", \
                                    "pptv+pptv", "mol ppt"/))) then
                if (output@units.eq."1e-9") then
                    output(:, 1:) = output(:, 1:) / 1.e3
                else
                    output(:, 1:) = output(:, 1:) / 1.e12
                end if
            else
                info_output("   --> Unrecognized input units " + units, \
                            verbosity, 1)
                info_output("   --> Skipping", verbosity, 1)
                delete(output)
                delete(alt)
                delete(data)
                continue
            end if
            end if
            end if

            ;; Append global attributes
            gAtt = True
            gAtt@comment       = "Altitude conversion based on " + \
                "standard pressures"
            gAtt@history       = "Created on " + systemfunc("date")
            gAtt@host          = systemfunc("echo $HOST")
            gAtt@user          = systemfunc("echo $USER")
            gAtt@lonmax        = where(maxlon.lt.0., maxlon + 360., maxlon)
            gAtt@lonmin        = where(minlon.lt.0., minlon + 360., minlon)
            gAtt@latmax        = maxlat
            gAtt@latmin        = minlat
            gAtt@location      = location
            gAtt@campaign      = campaign
            gAtt@period        = timerange
            gAtt@field         = FIELD
            gAtt@tier          = TIER
            gAtt@source        = "Not available"
            gAtt@reference     = "Emmons, L. K. et al., J. Geophys. Res., " + \
                "105, 20497-20538, doi:10.1029/2000JD900232, 2000."
            gAtt@title         = OBSNAME + \
                " vertical profiles reformatted for the ESMValTool"
            gAtt@conventions   = "CF/CMOR"

            ;; Outfile
            d1 = str_get_field(timerange, 1, "-")
            year1 = str_get_field(d1, 1, "/") + str_get_field(d1, 2, "/")
            d2 = str_get_field(timerange, 2, "-")
            year2 = str_get_field(d1, 1, "/") + str_get_field(d2, 2, "/")
            fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_" + \
                str_sub_str(campaign, "_", "-") + "-" + \
                str_sub_str(location, "_", "-") + "_" + FIELD + \
                "_" + NAMES(vID) + "_" + year1 + "-" + year2 + ".nc"

            ;; Write variable
            write_nc_profile(fout, NAMES(vID), output, gAtt)

            delete(output)
            delete(alt)
            delete(data)

        end do
        delete(flist)

    end do

    leave_msg(diag_script, "", 1)

end
