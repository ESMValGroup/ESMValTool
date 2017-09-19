;;#############################################################################
;; REFORMAT SCRIPT FOR CIRRUS OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier:
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Baumgardner et al., Environ. Res. Lett., 3, 025010,
;;    doi:10.1088/1748-9326/3/2/025010, 2008.
;;
;; Last access
;;    Original data from D. Baumgardner.
;;
;; Download and processing instructions
;;    Data has been converted with dos2unix and untabified.
;;
;; Caveats
;;    Data are binned into 1-km altitude slices.
;;    Data from all flights are collected in a sort of climatology.
;;    Levels with less than 10 valid points are not considered.
;;    Number concentrations are not included since information on the size cuts
;;    is not available.
;;
;; Modification history
;;    20140911-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "CIRRUS"

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

    ;; Selected variables
    VARIABLES = (/"MassNoBC", "MassBCtotal"/)

    ;; Corresponding standard names
    NAMES = (/"mmrbcfree", "mmrbc"/)

    ;; PERIOD
    PERIOD = "2006/11/24-2006/11/29"

end

load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_CIRRUS.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    files = systemfunc("ls " + INDIR + "merged*")

    ;; Check that data has consistent number of columns
    do ff = 0, dimsizes(files) - 1
        head = readAsciiHead(files(ff), 1)
        if (.not.isdefined("ncols")) then
            ncols = str_fields_count(head, " ")
        else
            if (ncols.ne.str_fields_count(head, " ")) then
                error_msg("f", diag_script, "", \
                          "inconsistent number of columns")
            end if
        end if
        delete(head)
    end do

    ;; Read and merge data in a single array
    nrows = 0
    do ff = 0, dimsizes(files) - 1
        ld = readAsciiTable(files(ff), 1, "string", 1)
        nrows = nrows + dimsizes(ld(:, 0))
        if (.not.isdefined("data")) then
            data = ld(:, 0)
        else
            tmp = array_append_record(data, ld(:, 0), 0)
            delete(data)
            data = tmp
            delete(tmp)
        end if
        delete(ld)
    end do
    if (dimsizes(data).ne.nrows) then
        error_msg("f", diag_script, "", "data_loss")
    end if

    ;; Read header (first file) and get column names
    head = readAsciiHead(files(0), 1)
    ncols = str_fields_count(head, " ")
    cols = new(ncols, "string")
    do ii = 1, ncols
        cols(ii - 1) = str_get_field(head, ii, " ")
    end do
    delete(ncols)
    delete(head)

    ;; Read latitude and longitude
    lat = tofloat(str_get_field(data, ind(cols.eq."LAT") + 1, " "))
    lon = tofloat(str_get_field(data, ind(cols.eq."LON") + 1, " "))
    lon = where(lon.lt.0., lon + 360., lon)

    ;; Read altitude
    alt = tofloat(str_get_field(data, ind(cols.eq."ALT") + 1, " "))
    alt = alt / 1000.

    ;; Read pressure and temperature
    press = tofloat(str_get_field(data, ind(cols.eq."PSTAT") + 1, " "))
    press = press * 100.  ; [hPa] --> [Pa]
    temp = tofloat(str_get_field(data, ind(cols.eq."SAT") + 1, " "))
    temp = temp + T0stp  ; [C] --> [K]

    ;; Loop over variables
    do vID = 0, dimsizes(VARIABLES) - 1

        ;; Read and filter data
        var = tofloat(str_get_field(data, \
                                    ind(cols.eq.VARIABLES(vID)) + 1, " "))

        ;; Filter data: exclude non-positive concentrations
        validx = ind(var.gt.0)
        var_f = var(validx)
        lat_f = lat(validx)
        lon_f = lon(validx)
        alt_f = alt(validx)
        press_f = press(validx)
        temp_f = temp(validx)
        delete(var)

        ;; Set region
        MINLAT = min(lat_f)
        MAXLAT = max(lat_f)
        MINLON = min(lon_f)
        MAXLON = max(lon_f)

        ;; Convert units
        var_f = var_f / 1.e12  ; [ng/m3] --> [kg/m3]
        var_f = var_f * Rgas * temp_f / \
            MWair / press_f * 1.e3  ; [kg/m3] --> [kg/kg]

        ;; Create altitude bins
        maxalt = toint(max(alt_f) + 1)
        ilev = tofloat(ispan(0, maxalt, 1))
        mlev = ilev(0:dimsizes(ilev) - 2) + 0.5

        ;; Define output array
        NLEVS = dimsizes(mlev)
        NCOLS = 6
        OUTPUT = new((/NLEVS, NCOLS/), float)
        OUTPUT!0 = "plev"
        OUTPUT!1 = "column"
        tmp = stdatmus_z2tdp(mlev * 1000.)
        OUTPUT&plev = todouble(100. * tmp(2, :))
        delete(tmp)
        OUTPUT&column = (/"N", "mean", "stddev", "25%", "median", "75%"/)
        OUTPUT&plev@long_name     = "pressure"
        OUTPUT&plev@positive      = "down"
        OUTPUT&plev@axis          = "Z"
        OUTPUT&plev@units         = "Pa"
        OUTPUT&plev@standard_name = "air_pressure"

        ;; Bin data
        ntot = 0
        do kk = 0, dimsizes(ilev) - 2
            kidx = ind(alt_f.ge.ilev(kk) .and. alt_f.lt.ilev(kk + 1))
            nidx = dimsizes(kidx)
            if (.not.all(ismissing(kidx))) then
                ntot = ntot + dimsizes(kidx)
                ;; Exclude levels with less than 10 points
                if (nidx.lt.10) then
                    delete(kidx)
                    continue
                end if
                slice = var_f(kidx)
                qsort(slice)
                OUTPUT(kk, 0) = (/dimsizes(kidx)/)
                OUTPUT(kk, 1) = (/avg(slice)/)
                OUTPUT(kk, 2) = (/stddev(slice)/)
                OUTPUT(kk, 3) = (/slice(round(0.25 * nidx, 3))/)
                OUTPUT(kk, 4) = (/dim_median(slice)/)
                OUTPUT(kk, 5) = (/slice(round(0.75 * nidx, 3))/)
                delete(slice)
            end if
            delete(kidx)
        end do

        if (dimsizes(validx).ne.ntot) then
            error_msg("f", diag_script, "", "data_loss")
        end if

        ;; Set variable attributes (CMOR)
        tmp = var_attrib(OUTPUT, NAMES(vID))
        delete(OUTPUT)
        OUTPUT = tmp
        delete(tmp)

        ;; Append global attributes
        gAtt = True
        gAtt@comment       = "Data binned in 1-km altitude slices, bins " + \
            "with less than 10 valid points excluded"
        gAtt@history       = "Created on " + systemfunc("date")
        gAtt@host          = systemfunc("echo $HOST")
        gAtt@user          = systemfunc("echo $USER")
        gAtt@lonmax        = where(MAXLON.lt.0., \
                                   MAXLON + 360., \
                                   MAXLON)
        gAtt@lonmin        = where(MINLON.lt.0., \
                                   MINLON + 360., \
                                   MINLON)
        gAtt@latmax        = MAXLAT
        gAtt@latmin        = MINLAT
        gAtt@info          = ""
        gAtt@period        = PERIOD
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = "Data provided by Darrel Baumgardner (UNAM)"
        gAtt@reference     = "Baumgardner et al., Environ. Res. Lett., 3, " + \
            "025010, doi:10.1088/1748-9326/3/2/025010, 2008."
        gAtt@title         = OBSNAME + \
            " vertical profiles reformatted for the ESMValTool"
        gAtt@conventions   = "CF/CMOR"

        ;; Outfile
        d1 = str_get_field(str_get_field(PERIOD, 1, "-"), 1, "/") + \
          str_get_field(str_get_field(PERIOD, 1, "-"), 2, "/")
        d2 = str_get_field(str_get_field(PERIOD, 2, "-"), 1, "/") + \
          str_get_field(str_get_field(PERIOD, 2, "-"), 2, "/")
        fout = OUTDIR + "OBS_" + OBSNAME + "_campaign_1_" + \
            FIELD + "_" + NAMES(vID) + "_" + d1 + "-" + d2 + ".nc"

        ;; Write variable
        write_nc_profile(fout, NAMES(vID), OUTPUT, gAtt)

        delete(validx)
        delete(var_f)
        delete(lat_f)
        delete(lon_f)
        delete(alt_f)
        delete(press_f)
        delete(temp_f)
        delete(OUTPUT)

    end do
    delete(lat)
    delete(lon)
    delete(alt)
    delete(press)
    delete(temp)

    leave_msg(diag_script, "", 1)

end
