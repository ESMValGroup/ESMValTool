;;#############################################################################
;; REFORMAT SCRIPT FOR GLOBALVIEW (CO) OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    ftp://aftp.cmdl.noaa.gov/products/globalview/co/gv_co.tar.gz
;;
;; Last access
;;    20150720
;;
;; Download and processing instructions
;;
;; Caveats
;;    Hardcoded for *_ext* files
;;    Ignores weighting (*_wts* files)
;;
;; Modification history
;;    20140623-A_gott_kl: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "GLOBALVIEW"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Input file selectors:
    ;; http://www.esrl.noaa.gov/gmd/ccgg/globalview/gv_filenames.html
    ;; NOTE: Wildcard "*" will match patterns of any length,
    ;;       while each "?" represents exactly one character
    SITE = (/"???"/)  ; sampling site
    GROUP = (/""/)  ; data grouping
    LAB = (/"??"/)  ; Measurement Laboratory ID Number
    STRAT = (/"?"/)  ; Sampling Strategy:
                     ; C=quasi-continuous, D=discrete_flask, P=flask_package
    PLATF = (/"0"/)  ; Sampling Platform:
                     ;0=surface, 1=ship, 2=aircraft, 3=tower
    GAS = (/"co"/)  ; Gas Name

    ;; Selected variable (standard name): one per corresponding element of GAS
    VAR = (/"vmrco"/)

    ;; Units conversion
    UNITS = (/"mole mole-1"/)  ; target units
    CONVERSION = (/1.e-9/)  ; conversion factor to be multiplied to raw data

    ;; number of header lines
    n_head = 17

    ;; Create soft links as workaround for file name detection
    soft_links = True

end

load "./interface_scripts/messaging.ncl"
load "./interface_scripts/auxiliary.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_GLOBALVIEW.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Input files ############################################################
    ;; Loop over input file selectors to create lists of input files
    do i_GAS = 0, dimsizes(GAS) - 1
    do i_PLATF = 0, dimsizes(PLATF) - 1
    do i_STRAT = 0, dimsizes(STRAT) - 1
    do i_LAB = 0, dimsizes(LAB) - 1
    do i_GROUP = 0, dimsizes(GROUP) - 1
    do i_SITE = 0, dimsizes(SITE) - 1

        ;; extended files
        i_ext = systemfunc("ls " + INDIR + \
                           SITE(i_SITE) + GROUP(i_GROUP) + "_" + \
                           LAB(i_LAB) + STRAT(i_STRAT) + \
                           PLATF(i_PLATF) + "_ext." + GAS(i_GAS))
        if (.not.all(ismissing(i_ext))) then
            if (.not.isvar("ext_files")) then
                ext_files = i_ext
            else
                temp = array_append_record(i_ext, ext_files, 0)
                delete(ext_files)
                ext_files = temp
                delete(temp)
            end if
        end if

    end do  ; SITE
    end do  ; GROUP
    end do  ; LAB
    end do  ; STRAT
    end do  ; PLATF
    end do  ; GAS

    if(.not.isvar("ext_files")) then
        error_msg("f", diag_script, "", "no input files found in " + INDIR)
    end if

    ;; make sure each file occurs only once
    i = UNIQ(ext_files)
    temp = ext_files(i)
    delete(ext_files)
    ext_files = temp
    delete(temp)
    delete(i)

    ;; Output files ###########################################################
    ;; One output file per input file
    out_files = new(dimsizes(ext_files), string)
    varnames = new(dimsizes(ext_files), string)
    fields = new(dimsizes(ext_files), string)
    stations = new(dimsizes(ext_files), string)
    do i = 0, dimsizes(ext_files) - 1

        temp1 = str_split(ext_files(i), "/")  ; separate from path
        temp2 = temp1(dimsizes(temp1) - 1)  ; file name
        temp3 = str_split(temp2, ".")
        temp4 = temp3(0)  ; base name
        temp5 = temp3(1)  ; extension
        temp6 = str_split(temp4, "_")
        platform = str_get_cols(temp6(1), -1, -1)

        j = ind(GAS.eq.temp5)
        varname = VAR(j(0))

        if(platform.eq."2") then  ; aircraft
            field = "T0M"
            i_out = OUTDIR + "OBS_" + OBSNAME + "_insitu_" + \
                temp6(0) + "-" + temp6(1) + "_" + field + "_" + varname + ".nc"
        else  ; something at or near the ground
            field = "T0Ms"
            i_out = OUTDIR + "OBS_" + OBSNAME + "_ground_" + \
                temp6(0) + "-" + temp6(1) + "_" + field + "_" + varname + ".nc"
        end if

        fields(i) = field
        out_files(i) = i_out
        varnames(i) = varname
        stations(i) = temp6(0)
        delete([/field, temp1, temp2, temp3, temp4, temp5, temp6, varname, j/])

    end do

    ;; Process files #######################################################
    do ff = 0, dimsizes(ext_files) - 1
        info_output("Reading " + ext_files(ff), verbosity, 1)
        varname = varnames(ff)
        out_file = out_files(ff)
        station = stations(ff)
        FIELD = fields(ff)

        ;: Read header information
        head = readAsciiHead(ext_files(ff), n_head - 1)
        specifier = str_squeeze(head(0))
        sampling_strategy = str_squeeze(head(3))
        sampling_platform = str_squeeze(head(5))
        location = str_squeeze(head(6))
        raw_file_date = str_squeeze(head(12))
        longitude = str_get_field(str_squeeze(head(10)), 2, " ")
        latitude = str_get_field(str_squeeze(head(10)), 1, " ")
        altitude = str_get_field(str_squeeze(head(10)), 3, " ") + \
            " masl"  ; NOTE: str_get_field numbering starts with 1
        provider = str_squeeze(head(2))

        ;; Read table (ext): UTC        S(t)      REF(t)        diff
        n_data = toint(str_squeeze(str_get_field(head(13), 2, ":")))
        data = readAsciiTable(ext_files(ff), 4, "float", n_head - 1)
        data@_FillValue = -999.9990
        temp = dimsizes(data)
        if(n_data.ne.temp(0)) then
            error_msg("f", diag_script, "", "expected # of lines (" + \
                      tostring(n_data) + ") does not match # of read " + \
                      "lines (" + tostring(temp(0)) + ") in " + ext_files(ff))
        end if
        delete(temp)

        ;; Create NCL compatible time coordinate from input data
        ut = yyyyfrac_to_ut(data(:, 0), 0, 0)  ; see reformat_obs_func.ncl
        uti = toint(ut)
        YEAR = tostring(min(uti(:, 0)))
        units = "days since " + YEAR + "-01-01 00:00:0.0"
        option = 0
        option@calendar = "standard"
        time = cd_inv_calendar(uti(:, 0), uti(:, 1), uti(:, 2), \
                               uti(:, 3), uti(:, 4), ut(:, 5), \
                               units, option)
        delete([/units, YEAR, option/])

        ;; Find start and end year
        ;; Consider whole years only: one year in the raw data
        ;;     is represented by 48 "weekly" values
        n = 0  ; init
        i = 0
        do while(n.lt.48)
            YEAR1 = uti(i, 0)
            n = dimsizes(ind(uti(:, 0).eq.YEAR1))
            i = i + 1  ; counter
        end do
        n = 0  ; reset
        i = n_data - 1
        do while(n.lt.48)
            YEAR2 = uti(i, 0)
            n = dimsizes(ind(uti(:, 0).eq.YEAR2))
            i = i - 1  ; counter
        end do
        delete([/i, n/])
        suti = tostring(uti)
        info_output(" Raw data cover " + suti(0, 0) + "/" + suti(0, 1) + "/" +\
                    suti(0, 2) + " to " + suti(n_data - 1, 0) + "/" + \
                    suti(n_data - 1, 1) + "/" + suti(n_data - 1, 2) + \
                    ". Using " + tostring(YEAR1) + " to " + tostring(YEAR2) + \
                    " (12 months per year).", verbosity, 1)

        ;; Create variable with time coordinate
        ivar = ind(VAR.eq.varname)
        var = (data(:, 2) + data(:, 3)) * CONVERSION(ivar)
        var!0 = "time"
        var&time = time
        var@units = UNITS(ivar)

        ;; Extract selected period
        Y1Y2 = time_operations(var, YEAR1, YEAR2, "extract", "", 0)

        ;; Create time coordinate
        ;; Pass some parameters needed by function time_attrib as attributes
        timec = create_timec(YEAR1, YEAR2)
        timec@YEAR1 = YEAR1
        timec@YEAR2 = YEAR2
        timec@FIELD = FIELD
        datec = ut_calendar(timec, 1)

        ;; Calculate mymm
        mymm = time_operations(Y1Y2, YEAR1, YEAR2, "average", "mymm", False)
        mymm!0 = "time"
        mymm&time = timec
        delete(Y1Y2)

        ;; Write data ########################################################
        output = var_attrib(mymm, VAR)

        ;; Reformat time coordinate
        ctime = time_attrib(output&time)
        delete(output&time)
        output&time = ctime
        delete(ctime)

        ;; Set global attributes
        gAtt = True
        gAtt@history = "Created on " + systemfunc("date")
        gAtt@host = systemfunc("echo $HOST")
        gAtt@user = systemfunc("echo $USER")
        gAtt@altitude = altitude
        gAtt@longitude = longitude
        gAtt@latitude = latitude
        gAtt@location = location
        gAtt@specifier = specifier
        gAtt@sampling_strategy = sampling_strategy
        gAtt@sampling_platform = sampling_platform
        gAtt@provider = provider
        gAtt@raw_file_date = raw_file_date
        gAtt@raw_file = ext_files(ff)
        gAtt@period = tostring(YEAR1) + "-" + tostring(YEAR2)
        gAtt@field = FIELD
        gAtt@tier = TIER
        gAtt@source = "http://www.esrl.noaa.gov/gmd/ccgg/globalview"
        gAtt@reference = "GLOBALVIEW-CO2: Cooperative Atmospheric Data " + \
            "Integration Project - Carbon Dioxide, CD-ROM, NOAA ESRL, " + \
            "Boulder, Colorado, 2010."
        gAtt@title = OBSNAME + " station data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Write variable
        out_file = str_sub_str(out_file, ".nc", "_" + \
                               YEAR1 + "01-" + YEAR2 + "12.nc")
        write_nc(out_file, VAR, output, gAtt)
        info_output("  Wrote " + out_file, verbosity, 1)
        delete(gAtt)
        delete(output)
        delete(data)

        ;; Create soft links (workaround to deal with naming conventions)
        if(soft_links) then
            linked_file = str_sub_str(out_file, VAR, VAR + "_" + station)
            system("ln -s " + out_file + " " + linked_file)
        end if

    end do  ; loop over input files

    leave_msg(diag_script, "", 1)

end
