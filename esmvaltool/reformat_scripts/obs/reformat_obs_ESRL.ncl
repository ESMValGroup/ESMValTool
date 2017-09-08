;;#############################################################################
;; REFORMAT SCRIPT FOR ESRL OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.esrl.noaa.gov/gmd/dv/data/index.php
;;
;; Last access
;;    20140115
;;
;; Download and processing instructions
;;    Select the following filters:
;;        Category: Greenhouse Gases.
;;        Parameter name: Carbon Dioxide.
;;        Type: Insitu.
;;        Frequency: Monthly Averages.
;;        Site: BRW, MLO, SMO, SPO.
;;    Convert the data to unix format using dos2unix.
;;    Get the list of stations at:
;;        http://www.esrl.noaa.gov/gmd/obop/
;;    This list must be manually formatted in a comma-separated table, given as
;;    FILE_STAT below.
;;    Check the period covered by the data and give it below as YEAR1 and
;;    YEAR2.
;;
;; Modification history
;;    20140122-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ESRL"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1973
    YEAR2 = 2012

    ;; Selected variable (standard name)
    VAR = "co2"

    ;; Field
    FIELD = "T0M"

    ;; Stations list
    FILE_STAT = INDIR + OBSNAME + "_stations.csv"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_ESRL.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create time coordinate
    timec = create_timec(YEAR1, YEAR2)
    datec = ut_calendar(timec, 1)

    ;; Read station info
    info_output("Reading station information", verbosity, 1)
    stat = readAsciiTable(FILE_STAT, 1, "string", 0)
    st_code = str_squeeze(str_get_field(stat(:, 0), 1, ","))
    st_name = str_squeeze(str_get_field(stat(:, 0), 2, ","))
    st_lat  = tofloat(str_get_field(stat(:, 0), 3, ","))
    st_lon  = tofloat(str_get_field(stat(:, 0), 4, ","))
    st_alt  = toint(str_get_field(stat(:, 0), 5, ","))
    nstats = dimsizes(st_code)
    delete(stat)

    ;; Loop over stations
    do stID = 0, nstats - 1

        infile = INDIR + VAR + "_" + str_lower(st_code(stID)) + \
            "_surface-insitu_1_ccgg_month.txt"

        if (isfilepresent_esmval(infile)) then

            head = readAsciiHead(infile, "# data_field")
            data = readAsciiTable(infile, 1, "string", dimsizes(head))
            delete(head)

            ;; Read site name
            tmp = str_get_field(data(:, 0), 1, " ")
            if (.not.all(tmp(0).eq.tmp)) then
                info_output("Inconsistent site name in " + files(ii), \
                            verbosity, 1)
            else
                site = tmp(0)
            end if
            delete(tmp)
            info_output("Processing site " + site, verbosity, 1)

            ;; Read time and extract full years only
            year = str_squeeze(str_get_field(data(:, 0), 2, " "))
            month = str_squeeze(str_get_field(data(:, 0), 3, " "))
            date = toint(year + month)
            var = tofloat(str_get_field(data(:, 0), 4, " "))
            delete(data)

            ;; Check time range
            if (min(toint(year)).lt.YEAR1 .or. max(toint(year)).gt.YEAR2) then
                error_msg("f", diag_script, "", "given time range (" + \
                          YEAR1 + "-" + YEAR2 + ")" + \
                          " not consistent with data (" + \
                          min(toint(year)) + "-" + max(toint(year)) + ")")
            end if
            delete(year)
            delete(month)

            ;; Create time coordinate
            timec = create_timec(YEAR1, YEAR2)
            datec = cd_calendar(timec, 1)

            ;; Create output variable
            output = new(dimsizes(timec), float)
            output!0 = "time"
            output&time = timec
            output@_FillValue = FILL

            do ii = 0, dimsizes(var) - 1
                output(ind(datec.eq.date(ii))) = (/var(ii)/)
            end do
            delete(var)
            delete(date)
            output = where(output.lt.0., output@_FillValue, output)

            ;; Set variable attributes (CMOR)
            tmp = var_attrib(output, VAR)
            delete(output)
            output = tmp
            delete(tmp)

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
            gAtt@altitude = st_alt(stID)
            gAtt@longitude = where(st_lon(stID) .lt. 0, \
                                   st_lon(stID) + 360., st_lon(stID))
            gAtt@latitude = st_lat(stID)
            gAtt@station = st_name(stID)
            gAtt@period = YEAR1 + "-" + YEAR2
            gAtt@field = FIELD
            gAtt@tier = TIER
            gAtt@source = "http://www.esrl.noaa.gov/gmd/dv/data/index.php"
            gAtt@reference = "Dr. Pieter Tans, NOAA/ESRL " + \
                "(www.esrl.noaa.gov/gmd/ccgg/trends/) and " + \
                "Dr. Ralph Keeling, Scripps Institution of Oceanography " + \
                "(scrippsco2.ucsd.edu/)"
            gAtt@title = OBSNAME + \
                " station data reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Outfile
            fout = OUTDIR + "OBS_" + OBSNAME + "_ground_" + \
                st_code(stID) + "_" + FIELD + "_" + VAR + "_" + \
                YEAR1 + "01-" + YEAR2 + "12.nc"

            ;; Write variable
            write_nc(fout, VAR, output, gAtt)
            delete(gAtt)

        end if

    end do

    leave_msg(diag_script, "", 1)

end
