;;#############################################################################
;; REFORMAT SCRIPT FOR CASTNET OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://java.epa.gov/castnet/clearsession.do
;;
;; Last access
;;    20130116
;;
;; Download and processing instructions
;;    Select the Prepackaged data and DRYCHEM.
;;    Extract the zip archive and specify the path as INDIR below.
;;    Convert the data to unix format using dos2unix.
;;    Get the list of stations at:
;;        http://java.epa.gov/castnet/epa_jsp/sites.jsp
;;    This list must be manually formatted in a comma-separated table, given as
;;    FILE_STAT below.
;;    Discontinued sites shall not be included in this table.
;;    Check the period covered by the data and give it below as YEAR1 and YEAR2
;;
;; Caveats
;;
;; Modification history
;;    20150409-A_righ_ma: added Na and Cl.
;;    20140123-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "CASTNET"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1987
    YEAR2 = 2012

    ;; Field
    FIELD = "T0M"

    ;; Selected variables (standard name)
    VARS   = (/"sconcso4", "sconcno3", "sconcnh4", "sconcna", "sconccl"/)
    NVARS = dimsizes(VARS)

    ;; Name in the raw data
    COMPONENT = (/"TSO4", "TNO3", "TNH4", "NA", "CL"/)

    ;; Unit conversion factor [ug/m3] --> [kg/m3]
    CONVERSION = (/ 1.e-9, 1.e-9, 1.e-9, 1.e-9, 1.e-9/)

    ;; Data filename
    INFILE = INDIR + "drychem.csv"

    ;; Stations list
    FILE_STAT = INDIR + OBSNAME + "_stations.csv"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_CASTNET.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create time coordinate
    timec = create_timec(YEAR1, YEAR2)
    datec = cd_calendar(timec, 1)

    ;; Read station list
    data = readAsciiTable(FILE_STAT, 1, "string", 0)
    st_code = str_squeeze(str_get_field(data(:, 0), 2, ","))
    st_name = str_squeeze(str_get_field(data(:, 0), 3, ","))
    st_lat  = tofloat(str_get_field(data(:, 0), 4, ","))
    st_lon  = tofloat(str_get_field(data(:, 0), 5, ","))
    st_alt  = toint(str_get_field(data(:, 0), 6, ","))
    nstats = dimsizes(st_code)
    delete(data)

    info_output(nstats + " stations found", verbosity, 1)

    head = readAsciiHead(INFILE, 1)

    ; Find column names
    colnames = new(str_fields_count(head, ","), string)
    do ii = 0, dimsizes(colnames) - 1
        colnames(ii) = str_get_field(head, ii + 1, ",")
    end do

    ; Loop over variables
    do vID = 0, NVARS - 1

        info_output("Processing " + VARS(vID), verbosity, 1)

        ; Find column number for the variable
        col = ind(colnames .eq. COMPONENT(vID)) + 1

        ; Flag is the next column
        flag = col + 1

        info_output("  FIELD  = " + colnames(col - 1), verbosity, 1)
        info_output("  FLAG   = " + colnames(flag - 1), verbosity, 1)

        ; Loop over stations
        do stID = 0, nstats - 1

            info_output("  Station " + st_code(stID), verbosity, 1)

            ; Define output array
            output = new(dimsizes(timec), float)
            output!0 = "time"
            output&time = timec
            output@_FillValue = FILL

            ; Define number of actual observation days (for weighted average)
            ndays = output

            ; Sample type
            sample = str_sub_str(systemfunc("more " + INFILE + " | grep " + \
                                 st_code(stID) + " | awk -F " + \
                                 str_get_dq + "," + str_get_dq + \
                                 " '{print $3}'"), str_get_dq, "")

            ; Time interval of each input
            date1 = str_sub_str(systemfunc("more " + INFILE + " | grep " + \
                                st_code(stID) + " | awk -F " + \
                                str_get_dq + "," + str_get_dq + \
                                " '{print $4}'"), str_get_dq, "")

            date2 = str_sub_str(systemfunc("more " + INFILE + " | grep " + \
                                st_code(stID) + " | awk -F " + \
                                str_get_dq + "," + str_get_dq + \
                                " '{print $5}'"), str_get_dq, "")

            ; Variable values
            tmp = systemfunc("more " + INFILE + " | grep " + st_code(stID) + \
                             " | awk -F " + str_get_dq + "," + str_get_dq + \
                             " '{print $" + col + "}'")

            ;; Search for duplicates
            ;; (valid types: A,B,C,H; multiple types valid if B,D,N)
            values = new(dimsizes(tmp), float)
            values@_FillValue = FILL
            if (.not.all(ismissing(tmp))) then
                do ii = 0, dimsizes(tmp) - 1
                    ; multiple values, B,D,N: choose B (combined N+D)
                    if (num(date1(ii).eq.date1).gt.1) then
                        if (sample(ii).eq."B") then
                            if (tmp(ii).ne."") then
                                values(ii) = tofloat(tmp(ii))
                            end if
                        else
                            if (sample(ii).ne."D" .and. sample(ii).ne."N") then
                                error_msg("w", diag_script, "", \
                                          "multiple sample " + sample(ii))
                            end if
                        end if
                    ; single value
                    else
                        if (any(sample(ii) .eq. (/"A", "B", "C", "H"/))) then
                            if (tmp(ii) .ne. "") then
                                values(ii) = tofloat(tmp(ii))
                            end if
                        else
                            error_msg("w", diag_script, "", \
                                      "skipping single sample type " + \
                                      sample(ii) + "  " + date1(ii))
                        end if
                    end if
                end do
            end if
            delete(tmp)

            ; Variable flags
            flags = systemfunc("more " + INFILE + " | grep " + \
                               st_code(stID) + " | awk -F " + \
                               str_get_dq + "," + str_get_dq + \
                               " '{print $" + flag + "}'")

            ; Filter
            values = where(values .lt. 0., values@_FillValue, values)
            values = where(flags.eq."I" .or. flags .eq. "M" .or. \
                           flags .eq. "N", values@_FillValue, values)

            if (all(ismissing(values))) then
                info_output("No valid data for this station", verbosity, 1)
                delete(date1)
                delete(date2)
                delete(values)
                delete(flags)
                delete(sample)
                continue
            end if

            ;; Calculate monthly means
            ;; Loop over all values and find position in the time array
            do ii = 0, dimsizes(values) - 1

                if (.not.ismissing(values(ii))) then

                    ; Date 1
                    yy1 = toint(str_get_cols(date1(ii), 0, 3))
                    mm1 = toint(str_get_cols(date1(ii), 5, 6))
                    dd1 = toint(str_get_cols(date1(ii), 8, 9))
                    hh1 = toint(str_get_cols(date1(ii), 11, 12))
                    pp1 = toint(str_get_cols(date1(ii), 14, 15))
                    ss1 = toint(str_get_cols(date1(ii), 17, 18))
                    tt1 = cd_inv_calendar(yy1, mm1, dd1, \
                                          hh1, pp1, ss1, TUNITS, 0)

                    ; Date 2
                    yy2 = toint(str_get_cols(date2(ii), 0, 3))
                    mm2 = toint(str_get_cols(date2(ii), 5, 6))
                    dd2 = toint(str_get_cols(date2(ii), 8, 9))
                    hh2 = toint(str_get_cols(date2(ii), 11, 12))
                    pp2 = toint(str_get_cols(date2(ii), 14, 15))
                    ss2 = toint(str_get_cols(date2(ii), 17, 18))
                    tt2 = cd_inv_calendar(yy2, mm2, dd2, \
                                          hh2, pp2, ss2, TUNITS, 0)

                    ; Positions in the time/date array
                    idx_1 = ind(datec .eq. (100 * yy1 + mm1))
                    idx_2 = ind(datec .eq. (100 * yy2 + mm2))

                    ; Loop over selected time intervals
                    do idx = idx_1, idx_2
                        cy = toint(datec(idx) / 100)
                        cm = toint(datec(idx) - 100 * cy)
                        pos_sta = cd_inv_calendar(cy, cm, 1, \
                                                  0, 0, 0, TUNITS, 0)
                        if (cm.eq.12) then
                            pos_end = cd_inv_calendar(cy + 1, 1, 1, \
                                                      0, 0, 0, TUNITS, 0)
                        else
                            pos_end = cd_inv_calendar(cy, cm + 1, 1, \
                                                      0, 0, 0, TUNITS, 0)
                        end if
                        pos_1 = max((/pos_sta, tt1/))
                        pos_2 = min((/pos_end, tt2/))
                        nd = tofloat(pos_2 - pos_1)
                        if (ismissing(output(idx))) then
                            output(idx) = values(ii) * nd
                            ndays(idx)  = nd
                        else
                            output(idx) = output(idx) + values(ii) * nd
                            ndays(idx) = ndays(idx) + nd
                        end if
                        delete(cy)
                        delete(cm)
                        delete(pos_sta)
                        delete(pos_end)
                        delete(pos_1)
                        delete(pos_2)
                        delete(nd)
                    end do
                end if
            end do

            ; Check number of days
            do ii = 0, dimsizes(datec) - 1
                cy = toint(datec(ii) / 100)
                cm = toint(datec(ii) - 100 * cy)
                maxnd = days_in_month(cy, cm)
                if (.not.ismissing(ndays(ii))) then
                    if ((ndays(ii) - maxnd) .gt. 1.e-3) then
                        error_msg("f", diag_script, "", \
                                  "total number of days > " + \
                                  "days in month for " + datec(ii) + " (" + \
                                  ndays(ii) + " > " + maxnd + ")")
                    end if
                else
                    if (.not.ismissing(output(ii))) then
                        error_msg("f", diag_script, "", \
                                  "missing ndays for non-missing output")
                    end if
                end if
            end do

            ; Weighted average
            output = output / ndays

            ; Unit conversion
            output = output * CONVERSION(vID)  ; unit conversion

            ;; Write output
            if (.not.all(ismissing(output))) then

                ;; Set variable attributes (CMOR)
                tmp = var_attrib(output, VARS(vID))
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
                gAtt@source = "http://epa.gov/castnet/javaweb/index.html"
                gAtt@reference = "Edgerton, E. et al., National dry " + \ \
                    "deposition network: Second annual progress report, " + \
                    "Environmental Protection Agency, 1990."
                gAtt@title = OBSNAME + \
                    " station data reformatted for the ESMValTool"
                gAtt@conventions = "CF/CMOR"

                ;; Outfile
                fout = OUTDIR + "OBS_" + OBSNAME + "_ground_" + \
                    st_code(stID) + "_" + FIELD + "_" + VARS(vID) + \
                    "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

                ;; Write variable
                write_nc(fout, VARS(vID), output, gAtt)
                delete(gAtt)
            end if

            delete(date1)
            delete(date2)
            delete(values)
            delete(flags)
            delete(sample)
            delete(output)
            delete(ndays)

        end do  ; stations

    end do  ; variables

    leave_msg(diag_script, "", 1)

end
