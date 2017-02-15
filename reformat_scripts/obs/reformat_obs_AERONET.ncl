;;#############################################################################
;; REFORMAT SCRIPT FOR AERONET OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://aeronet.gsfc.nasa.gov/cgi-bin/combined_data_access_new
;;
;; Last access
;;    20150327
;;
;; Download and processing instructions
;;    Select Level 2.0 AOD, Version 2, Monthly Averages.
;;    Extract the tar.gz archive and specify the path as INDIR below.
;;
;; Caveats
;;    In the (numerous) cases where AOT 550 is not available in the data, this
;;    script calculates it from the AOT at other wavelenghts using the Angstrom
;;    exponent (de Meji & Lelieveld, Atm. Res. 99, 415-433, 2011, Eq. 2):
;;
;;       AOT(550) = AOT(440) * (550/440)^(-ang)
;;
;;    where ang is the Angstrom exponent 440/870.
;;
;; Modification history
;;    20140122-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "AERONET"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/AOT/LEV20/MONTHLY/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1992
    YEAR2 = 2015

    ;; Selected variable (standard name)
    VAR = "od550aer"

    ;; Field
    FIELD = "T0M"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_AERONET.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create time coordinate
    timec = create_timec(YEAR1, YEAR2)
    datec = ut_calendar(timec, 1)

    ;; List of input files
    raw_files = systemfunc("ls " + INDIR + "*.lev20")

    ;; Month names
    monthnames = (/"JAN", "FEB", "MAR", "APR", "MAY", "JUN",\
                                 "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"/)

    ;; Loop over input files
    do ff = 0, dimsizes(raw_files) - 1

        info_output(ff + " - " + systemfunc("basename " + raw_files(ff)), \
                    verbosity, 1)

        ;: Read header information
        head = readAsciiHead(raw_files(ff), "Month")
        st_name = str_squeeze(str_sub_str(str_get_field( \
            head(2), 1, ","), "Location=", ""))
        st_lon = tofloat(str_sub_str(str_get_field( \
            head(2), 2, ","), "long=", ""))
        st_lat = tofloat(str_sub_str(str_get_field( \
            head(2), 3, ","), "lat=", ""))
        st_alt = toint(str_sub_str(str_get_field( \
            head(2), 4, ","), "elev=", ""))
        ncols = str_fields_count(head(4), ",")

        ;; Read column names
        colnames = new(ncols, string)
        do ii = 1, ncols
            colnames(ii - 1) = str_squeeze(str_get_field(head(4), ii, ","))
        end do

        ;; Find required indexes
        month_ind = ind(colnames .eq. "Month")
        aot550_ind = ind(colnames .eq. "AOT_551 of Daily Avg")
        aot440_ind = ind(colnames .eq. "AOT_440 of Daily Avg")
        angstr_ind = ind(colnames .eq. "440-870 Angstrom of Daily Avg")

        ;; Read data
        data = readAsciiTable(raw_files(ff), 1, "string", dimsizes(head))
        month = str_get_field(data(:, 0), month_ind + 1, ",")
        date     = new(dimsizes(month), integer)
        aot550 = new(dimsizes(month), float)
        aot440 = new(dimsizes(month), float)
        angstr = new(dimsizes(month), float)

        aot550 = tofloat(str_sub_str(str_get_field(\
                         data(:, 0), aot550_ind + 1, ","), "N/A", -999))
        aot440 = tofloat(str_sub_str(str_get_field(\
                         data(:, 0), aot440_ind + 1, ","), "N/A", -999))
        angstr = tofloat(str_sub_str(str_get_field(\
                         data(:, 0), angstr_ind + 1, ","), "N/A", -999))

        aot550 = where(aot550 .eq. -999, aot550@_FillValue, aot550)
        aot440 = where(aot440 .eq. -999, aot440@_FillValue, aot440)
        angstr = where(angstr .eq. -999, angstr@_FillValue, angstr)

        ;; Define date
        do ii = 0, dimsizes(month) - 1
            date(ii) = 100 * toint(str_get_field(month(ii), 1, "-")) + \
                ind(str_get_field(month(ii), 2, "-") .eq. monthnames) + 1
        end do

        ;; Calculate alternative aot(550) from aot(440) and Angstrom exponent
        aot550_alt = aot440
        aot550_alt = aot440 * (550/440.) ^ (-1. * angstr)

        ;; Choose the most complete series
        if (num(ismissing(aot550)) .gt. num(ismissing(aot550_alt))) then
            aot550 = aot550_alt
        end if

        ;; Define output array
        output = new(dimsizes(timec), float)
        output!0 = "time"
        output&time = timec
        output@_FillValue = FILL
        do ii = 0, dimsizes(aot550) - 1
            didx = ind(date(ii) .eq. datec)
            output(didx) = (/aot550(ii)/)
        end do

        ;; Write data
        if (.not.all(ismissing(aot550))) then

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
            gAtt@altitude = st_alt
            gAtt@longitude = where(st_lon .lt. 0, st_lon + 360., st_lon)
            gAtt@latitude = st_lat
            gAtt@station = st_name
            gAtt@period = YEAR1 + "-" + YEAR2
            gAtt@field = FIELD
            gAtt@tier = TIER
            gAtt@source = "http://aeronet.gsfc.nasa.gov/cgi-bin/" + \
                "combined_data_access_new"
            gAtt@reference = "Holben, B. N. et al., Rem. Sens. Environ., " + \
                "16, 1-16, doi:10.1016/S0034-4257(98)00031-5, 1998."
            gAtt@title = OBSNAME + \
                " station data reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Outfile
            fout = OUTDIR + "OBS_" + OBSNAME + "_ground_" + \
                str_sub_str(st_name, "_", "-") + "_" + FIELD + \
                "_" + VAR + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

            ;; Write variable
            write_nc(fout, VAR, output, gAtt)
            delete(gAtt)
            delete(output)
        end if

        delete(head)
        delete(data)
        delete(st_name)
        delete(st_lat)
        delete(st_lon)
        delete(st_alt)
        delete(ncols)
        delete(colnames)
        delete(month_ind)
        delete(aot550_ind)
        delete(aot440_ind)
        delete(angstr_ind)
        delete(month)
        delete(date)
        delete(aot550)
        delete(aot550_alt)
        delete(aot440)
        delete(angstr)

    end do

    leave_msg(diag_script, "", 1)

end
