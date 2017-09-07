;;#############################################################################
;; REFORMAT SCRIPT FOR AERONET OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    Che et al., Atmos. Chem. Phys., 15, 7619-7652,
;;    doi:10.5194/acp-15-7619-2015, 2015.
;;
;; Last access
;;    Data as published in the source.
;;
;; Download and processing instructions
;;
;; Caveats
;;    AOT(550) is calculated from AOT(440) using the Angstrom exponent
;;
;;       AOT(550) = AOT(440) * (550/440)^(-ang)
;;
;;    where ang is the Angstrom exponent 440/870.
;;
;; Modification history
;;    20150717-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "CARSNET"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period (climatology 2002-2013)
    YEAR1 = 2007
    YEAR2 = 2007

    ;; Selected variable (standard name)
    VAR = "od550aer"

    ;; Field
    FIELD = "T0A"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_CARSNET.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create time coordinate
    timec = create_timec(2007, 2007)
    datec = ut_calendar(timec, 1)

    ;; List of input files
    raw_file = "Che2015.asc"

    ;; Month names
    monthnames = (/"jan", "feb", "mar", "apr", "may", "jun",\
                   "jul", "aug", "sep", "oct", "nov", "dec"/)

    ;: Read header information
    head = readAsciiHead(INDIR + raw_file, 11)

    ;; Columns
    nc = str_fields_count(head(10), ",")
    cols = new(nc, string)
    do ii = 0, nc - 1
        cols(ii) = str_get_field(head(10), ii + 1, ",")
    end do

    ;; Read data
    data = readAsciiTable(INDIR + raw_file, 1, "string", 12)
    nrows = dimsizes(data(:, 0))

    ;; Station list
    st_name = str_get_field(data(:, 0), ind(cols.eq."site_name") + 1, ",")
    st_lon = tofloat(str_get_field(data(:, 0), ind(cols.eq."lon") + 1, ","))
    st_lat = tofloat(str_get_field(data(:, 0), ind(cols.eq."lat") + 1, ","))
    st_alt = toint(str_get_field(data(:, 0), ind(cols.eq."alt_m") + 1, ","))
    st_type = str_get_field(data(:, 0), ind(cols.eq."site_type") + 1, ",")

    ;; Loop over stations
    do ii = 0, dimsizes(st_name) - 1

        info_output("Processing station " + st_name(ii), verbosity, 1)

        ;; Define output array
        output = new(dimsizes(timec), float)
        output!0 = "time"
        output&time = timec
        output@_FillValue = FILL

        ;; Temporary arrays
        aod440 = new(12, float)
        aod440@_FillValue = -999.

        ang440870 = new(12, float)
        ang440870@_FillValue = -999.

        ;; Loop over months
        do mm = 0, 11
            idc = ind(cols.eq."odaer440_mean_" + monthnames(mm)) + 1
            aod440(mm) = tofloat(str_get_field(data(ii, 0), idc, ","))
            idc = ind(cols.eq."angstrom440-870_mean_" + monthnames(mm)) + 1
            ang440870(mm) = tofloat(str_get_field(data(ii, 0), idc, ","))
        end do

        output = aod440 * (550/440.) ^ (-1. * ang440870)
        delete(aod440)
        delete(ang440870)

        ;; Write data
        if (.not.all(ismissing(output))) then

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
            gAtt@type = st_type(ii)
            gAtt@altitude = st_alt(ii)
            gAtt@longitude = \
                where(st_lon(ii).lt.0, st_lon(ii) + 360., st_lon(ii))
            gAtt@latitude = st_lat(ii)
            gAtt@station = st_name(ii)
            gAtt@period = "2002-2013"
            gAtt@field = FIELD
            gAtt@tier = TIER
            gAtt@source = "Data as published in the source."
            gAtt@reference = "Che et al., Atmos. Chem. Phys., " + \
                "15, 7619-7652, doi:10.5194/acp-15-7619-2015, 2015."
            gAtt@title = OBSNAME + \
                " station data reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Outfile
            fout = OUTDIR + "OBS_" + OBSNAME + "_ground_" + \
                str_sub_str(st_name(ii), "_", "-") + "_" + FIELD + \
                "_" + VAR + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

            ;; Write variable
            write_nc(fout, VAR, output, gAtt)
            delete(gAtt)
            delete(output)

        end if

    end do

    leave_msg(diag_script, "", 1)

end
