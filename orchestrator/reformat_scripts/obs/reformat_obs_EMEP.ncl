;;#############################################################################
;; REFORMAT SCRIPT FOR EMEP OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://ebas.nilu.no
;;
;; Last access
;;    20150410
;;
;; Download and processing instructions
;;    Select the following (one at a time!):
;;      framework    matrix      component
;;        EMEP         aerosol     ammonium
;;        EMEP         aerosol     chloride
;;        EMEP         aerosol     nitrate
;;        EMEP         aerosol     sodium
;;        EMEP         aerosol     sulphate_total
;;        EMEP         pm10        pm10_mass
;;        EMEP         pm25        pm25_mass
;;    Select the time period 1970-2014 (or later).
;;    Convert the file format using dos2unix.
;;
;; Caveats
;;    Some data files have duplicated time coordinate in the data which must be
;;    removed before processing (duplicated data are spotted by the
;;    "process_EBAS_data" function called by this script).
;;
;; Modification history
;;    20150413-A_righ_ma: revised and improved.
;;    20140127-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "EMEP"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1970
    YEAR2 = 2014

    ;; Field
    FIELD = "T0M"

    ;; Selected variables (standard name)
    VARS = (/"sconcnh4", "sconccl", "sconcno3", "sconcpm10", \
             "sconcpm2p5", "sconcna", "sconcso4"/)
    NVARS = dimsizes(VARS)

    ;; Matrix
    MATRIX = (/"aerosol", "aerosol", "aerosol", "pm10", \
               "pm25", "aerosol", "aerosol"/)

    ;; Name in the raw data
    COMPONENT = (/"ammonium", "chloride", "nitrate", "pm10_mass", \
                  "pm25_mass", "sodium", "sulphate_total"/)

    ;; Expected units in the input data
    RAWUNITS = (/"ug N/m3", "ug/m3", "ug N/m3", "ug/m3", \
                 "ug/m3", "ug/m3", "ug S/m3"/)

    ;; Corresponding conversion factor to CMOR standard units [kg(species)/m3]
    CONVERSION = (/1.287955e-9, 1.e-9, 4.430737e-9, 1.e-9, \
                  1.e-9, 1.e-9, 2.99579e-9/)

    ;; Name in the column header
    COLUMN = (/"NH4+", "Cl-", "NO3-", "PM10", \
               "PM2.5", "Na+", "SO4--"/)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_EMEP.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Reading station list
    info_output("Extracting stations list", verbosity, 1)
    tmp = systemfunc("find " + INDIR + " -name '*.nas' " + \
                     "-exec basename {} \; | awk -F '.' '{print $1}'")
    stations = tmp(UNIQ(tmp))
    delete(tmp)
    sqsort(stations)
    info_output(dimsizes(stations) + " stations found", verbosity, 1)

    do st = 0, dimsizes(stations) - 1

        ST_CODE = stations(st)

        cnt = st + 1
        info_output("", verbosity, 1)
        info_output("Processing station " + ST_CODE + " (" + cnt + " of " + \
                    dimsizes(stations) + ")", verbosity, 1)

        ;; List of input files for this station
        FILES = systemfunc("ls " + INDIR + ST_CODE + "*.nas")
        if (all(ismissing(FILES))) then
            info_output("  --> No data for this station", verbosity, 1)
            delete(FILES)
            continue
        end if

        ;; Process station
        out = process_EBAS_data(VARS, RAWUNITS, MATRIX, COMPONENT, COLUMN, \
                                INDIR, ST_CODE, YEAR1, YEAR2)

        ;; Get station info from the first file in the list
        head = readAsciiHead(FILES(0), "starttime")
        do hh = 0, dimsizes(head) - 1
            cline = head(hh)

            ;; Station name
            if (isStrSubset(cline, "Station name:")) then
                st_name = str_squeeze(str_get_field(cline, 2, ":"))
            end if

            ;; Latitude
            if (isStrSubset(cline, "Station latitude:")) then
                st_lat = tofloat(str_get_field(cline, 2, ":"))
            end if

            ;; Longitude
            if (isStrSubset(cline, "Station longitude:")) then
                st_lon = tofloat(str_get_field(cline, 2, ":"))
            end if

            ;; Altitude
            if (isStrSubset(cline, "Station latitude:")) then
                st_alt = toint(str_sub_str( \
                    str_get_field(cline, 2, ":"), "m", ""))
            end if

            ;; Type
            if (isStrSubset(cline, "Station setting:")) then
                st_type = str_squeeze(str_get_field(cline, 2, ":"))
            end if

        end do
        delete(head)
        delete(FILES)

        ;; Write output
        do vID = 0, NVARS - 1
            output = out(:, vID)
            if (all(ismissing(output))) then
                delete(output)
                continue
            end if

            ;; Convert units
            output = output * CONVERSION(vID)

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
            gAtt@altitude = st_alt
            gAtt@longitude = where(st_lon .lt. 0, st_lon + 360., st_lon)
            gAtt@latitude = st_lat
            gAtt@station = st_name
            if (isdefined("st_type")) then
                gAtt@type = st_type
            end if
            gAtt@period = YEAR1 + "-" + YEAR2
            gAtt@field = FIELD
            gAtt@tier = TIER
            gAtt@source = "http://ebas.nilu.no/"
            gAtt@reference = "Hjellbrekke, A.-G. and Fjaeraa, A. M., " + \
                "Acidifying and eutrophying compounds and particulate " + \
                "matter, Norwegian Meteorological Institute, 2011."
            gAtt@title = OBSNAME + \
                " station data reformatted for the ESMValTool"
            gAtt@conventions = "CF/CMOR"

            ;; Outfile
            fout = OUTDIR + "OBS_" + OBSNAME + "_ground_" + ST_CODE + \
                "_" + FIELD + "_" + VARS(vID) + "_" + YEAR1 + "01-" + \
                YEAR2 + "12.nc"

            ;; Write variable
            write_nc(fout, VARS(vID), output, gAtt)
            delete(gAtt)
            delete(output)
        end do
        delete(out)

        info_output("", verbosity, 1)
        info_output("Completed!", verbosity, 1)

    end do

    leave_msg(diag_script, "", 1)

end
