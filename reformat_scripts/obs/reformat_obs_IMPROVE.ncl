;;#############################################################################
;; REFORMAT SCRIPT FOR IMPROVE OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://views.cira.colostate.edu/web/DataWizard/
;;
;; Last access
;;    20150904
;;
;; Download and processing instructions
;;    Select the following:
;;        Reports: RAW-DATA.
;;        Dataset: IMPROVE Aerosol.
;;        Sites: All.
;;        Parameters: NH4f, ECf, OCf, CHLf, MT, MF, RCFM, NO3f, SeaSaltf, NAf,
;;                    SO4f.
;;        Dates: 1978-2012 (or later).
;;        Fields: Dataset, Site, Date, Parameter, POC, Data Value,
;;                Unit, Status Flag.
;;        Filters: all valid data (V*).
;;        Options: text file, fixed width, standard.
;;    Convert the file format using dos2unix.
;;
;; Caveats
;;
;; Modification history
;;    20140127-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "IMPROVE"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = "/export/pa_data01/ESMVal/obs/RAW/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = "/export/pa_data02/ESMVal/obs/Tier" + TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1988
    YEAR2 = 2014

    ;; Selected variables (standard name)
    VARS   = (/"sconcnh4", "sconcbc", "sconcoa", "sconccl", "sconcpm10", \
               "sconcpm2p5", "sconcno3", "sconcss", "sconcna", "sconcso4"/)

    NVARS = dimsizes(VARS)

    ;; Name in the raw data
    COMPONENT = (/"NH4f", "ECf", "OCf", "CHLf", "MT", \
                  "MF", "NO3f", "SeaSaltf", "NAf", "SO4f"/)

    ;; Factors for unit conversion to the CMOR standard
    CONVERSION = (/ 1.e-9, 1.e-9, 1.4e-9, 1.e-9, 1.e-9, \
                    1.e-9, 1.e-9, 1.e-9, 1.e-9, 1.e-9/)

    ;; General comments (appended as global attribute to the output file)
    GATT_COMM = \
        (/"Mass of ammonium particles < 2.5 um in diameter", \  ; NH4f
          "From TOR carbon fractions (E1+E2+E3-OP)", \  ; ECf
          "From TOR carbon fractions (OC1+OC2+OC3+OC4+OP), " + \
          "multiplied by " + CONVERSION(2)*1.e9, \  ; OCf
          "Chloride Elemental Concentration FINE Size Fraction", \  ; CHLf
          "Gravimetric mass < 10 um in diameter", \  ; MT
          "Gravimetric fine mass", \  ; MF
          "Mass of nitrate particles < 2.5 um in diameter", \  ; NO3f
          "1.8 x [Chloride], or 1.8 x [Chlorine] if the chloride " + \
          "measurement is below detection " + \
          "limits, missing or invalid", \  ; SeaSaltf
          "Mass of sodium particles < 2.5 um in diameter", \  ; NAf
          "Mass of sulfate particles < 2.5 um in diameter"/)  ; SO4f

    ;; Field
    FIELD = "T0M"

    ;; Data filename
    INFILE = INDIR + "IMPROVE_19880302_20140530.dat"

end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = 1
    diag_script = "reformat_obs_IMPROVE.ncl"

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create time coordinate
    timec = create_timec(YEAR1, YEAR2)
    datec = ut_calendar(timec, 1)

    ;; Find first and last line with the station info in the input
    ls1 = toint(systemfunc("head -n 500 " + INFILE + \
                " | nl -ba | grep 'Site        ' | awk '{print $1}'")) + 1
    ls2 = toint(systemfunc("head -n 500 " + INFILE + \
                " | nl -ba | grep Parameters | awk '{print $1}'")) - 2

    dd = ls2 - ls1 + 1

    stations = systemfunc("head -n " + ls2 + " " + INFILE + " | tail -" + dd)

    st_name = str_squeeze(str_get_cols(stations, 0, 43))
    st_code = str_squeeze(str_get_cols(stations, 44, 50))
    st_lat  = tofloat(str_squeeze(str_get_cols(stations, 128, 137)))
    st_lon  = tofloat(str_squeeze(str_get_cols(stations, 138, 149)))
    st_alt  = toint(str_squeeze(str_get_cols(stations, 150, 160)))

    st_lat@_FillValue = -999.
    st_lon@_FillValue = -999.
    st_alt@_FillValue = -999

    st_lat = where(st_lat.eq.-999, st_lat@_FillValue, st_lat)
    st_lon = where(st_lon.eq.-999, st_lon@_FillValue, st_lon)
    st_alt = where(st_alt.eq.-999, st_alt@_FillValue, st_alt)

    nstats = dimsizes(st_name)

    ;; Remove double entries and invalid stations
    do stID = 0, nstats - 1
        if (ismissing(st_lat(stID)) .or. ismissing(st_lon(stID))) then
            st_code(stID) = -999
        end if
        if (num(st_code(stID) .eq. st_code) .ne. 1) then
            st_code(stID) = -999
        end if
    end do

    qq = ind(st_code .ne. -999)

    tmp = st_code(qq)
    delete(st_code)
    st_code = tmp
    delete(tmp)

    tmp = st_name(qq)
    delete(st_name)
    st_name = tmp
    delete(tmp)

    tmp = st_lat(qq)
    delete(st_lat)
    st_lat = tmp
    delete(tmp)

    tmp = st_lon(qq)
    delete(st_lon)
    st_lon = tmp
    delete(tmp)

    tmp = st_alt(qq)
    delete(st_alt)
    st_alt = tmp
    delete(tmp)

    nstats = dimsizes(st_name)

    info_output(nstats + " stations found", verbosity, 1)

    ;; Find first and last line with data info in the input
    ld1 = toint(systemfunc("head -n 12000 " + INFILE + " | nl -ba | grep " + \
                "Dataset | grep    SiteCode | awk '{print $1}'"))
    ld2 = toint(systemfunc("more " + INFILE + " | wc -l"))
    ldd = ld2 - ld1

    ;; Header
    head = systemfunc("head -n " + ld1 + " " + INFILE + " | tail -1")

    ;; Find column indexes for the date
    d1 = str_index_of_substr(head, "Date", 0)
    d2 = d1 + 15

    ;; Find column indices for site code
    s1 = str_index_of_substr(head, "SiteCode", 0)
    s2 = s1 + 13

    ;; Loop over variables
    do vID = 0, NVARS - 1

        info_output("Processing " + VARS(vID), verbosity, 1)

        ;; Find column indeces for the variable value
        c1 = str_index_of_substr(head, COMPONENT(vID) + ":Value", 0)
        c2 = c1 + 20

        ;; Find column indeces for the variable flag
        f1 = str_index_of_substr(head, COMPONENT(vID) + ":StatusFlag", 0)
        f2 = f1 + 20

        info_output("  FIELD  = " + str_get_cols(head, c1, c2), verbosity, 1)
        info_output("  FLAG   = " + str_get_cols(head, f1, f2), verbosity, 1)

        ;; Loop over stations
        do stID = 0, nstats - 1

            info_output("  Station " + st_code(stID), verbosity, 1)

            ;; Get data for the current stations
            data = systemfunc("more " + INFILE + " | tail -" + ldd + \
                              " | grep " + st_code(stID))

            ;; Variable values
            values = tofloat(str_squeeze(str_get_cols(data, c1, c2)))
            values@_FillValue = FILL

            ;; Flag
            flags = str_squeeze(str_get_cols(data, f1, f2))

            ;; Filter
            values = where(values.lt.0., values@_FillValue, values)  ; missing
            values = where(flags.ne."V0" .and. \
                           flags.ne."V1" .and. \
                           flags.ne."V2" .and. \
                           flags.ne."V3" .and. \
                           flags.ne."V4" .and. \
                           flags.ne."V5" .and. \
                           flags.ne."V6" .and. \
                           flags.ne."V7" .and. \
                           flags.ne."VM", \
                           values@_FillValue, values)  ; invalid flag

            if (all(ismissing(values))) then
                delete(data)
                delete(values)
                delete(flags)
                continue
            end if

            ;; Get date
            year  = toint(str_get_cols(str_squeeze( \
                          str_get_cols(data, d1, d2)), 0, 3))
            month = toint(str_get_cols(str_squeeze( \
                          str_get_cols(data, d1, d2)), 4, 5))

            ;; Define output array
            output = new(dimsizes(timec), float)
            output!0 = "time"
            output&time = timec
            output@_FillValue = FILL

            ;; Time loop
            do yy = YEAR1, YEAR2
                do mm = 1, 12
                    date = 100 * yy + mm
                    didx = ind(datec.eq.date)
                    vidx = ind(yy.eq.year .and. mm.eq.month)
                    if (.not.all(ismissing(vidx))) then
                        output(didx) = dim_avg(values(vidx))
                    end if
                    delete(date)
                    delete(didx)
                    delete(vidx)
                end do
            end do

            if (.not.all(ismissing(output))) then

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
                gAtt@comment = GATT_COMM(vID)
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
                gAtt@source = "http://views.cira.colostate.edu/web/DataWizard/"
                gAtt@reference = "Hand, J. L. et al., Spatial and " + \
                    "seasonal patterns and temporal variability of haze " + \
                    "and its constituents in the United States, " + \
                    "Interagency Monitoring of Protected Visual " + \
                    "Environments, 2011."
                gAtt@title = OBSNAME + \
                    " station data reformatted for the ESMValTool"
                gAtt@conventions = "CF/CMOR"

                ;; Outfile
                fout = OUTDIR + "OBS_" + OBSNAME + "_ground_" + \
                    st_code(stID) + "_" + FIELD + "_" + VARS(vID) + ".nc"

                ;; Write variable
                write_nc(fout, VARS(vID), output, gAtt)
                delete(gAtt)

            end if

            delete(data)
            delete(values)
            delete(flags)
            delete(year)
            delete(month)
            delete(output)

        end do  ; stations
    end do  ; variables

end
