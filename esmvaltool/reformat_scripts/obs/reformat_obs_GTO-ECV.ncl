;;#############################################################################
;; REFORMAT SCRIPT FOR GTO-ECV OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    Loyola et al., Int. J. Remote Sens., 30, 4295-4318,
;;    doi:10.1080/01431160902825016, 2009.
;;
;; Last access
;;    Data provided by D. Loyola.
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20140707-A_gott_kl: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "GTO-ECV"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period: RAW data cover 199506 - 201012
    YEAR1 = 1996  ; yyyy
    YEAR2 = 2010  ; yyyy

    ;; Selected variables
    rawVARS = (/"GTO"/)  ; raw name
    VARS = (/"toz"/)  ; standard name

    ;; Units in the input data and corresponding conversion factor
    ;; to CMOR standard units
    RAWUNITS = (/"Dobson Units"/)
    CONVERSION = (/1./)

    ;; Field
    FIELD = "T2Ms"

    ;; Gridtype
    GRIDTYPE = "Scanning Nadir"

    ;; missing values
    FillValue = -999  ; just something impossible

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_GTO-ECV.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create coordinates
    temp = create_timec(YEAR1, YEAR2)
    timec = time_attrib(temp)
    delete(temp)
    datec = ut_calendar(timec, 0)
    temp = fspan(89.5, -89.5, 180)
    lat = lat_attrib(temp)
    delete(temp)
    temp = fspan(0.5, 359.5, 360)
    lon = lon_attrib(temp)
    delete(temp)

    ;; Loop over variables
    do ivar = 0, dimsizes(VARS) - 1

        ;; Target fields
        data_var = new((/dimsizes(timec), dimsizes(lat), \
                        dimsizes(lon)/), float)
        data_stddev = new((/dimsizes(timec), dimsizes(lat), \
                           dimsizes(lon)/), float)
        original_files = new(dimsizes(timec), string)

        ;; Time loop = loop over (monthly) raw files
        itime = 0
        do iyear = YEAR1, YEAR2, 1
            do imonth = 1, 12, 1

                ;; Raw file name
                year = tostring(iyear)
                if (imonth .lt. 10)then
                    month = "0" + tostring(imonth)
                else
                    month = tostring(imonth)
                end if
                INFILE = rawVARS(ivar) + "-ECV_V1_" + year + \
                    "_" + month + ".dat"
                info_output("Processing " + INDIR + INFILE, verbosity, 1)

                ;; collect file names -> global attribute in outfile
                original_files(itime) = INFILE

                ;; Read raw data
                ;; Data structure of raw files: lat, lon, mean ozone, stddev
                data  = readAsciiTable(INDIR + INFILE, 4, "float", 0)
                    ;; dimsizes(data) = (/number of lines, 4/)

                do ilon = 0, dimsizes(lon) - 1
                    index1 = ind(data(:, 1).eq.lon(ilon))
                    data_tmp = data(index1, :)

                    do ilat = 0, dimsizes(lat) - 1
                        index2 = ind(data_tmp(:, 0).eq.lat(ilat))

                        if (.not.ismissing(index2)) then
                            data_var(itime, ilat, ilon) = data_tmp(index2, 2)
                            data_stddev(itime, ilat, ilon) = \
                                data_tmp(index2, 3)
                        else
                            data_var(itime, ilat, ilon) = FillValue
                            data_stddev(itime, ilat, ilon) = FillValue
                        end if

                        delete(index2)
                    end do

                    delete(data_tmp)
                    delete(index1)
                end do

                delete(data)
                itime = itime + 1
            end do
        end do
        delete(itime)

        ;; Scale
        data_var = data_var * CONVERSION(ivar)
        data_var@_FillValue = FillValue
        data_stddev = data_stddev * CONVERSION(ivar)
        data_stddev@_FillValue = FillValue

        ;; Coordinates
        data_var!0 = "time"
        data_var!1 = "lat"
        data_var!2 = "lon"
        data_var&lon = lon
        data_var&lat = lat
        data_var&time = timec
        copy_VarCoords(data_var, data_stddev)

        ;; Set variable attributes
        output = var_attrib(data_var, VARS(ivar))
        delete(data_var)
        output_stddev = var_attrib(data_stddev, VARS(ivar))
        delete(data_stddev)

        ;; Set global attributes
        gAtt = True
        gAtt@original_files = str_concat(original_files + ", ")
        gAtt@gridtype      = GRIDTYPE
        gAtt@history       = "Created on " + systemfunc("date")
        gAtt@host          = systemfunc("echo $HOST")
        gAtt@user          = systemfunc("echo $USER")
        gAtt@period        = tostring(YEAR1) + "-" + tostring(YEAR2)
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = "Data provided by Diego Loyola (DLR)"
        gAtt@reference     = "Loyola et al., Int. J. Remote Sens., 30, " + \
            "4295-4318, doi:10.1080/01431160902825016, 2009."
        gAtt@title = OBSNAME + " satellite data reformatted for the " + \
            "ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Write variable
        fout = OUTDIR + "OBS_" + OBSNAME + "_sat_1_" + FIELD + \
            "_" + VARS(ivar) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"
        write_nc(fout, VARS(ivar), output, gAtt)
        delete(output)
        info_output("Wrote " + fout, verbosity, 1)

        ;; Write stddev
;;         fout = OUTDIR + "OBS_" + OBSNAME + "_sat_1_" + FIELD + \
;;             "_" + VARS(ivar) + "_stdddev.nc"
;;         write_nc(fout, VARS(ivar) + "_stddev", output_stddev, gAtt)
;;         delete(output_stddev)
;;         info_output("Wrote " + fout, verbosity, 1)

        delete(gAtt)
    end do  ; ivar

    leave_msg(diag_script, "", 1)

end
