;;#############################################################################
;; REFORMAT SCRIPT FOR ERA-Interim PRECIPITATION DATA (FORECAST)
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    http://apps.ecmwf.int/datasets/data/interim-full-moda/
;;
;; Last access
;;    20150915
;;
;; Download and processing instructions
;;    Select "Synoptic Monthly Means". Download "Total Prepcipitation" for
;;    time 00, step 12 and time 12, step 12. These two datasets have to be
;;    added to get monthly means for precipitation, see:
;;
;;    http://www.ecmwf.int/en/how-can-i-get-monthly-means-daily-forecast-accumulations-eg-precipitation-era-interim-data-server
;;
;;    [...] you can calculate the monthly means of daily forecast accumulations
;;    (for step 12) by adding together the synoptic monthly means at time 00,
;;    step 12 and time 12, step 12.
;;
;;    A registration is required for downloading the data.
;;
;; Caveats
;;
;; Modification history
;;    20150915-A_laue_ax: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "ERA-Interim"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = "/export/pa_data01/ESMVal/obs/RAW/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = "/export/pa_data02/ESMVal/obs/Tier" + TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1979
    YEAR2 = 2014

    ;; Selected variable (standard name)
    VARS = (/"pr"/)

    ;; Name in the raw data
    NAME = (/"tp"/)

    ;; Units in the input  data and corresponding conversion factor
    ;; to CMOR standard units
    RAWUNITS = (/"m"/)

    CONVERSION = (/0.011574/)

    ;; Fields
    FLDARR = (/"T2Ms"/)

    ;; Initialize global variable
    FIELD = ""

end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = 1
    diag_script = "reformat_obs_ERA-Interim_precip.ncl"

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        ;; Set field
        FIELD = FLDARR(vID)

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ", \
                    verbosity, 1)

        fill_start = 0

        ;; Read file
        fname1 = INDIR + "ERA-Interim_" + NAME(vID) + "_time00_step12.nc"
        fname2 = INDIR + "ERA-Interim_" + NAME(vID) + "_time12_step12.nc"

        f1 = addfile(fname1, "r")
        f2 = addfile(fname2, "r")

        ;; Check time coordinate
        date = cd_calendar(f1->time, 0)
        if (any(date(:, 0).ne.yy)) then
            error_msg("f", diag_script, "", "time coordinate not " + \
                      "corresponding to selected year")
        end if
        delete(date)

        ;; Read variable
        xx1 = f1->$NAME(vID)$
        xx2 = f2->$NAME(vID)$

        ;; Add to missing values count
        fill_start = fill_start + num(ismissing(xx1))
        fill_start = fill_start + num(ismissing(xx2))

        ;; Unpack variable according to metadata information
        tmp1 = todouble(xx1)
        tmp1 = tmp1 * xx1@scale_factor + xx1@add_offset
        tmp2a = tofloat(tmp1)
        copy_VarMeta(xx1, tmp2a)
        delete(xx1)
        delete(tmp1)

        tmp1 = todouble(xx2)
        tmp1 = tmp1 * xx2@scale_factor + xx2@add_offset
        tmp2b = tofloat(tmp1)
        copy_VarMeta(xx2, tmp2b)
        delete(xx2)
        delete(tmp1)

        ; add time 00 and time 12

        output = tmp2a + tmp2b
        copy_VarMeta(tmp2a, output)

        output = where(output.ge.0, output, 0.0)

        ; -----------------------------------------------------------------

        ;; Convert units
        output = output * CONVERSION(vID)

        ;; Set variable attributes
        tmp = var_attrib(output, VARS(vID))
        delete(output)
        output = tmp
        delete(tmp)

        ;; Format time coordinate
        ctime = time_attrib(output&time)
        delete(output&time)
        output&time = ctime
        delete(ctime)

        ;; Format latitude coordinate
        output!1 = "lat"
        if (isMonotonic(output&lat) .eq. 0) then
            error_msg("f", diag_script, "", \
                      "non-monotonic latitude coordinate")
        end if
        if (isMonotonic(output&lat) .eq. -1) then
            output = output(:, ::-1, :)
        end if
        clat = lat_attrib(output&lat)
        delete(output&lat)
        output&lat = clat
        delete(clat)

        output!2 = "lon"
        if (isMonotonic(output&lon) .eq. 0) then
            error_msg("f", diag_script, "", \
                      "non-monotonic longitude coordinate")
        end if
        if (any(output&lon.lt.0.)) then
            output = lonFlip(output)
        end if
        clon = lon_attrib(output&lon)
        delete(output&lon)
        output&lon = clon
        delete(clon)

        info_output("  Range: " + min(output) + "; " + \
                    max(output), verbosity, 1)

        ;; Set global attributes
        gAtt = True
        gAtt@history       = "Created on " + systemfunc("date")
        gAtt@host          = systemfunc("echo $HOST")
        gAtt@user          = systemfunc("echo $USER")
        gAtt@period        = YEAR1 + "-" + YEAR2
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = "http://apps.ecmwf.int/datasets/data/" + \
            "interim_full_moda/"
        gAtt@reference     = "Dee, D. P. et al., Q. J. Roy. Meteor. Soc., " + \
            "137, 553-597, doi:10.1002/qj.828, 2011."
        gAtt@title = OBSNAME + \
            " reanalysis data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Check fill values
        fill_end = num(ismissing(output))
        if (fill_start.ne.fill_end) then
            error_msg("f", diag_script, "", "missing values lost " + \
                      "during conversion")
        end if
        delete(fill_start)
        delete(fill_end)

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
            "_" + VARS(vID) + ".nc"

        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)
        delete(gAtt)
        delete(output)

    end do

end
