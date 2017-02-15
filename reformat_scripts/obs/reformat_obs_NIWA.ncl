;;#############################################################################
;; REFORMAT SCRIPT FOR NIWA OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    http://www.bodekerscientific.com/data/total-column-ozone
;;
;; Last access
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20140528-A_gott_kl: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "NIWA"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1980
    YEAR2 = 2010

    ;; Selected variable (standard name)
    VARS = (/"toz"/)

    ;; Name in the raw data
    NAME = (/"TCO"/)

    ;; Units in the input  data and corresponding conversion factor
    ;; to CMOR standard units
    RAWUNITS = (/"Dobson Units (DU); 1DU=2.69×10^16 molecules/cm2"/)
    CONVERSION = (/1./)

    ;; Field
    FIELD = "T2Ms"

    ;; missing values in the input data
    FillValue = -999

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_NIWA.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    original_files = ""

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        info_output("Processing " + VARS(vID), verbosity, 1)

        fill_start = 0

        ;; Create timeseries
        do yy = YEAR1, YEAR2

            info_output("  " + yy, verbosity, 1)

            ;; Read file
            fname = "NIWAAssimOzoneV2.7_" + yy + "_Patched.nc"
            f = addfile(INDIR + fname, "r")

            ;; collect file names -> global attribute in outfile
            original_files = fname + ", " + original_files

            ;; Check time coordinate
            date = cd_calendar(f->time, 0)
            if (any(date(:, 0).ne.yy)) then
                error_msg("f", diag_script, "", \
                          "Time coordinate not corresponding to selected year")
            end if

            ;; Read variable
            xx = f->$NAME(vID)$
            if(.not.isatt(xx, "_FillValue")) then
                xx@_FillValue = FillValue
            end if

            ;; Add to missing values count
            fill_start = fill_start + num(ismissing(xx))

            ;; Add to time series
            if (.not.isdefined("output")) then
                output = xx
            else
                tmp3 = array_append_record(output, xx, 0)
                delete(output)
                output = tmp3
                delete(tmp3)
            end if
            delete(xx)

        end do

        ;; Convert units
        if (output@Units.ne.RAWUNITS(vID)) then
            error_msg("f", diag_script, "", \
                      "Unexpected units for " + NAME(vID))
        else
            output = output * CONVERSION(vID)
        end if

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
                      "Non-monotonic latitude coordinate")
        end if
        if (isMonotonic(output&lat) .eq. -1) then
            output = output(:, :, ::-1, :)
        end if
        clat = lat_attrib(output&lat)
        delete(output&lat)
        output&lat = clat
        delete(clat)

        ;; Format longitude coordinate
        output!2 = "lon"
        if (isMonotonic(output&lon) .eq. 0) then
            error_msg("f", diag_script, "", \
                      "Non-monotonic longitude coordinate")
        end if
        if (any(output&lon.lt.0.)) then
            output = lonFlip(output)
        end if
        clon = lon_attrib(output&lon)
        delete(output&lon)
        output&lon = clon
        delete(clon)

        info_output("  Range: " + min(output) + "; " + max(output), \
                    verbosity, 1)

        ;; Set global attributes
        gAtt = True
        gAtt@history       = "Created on " + systemfunc("date")
        gAtt@host          = systemfunc("echo $HOST")
        gAtt@user          = systemfunc("echo $USER")
        gAtt@period        = YEAR1 + "-" + YEAR2
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = "http://www.bodekerscientific.com/data/" + \
                             "total-column-ozone"
        gAtt@reference     = "Bodeker, G. E. et al., Atmos. Chem. Phys., " + \
            "5, 2603-2615, doi:10.5194/acp-5-2603-2005, 2005."
        gAtt@title = OBSNAME + \
            " reanalysis data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Check fill values
        fill_end = num(ismissing(output))
        if (fill_start.ne.fill_end) then
            error_msg("f", diag_script, "", \
                      "Missing values lost during conversion")
        end if
        delete(fill_start)
        delete(fill_end)

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)
        info_output("Wrote " + fout, verbosity, 1)
        delete(gAtt)
        delete(output)

    end do

    leave_msg(diag_script, "", 1)

end
