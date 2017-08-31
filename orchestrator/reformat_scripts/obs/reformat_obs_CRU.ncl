;;#############################################################################
;; REFORMAT SCRIPT FOR CRU REANALYISIS DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    http://browse.ceda.ac.uk/browse/badc/cru/data/cru_ts/cru_ts_3.00/data
;;
;; Last access
;;    20151027
;;
;; Download and processing instructions
;;
;; Caveats
;;
;; Modification history
;;    20151027-A_righ_ma: fixed units conversion for pr (mm/month not mm/day).
;;    20151014-A_wenz_sa: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "CRU"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1901
    YEAR2 = 2006

    ;; Selected variable (standard name)
    VARS = (/"pr", "tas"/)

    ;; Name in the raw data
    NAME = (/"pre", "tmp"/)

    ;; Units in the input  data and corresponding conversion factor
    ;; to CMOR standard units
    RAWUNITS = (/ "mm", "degrees Celsius"/)

    ;; Fields
    FLDARR = (/ "T2Ms", "T2Ms"/)

    ;; Initialize global variable
    FIELD = ""

end

load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_CRU.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent_esmval(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        ;; Set field
        FIELD = FLDARR(vID)

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ", \
                    verbosity, 1)

        ;; Read file
        fname = INDIR + "cru_ts_3_00." + YEAR1 + "." + YEAR2 + "." + \
            NAME(vID) + ".nc"
        f = addfile(fname, "r")

        ;; Read variable
        xx = f->$NAME(vID)$

        ;; Unpack variable according to metadata information
        tmp1 = todouble(xx)
        tmp1 = tmp1 * xx@scale_factor
        output = tofloat(tmp1)
        copy_VarMeta(xx, output)
        delete(xx)
        delete(tmp1)
        output@_FillValue = -99.9
        output@missing_Value = output@_FillValue

        ;; Convert units
        if (output@units.ne.RAWUNITS(vID)) then
            error_msg("f", diag_script, "", \
                      "unexpected units for " + NAME(vID))
        else
            conversion = False

            ;; [degC] --> [K]
            if (NAME(vID).eq."tmp") then
                output = output + T0stp
                conversion = True
            end if

            ;; [mm/month] --> [kg/m2/s]
            if (NAME(vID).eq."pre") then
                date = cd_calendar(output&time, 0)
                dm = days_in_month(toint(date(:, 0)), toint(date(:, 1)))
                output = output / 24. / 3600. / conform(output, dm, 0)
                conversion = True
                delete(date)
                delete(dm)
            end if

            ;; Check conversion
            if (.not.conversion) then
                error_msg("f", diag_script, "", \
                          "no unit conversion applied")
            end if
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
            error_msg("fatal", diag_script, "",\
                      "non-monotonic latitude coordinate")
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
            error_msg("fatal", diag_script, "",\
                      "non-monotonic longitude coordinate")
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
        gAtt@source        = "http://browse.ceda.ac.uk/browse/badc/cru/" + \
            "data/cru_ts/cru_ts_3.00/data"
        gAtt@reference     = "Not available"
        gAtt@title = OBSNAME + \
            " reanalysis data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)
        delete(gAtt)
        delete(output)

    end do

    leave_msg(diag_script, "", 1)

end
