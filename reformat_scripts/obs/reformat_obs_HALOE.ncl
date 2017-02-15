;;#############################################################################
;; REFORMAT SCRIPT FOR HALOE CLIMATOLOGY DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    Grooss, J.-U. and Russell III, J. M., Atmos. Chem. Phys., 5, 2797-2807,
;;    doi:10.5194/acp-5-2797-2005, 2005.
;;
;; Last access
;;    20140707
;;
;; Download and processing instructions
;;    Supplementary material of the given reference.
;;
;; Caveats
;;    The time period includes incomplete years. According to the given
;;    reference (Section 3): "Between October 1991 and August 2002, a total of
;;    78600 HALOE profiles of O3, CH4, H2O, HCl, and HF mixing ratio were then
;;    combined into bins of 5 degree equivalent latitude for 22 pressure levels
;;    between 316 and 0.1 hPa and 12 months".
;;    In order to create T3M fields a fake longitude coordinate is added.
;;
;; Modification history
;;    20140707-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "HALOE"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1991
    YEAR2 = 2002
    NYEARS = YEAR2 - YEAR1 + 1

    ;; Selected variable (standard name)
    VARS = (/"vmrh2o"/)

    ;; Name in the raw data
    NAME = (/"H2O"/)

    ;; Units in the input  data and corresponding conversion factor
    ;; to CMOR standard units
    RAWUNITS = (/"m^3/m^3"/)
    CONVERSION = (/1./)

    ;; Field
    FIELD = "T3M"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity =  stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_HALOE.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Create time coordinate
    TIME = new(12 * NYEARS, double)
    do yy = YEAR1, YEAR2
        do mm = 1, 12
            TIME(12 * (yy - YEAR1) + mm - 1) = \
                cd_inv_calendar(yy, mm, 15, 0, 0, 0, TUNITS, 0)
        end do
    end do

    ;; Open input file
    f = addfile(INDIR + "haloe_eqlat_monthly.nc", "r")

    ;; Loop over variables
    do vID = 0, dimsizes(VARS) - 1

        info_output("Processing " + VARS(vID), verbosity, 1)

        xx = f->$NAME(vID)$

        tmp = new((/12 * NYEARS, dimsizes(xx&press), dimsizes(xx&eqlat)/), \
                  float)
        tmp!0 = "time"
        tmp&time = TIME

        ;; Extract indexes for valid time range (see caveats)
        idx1 = ind(cd_calendar(TIME, -1).eq.199110)
        idx2 = ind(cd_calendar(TIME, -1).eq.200208)

        tmp(idx1:idx2, :, :) = (/xx/)

        ;; Create a fake longitude coordinate
        nlon = 50
        dims = array_append_record(dimsizes(tmp), nlon, 0)
        output = new(dims, float)
        output!0 = "time"
        output!1 = "plev"
        output!2 = "lat"
        output!3 = "lon"
        output&time = TIME
        output&plev = f->press
        output&lat = f->eqlat
        output&lon = lonGlobeF(nlon, "lon", "latitude", "degrees_East")
        do ii = 0, dimsizes(output&lon) - 1
            output(:, :, :, ii) = (/tmp/)
        end do
        copy_VarAtts(xx, output)
        delete(tmp)
        delete(xx)

        ;; Convert units
        if (output@units.ne.RAWUNITS(vID)) then
            error_msg("f", diag_script, "", \
                      "unexpected units for " + NAME(vID))
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

        ;; Format vertical coordinate
        output!1 = "plev"
        if (isMonotonic(output&plev) .eq. 0) then
            error_msg("f", diag_script, "", \
                      "non-monotonic vertical coordinate")
        end if
        if (isMonotonic(output&plev).eq.1) then
            output = output(:, ::-1, :, :)
        end if
        if (output&plev@units.eq."hPa") then
            output&plev = output&plev * 100
        else
            error_msg("f", diag_script, "", \
                      "unexpected units for plev " + output&plev@units)
        end if
        cplev = plev_attrib(output&plev)
        delete(output&plev)
        output&plev = cplev
        delete(cplev)

        ;; Format latitude coordinate
        output!2 = "lat"
        if (isMonotonic(output&lat) .eq. 0) then
            error_msg("f", diag_script, "", \
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
        output!3 = "lon"
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

        ;; Set global attributes
        gAtt = True
        gAtt@history       = "Created on " + systemfunc("date")
        gAtt@host          = systemfunc("echo $HOST")
        gAtt@user          = systemfunc("echo $USER")
        gAtt@period        = YEAR1 + "-" + YEAR2
        gAtt@field         = FIELD
        gAtt@tier          = TIER
        gAtt@source        = "Supplmentary material of the given reference."
        gAtt@reference     = "Grooss, J.-U. and Russell III, J. M., " + \
            "Atmos. Chem. Phys., 5, 2797-2807, " + \
            "doi:10.5194/acp-5-2797-2005, 2005"
        gAtt@title = OBSNAME + \
            " satellite data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_sat_1_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), output, gAtt)
        delete(gAtt)
        delete(output)

    end do

    leave_msg(diag_script, "", 1)

end
