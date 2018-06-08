;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-AEROSOL SATELLITE DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    ftp://ftp.icare.univ-lille1.fr
;;
;; Last access
;;    20180110
;;
;; Download and processing instructions
;;    Download the data from:
;;       CCI-Aerosols/IASI_ULB_V8/L3_monthly/
;;    Consider only the complete years (2008-2016).
;;    This script uses the ALL files only, which combine the ascending (AN,
;;    night observations) and descending (DN, day observations) part of the
;;    orbit.
;;
;; Caveats
;;
;; Modification history
;;    20180110-A_righ_ma: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "IASI"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 2008
    YEAR2 = 2016

    ;; Selected variable (standard name)
    VARS = (/"od550dust"/)

    ;; Name in the raw data
    NAME = (/"D_AOD550_mean"/)

    ;; Field
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_IASI.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Set input files
    files = systemfunc("ls " + INDIR + "*ALL*.nc")
    nfiles = dimsizes(files)
    ntime = 12 * (YEAR2 - YEAR1 + 1)
    if (nfiles .ne. ntime) then
        error_msg("f", diag_script, "", "number of input files not " + \
                  "consistent with the specified time period: " + \
                  nfiles + " != " + ntime)
    end if

    f = addfiles(files, "r")
    ListSetType(f, "join")

    do vID = 0, dimsizes(VARS) - 1

        info_output("", verbosity, 1)
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ", \
                    verbosity, 1)

        OUTPUT = f[:]->$NAME(vID)$
        OUTPUT!0 = "time"
        OUTPUT&time = create_timec(YEAR1, YEAR2)

        ;; Set fill value
        OUTPUT = where(OUTPUT.eq.-999, OUTPUT@_FillValue, OUTPUT)

        ;; Format time coordinate
        ctime = time_attrib(OUTPUT&time)
        delete(OUTPUT&time)
        OUTPUT&time = ctime
        delete(ctime)

        ;; Format latitude coordinate
        OUTPUT!1 = "lat"
        if (isMonotonic(OUTPUT&lat) .eq. 0) then
            error_msg("f", diag_script, \
                      "", "non-monotonic latitude coordinate")
        end if
        if (isMonotonic(OUTPUT&lat) .eq. -1) then
            OUTPUT = OUTPUT(:, ::-1, :)
        end if
        clat = lat_attrib(OUTPUT&lat)
        delete(OUTPUT&lat)
        OUTPUT&lat = clat
        delete(clat)

        ;; Format longitude coordinate
        OUTPUT!2 = "lon"
        if (isMonotonic(OUTPUT&lon) .eq. 0) then
            error_msg("f", diag_script, \
                      "", "non-monotonic longitude coordinate")
        end if
        if (any(OUTPUT&lon.lt.0.)) then
            OUTPUT = lonFlip(OUTPUT)
        end if
        clon = lon_attrib(OUTPUT&lon)
        delete(OUTPUT&lon)
        OUTPUT&lon = clon
        delete(clon)

        ;; Set variable attributes
        tmp = var_attrib(OUTPUT, VARS(vID))
        delete(OUTPUT)
        OUTPUT = tmp
        delete(tmp)

        info_output("Range: " + min(OUTPUT) + "; " + \
                    max(OUTPUT), verbosity, 1)

        ;; Set global attributes
        gAtt = True
        gAtt@comment = "IASI dataset based on the Universite' libre de " + \
            "Bruxelles (ULB) algorithm (monthly mean L3 data)"
        gAtt@history = "Created on " + systemfunc("date")
        gAtt@host    = systemfunc("echo $HOST")
        gAtt@user    = systemfunc("echo $USER")
        gAtt@period  = YEAR1 + "-" + YEAR2
        gAtt@field   = FIELD
        gAtt@tier    = TIER
        gAtt@source  = "ftp://ftp.icare.univ-lille1.fr/" + \
            "CCI-Aerosols/IASI_ULB_V8/L3_monthly/"
        gAtt@reference = "Contact Lieven Clarisse (lclariss@ulb.ac.be)"
        gAtt@title = OBSNAME + " satellite data reformatted for the ESMValTool"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_sat_ULB-v8_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), OUTPUT, gAtt)
        delete(gAtt)
        delete(OUTPUT)

    end do

    leave_msg(diag_script, "", 1)

end
