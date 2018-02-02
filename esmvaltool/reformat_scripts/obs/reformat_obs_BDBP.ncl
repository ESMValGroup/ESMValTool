;;#############################################################################
;; REFORMAT SCRIPT FOR ESACCI-AEROSOL SATELLITE DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: on request
;;
;; Source
;;    greg@bodekerscientific.com
;;
;; Last access
;;    20150119
;;
;; Download and processing instructions
;;    Email to greg@bodekerscientific.com
;;
;; Caveats
;;
;; Modification history
;;    20160801-A_wenz_sa: written based on reformat_obs_ESACCI-AEROSOL.ncl.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "BDBP"

    ;; Tier
    TIER = 3

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1979
    YEAR2 = 2007

    ;; Selected variable (standard name)
    VARS = (/"tro3prof"/)

    ;; Name in the raw data
    NAME = (/"O3"/)

    ;; Field
    FIELD = "T2Mz"

end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/logging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    diag_script = "reformat_obs_BDBP.ncl"

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    do vID = 0, dimsizes(VARS) - 1

        info_output("")
        info_output("Processing " + VARS(vID) + " (" + FIELD + ") ")

        time = create_timec(YEAR1, YEAR2)
        date = cd_calendar(time, 1)

        ;; Read file
        str = "find " + INDIR + " -type f -name 'CCMVal2_" \
              + "REF-B1_BSOzone-MR-PRS_Tier0_T2Mz_O3.nc'"
        fname = systemfunc(str)

        ;; No files found
        if (all(ismissing(fname))) then
            continue
        end if

        ;; Extract data
        f = addfile(fname(0), "r")
        xx = f->$NAME(vID)$

        ;; Assign to global array
        OUTPUT = xx

        delete(fname)

        ;; Set fill value
        ;OUTPUT = where(OUTPUT.eq.-9.9, OUTPUT@_FillValue, OUTPUT)

        ;; Set variable attributes
        tmp = var_attrib(OUTPUT, VARS(vID))
        delete(OUTPUT)
        OUTPUT = tmp
        delete(tmp)

        ;; Format time coordinate
        ctime = time_attrib(OUTPUT&time)
        delete(OUTPUT&time)
        OUTPUT&time = ctime
        delete(ctime)

        ;; Format longitudepressure level coordinate
        OUTPUT!1 = "plev"
        if (isMonotonic(OUTPUT&plev) .eq. 0) then
            error_msg("f", diag_script, \
                      "", "non-monotonic pressure levels coordinate")
        end if
        cplev = plev_attrib(OUTPUT&plev)
        delete(OUTPUT&plev)
        OUTPUT&plev = cplev
        delete(cplev)

        ;; Format latitude coordinate
        OUTPUT!2 = "lat"
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

        info_output("  Range: " + min(OUTPUT) + "; " + max(OUTPUT))

        ;; Set global attributes
        gAtt = True
        gAtt@comment = "Bodeker Scientific vertical ozone profile database" + \
                       "This is a preliminary version of the data base"
        gAtt@history = "Created on " + systemfunc("date")
        gAtt@host    = systemfunc("echo $HOST")
        gAtt@user    = systemfunc("echo $USER")
        gAtt@period  = YEAR1 + "-" + YEAR2
        gAtt@field   = FIELD
        gAtt@tier    = TIER
        gAtt@source  = "Greg Bodeker (greg@bodekerscientific.com)" + \
                       "and Birgit Hassler (birgit.hassler@googlemail.com)"
        gAtt@reference   = "NaN"
        gAtt@conventions = "CF/CMOR"

        ;; Outfile
        fout = OUTDIR + "OBS_" + OBSNAME + "_reanaly_1_" + FIELD + \
            "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

        ;; Write variable
        write_nc(fout, VARS(vID), OUTPUT, gAtt)
        delete(gAtt)
        delete(OUTPUT)

    end do

end
