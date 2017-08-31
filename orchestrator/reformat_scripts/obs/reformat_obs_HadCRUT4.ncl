;;###############################################################################
;; REFORMAT SCRIPT FOR THE HadCRUT4 TEMPERATURE (tas) OBSERVATIONAL DATA
;;
;; Note: In contrast to the HadCRUT3 reformat script which produces temperature
;;       anomalies (relative to the 1961-90 climatology), this script calculates
;;       absolute tas by adding the climatology ("absolute.nc") to the anomalies
;;       ("HadCRUT.4.4.0.0.median.nc").
;;###############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    https://crudata.uea.ac.uk/cru/data/temperature/
;;
;; Last access
;;    20160203
;;
;; Download and processing instructions
;;    Download the dataset "HadCRUT4" (median temperature anomalies) and
;;    the dataset "Absolute" (absolute temperatures for the base period
;;    1961-90 on a 5x5 grid (Jones et al., 1999)
;;
;; References
;;    Jones, P.D., New, M., Parker, D.E., Martin, S. and Rigor, I.G., 1999:
;;       Surface air temperature and its variations over the last 150 years.
;;       Reviews of Geophysics 37, 173-199, doi:10.1029/1999RG900002.
;;    Morice, C.P., Kennedy, J.J., Rayner, N.A. and Jones, P.D., 2012:
;;       Quantifying uncertainties in global and regional temperature change
;;       using an ensemble of observational estimates: the HadCRUT4 dataset.
;;       J. Geophys. Res., 117, D08101, doi:10.1029/2011JD017187.
;; Caveats
;;    none
;;
;; Modification history
;;    20160203-A_laue_ax: written.
;;
;;###############################################################################

begin

    ;; Source name
    OBSNAME = "HadCRUT4"

    ;; Tier
    TIER = 2

    ;; Input dir (raw data)
    INDIR = getenv("ESMValTool_RAWOBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Output dir (CMOR-ized data)
    OUTDIR = getenv("ESMValTool_OBSPATH") + "/Tier" + \
        TIER + "/" + OBSNAME + "/"

    ;; Period
    YEAR1 = 1850
    YEAR2 = 2015

    ;; Selected variable (standard name)
    VARS = (/"tas"/)

    ;; Field
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_HadCRUT4.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Read file
    fname  = INDIR + "HadCRUT.4.4.0.0.median.nc"
    fname2 = INDIR + "absolute.nc"
    f  = addfile(fname, "r")
    f2 = addfile(fname2, "r")

    ;; Read variables
    anomaly = f->temperature_anomaly
    climatology = tofloat(where(.not.ismissing(f2->tem), \
        f2->tem * f2->tem@scale_factor + 273.15, anomaly@_FillValue))
    climatology@_FillValue = anomaly@_FillValue
    climatology@missing_value = anomaly@_FillValue
    info_output("  Climatology range: " + min(climatology) + \
                " K to " + max(climatology) + " K", verbosity, 1)

    output = anomaly
;    printVarSummary(output)
;    printVarSummary(climatology)
    dims = dimsizes(output)
    do yr = 0, dims(0) / 12 - 1
        m1 = yr * 12
        m2 = m1 + 11
        output(m1:m2, :, :) = where(.not.ismissing(climatology), \
            output(m1:m2, :, :) + climatology, anomaly@_FillValue)
    end do

    info_output("", verbosity, 1)
    info_output("Processing " + VARS(0), verbosity, 1)
 
    ;; Set variable attributes
    tmp = var_attrib(output, VARS(0))
    output = tmp
    delete(tmp)

;    ;; Delete level coordinate (dimension 1)
;    tmp = rm_single_dims(output)
;    delete(output)
;    output = tmp
;    delete(tmp)

    ;; Format time coordinate
;    output_new = output(0:dimsizes(output&time)-3,:,:)
;    delete(output)
;    output = output_new
;    delete(output_new)
;    output!0 = "time"
    ctime = time_attrib(output&time)
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ;; Format latitude coordinate
    output!1 = "lat"
    if (isMonotonic(output&lat) .eq. 0) then
        error_msg("f", diag_script, \
                  "", "non-monotonic latitude coordinate")
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
        error_msg("f", diag_script, \
                  "", "non-monotonic longitude coordinate")
    end if
    clon = lon_attrib(output&lon)
    delete(output&lon)
    output&lon = clon
    delete(clon)

    info_output("  Range: " + min(output) + "; " + max(output), verbosity, 1)

    ;; Set global attributes
    gAtt = True
    gAtt@history       = "Created on " + systemfunc("date")
    gAtt@host          = systemfunc("echo $HOST")
    gAtt@user          = systemfunc("echo $USER")
    gAtt@period        = YEAR1 + "-" + YEAR2
    gAtt@field         = FIELD
    gAtt@tier          = TIER
    gAtt@source        = "https://crudata.uea.ac.uk/cru/data/temperature/"
    gAtt@reference     = "Morice et al., J. Geophys. Res. (2012); " + \
        "Jones et al., J. Reviews of Geophysics (1999)"
    gAtt@title         = "HadCRUT ground data " + \
        "reformatted for the ESMValTool"
    gAtt@conventions   = "CF/CMOR"

    ;; Outfile
    fout = OUTDIR + "OBS_HadCRUT4_ground_1_T2Ms_tas_" + \
        YEAR1 + "01-" + YEAR2 + "12.nc"

    ;; Write variable
    write_nc(fout, VARS(0), output, gAtt)
    delete(gAtt)
    delete(output)

    leave_msg(diag_script, "", 1)

end

