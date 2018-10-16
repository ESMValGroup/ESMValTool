;;#############################################################################
;; REFORMAT SCRIPT FOR THE HadCRUT3v TEMPERATURE (tas) OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.metoffice.gov.uk/hadobs/hadcrut3/data/download.html
;;
;; Last access
;;    20140311
;;
;; Download and processing instructions
;;    The original data file is compressed to netCDF4 ZIP format.
;;    To decompress use "cdo -f nc -copy <original-file> <unzipped-file>"
;;
;; Caveats
;;    No caveats, but note that HadCRUT3v variance adjusted tas data is
;;    provided as anomalies with respect to the period 1958-2001.
;;
;; Modification history
;;    20150330-A_righ_ma: updated paths and global attributes.
;;    20140311-A_senf_da: written.
;;
;;#############################################################################

begin

    ;; Source name
    OBSNAME = "HadCRUT"

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
    YEAR2 = 2013

    ;; Selected variable (standard name)
    VARS = (/"tas"/)

    ;; Name in the raw data
    NAME = (/"temp"/)

    ;; Field
    FIELD = "T2Ms"

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    diag_script = "reformat_obs_HadCRUT.ncl"
    enter_msg(diag_script, "", 1)

    ;; Create output dir
    if (.not.isfilepresent(OUTDIR)) then
        system("mkdir -p " + OUTDIR)
    end if

    ;; Read file
    fname = INDIR + "HadCRUT3v.nc"
    f = addfile(fname, "r")

    ;; Read variable
    output = f->temp

    info_output("", verbosity, 1)
    info_output("Processing " + VARS(0), verbosity, 1)

    ;; Set variable attributes
    tmp = var_attrib(output, VARS(0))
    output = tmp
    delete(tmp)

    ;; Delete level coordinate (dimension 1)
    tmp = rm_single_dims(output)
    delete(output)
    output = tmp
    delete(tmp)

    ;; Format time coordinate
    output_new = output(0:dimsizes(output&t) - 3, :, :)
    delete(output)
    output = output_new
    delete(output_new)
    output!0 = "time"
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
    gAtt@source        = \
        "http://www.metoffice.gov.uk/hadobs/hadcrut3/data/download.html"
    gAtt@reference     = "Not available"
    gAtt@title         = "HadCRUT ground data " + \
        "reformatted for the ESMValTool"
    gAtt@conventions   = "CF/CMOR"

    ;; Outfile
    fout = OUTDIR + "OBS_HadCRUT_ground_3v_T2Ms_tas_" + \
        YEAR1 + "01-" + YEAR2 + "12.nc"

    ;; Write variable
    write_nc(fout, VARS(0), output, gAtt)
    delete(gAtt)
    delete(output)

    leave_msg(diag_script, "", 1)

end
