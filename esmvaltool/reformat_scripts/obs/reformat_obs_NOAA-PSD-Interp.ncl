;;#############################################################################
;; REFORMAT SCRIPT FOR Interpolated NOAA rlut OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://www.esrl.noaa.gov/psd/data/gridded/data.interp_OLR.html
;;
;; Last access
;;    20150810
;;
;; Download and processing instructions
;;    Download the file "olr.day.mean.nc" from the above link,
;;    update paths in script and execute
;;
;; Caveats
;;    Requires NCL v6.1.1
;;    If you use the Interpolated OLR data in a publication, please cite
;;    Liebmann  and Smith, Bull. Am. Meteorol. Soc. 77, 1275-1277 (1996).
;;
;;    Please note: If you acquire Interpolated OLR data products from PSD, we
;;    ask that you acknowledge us in your use of the data. This may be done by
;;    including text such as Interpolated OLR data provided by the
;;    NOAA/OAR/ESRL PSD, Boulder, Colorado, USA, from their Web site at
;;    http://www.esrl.noaa.gov/psd/ in any documents or publications using
;;    these data. We would also appreciate receiving a copy of the relevant
;;    publications. This will help PSD to justify keeping the Interpolated OLR
;;    data set freely available online in the future. Thank you!
;;
;; Modification history
;;    20151116-A_laue_ax: changed field code in output filename from "T2Ms"
;;                        to "T2Ds"
;;    20150810-A_eval_ma: written.
;;
;;#############################################################################

function get_settings(res[1]:logical)
;; Settings for the various variables availabe in the WHOI data set
begin

    ;; Default settings
    res@year1 = 1975
    res@year2 = 2013
    res@orig_missing_value = 32766.
    res@variable_offset = 327.65
    res@scaling_factor = 0.01

    res@attrCopy_err2err = (/"comments"/)
    res@attrCopy_var2err = (/"units"/)

    res@raw_data_folder = getenv("ESMValTool_RAWOBSPATH") + "/Tier2/NOAA-PSD-Interp/"
    res@processed_data_folder = getenv("ESMValTool_OBSPATH") + "/Tier2/NOAA-PSD-Interp/"

    ;; Variable specific settings
    ;; Monthly mean surface latent heat flux, positive upward"
    res@input_filename = "olr.day.mean.nc"
    res@orig_vars = "olr"
    res@vars = "rlut"

    return(res)
end

function append_history(infile[1]:file,
                        to_add[1]:string)
;; Creates (or appends) a history string
local newline, f, datestamp, history
begin

    newline = str_get_nl()
    datestamp = systemfunc("date -u")

    ;; Check existing history attribute
    if (isatt(infile, "history")) then
        history = infile@history
    else
        history = ""
    end if

    ;; Create or append history
    if (strlen(history) .gt. 0) then
        history = datestamp + ": " + to_add + newline + history
    else
        history = datestamp + ": " + to_add
    end if
    return(history)

end

function get_output_filename(res[1]:logical,
                             year1[1]:integer,
                             year2[1]:integer)
;; Define the CMIP5 DRS output filename
local outfile_string
begin

    outfile_string = res@processed_data_folder + "/" + \
        "OBS_NOAA-PSD-Interp_ground_1_T2Ds_" + res@vars + "_" + \ 
        tostring(year1) + "01-" + tostring(year2) + "12.nc"

    if (isfilepresent(outfile_string)) then
        system("rm -f " + outfile_string)
    end if

    return(outfile_string)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/constants.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin
    res = True
    res = get_settings(res)

    diag_script = "reformat_obs_NOAA-PSD-interp-rlut.ncl"

    ;; Create output dir
    if (.not. isfilepresent_esmval(res@processed_data_folder)) then
        system("mkdir -p " + res@processed_data_folder)
    end if

    ;; Loop over all years (WHOI data is provided in yearly files)

    ;; Read input variable and fix _FillValue
    ;; (see caveat abover regarding the warning)
    infile = addfile(res@raw_data_folder + "/"\
                     + res@input_filename, "r")
    invar = infile->$res@orig_vars$
    invar_coords = invar
    invar := tofloat(invar)
    invar := where(invar .eq. 32766., default_fillvalue("float"), invar)
    invar@_FillValue = default_fillvalue("float")
    invar = invar * res@scaling_factor
    invar = invar + res@variable_offset
    copy_VarCoords(invar_coords, invar)

    ;; Reformat the time coordinate
    time_conv = toint(cd_calendar(invar&time, 0))
    time_conv(:, 3) = 12
    time_conv(:, 4) = 0
    time_conv(:, 5) = 0
    time_ind = ind((time_conv(:, 0) .ge. res@year1) .and. (time_conv(:, 0) .le. res@year2))
    time_subset = time_conv(time_ind, :)

    new_units = "days since 1950"
    cal_res = 0
    cal_res@calendar = "proleptic_gregorian"
    time     = cd_inv_calendar(time_subset(:, 0),\
                               time_subset(:, 1),\
                               time_subset(:, 2),\
                               time_subset(:, 3),\
                               time_subset(:, 4),\
                               time_subset(:, 5),\
                                     new_units,\
                                       cal_res)

    FIELD = "T2Ds"
    time@YEAR1 = time_subset(0, 0)
    time@YEAR2 = time_subset(dimsizes(time_ind) - 1, 0)
    time@FIELD = FIELD

    time := time_attrib(time)
    time_size = dimsizes(time)

    ;; Create output variable
    outvar = new((/time_size, dimsizes(invar&lat), dimsizes(invar&lon)/), float)

    ;; Set variable attributes
    outvar := var_attrib(invar(time_ind, :, :), res@vars)

    delete(outvar&time)
    outvar!0 = "time"
    outvar&time = time

    ;; Reformat latitude coordinate
    outvar!1 = "lat"
    outvar&lat = invar&lat
    if (isMonotonic(outvar&lat) .eq. 0) then
        error_msg("f", diag_script, "",\
                  "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(outvar&lat) .eq. -1) then
        outvar = outvar(:, ::-1, :)
    end if
    clat := lat_attrib(outvar&lat)
    outvar&lat := clat

    ;; Reformat longitude coordinate
    outvar!2 = "lon"
    outvar&lon = invar&lon
    if (isMonotonic(outvar&lon) .eq. 0) then
        error_msg("f", diag_script, "",\
                  "non-monotonic longitude coordinate")
    end if
    if (any(outvar&lon .lt. 0.)) then
        outvar = lonFlip(outvar)
    end if
    clon := lon_attrib(outvar&lon)
    outvar&lon := clon

    ;; Write reformatted variable/error estimate to file
    outfile_string = get_output_filename(res, res@year1, res@year2)
    outfile = addfile(outfile_string, "c")
    filedimdef(outfile, "time", -1, True)

    outfile->$res@vars$ = outvar
    delete(time)
    delete(invar)
    delete(invar_coords)
    delete(outvar)

    ;; Copy along global attributes
    global_attnames = getvaratts(infile)
    do i = 0, dimsizes(global_attnames) - 1
        old_attr_prefix = "original-attr_"
        attr = old_attr_prefix + global_attnames(i)
        ;; 'history' attribute is handled below
        if (.not. (global_attnames(i) .eq. "history")) then
            outfile@$attr$ =  infile@$global_attnames(i)$
        end if
    end do

    ;; Set global attributes
    outfile@history = \
        append_history(infile, "data reformatted to be used in ESMValTool")
    outfile@host    = systemfunc("hostname")
    outfile@user    = systemfunc("whoami")
    outfile@field   = FIELD
    outfile@tier    = 2
    outfile@source  = \
        "http://www.esrl.noaa.gov/psd/data/gridded/data.interp_OLR.html"
    outfile@reference = \
        "Liebmann and Smith (Bulletin of the American Meteorological " + \
        "Society, 77, 1275-1277, June 1996)"
    outfile@title   = "NOAA Interpolated Outgoing Longwave Radiation (OLR)"
    outfile@conventions = "CF/CMOR"

end
