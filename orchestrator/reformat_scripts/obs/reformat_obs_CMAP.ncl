;;#############################################################################
;; REFORMAT SCRIPT FOR CPC Merged Analysis of Precipitation (CMAP)
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    https://climatedataguide.ucar.edu/climate-data/cmap-cpc-merged-analysis-precipitation
;;      and actual download via
;;    http://www.esrl.noaa.gov/psd/data/gridded/data.cmap.html
;;
;;    The data set fetched for this script is the Precipitation Mean (Standard
;;         Monthly), the file "precip.mon.mean.nc"
;;
;; Required standard acknowledgment:
;; --------------------------------
;;    Xie, P., and P.A. Arkin, 1997: Global precipitation: A 17-year monthly
;;    analysis based on gauge observations, satellite estimates, and numerical
;;    model outputs. Bull. Amer. Meteor. Soc., 78, 2539 - 2558.
;;
;;    """Please note: If you acquire CMAP Precipitation data products from PSD,
;;    we ask thatyou acknowledge us in your use of the data. This may be done by
;;    includingtext such as CMAP Precipitation data provided by the
;;    NOAA/OAR/ESRL PSD, Boulder,Colorado, USA, from their Web site at
;;    http://www.esrl.noaa.gov/psd/ in any documents or publications using these
;;    data. We would also appreciate receiving a copy of the relevant
;;    publications. This will help PSD to justify keeping the CMAP Precipitation
;;    data set freely available online in the future. Thank you!"""
;;
;; Other references:
;; ------------------
;;    Huffman, G. J. and co-authors, 1997: The Global Precipitation Climatology
;;    Project (GPCP) combined data set. Bull. Amer. Meteor. Soc., 78, 5-20.
;;
;;    Reynolds, R. W., 1988: A real-time global sea surface temperature
;;    analysis. J. Climate, 1, 75-86.
;;
;;    Spencer, R. W., 1993: Global oceanic precipitation from the MSU during
;;    1979-91 and comparisons to other climatologies. J. Climate, 6, 1301-1326.
;;
;;    Xie P., and P. A. Arkin, 1997: Global precipitation: a 17-year monthly
;;    analysis based on gauge observations, satellite estimates, and numerical
;;    model outputs. Bull. Amer. Meteor. Soc., 78, 2539-2558.
;;
;; Last access
;;    20150828
;;
;; Download and processing instructions
;;    Download data follow instructions on the above link
;;    Edit the settings available in the
;;           'get_settings(...)'
;;    function below to be consistent with your system.
;;    Run with:
;;       ncl reformat_scripts/obs/reformat_CMAP.ncl
;;
;; Caveats
;;    - Requires NCL v6.1.1 (due to the ':='-operator
;;
;; Modification history
;;    20150828-eval_ma: fixed header
;;    20150430-eval_ma: written.
;;
;;#############################################################################

load "./reformat_scripts/obs/reformat_obs_func-local.ncl"

function get_settings(res[1]:logical)
; Settings for the various variables availabe in the WHOI data set
begin
    ;; Input/ouput/working foldesr
;    res@raw_data_inputfile = "/nobackup/rossby17/sm_maeva/tmp/tmp/in/precip.mon.mean.nc"
    res@raw_data_inputfile = getenv("ESMValTool_RAWOBSPATH") + "/Tier2/CMAP/precip.mon.mean.nc"
;    res@processed_data_folder = "/nobackup/rossby17/sm_maeva/tmp/tmp/out/"
    res@processed_data_folder = getenv("ESMValTool_OBSPATH") + "/Tier2/CMAP/"
;    res@start_year = 1980
    res@start_year = 1979
    res@end_year = 2013

;-;    res@variables_to_extract = (/"precip"/)
;-;    res@varname_mapped2cmor = (/(/"precip", "pr"/)/)
    res@variables_to_extract = (/"precip"/)
    res@varname_mapped2cmor = (/(/"precip", "pr"/)/)

    return res
end

procedure set_time_attributes(time[1]:double)
; Set CMIP5 time attributes
begin
    time@standard_name = "time"
    time@long_name = "time"
    time@axis = "T"
end

function append_history(infile[1]:file,
                        to_add[1]:string)
; Creates (or appends) a history string
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
    return history
end

function get_output_filename(res[1]:logical,
                        var_name[1]:string,
                             mip[1]:string,
                           year1[1]:integer,
                           year2[1]:integer)
; Define the CMIP5 DRS output filename
local outfile_string, month_str
begin
    outfile_string = res@processed_data_folder + "/" + \
        "OBS_CMAP_mergedanaly_1_T2Ms_pr_" + \
        tostring(year1) + "01-" + tostring(year2) + "12.nc"
    if (isfilepresent(outfile_string)) then
        system("rm -f " + outfile_string)
    end if
    return outfile_string
end

function get_cmip_name(resource : logical,
                       index[1]:integer)
; Fetch item from a list or a list of lists
; transparently
local ret_val, suff_sizes, var2cmor_array
begin
    var2cmor_array = resource@varname_mapped2cmor
    lf_vars = resource@variables_to_extract

    str_sizes = dimsizes(dimsizes(var2cmor_array))
    if (str_sizes .eq. 1) then
        ret_val = var2cmor_array(1)
        if (lf_vars(0) .ne. var2cmor_array(0)) then
            print("fatal: in get_cmip_name (reformat_obs_CMAP.ncl), "\
                        + "requeste variable mismatch: " + lf_vars(0)\
                                                         + " .ne. "\
                                                         + var2cmor_array(0))
            status_exit(1)
        end if
    else if(str_sizes .eq. 2) then
        ret_val = var2cmor_array(index, 1)
        if (lf_vars(index) .ne. var2cmor_array(index, 0)) then
            print("fatal: in get_cmip_name (reformat_obs_CMAP.ncl), "\
                    + "requeste variable mismatch: " + lf_vars(index)\
                                                     + " .ne. "\
                                                     + var2cmor_array(index, 0))
            status_exit(1)
        end if
    else
        print("fatal: in get_cmip_name (reformat_obs_CMAP.ncl), "\
              + "invalid array size")
        status_exit(1)
    end if
    end if
    return ret_val
end

begin
    res = True
    res = get_settings(res)

    diag_script = "reformat_obs_CMAP.ncl"

    ;; Create output dir
;    if (.not. fileexists(res@processed_data_folder)) then    ; NCL v6.2.1
    if (.not. isfilepresent(res@processed_data_folder)) then  ; NCL v6.2.0
        system("mkdir -p " + res@processed_data_folder)
    end if

    fin = addfile(res@raw_data_inputfile, "r")
    requested_vars = res@variables_to_extract

    ;; Fix the time dimension
    time = fin->time

    ;; Reformat the time coordinate
    cd_time = cd_calendar(time, 0)
    sy_idx = min(ind(cd_time(:, 0) .eq. res@start_year))
    ey_idx = max(ind(cd_time(:, 0) .eq. res@end_year))

    cal_res = 0
    cal_res@calendar = "proleptic_gregorian"
    day = ispan(sy_idx, ey_idx, 1)
    hour = ispan(sy_idx, ey_idx, 1)
    minute = ispan(sy_idx, ey_idx, 1)
    second = ispan(sy_idx, ey_idx, 1)
    day = 15
    hour = 0
    minute = 0
    second = 0
    time_out = cd_inv_calendar(toint(cd_time(sy_idx:ey_idx, 0)),\
                               toint(cd_time(sy_idx:ey_idx, 1)),\
                               day,\
                               hour,\
                               minute,\
                               second,\
                               TUNITS,\
                               cal_res)

    time_out@YEAR1 = res@start_year
    time_out@YEAR2 = res@end_year
    time_out@FIELD = "T2Ms"
    time_out := time_attrib(time_out)

    ;; Copy along time-coordinate attributes
    time_attnames = getvaratts(time)
    do i = 0, dimsizes(time_attnames) - 1
        old_attr_prefix = "original-attr_"
        attr = old_attr_prefix + time_attnames(i)
        time_out@$attr$ = time@$time_attnames(i)$
    end do

    ;; Loop over all indicated variables
    do idx_var = 0, dimsizes(requested_vars) - 1

        invar_str = requested_vars(idx_var)
        invar = fin->$invar_str$

        invar_fv = invar@_FillValue
        invar_coords = invar
        invar := tofloat(invar)
        invar := where(invar .eq. tofloat(invar_fv),\
                       default_fillvalue("float"),\
                       invar)
        invar@_FillValue = default_fillvalue("float")
        copy_VarCoords(invar_coords, invar)
        invar := invar(sy_idx:ey_idx, :, :)
        do i = 0, dimsizes(time_attnames) - 1
            delete(invar&time@$time_attnames(i)$)
        end do

        ;; Output variable dimensions
        lat = fin->lat
        lon = fin->lon

        outvar = new((/dimsizes(time_out),\
                       dimsizes(lat),\
                       dimsizes(lon)/),\
                       float)
        outvar_str = get_cmip_name(res, idx_var)

        outvar_attr = True
        outvar_attr@standard_name = "precipitation_flux"
        outvar_attr@units = "kg m-2 s-1"
        outvar_attr@cell_methods = "time: mean"
        outvar_attr@cell_measures = "area: areacella"
        outvar_attr@long_name = "Precipitation"
        outvar_attr@comment = "at surface; includes both liquid and solid phases from all types of clouds (both large-scale and convective)"

        outvar(:, :, :) = var_attrib(invar, outvar_attr)
; laue_ax+
        ; convert from mm/day to kg m-2 s-1
        outvar(:, :, :) = outvar(:, :, :) / (24.0 * 3600.0)
; laue_ax-

        outvar!0 = "time"
        outvar&time = time_out

        ;; Set reformatted latitude dimension
        outvar!1 = "lat"
        outvar&lat = lat
        if (isMonotonic(outvar&lat) .eq. 0) then
            print("fatal: in "\
                  + diag_script\
                  + ", non-monotonic latitude coordinate")
            status_exit(1)
        end if
        if (isMonotonic(outvar&lat) .eq. -1) then
            outvar = outvar(:, ::-1, :)
        end if
        clat := lat_attrib(outvar&lat)
        outvar&lat := clat

        ;; Set reformatted longitude dimension
        outvar!2 = "lon"
        outvar&lon = lon
        if (isMonotonic(outvar&lon) .eq. 0) then
            print("fatal: in "\
                  + diag_script\
                  + ", non-monotonic longitude coordinate")
            status_exit(1)
        end if
        if (any(outvar&lon .lt. 0.)) then
            outvar = lonFlip(outvar)
        end if
        clon := lon_attrib(outvar&lon)
        outvar&lon := clon

        ;; Define output file with unlimited time dimension
        outfile = get_output_filename(res,\
                                      outvar_str,\
                                      "Amon",\
                                      toint(res@start_year),\
                                      toint(res@end_year))
        fout = addfile(outfile, "c")
        filedimdef(fout, "time", -1, True)

        ;; Copy along global attributes
        global_attnames = getvaratts(fin)
        do i = 0, dimsizes(global_attnames) - 1
            old_attr_prefix = "original-attr_"
            attr = old_attr_prefix + global_attnames(i)
            ;; 'history' attribute is handled below
            if (ismissing(str_match(global_attnames(i), "Metadata_"))) then
                fout@$attr$ =  fin@$global_attnames(i)$
            end if
        end do

        ;; Copy along variable attributes
        variable_attnames = getvaratts(invar_coords)
        do i = 0, dimsizes(variable_attnames) - 1
            old_attr_prefix = "original-attr_"
            attr = old_attr_prefix + variable_attnames(i)
            ;; 'history' attribute is handled below
            if (.not. (variable_attnames(i) .eq. "history")) then
                outvar@$attr$ = invar_coords@$variable_attnames(i)$
            end if
        end do

        ;; Set global attributes sprecific for this script
        fout@history = append_history(fin, "data reformatted to be used in ESMValTool")
        fout@host    = systemfunc("hostname")
        fout@user    = systemfunc("whoami")
        fout@field   = "T2Ms"
        fout@tier    = 2
        fout@source  = "http://www.esrl.noaa.gov/psd/data/gridded/data.cmap.html"
        fout@title   = "CPC Merged Analysis of Precipitation (CMAP)"
        fout->$outvar_str$ = outvar

        delete(invar)
        delete(invar_coords)
        delete(outvar)
        delete(variable_attnames)
    end do
end
