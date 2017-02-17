;;#############################################################################
;; REFORMAT SCRIPT FOR MODIS Level 3 OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 2: other freely-available dataset.
;;
;; Source
;;    http://modis-atmos.gsfc.nasa.gov/index.html
;;    The data set fetched is the,
;;         "MOD 08 - Gridded Atmospheric Product"
;;    i.e., monthly level 3 files from collection 6,
;;    these files are named,
;;         'MYD08_M3.AYYYYDDD.006*'
;;
;; Recommended acknowledgment for MODIS data (fetched from MODAPS)
;; Q: How do I acknowledge the use of MODAPS Data?
;;
;; A: The distribution of the MODAPS data sets are funded by NASA's
;; Earth-Sun System Division (ESSD). The data are not copyrighted;
;; however, in the event that you publish data or results using
;; these data, we request that you include the following acknowledgment:
;;
;; "The data used in this study were acquired as part of the NASA's
;; Earth-Sun System Division and archived and distributed by the
;; MODIS Adaptive Processing System (MODAPS)."
;;
;; We would appreciate receiving a copy of your publication, which
;; can be forwarded to MODAPSUSO@lists.nasa.gov.
;;
;; Last access
;;    20150504
;;
;; Download and processing instructions
;;    Download data follow instructions on the above link
;;    Edit the settings available in the
;;           'get_settings(...)'
;;    function below to be consistent with your system.
;;    Run with:
;;       ncl reformat_scripts/obs/reformat_obs_modis.ncl
;;
;; Caveats
;;    - Requires NCL v6.1.1 (due to the ':='-operator
;;    - Data is converted into separate monthly files, not yearly
;;
;; Modification history
;;    20151118-A-laue_ax: bugfix: added unit conversion
;;                        clivi, clwvi: g/m2 --> kg/m2
;;                        clt: frac --> %
;;    20150430-eval_ma: written.
;;
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./reformat_scripts/constants.ncl"

function get_settings(res [1]:logical)
; Settings for the various variables availabe in the WHOI data set
begin
    ;; Default settings
    res@year1 = 2002  ; start year
    res@year2 = 2014  ; end year
    res@modis_filename_prefix = "MYD08_M3.A"

    ;; Dimension names in the HDF file (lat/lon)
    res@lat_var = "YDim"
    res@lon_var = "XDim"
    res@lat_name = "YDim_mod08"
    res@lon_name = "XDim_mod08"

    ;; Variables to extract from HDF
    ;;                        MODIS variable name ("prefix")             CMIP5 name (prefix)
    res@vars2copy_root = (/(/"Cloud_Fraction_Mean",                     "clt"             /),\
                           (/"Cloud_Water_Path_Liquid_Mean",            "clwvi"           /),\
                           (/"Cloud_Water_Path_Ice_Mean",               "clivi"           /)/)
;                          (/"Cloud_Optical_Thickness_1L_Liquid_Mean",  ""                /),\
;                          (/"Cloud_Optical_Thickness_1L_Ice_Mean"   ,  ""                /)/)

    ;;                       MODIS variable name ("suffix")              CMIP5 name (suffix)
    res@vars2copy_var = (/(/"Mean",                                     ""                /)/)
;                          (/"Std",                                     "-std"            /),\
;                          (/"Min",                                     "-min"            /),\
;                          (/"Max",                                     "-max"            /)/)

    ;; Input/ouput/working foldesr
;    res@raw_data_folder = "/nobackup/smhid11/sm_abhay/MODIS_L3_C6"
;    res@processed_data_folder = "/nobackup/rossby17/sm_maeva/requests/modis-collection6/processed_data"
    res@raw_data_folder = "/data/ESMValTool/obs/RAW/Tier2/MODIS_L3_C6"
    res@processed_data_folder = "/data/ESMValTool/obs/Tier2/MODIS_L3_C6"
    res@scratch_dir = "/scratch/local/"

    return res
end

function get_time_from_filename(infile [1]:string)
; Compute the current time using the MODIS filename structure
local infile_str, yyyyddd, yyyymmdd, year, month,
    day, cal_res, time
begin

    yyyyddd  = toint(str_get_cols(infile, 10, 16))
    yyyymmdd = yyyyddd_to_yyyymmdd(yyyyddd)
    year     = toint(str_get_cols(tostring(yyyymmdd), 0, 3))
    month    = toint(str_get_cols(tostring(yyyymmdd), 4, 5))
    day      = toint(str_get_cols(tostring(yyyymmdd), 6, 7))
    cal_res = 0
    cal_res@calendar = "proleptic_gregorian"
    time = cd_inv_calendar(year, month, day, 0, 0, 0, TUNITS, cal_res)
    return time
end

function get_year_from_filename(infile [1]:string)
; Wrapper 'function to get_time_from_filename(...)' to
; retrieve current year
local time, times, year
begin
    time = get_time_from_filename(infile)
    times = cd_calendar(time, 0)
    year = toint(times(:, 0))
    return year
end

function get_month_from_filename(infile [1]:string)
; Wrapper 'function to get_time_from_filename(...)' to
; retrieve current month
local time, times, month
begin
    time = get_time_from_filename(infile)
    times = cd_calendar(time, 0)
    month = toint(times(:, 1))
    return month
end

procedure set_time_attributes(time [1]:double)
; Set CMIP5 time attributes
begin
    time@standard_name = "time"
    time@long_name = "time"
    time@axis = "T"
end

function append_history(infile [1]:file,
                        to_add [1]:string)
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

function get_output_filename(res [1]:logical,
                        var_name [1]:string,
                             mip [1]:string,
                            year [1]:integer,
                           month [1]:integer)
; Define the CMIP5 DRS output filename
local outfile_string, month_str
begin
    month_str = sprinti("%02i", month)
    outfile_string = res@processed_data_folder + "/"\
                                               + var_name + "_"\
                                               + mip + "_"\
                                               + "MODIS-L3_"\
                                               + "obs_"\
                                               + "C6_"\
                                               + tostring(year) + month_str + "-"\
                                               + tostring(year) + month_str + ".nc"
    if (isfilepresent(outfile_string)) then
        system("rm -f " + outfile_string)
    end if
    return outfile_string
end

function get_modis_suffix(sufficies [*]: string,
                          index [1]:integer)
; Fetch the last part of the Modis variable. A wrapper
; is needed as the array might be both 1D or 2D
local ret_val, suff_sizes
begin
    suff_sizes = dimsizes(dimsizes(sufficies))
    if (suff_sizes .eq. 1) then
        ret_val = sufficies(0)
    else if(suff_sizes .eq. 2) then
        ret_val = sufficies(index, 0)
    else
        print("fatal: in get_modis_suffix (reformat_obs_modis.ncl) invalid array size")
    end if
    end if
    return ret_val
end

function get_cmip_suffix(sufficies [*]: string,
                         index [1]:integer)
; Fetch the last part of the Modis variable. A wrapper
; is needed as the array might be both 1D or 2D
local ret_val, suff_sizes
begin
    suff_sizes = dimsizes(dimsizes(sufficies))
    if (suff_sizes .eq. 1) then
        ret_val = sufficies(1)
    else if(suff_sizes .eq. 2) then
        ret_val = sufficies(index, 1)
    else
        print("fatal: in get_cmip_suffix (reformat_obs_modis.ncl) invalid array size")
    end if
    end if
    if (strlen(ret_val) .gt. 0) then
        ret_val = "-" + ret_val
    end if
    return ret_val
end

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin
    res = True
    res = get_settings(res)

    diag_script = "reformat_obs_modis.ncl"

    ;; Create output dir
    if (.not. isfilepresent_esmval(res@processed_data_folder)) then
        system("mkdir -p " + res@processed_data_folder)
    end if

    ;; Loop over all indicated years
    do yy = res@year1, res@year2
        timestamp = systemfunc("date")
        print(timestamp + "  --  year: " + yy\
                        + " (in range " + res@year1\
                        + " to " + res@year2 + ")")
        files := systemfunc("ls -1 " + res@raw_data_folder + "/"\
                                     + res@modis_filename_prefix\
                                     + yy + "*")

        ;; Loop over all files in indicated year
        ;; (MODIS is provided in monthly files)
        do file_idx = 0, dimsizes(files) - 1
            inpath_str = files(file_idx)
            infile_str = systemfunc("basename " + inpath_str)
            fin = addfile(inpath_str, "r")

            ;; Current date is encoded in filename
            time = get_time_from_filename(infile_str)
            year = get_year_from_filename(infile_str)
            month = get_month_from_filename(infile_str)
            set_time_attributes(time)

            ;; Loop over variables to fecth from input file
            no_vars_root = dimsizes(res@vars2copy_root)
            ;; The "root" part of the MODIS variable filename
            do var_idx_root = 0, no_vars_root(0) - 1
                vars_root = res@vars2copy_root
                vars = res@vars2copy_var
                no_vars = dimsizes(vars)
                ;; The "suffix" (last) part of the MODIS variable filename
                do var_idx = 0, no_vars(0) - 1
                    invar_str = vars_root(var_idx_root, 0) + "_" + get_modis_suffix(vars, var_idx)
                    outvar_str = vars_root(var_idx_root, 1)
                    outvar_str = outvar_str + get_cmip_suffix(vars, var_idx)

                    ;; Read file, scale data, and copy along meta data
                    invar = fin->$invar_str$
                    invar_fv = invar@_FillValue
                    invar_coords = invar
                    invar := tofloat(invar)
                    invar := where(invar .eq. tofloat(invar_fv), default_fillvalue("float"), invar)
; A-laue_ax+
                    ; convert to correct units
                    if (outvar_str.eq."clt") then
                         invar = 100.0 * invar  ; frac --> %
                    end if
                    if ((outvar_str.eq."clivi").or.(outvar_str.eq."clwvi")) then
                         invar = 0.001 * invar  ; g/m2 --> kg/m2
                    end if
; A-laue_ax-
                    invar@_FillValue = default_fillvalue("float")
                    copy_VarCoords(invar_coords, invar)
                    if (isatt(invar_coords, "scale_factor")) then
                        invar = invar * tofloat(invar_coords@scale_factor)
                    end if
                    if (isatt(invar_coords, "add_offset")) then
                        invar = invar + tofloat(invar_coords@add_offset)
                    end if

                    ;; Create output variable
                    lat = fin->$res@lat_var$
                    lon = fin->$res@lon_var$
                    outvar = new((/1, dimsizes(lat), dimsizes(lon)/), float)

                    ;; Set variable attributes
                    if (get_modis_suffix(vars, var_idx) .eq. "Mean") then
                        outvar(0, :, :) = var_attrib(invar, outvar_str)
                    end if

                    ;; Set reformatted time dimension
                    outvar!0 = "time"
                    outvar&time = time

                    ;; Set reformatted latitude dimension
                    outvar!1 = "lat"
                    outvar&lat = lat
                    if (isMonotonic(outvar&lat) .eq. 0) then
                        error_msg("fatal",\
                                  diag_script,\
                                  "",\
                                  "non-monotonic latitude coordinate")
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
                        error_msg("fatal",\
                                  diag_script,\
                                  "",\
                                  "non-monotonic longitude coordinate")
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
                                                  year,\
                                                  month)
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
                    fout@title   = "MOD 08 - Gridded Atmospheric Product reformatted for ESMValTool"
                    fout->$outvar_str$ = outvar

                    delete(invar)
                    delete(invar_coords)
                    delete(outvar)
                    delete(variable_attnames)
                end do
            end do
        end do
    end do
end
