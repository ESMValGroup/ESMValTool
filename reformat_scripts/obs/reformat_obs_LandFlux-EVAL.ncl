;;#############################################################################
;; REFORMAT SCRIPT FOR LandFlux-EVAL OBSERVATIONAL DATA
;;#############################################################################
;;
;; Tier
;;    Tier 3: restricted dataset.
;;
;; Source
;;    http://www.iac.ethz.ch/groups/seneviratne/research/LandFlux-EVAL
;;    https://data.iac.ethz.ch/landflux/  (login required)
;;
;; Last access
;;    20150722
;;
;; Download and processing instructions
;;    Register to get user name and password
;;    Download the file LandFluxEVAL.merged.89-05.monthly.all.nc
;;    Edit the settings available in the get_settings(...) function below to be
;;    consistent with your system.
;;
;; Caveats
;;    Requires NCL v6.1.1 (due to the ':='-operator)
;;    See policy info at:
;;    http://www.iac.ethz.ch/groups/seneviratne/research/LandFlux-EVAL/DataPolicy
;;
;; Modification history
;;    20150722-righ_ma: outfile name changed according to OBS class and pep8.
;;    20150430-eval_ma: written.
;;
;;#############################################################################

load "./reformat_scripts/constants.ncl"

function get_settings(res[1]:logical)
; Settings for the various variables availabe in the WHOI data set
begin

    ;; Input/ouput/working foldesr
    res@raw_data_inputfile = getenv("ESMValTool_RAWOBSPATH") + "/Tier3/" + \
        "LandFlux-EVAL/LandFluxEVAL.merged.89-05.monthly.all.nc"
    res@processed_data_folder = getenv("ESMValTool_OBSPATH") + "/Tier3/" + \
        "LandFlux-EVAL/"

    res@variables_to_extract = (/"ET_mean",\
                                 "ET_sd"/)
    res@varname_mapped2cmor = (/(/"ET_mean", "et"/),\
                                (/"ET_sd", "et-sd"/)/)

    return(res)
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
    return(history)

end

function get_output_filename(res[1]:logical,
                             var_name[1]:string,
                             mip[1]:string,
                             year1[1]:integer,
                             year2[1]:integer)
; Define the output filename
local outfile_string, month_str
begin

    outfile_string = res@processed_data_folder + "/" + \
        "OBS_LandFlux-EVAL_ground_1_T2Ms_" + var_name + "_" + \
        tostring(year1) + "01-" + tostring(year2) + "12.nc"

    if (isfilepresent(outfile_string)) then
        system("rm -f " + outfile_string)
    end if
    return(outfile_string)

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
    if (str_sizes.eq.1) then
        ret_val = var2cmor_array(1)
        if (lf_vars(0) .ne. var2cmor_array(0)) then
            error_msg("f", "get_cmip_name", "reformat_obs_LandFlux-EVAL.ncl",\
                      "requeste variable mismatch: " + lf_vars(0) + \
                      " .ne. " + var2cmor_array(0))
        end if
    else if(str_sizes .eq. 2) then
        ret_val = var2cmor_array(index, 1)
        if (lf_vars(index) .ne. var2cmor_array(index, 0)) then
            error_msg("f", "get_cmip_name", "reformat_obs_LandFlux-EVAL.ncl", \
                      "requeste variable mismatch: " + lf_vars(index) + \
                      " .ne. " + var2cmor_array(index, 0))
        end if
    else
        error_msg("f", "get_cmip_name", "reformat_obs_LandFlux-EVAL.ncl",\
                  "invalid array size")
    end if
    end if
    return(ret_val)

end

load "./interface_scripts/messaging.ncl"
load "./reformat_scripts/obs/reformat_obs_func.ncl"

begin

    res = True
    res = get_settings(res)

    diag_script = "reformat_obs_LandFlux-EVAL.ncl"

    ;; Create output dir
    if (.not. isfilepresent_esmval(res@processed_data_folder)) then
        system("mkdir -p " + res@processed_data_folder)
    end if

    fin = addfile(res@raw_data_inputfile, "r")
    landflux_vars = res@variables_to_extract

    ;; Loop over all indicated variables
    do idx_var = 0, dimsizes(landflux_vars) - 1
        invar_str = landflux_vars(idx_var)
        invar = fin->$invar_str$

        invar_fv = invar@_FillValue
        invar_coords = invar
        invar := tofloat(invar)
        invar := where(invar.eq.tofloat(invar_fv), \
                       default_fillvalue("float"), invar)
        invar@_FillValue = default_fillvalue("float")
        copy_VarCoords(invar_coords, invar)

        ;; Create output variable
        lat = fin->lat
        lon = fin->lon
        time = fin->time
        outvar = new((/dimsizes(time), dimsizes(lat), dimsizes(lon)/), float)
        outvar_str = get_cmip_name(res, idx_var)

        ;; Set variable attributes (but only for ET_mean
        if (invar_str .eq. "ET_mean") then
            outvar(:, :, :) = var_attrib(invar, outvar_str)
        end if

        ;; Reformat the time coordinate
        time_out = new(dimsizes(time), double)
        cal_res = 0
        cal_res@calendar = "proleptic_gregorian"
        do idx_time = 0, dimsizes(time) - 1
            yyyymm = tostring(time(idx_time))
            yyyy = toint(str_get_cols(yyyymm, 0, 3))
            mm = toint(str_get_cols(yyyymm, 4, 5))
            time_out(idx_time) = \
                cd_inv_calendar(yyyy, mm, 15, 0, 0, 0, TUNITS, cal_res)
        end do
        yyyy1 = floor(time(0) / 100.)
        yyyy2 = floor(time(dimsizes(time) - 1) / 100.)
        time_out@YEAR1 = yyyy1
        time_out@YEAR2 = yyyy2
        time_out@FIELD = "T2Ms"
        time_out := time_attrib(time_out)
        outvar!0 = "time"
        outvar&time = time_out

        ;; Set reformatted latitude dimension
        outvar!1 = "lat"
        outvar&lat = lat
        if (isMonotonic(outvar&lat) .eq. 0) then
            error_msg("f", diag_script, "", \
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
            error_msg("fatal", diag_script, "", \
                      "non-monotonic longitude coordinate")
        end if
        if (any(outvar&lon .lt. 0.)) then
            outvar = lonFlip(outvar)
        end if
        clon := lon_attrib(outvar&lon)
        outvar&lon := clon

        ;; Define output file with unlimited time dimension
        outfile = get_output_filename(res, outvar_str, "Amon", toint(yyyy1), \
                                      toint(yyyy2))
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
        fout@history = \
            append_history(fin, "data reformatted to be used in ESMValTool")
        fout@host    = systemfunc("hostname")
        fout@user    = systemfunc("whoami")
        fout@field   = "T2Ms"
        fout@period  = toint(yyyy1) + "-" + toint(yyyy2)
        fout@tier    = 3
        fout@source  = \
            "http://www.iac.ethz.ch/groups/seneviratne/research/LandFlux-EVAL"
        fout@reference = "Mueller, B. et al., Hydrol. Earth Syst. Sci., " + \
            "17, 3707-3720, doi:10.5194/hess-17-3707-2013, 2013."
        fout@title   = "LandFlux-EVAL benchmark synthesis products " + \
            "reformatted for the ESVMalTool"
        fout@conventions = "CF/CMOR"
        fout->$outvar_str$ = outvar

        delete(invar)
        delete(invar_coords)
        delete(outvar)
        delete(variable_attnames)
    end do
end
