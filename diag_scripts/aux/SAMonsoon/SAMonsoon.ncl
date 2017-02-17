;;#############################################################################
;; GENERAL ROUTINES FOR MONSOON DIAGOSTICS
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;     function compute_stat
;;     function fetch_years
;;     function fetch_unique_years
;;     function remove_figure_borders
;;     function set_string_once
;;     function get_plot_main_string
;;     function get_diffplot_right_string
;;     function getFillValue_with_default
;;     function get_diffplot_left_string
;;     function set_default_switches
;;     function single_contour
;;     function single_vector
;;     function map_zoom
;;     function apply_transforms
;;     function multipanel_get_no_cols
;;     function multipanel_get_no_rows
;;     function get_figures_per_page
;;     function get_latlon_string
;;     function fetch_first_last_year
;;     function remove_entry_from_set
;;     function list_attributes
;;     function attach_attributes
;;     function count_number_of_entries
;;     function get_entry_indicies
;;     function retrieve_storage_record
;;     function prefix_filter
;;     function get_stored_attribute_names
;;     function find_dim
;;     function flatten_to_string
;;     function string_to_array
;;     function retrieve_from_vault
;;     function netcdf_to_vault
;;     function band_pass_time_lat
;;     function time_xcor_lag
;;     function mjo_xcor_lag_ovly_local
;;     function clmDayTLL_calendar
;;     function smthClmDayTLL_calendar
;;     function calcDayAnomTLL_calendar
;;     function write_tuple_with_dash
;;     function no_unique_vault_entries
;;     function empty_string_if_missing_attribute
;;     function keep_only_n_cols
;;     procedure set_res_to_false
;;     procedure inset_top_text
;;     procedure rebalance_wrt_min
;;     procedure place_models_on_pages
;;     procedure add_storage_record
;;     procedure vault_to_netcdf
;;     procedure store_variable
;;     procedure set_missing_years
;;
;;#############################################################################
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
;load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"

undef("compute_stat")
function compute_stat(requested_stats[*] : string,
                      valid_stats[*] : string,
                      data1 : numeric,
                      data2 : numeric)
;; Arguments
;;    requested_stats: Request statistics
;;        valid_stats: Valid statistics
;;              data1: First data set
;;                        @yrs_intersection: intersection of years
;;                        @idx: model index
;;                        @refA: reference area in lat/lon
;;                        @latB: latitude boundary
;;                        @lonB: longitude boundary
;;                        @bob: Bay of Bengal area avearge
;;                        @eio: eastern equatorial Indian Ocean
;;                        @sa34: India region area average
;;                        @en34: enso34 area average
;;              data2: Second data set
;;
;; Return value
;;    statistics: a logical with all the requested (and valid) statistics
;;                attached
;;
;; Caveats
;;
;; References
;;
;; Modifiation history
;;    2015-07-02-A_eval_ma: written.
;;
local all_dates,    all_dates_dim,  compute_these_stats,\
    first_year,   idx,            last_year,\
    stat_format,  statistics,     verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering compute_stat", verbosity, 8)
    stat_format = "%4.2f"
    statistics = True
    compute_these_stats = intersection(requested_stats, valid_stats)

    ;; Check if compute stat is empty (logical instead of intersection array)
    if(.not. islogical(compute_these_stats)) then
        do idx = 0, dimsizes(compute_these_stats) - 1
            if (compute_these_stats(idx) .eq. "mean") then
                statistics@mean = sprintf(stat_format, coswgt_areaave(data1))
            else if (compute_these_stats(idx) .eq. "yrs") then
                if (isatt(data1, "yrs_intersection")) then
                    first_year = data1@yrs_intersection(0)
                    last_year = data1@yrs_intersection(dimsizes(data1@yrs_intersection) - 1)
                else
                    first_year = models@start_year(data1@idx)
                    last_year = models@end_year(data1@idx)
                end if
                statistics@yrs = first_year + "-" + last_year
            else if (compute_these_stats(idx) .eq. "refA") then
                statistics@refA = data1@refA
            else if (compute_these_stats(idx) .eq. "latB") then
                statistics@latB = data1@latB
            else if (compute_these_stats(idx) .eq. "lonB") then
                statistics@lonB = data1@lonB
            else if (compute_these_stats(idx) .eq. "rmse") then
                statistics@rmse = sprintf(stat_format, coswgt_arearmse(data1, data2))
            else if (compute_these_stats(idx) .eq. "corr") then
                statistics@corr = sprintf(stat_format, coswgt_pattern_cor(data1, data2))
            else if (compute_these_stats(idx) .eq. "bob") then
                statistics@bob = sprintf(stat_format, data1@bob)
            else if (compute_these_stats(idx) .eq. "eio") then
                statistics@eio = sprintf(stat_format, data1@eio)
            else if (compute_these_stats(idx) .eq. "sahel") then
                statistics@sahel = sprintf(stat_format, data1@sahel)
            else if (compute_these_stats(idx) .eq. "sa34") then
                statistics@sa34 = sprintf(stat_format, data1@sa34)
            else if (compute_these_stats(idx) .eq. "en34") then
                statistics@en34 = sprintf(stat_format, data1@en34)
            end if
            end if
            end if
            end if
            end if
            end if
            end if
            end if
            end if
            end if
            end if
            end if
        end do
    end if
    info_output(">>>>>>>> Leaving compute_stat", verbosity, 8)
    return statistics
end

undef("fetch_years")
function fetch_years(dataset : numeric)
;;                  return_val[*] : numeric
;; Arguments:
;;     dataset   -- A dataset with the "time"-dimension
;; Description:
;;     Returns an array where each entry in the time-dimension is represented
;;     by its year
local dataset_time, dataset_years, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering fetch_years", verbosity, 8)

    dataset_time = cd_calendar(dataset&time, 0)
    dataset_years = dataset_time(:, 0)

    info_output(">>>>>>>> Leaving fetch_years", verbosity, 8)
    return dataset_years
end

undef("fetch_unique_years")
function fetch_unique_years(dataset : numeric)
;;                  return_val[*] : numeric
;; Arguments:
;;     dataset   -- A dataset with the "time"-dimension
;; Description:
;;     Returns an array with the years (doublets removed)
local dataset_years, verbosity, uniq_ind, years
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering fetch_unique_years", verbosity, 8)

    dataset_years = fetch_years(dataset)
    uniq_ind = UNIQ(dataset_years)
    years = dataset_years(uniq_ind)

    info_output(">>>>>>>> Leaving fetch_unique_years", verbosity, 8)
    return years
end

undef("set_res_to_false")
procedure set_res_to_false(res[1] : logical, \
                           res_list[*] : string)
;; Arguments:
;;    res - Plot resource variable
;;    res_list - Specific resources to update
;; Description:
;;    Set all res-resources in res_list to False
local verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering set_res_to_false", verbosity, 8)
    dim_RES = dimsizes(res_list)
    do idx = 0, dim_RES - 1
        res@$res_list(idx)$ = False
    end do
    info_output(">>>>>>>> Leaving set_res_to_false", verbosity, 8)
end

undef("remove_figure_borders")
function remove_figure_borders(res[1] : logical)
;; Arguments:
;;    res - Plot resource variable
;; Returns:
;;    res - Updated resource
;; Description:
;;    Removes borders and tickmarks from a plot resources
;;    their difference on a single page
local verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering remove_figure_borders", verbosity, 6)
    set_res_to_false(res, (/"tmXBLabelsOn", "tmXBOn", "tmYLOn", \
                            "tmXTOn", "tmYROn", "tmXTBorderOn", \
                            "tmXBBorderOn", "tmYLBorderOn", \
                            "tmYRBorderOn"/))
    info_output(">>>>>>>> Leaving remove_figure_borders", verbosity, 6)
    return res
end

undef("set_string_once")
function set_string_once(title[1] : string,
                         idx_mod[1] : integer)
;;                       return_val[1] : string
;; Arguments:
;;    title - Input title
;;    idx_mod - Current model index
;; Returns:
;;    use_title - Title to use
;; Description:
;;    Uses the input title for model index 0,
;;    otherwise empty string
local title, idx_mod, use_title, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering set_string_once", verbosity, 6)
    use_title = ""
    if (idx_mod .eq. 0) then
        use_title = title
    end if
    info_output(">>>>>>>> Leaving set_string_once", verbosity, 6)
    return use_title
end

undef("get_plot_main_string")
function get_plot_main_string(title_prefix[1] : string,
                                  title_suffix[1] : string,
                                  idx[*] : integer,
                                  switches[1] : logical)
local title_string
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering get_plot_main_string", verbosity, 8)
    if (switches@skip_default_title) then
        title_string = title_prefix
    else
        if (dimsizes(idx) .gt. 1) then
            if (idx(0) .eq. idx(1)) then
                title_string = title_prefix + models@name(idx(0)) + " "\
                                            + models@experiment(idx(0))\
                                            + title_suffix
            else
                title_string = title_prefix + models@name(idx(0)) + "/"\
                                            + models@name(idx(1))\
                                            + title_suffix
            end if
        else
            title_string = title_prefix + models@name(idx(0)) + " "\
                                        + models@experiment(idx(0))\
                                        + title_suffix
        end if
    end if
    info_output(">>>>>>>> Leaving get_plot_main_string", verbosity, 8)
    return title_string
end

undef("get_diffplot_right_string")
function get_diffplot_right_string(switches[1] : logical)
local gsn_right_string
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering get_diffplot_right_string", verbosity, 8)
    gsn_right_string = ""
    if (isatt(switches, "gsn_right_string")) then
        if(islogical(switches@gsn_right_string)) then
            gsn_right_string = ""
        else  ; Assume string
            if (switches@gsn_right_string .eq. "False") then
                gsn_right_string = ""
            else
                gsn_right_string = switches@gsn_right_string
            end if
        end if
    end if
    info_output(">>>>>>>> Leaving get_diffplot_right_string", verbosity, 8)
    return gsn_right_string
end

undef("getFillValue_with_default")
function getFillValue_with_default(data : numeric)
; Description:
;    Fetches FillValue from 'data', if not present, set it to the
;    default value
local data_FillValue, data, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering getFillValue_with_default", verbosity, 8)
    data_FillValue = getFillValue(data)
    if (data_FillValue .eq. "No_FillValue") then
        if (isinteger(data)) then
            data_FillValue = default_fillvalue("integer")
        else if (isfloat(data)) then
            data_FillValue = default_fillvalue("float")
        else if (isdouble(data)) then
            data_FillValue = default_fillvalue("double")
        else if (isshort(data)) then
            data_FillValue = default_fillvalue("short")
        end if
        end if
        end if
        end if
        data@_FillValue = data_FillValue
    end if
    info_output(">>>>>>>> Leaving getFillValue_with_default", verbosity, 8)
    return data_FillValue
end

undef("get_diffplot_left_string")
function get_diffplot_left_string(statistics[1] : logical)
;                           retrun_val[1] : string
; Arguments
;    statistics -- switch controling which statistics to display
; Returns
;    gsnLeftString -- text string with the req. stats concatenated
local gsnLeftString, present_statsitics, idx, stat_order, stat_idx
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering get_diffplot_left_string", verbosity, 8)
    present_statistics = getvaratts(statistics)
    gsnLeftString = ""

    ;; Add requested statistics
    stat_order = (/"yrs", "refA", "latB", "lonB", "mean", "rmse", "corr", "bob", "eio", "sahel", "sa34", "en34"/)
    if (.not. any(ismissing(present_statistics))) then
        if (is_array_subset(present_statistics, stat_order)) then
            do idx = 0, dimsizes(stat_order) - 1
                stat_idx = ind(stat_order(idx) .eq. present_statistics)
                if (isatt(statistics, stat_order(idx))) then
                    gsnLeftString = gsnLeftString + "~F22~" + present_statistics(stat_idx) + ":~F21~ "\
                                                  + statistics@$present_statistics(stat_idx)$ + " "
                end if
            end do
        else
            ttt = relative_complement(stat_order, present_statistics)
            all_statistics_defined = (.not. ttt)
            if (all_statistics_defined) then
                do idx = 0, dimsizes(stat_order) - 1
                    stat_idx = ind(stat_order(idx) .eq. present_statistics)
                    if (isatt(statistics, stat_order(idx))) then
                        gsnLeftString = gsnLeftString + "~F22~" + present_statistics(stat_idx) + ":~F21~ "\
                                                      + statistics@$present_statistics(stat_idx)$ + " "
                    end if
                end do
            else
                error_msg("fatal",\
                          "SAMonsoon.ncl",\
                          "get_diffplot_left_string",\
                          "There are non-defined statistics: (=" + non_defined_statistics)
                status_exit(1)
            end if
        end if
    end if

    info_output(">>>>>>>> Leaving get_diffplot_left_string", verbosity, 8)
    return gsnLeftString
end

undef("set_default_switches")
function set_default_switches(possible_switches[*] : string,
                              present_switches[1] : logical)
;                           retrun_val[1] : logical
; Arguments:
;   possible_switches  -- array with possible logical switches
;   present_switches   -- logical listing actual switches

; Description:
;   Extracts resouces that are attached to 'present_switches'
;   These resources are moved to the return logical, default value is
;   given to all other possible switches
local idx, curr_value, res_return
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering set_default_switches", verbosity, 8)
    res_return = True
    do idx = 0, dimsizes(possible_switches) - 1
        curr_value = possible_switches(idx)
        res_return@$curr_value$ = False  ; Default value
        if (isatt(present_switches, curr_value)) then
            delete(res_return@$curr_value$)
            ;; If "False"/"True", convert to logcial
            if (isstring(present_switches@$curr_value$)) then
                if (present_switches@$curr_value$ .eq. "False") then
                    delete(present_switches@$curr_value$)
                    present_switches@$curr_value$ = False
                else if (present_switches@$curr_value$ .eq. "True") then
                    delete(present_switches@$curr_value$)
                    present_switches@$curr_value$ = True
                end if
                end if
            end if
            res_return@$curr_value$ = present_switches@$curr_value$
        end if
    end do
    info_output(">>>>>>>> Leaving set_default_switches", verbosity, 8)
    return res_return
end

undef("remove_if_default_axis_string")
procedure remove_if_default_axis_string(res [1]:logical)
begin
    if(isatt(res, "tiYAxisOn")) then
        if(.not. isatt(res, "tiYAxisString")) then
            res@tiYAxisString = ""
        end if
    end if
    if(isatt(res, "tiXAxisOn")) then
        if(.not. isatt(res, "tiXAxisString")) then
            res@tiXAxisString = ""
        end if
    end if
end

undef("single_contour")
function single_contour(wks[1] : graphic,
                        plot_data[*][*] : numeric,
                        title_prefix[1] : string,
                        title_suffix[1] : string,
                        plot_settings[1] : logical,
                        debugdisplay[1] : logical,
                        statistics[1] : logical,
                        res[1] : logical)
;; Arguments:
;;    wks          - Workstation for abs value plots
;;    plot_data   - The data set to plot
;;    title_prefix - String with main title prefix
;;    title_suffix - String with main title suffix
;;    plot_settings - Logical with plot settings
;;    degubdisplay - Degug information to add to the plot
;;    statistics - Statistical information to add to the plot headers
;;    res          - Plot resource variable
;; Description:
;;    Creates a contour plot complete with headers etc..
local idx, mean, rmse, possible_switches, switches
begin
    idx = plot_settings@idx
    possible_switches = (/"skip_statistics_in_plots",\
                          "skip_long_name_in_plots",\
                          "skip_default_title",\
                          "gsn_right_string"/)
    switches = set_default_switches(possible_switches, plot_settings)

    res@tiMainString = get_plot_main_string(title_prefix,\
                                            title_suffix,\
                                            idx,\
                                            switches)
    res@gsnLeftString = get_diffplot_left_string(statistics)
    res@gsnRightString = get_diffplot_right_string(switches)
    remove_if_default_axis_string(res)

    if (isatt(plot_settings, "xcor_lagged_plot")) then
        plot = mjo_xcor_lag_ovly_local(wks, plot_data, res)
    else if (isatt(plot_settings, "plot_pres_hgt")) then
        plot = gsn_csm_pres_hgt(wks, plot_data, res)
    else
        plot = gsn_csm_contour_map_ce(wks, plot_data, res)
    end if
    end if
    return plot
end

undef("single_vector")
function single_vector(wks[1] : graphic,
                       speed_ref_max[1] : numeric,
                       plot_data_speed[*][*] : numeric,
                       plot_data_ua[*][*] : numeric,
                       plot_data_va[*][*] : numeric,
                       title_prefix[1] : string,
                       title_suffix[1] : string,
                       plot_settings[1] : logical,
                       debugdisplay[1] : logical,
                       res[1] : logical)
;; Arguments:
;;    wks          - Workstation for abs value plots
;;    plot_data   - The data set to plot
;;    title_prefix - String with main title prefix
;;    title_suffix - String with main title suffix
;;    plot_settings - Logical with plot settings
;;    degubdisplay - Degug information to add to the plot
;;    res          - Plot resource variable
;; Description:
;;    Creates a contour plot complete with headers etc..
local idx, mean, rmse, possible_switches, switches, vres, plot_v
begin
    idx = plot_settings@idx

    vres = True
    vres@vcLineArrowColor      = "black"  ; Change vector color
    vres@vcLineArrowThicknessF = 1.0      ; Change vector thickness
    vres@vcMinDistanceF        = 0.04
    vres@vcRefMagnitudeF        = speed_ref_max  ; Define vector ref mag
    vres@vcRefLengthF           = 0.045    ; Define length of vec ref
    vres@vcGlyphStyle           = "LineArrow"
    vres@vcLineArrowHeadMaxSizeF = 0.015
    vres@vcRefLengthF          = 0.1
    vres@gsnDraw               = False
    vres@gsnFrame              = False
    vres@gsnLeftString         = " "
    vres@gsnRightString         = " "

    if (res@tmXBLabelsOn) then
        vres@vcRefAnnoOn = False
    else
        vres@vcRefAnnoOn = True
        vres@vcRefAnnoOrthogonalPosF   = -1.0            ; move ref vector
        vres@vcRefAnnoString2      = " "
    end if

    plot_v = gsn_csm_vector(wks, plot_data_ua, plot_data_va, vres)
    return plot_v
end

undef("map_zoom")
function map_zoom(res, zoom)
begin
    if(isatt(zoom, "latrange")) then
        res@mpLimitMode = "LatLon"  ; Zoom in on map
        res@mpMinLatF = zoom@$"latrange"$(0)
        res@mpMaxLatF = zoom@$"latrange"$(1)
    end if
    if(isatt(zoom, "lonrange")) then
        res@mpLimitMode = "LatLon"  ; Zoom in on map
        res@mpMinLonF = zoom@$"lonrange"$(0)
        res@mpMaxLonF = zoom@$"lonrange"$(1)
    end if
    return res
end

undef("apply_transforms")
function apply_transforms(transform[1] : logical,
                          data : numeric)
;; Arguments:
;;    transform - required_transforms, a logical with the following attributes,
;;                 @apply_these) string list, (/"name1", "name2", .../)
;;                 @transform1) first transform value, e.g., "JJAS" for a season
;;                 @transform2) second transform value, e.g., (/-30, 60/)
;;                              for min/max lat
;;    data      - the data set to transform
;; Description:
;;    Loop through the string list of transforms and applies them to the data
;;    set (second argument).
;;
local data_regrid,  data_tmp,   data_trans,  dims_3d,\
    idx_trans,    newdims,      ntimes,      Opt,\
    t,            test_regrid,  verbosity,   SrcMask2D,\
    weight,       source,       destination, lon, \
    nlon
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering apply_transforms", verbosity, 6)
    data_trans = data
    dim_trans = dimsizes(transform@apply_these)
    data_dims = dimsizes(dimsizes(data))

    do idx_trans = 0, dim_trans - 1

        ;; Run lonflip first if present
        lonflip_ind = ind("lonflip" .eq. transform@apply_these)

        if(.not. ismissing(lonflip_ind)) then
            lon = data_trans&lon
            nlon = dimsizes(lon)
            if (lon(nlon-1) .gt. 180) then
                data_tmp = lonFlip(data_trans)
                delete(data_trans)
                data_trans = data_tmp
                delete(data_tmp)
            end if
            delete(lon)
            delete(nlon)
            ; Make sure lonflip isn't run again
            transform@apply_these(lonflip_ind) = "lonflip_done"
        end if

        ;; Interpolate on a rectilinear grid
        interp_ind = ind("interp" .eq. transform@apply_these)
        ;info_output("<<<<<<<< Entering apply_transforms - interp on grid " + transform@destgrid, verbosity, 6)

        if (.not. ismissing(interp_ind)) then

         ;info_output("<<<<<<<< Entering apply_transforms - interp on grid " + transform@destgrid, verbosity, 6)
         
         if (transform@destgrid .eq. "1x1") then 
           mlon = 360          ; number of output longitudes
           mlat = 180          ; number of output latitudes

           klon0 = ispan(-180,mlon-1,1)*1.
           klat0 = ispan(-90,mlat-1,1)*1.
         else if (transform@destgrid .eq. "2.5x2.5") then
           mlon = toint(360/2.5)          ; number of output longitudes
           mlat = toint(180/2.5)          ; number of output latitudes

;CR CAREFUL: diurnal cycle modification: the grid definition for 2.5x2.5 seems wrong to me
;           klon0 = ispan(-180,mlon-1,1)*2.5
;           klat0 = ispan(-90,mlat-1,1)*2.5
           klon0 = ispan(-mlon/2+1,mlon/2-1,1)*2.5
           klat0 = ispan(-mlat/2+1,mlat/2-1,1)*2.5
;RC
         else
           print("fatal:interpolation unknown:" + transform@destgrid)
           status_exit(1)
         end if
         end if

         klat = tofloat(klat0)
         klon = tofloat(klon0)

         klat@units = "degrees_north"
         klon@units = "degrees_east"

         data_tmp = linint2_Wrap(data_trans&lon,data_trans&lat,data_trans,True,klon,klat,0)
         data_tmp!0 = "time"
         data_tmp!1 = "lat"
         data_tmp!2 = "lon"

         data_tmp&time = data_trans&time
         data_tmp&lat = klat
         data_tmp&lon = klon

         delete(data_trans)
         data_trans = data_tmp
         delete(data_tmp)

         

         delete(mlon)
         delete(mlat)
         delete(klon)
         delete(klat)
         ; Make sure interpolation isn't run again
         transform@apply_these(interp_ind) = "interp_done"
        end if


        ;; Restrict data set to a certain latitude range
        if(transform@apply_these(idx_trans) .eq. "latrange") then
            if (data_dims .eq. 3) then
                data_tmp = data_trans(time|:, \
                                      {lat|transform@$"latrange"$(0):transform@$"latrange"$(1)},\
                                      lon|:)
            else if (data_dims .eq. 4) then
                data_tmp = data_trans(time|:, \
                                      plev|:, \
                                      {lat|transform@$"latrange"$(0):transform@$"latrange"$(1)},\
                                      lon|:)
            end if
            end if
            delete(data_trans)
            data_trans = data_tmp
            delete(data_tmp)

        end if

        ;; Restrict data set to a certain longitude range
        if(transform@apply_these(idx_trans) .eq. "lonrange") then
            if (data_dims .eq. 3) then
                data_tmp = data_trans(time|:,\
                                      lat|:,\
                                      {lon|transform@$"lonrange"$(0):transform@$"lonrange"$(1)})
            else if (data_dims .eq. 4) then
                data_tmp = data_trans(time|:,\
                                      plev|:,\
                                      lat|:,\
                                      {lon|transform@$"lonrange"$(0):transform@$"lonrange"$(1)})
            end if
            end if
            delete(data_trans)
            data_trans = data_tmp
            delete(data_tmp)
        end if

        ;; Restrict data set to certain pressure levels
        if(transform@apply_these(idx_trans) .eq. "levrange") then
            if (data_dims .eq. 3) then
                noop()  ; Gracefully ignore request
            else if (data_dims .eq. 4) then
                levrange = transform@$"levrange"$
                if (levrange(0).gt.levrange(1)) then
                  levrange = (/levrange(1),levrange(0)/)
                end if

                data_tmp = data_trans(time|:,\
                                      {plev|levrange(0):levrange(1)},\
                                      lat|:,\
                                      lon|:)
            end if
            end if
            delete(data_trans)
            data_trans = data_tmp
            delete(data_tmp)
        end if

        ;; Restrict data set to a certain season
        if(transform@apply_these(idx_trans) .eq. "season") then
;CR: CAREFUL DIURNAL CYCLE MODIFICATION: extract 3hr data instead of monthy data
;            data_tmp = month_to_season_extended(data_trans, transform@season)
            data_tmp = extract_season(data_trans, transform@season)
;RC
            delete(data_trans)
            data_trans = data_tmp
            delete(data_tmp)
        end if

 
        ;; Restrict data set to certain years
        if(transform@apply_these(idx_trans) .eq. "years") then
            first_year = transform@years(0)
            last_year = transform@years(dimsizes(transform@years) - 1)
            data_tmp = extract_years(data_trans, first_year, last_year)
            delete(data_trans)
            data_trans = data_tmp
            delete(data_tmp)
        end if

        ;; Fill missing values using poisson
        if(transform@apply_these(idx_trans) .eq. "fill_missing_w_poisson") then
            if (isatt(transform, "is_cyclic")) then
                is_cyclic = transform@is_cyclic
            else
                info_output("fatal: transform 'fill_missing_w_poisson' doesn't know whether data is cyclic", verbosity, 1)
                status_exit(1)
            end if

            if any(ismissing(data_trans)) then
                guess     = 1                ; use zonal means
                nscan     = 1500             ; usually much less than this
                eps       = 1.e-2            ; variable dependent
                relc      = 0.6              ; relaxation coefficient
                opt       = 0                ; not used
                poisson_grid_fill(data_trans,\
                                  is_cyclic,\
                                  guess,\
                                  nscan,\
                                  eps,\
                                  relc,\
                                  opt)
                info_output("info: using 'poisson_grid_fill(..)' to remove missing values", verbosity, 1)
            end if
        end if

        ;; Interpolate to regular grid
        if(transform@apply_these(idx_trans) .eq. "interpolation") then
            Opt                   = True
            Opt@DstLLCorner       = (/transform@$"latrange"$(0), transform@$"lonrange"$(0)/)
            Opt@DstURCorner       = (/transform@$"latrange"$(1), transform@$"lonrange"$(1)/)
            Opt@DstGridType       = transform@interpolation  ; Destination grid
            Opt@ForceOverwrite    = True
            Opt@PrintTimings      = False
            Opt@Debug             = False

            ;; "bilinear" ; "patch"  ; "conserve" fails
            Opt@InterpMethod      = "patch"

            ;; Save regrid files for later use: careful with naming!
            weight_file = interface_get_figure_filename(transform@diag_script_base,\
                                                        transform@variable,\
                                                        transform@field_type,\
                                                        "weight",\
                                                        transform@idx)
            source_file = interface_get_figure_filename(transform@diag_script_base,\
                                                        transform@variable,\
                                                        transform@field_type,\
                                                        "source",\
                                                        transform@idx)
            destination_file = interface_get_figure_filename(transform@diag_script_base,\
                                                             transform@variable,\
                                                             transform@field_type,\
                                                             "destination",\
                                                             transform@idx)
            if (.not. isfilepresent_esmval(regridding_dir)) then
                system("mkdir -p " + regridding_dir)
            end if
            Opt@WgtFileName = regridding_dir + "/" + weight_file + ".nc"
            Opt@SrcFileName = regridding_dir + "/" + source_file + ".nc"
            Opt@DstFileName = regridding_dir + "/" + destination_file + ".nc"

            Opt@printTimings = False
            ;; ESMF_regrid may issue a warning
            ;; (e.g. if Dst is not covered Src grid)
            ;; -> set "exit_on_warning" to False (in main nml) to continue to
            ;; the next diag

            ;; Setting Opt@SrcRegional=True avoids artefacts
            ;; from incomplete Src grids

            ;; Activate if only part of the sphere is given
            Opt@SrcRegional = True

            ;; reuse regrid files if present
            if (isfilepresent_esmval(Opt@SrcFileName)) then
                Opt@SkipSrcGen = True
            end if
            if (isfilepresent_esmval(Opt@DstFileName)) then
                Opt@SkipDstGen = True
            end if

            ntimes = dimsizes(data(:, 0, 0))  ; get number of time steps

            Opt@SrcMask2D = where(.not. ismissing(data(0, :, :)), 1, 0)
            ;; Only to get dim sizes
            test_regrid = ESMF_regrid(data(0, :, :), Opt)
            newdims = dimsizes (test_regrid)

            dims_3d = (/ntimes, newdims(0), newdims(1)/)
            data_regrid = new(dims_3d, float, test_regrid@_FillValue)

            do t = 0, ntimes - 1
                SrcMask2D = where(.not. ismissing(data(t, :, :)), 1, 0)
                if (any(SrcMask2D - Opt@SrcMask2D .gt. 0)) then
                    Opt@SrcMask2D = SrcMask2D
                    data_regrid(t, :, :) = ESMF_regrid(data(t, :, :), Opt)
                else
                    data_regrid(t, :, :) = ESMF_regrid_with_weights(data(t, :, :), Opt@WgtFileName, Opt)
                end if
            end do
            delete(data_trans)
            data_trans = data_regrid
            copy_VarCoords_2(data, data_trans)
        end if
    end do

    ;; Reset lonflip in case it's need again
    lonflip_ind = ind("lonflip_done" .eq. transform@apply_these)
    if(.not. ismissing(lonflip_ind)) then
        transform@apply_these(lonflip_ind) = "lonflip"
    end if

    interp_ind = ind("interp_done" .eq. transform@apply_these)
    if(.not. ismissing(interp_ind)) then
        transform@apply_these(interp_ind) = "interp"
    end if

    info_output(">>>>>>>> Leaving apply_transforms", verbosity, 6)
    return data_trans
end

undef("inset_top_text")
procedure inset_top_text(wks[1] : graphic,
                         plot[1] : graphic,
                         inset_string[1] : string,
                         res[1] : logical)
;; Arguments:
;;     wks       -- current workstation
;;     plot       -- graphic object to draw text onto
;;     insert_string  -- string to attach
;;     res       -- resource settings for string box
;; Description:
;;     Places the text string as inset at the top of the figure
local i, txvpWidth, txvpHeight, tx, x, y, tickmark_height, res
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering inset_top_text", verbosity, 6)
    res@gsnDraw               = False
    res@gsnFrame              = False
    res@txFontHeightF         = 0.012
    res@txBackgroundFillColor = 0
    res@txFuncCode            = "~"
    res@txPerimOn             = True
    tx = gsn_create_text(wks, inset_string, res)
    getvalues tx
        "vpWidthF" : txvpWidth
        "vpHeightF" : txvpHeight
    end getvalues
    getvalues plot
        "vpWidthF" : parent_width
        "vpHeightF" : parent_height
;A-laue_ax        "tmXBLabelFontHeightF" : tickmark_height
    end getvalues
    mvxf = txvpWidth / parent_width
    mvyf = txvpHeight / parent_height

    amres = True
    amres@amZone           = 2
    amres@amSide           = "Top"
    amres@amJust           = "BottomCenter"
    amres@amParallelPosF   = 0.5
    amres@amOrthogonalPosF = -mvyf
    annoid = gsn_add_annotation(plot, tx, amres)
    info_output(">>>>>>>> Leaving inset_top_text", verbosity, 6)
end

undef("multipanel_get_no_cols")
function multipanel_get_no_cols(no_dims[1] : integer,
                                max_cols[1]: integer)
;; Arguments:
;;     no_dims    -- number of models to panel
;;     max_cols   -- do not allow more columns than this
;; Description:
;;     Based on the number of models, returns the number of columns
;;     that should be used for panelling
local cols, rows
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering multipanel_get_no_cols", verbosity, 6)
    if (max_cols .eq. 3) then
        cols = 2
        rows = no_dims / tofloat(cols)
        if (rows .gt. 3) then
            cols = 3
            rows = no_dims / tofloat(cols)
        end if
        if (no_dims .le. 3) then
            cols = 1
            rows = no_dims
        end if
    else if (max_cols .eq. 2) then
        cols = 1
        rows = no_dims / tofloat(cols)
        if (rows .gt. 3) then
            cols = 2
            rows = no_dims / tofloat(cols)
        end if
        if (no_dims .le. 3) then
            cols = 1
            rows = no_dims
        end if
    else
        print("fatal: multipanel_get_no_cols(...), illegal value for max_cols=" + max_cols)
        status_exit(1)
    end if
    end if
    info_output(">>>>>>>> Leaving multipanel_get_no_cols", verbosity, 6)
    return cols
end

undef("multipanel_get_no_rows")
function multipanel_get_no_rows(no_dims[1] : integer,
                                max_cols[1]: integer)
;; Arguments:
;;     no_dims    -- number of models to panel
;;     max_cols   -- do not allow more columns than this
;;     max_rows   -- do not allow more columns than this
;; Description:
;;     Based on the number of models, returns the number of columns
;;     that should be used for panelling
local cols, rows
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering multipanel_get_no_rows", verbosity, 6)
    if (max_cols .eq. 3) then
        cols = 2
        rows = no_dims / tofloat(cols)
        if (rows .gt. 3) then
            cols = 3
            rows = no_dims / tofloat(cols)
        end if
        if (no_dims .le. 3) then
            cols = 1
            rows = no_dims
        end if
    else if (max_cols .eq. 2) then
        cols = 1
        rows = no_dims / tofloat(cols)
        if (rows .gt. 3) then
            cols = 2
            rows = no_dims / tofloat(cols)
        end if
        if (no_dims .le. 3) then
            cols = 1
            rows = no_dims
        end if
    else
        print("fatal: multipanel_get_no_cols(...), illegal value for max_cols=" + max_cols)
        status_exit(1)
    end if
    end if
    info_output(">>>>>>>> Leaving multipanel_get_no_rows", verbosity, 6)
    return rows
end

undef("rebalance_wrt_min")
procedure rebalance_wrt_min(no_figures_on_page[1] : integer,
                            no_models[1] : integer,
                            min_pp[1] : integer)
local min_diff
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering rebalance_wrt_min", verbosity, 8)
    if (no_models .gt. 0) then
        if (no_models .lt. min_pp) then
            min_diff = min_pp - no_models
            no_figures_on_page = no_figures_on_page - min_diff
            no_models = no_models + min_diff
        end if
    end if
    info_output(">>>>>>>> Leaving rebalance_wrt_min", verbosity, 8)
end

undef("get_figures_per_page")
function get_figures_per_page(total_no_models[1] : integer,
                              max_pp[1] : integer,
                              min_pp[1] : integer)
;                           return_val[*] : integer
;; Arguments:
;; total_no_models  -- number of models to panel (including reference,
;;                                                counted once)
;; max_pp   -- max number of figures per page
;; min_pp   -- min number of figures per page
;; Description:
;;     Based on the number of models and the max number of figures
;;     on each page, panel the models (make sure the reference
;;     model is repeated on each page)
local verbosity,               no_refs,     no_models,\
    no_figures_first_page,   no_models,   models_pp,\
    no_figures,              no_models,   models_pp_tmp,\
    models_pp
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering get_figures_per_page", verbosity, 6)
    no_refs = 1
    no_models = total_no_models - no_refs

    no_figures_first_page = min((/total_no_models, max_pp/))
    no_models = no_models - (no_figures_first_page - no_refs)

    rebalance_wrt_min(no_figures_first_page, no_models, min_pp)

    models_pp = (/no_figures_first_page/)

    do while (no_models .gt. 0)
        no_figures = min((/no_models + no_refs, max_pp/))
        no_models = no_models - (no_figures - no_refs)
        rebalance_wrt_min(no_figures, no_models, min_pp)

        models_pp_tmp = array_append_record(models_pp, (/no_figures/), 0)
        delete(models_pp)
        models_pp = models_pp_tmp
        delete(models_pp_tmp)
    end do

    info_output(">>>>>>>> Leaving get_figures_per_page", verbosity, 6)
    return models_pp
end

undef("place_models_on_pages")
procedure place_models_on_pages(models[1] : logical,
                                idx_ref[*] : integer,
                                figures_per_page[*] : integer,
                                model_placement_by_idx[*] : integer,
                                figure_panel_idx[*] : integer)
;; Description:
;;   Maps the model indices onto each panel plot page
local curr, curr_mod
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering place_models_on_pages", verbosity, 6)
    curr = 0
    curr_mod = 0
    do j = 0, dimsizes(figures_per_page) - 1
        model_placement_by_idx(curr) = idx_ref(0)
        figure_panel_idx(curr) = 0
        curr = curr + 1
        do k = 1, figures_per_page(j) - 1
            figure_panel_idx(curr) = k
            ;; Assume there is only one ref-model
            if (any(curr_mod .eq. idx_ref)) then
                curr_mod = curr_mod + 1
                model_placement_by_idx(curr) = curr_mod
                curr = curr + 1
                curr_mod = curr_mod + 1
                continue
            end if
            model_placement_by_idx(curr) = curr_mod
            curr = curr + 1
            curr_mod = curr_mod + 1
        end do
    end do
    info_output(">>>>>>>> Leaving place_models_on_pages", verbosity, 6)
end

undef("get_latlon_string")
function get_latlon_string(min_lat[1] : numeric,
                           max_lat[1] : numeric,
                           min_lon[1] : numeric,
                           max_lon[1] : numeric)
;; Description:
;;   Converts min/max lat/lon to a string
local min_lat_str, max_lat_str, min_lon, max_lon_str, latlon_string
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering get_latlon_string", verbosity, 6)
    min_lat_str = sprinti("%i", min_lat)
    max_lat_str = sprinti("%i", max_lat)
    min_lon_str = sprinti("%i", min_lon)
    max_lon_str = sprinti("%i", max_lon)
    latlon_string = "~F22~lat:~F21~ " + min_lat_str + "-" + max_lat_str + ", "\
                  + "~F22~lon:~F21~ " + min_lon_str + "-" + max_lon_str
    info_output(">>>>>>>> Leaving get_latlon_string", verbosity, 6)
    return latlon_string
end

undef("fetch_first_last_year")
function fetch_first_last_year(dataset[*] : numeric)
;; Description:
local dataset_size, first_year, last_year
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering fetch_first_last_year", verbosity, 6)
    dataset_years = fetch_years(dataset)
    years_size = dimsizes(dataset_years)
    first_year = dataset_years(0)
    last_year = dataset_years(years_size - 1)
    info_output(">>>>>>>> Leaving fetch_first_last_year", verbosity, 6)
    return (/first_year, last_year/)
end

undef("remove_entry_from_set")
function remove_entry_from_set(idx_to_remove[1] : integer,
                               other_idxs[*] : integer)
;; Description:
;;   A set operator
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering remove_entry_from_set", verbosity, 6)
    info_output(">>>>>>>> Leaving remove_entry_from_set", verbosity, 6)
    return other_idxs(ind(other_idxs .ne. idx_to_remove))
end

undef("list_attributes")
function list_attributes(attribute_holder)
; Description:
;   Copies the attribute names from variable 'attribute_holder'
;   onto 'attribute_list'
local attributes, idx, attribute_list, attstring
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering list_attributes", verbosity, 8)
    attributes = getvaratts(attribute_holder)

    attribute_list = True
    do idx = 0, dimsizes(attributes) - 1
        attstring = attributes(idx)
        ;; Handle '_FillValue' problem due to automatic casting
        if (attstring .eq. "_FillValue") then
            attstring = "orig_FillValue"
        end if
        attribute_list@$attstring$ = attribute_holder@$attributes(idx)$
    end do
    info_output(">>>>>>>> Leaving list_attributes", verbosity, 8)
    return attribute_list
end

undef("attach_attributes")
function attach_attributes(variable_holder[1] : logical,
                           attachment_base[1] : string,
                           attachment[1] : logical)
local idx, attachment_atts
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering attach_attributes", verbosity, 8)
    attachment_atts = getvaratts(attachment)
    do idx = 0, dimsizes(attachment_atts) - 1
        attachment_name = attachment_base + "_" + attachment_atts(idx)
        variable_holder@$attachment_name$ = attachment@$attachment_atts(idx)$
    end do
    info_output(">>>>>>>> Leaving attach_attributes", verbosity, 8)
    return variable_holder
end

undef("count_number_of_entries")
function count_number_of_entries(entries[*],
                                 matching_substring[1] : string)
local idx, is_match, no_entries
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering count_number_of_entries", verbosity, 8)
    no_entries = 0
    ;; Count
    do idx = 0, dimsizes(entries) - 1
        is_match = str_index_of_substr(entries(idx), matching_substring, -1)
        if (.not. ismissing(is_match)) then
            no_entries = no_entries + 1
        end if
    end do
    info_output(">>>>>>>> Leaving count_number_of_entries", verbosity, 8)
    return no_entries
end

undef("get_entry_indicies")
function get_entry_indicies(no_entries[1] : integer,
                            variable_holder[1] : logical,
                            matching_substring[1] : string)
local varatts, entry_idx, idx, is_match, curr_entry,
    last_underscore, str_index, int_index, index_array
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering get_entry_indicies", verbosity, 8)
    varatts = getvaratts(variable_holder)
    index_array = new((/no_entries, 2/), integer)

    entry_idx = 0
    do idx = 0, dimsizes(varatts) - 1
        is_match = str_index_of_substr(varatts(idx), matching_substring, -1)
        if (.not. ismissing(is_match)) then
            curr_entry = varatts(idx)
            last_underscore = str_index_of_substr(curr_entry, "_", -1)
            str_index = str_get_cols(curr_entry, last_underscore + 1, -1)
            int_index = stringtointeger(str_index)
            index_array(entry_idx, 0) = idx
            index_array(entry_idx, 1) = int_index
            entry_idx = entry_idx + 1
        end if
    end do
    info_output(">>>>>>>> Leaving get_entry_indicies", verbosity, 8)
    return index_array
end

undef("retrieve_storage_record")
function retrieve_storage_record(variable_holder[1] : logical,
                                 storage_label[1] : string)
local varatts, idx, no_entries, storage_array, tmp, i_perm
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering retrieve_storage_record", verbosity, 8)
    varatts = getvaratts(variable_holder)
    no_entries = count_number_of_entries(varatts, storage_label)
    if (no_entries .gt. 0) then
        i_perm = get_entry_indicies(no_entries, variable_holder, storage_label)
        storage_array = variable_holder@$varatts(i_perm(i_perm(0, 1), 0))$

        do idx = 1, no_entries - 1
            tmp = array_append_record(storage_array,\
                                      variable_holder@$varatts(i_perm(i_perm(idx, 1), 0))$, 0)
            delete(storage_array)
            storage_array = tmp
            delete(tmp)
        end do
    else
        storage_array = 1e+20
        storage_array@_FillValue = 1e+20
    end if
    info_output(">>>>>>>> Leaving retrieve_storage_record", verbosity, 8)
    return storage_array
end

undef("add_storage_record")
procedure add_storage_record(variable_holder[1] : logical,
                             storage_label[1] : string,
                             storage[*])
local no_entries, idx, existing_entries
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering add_storage_record", verbosity, 6)
    existing_entries = retrieve_storage_record(variable_holder,\
                                               storage_label)
    if (.not. any(ismissing(existing_entries))) then
        extended_storage = array_append_record(existing_entries, storage, 0)
    else
        extended_storage = storage
    end if

    no_entries = dimsizes(extended_storage)
    do idx = 0, no_entries - 1
        storage_string = storage_label + "_" + sprinti("%i", idx)
        variable_holder@$storage_string$ = extended_storage(idx)
    end do
    info_output(">>>>>>>> Leaving add_storage_record", verbosity, 6)
end

undef("prefix_filter")
function prefix_filter(str_prefix[1] : string,
                       array_string[*] : string)
local idx, mathcing_attribute, return_array, tmp_array
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering prefix_filter", verbosity, 8)
    do idx = 0, dimsizes(array_string) - 1
        matching_string = str_index_of_substr(array_string(idx), str_prefix, -1)
        if (.not. ismissing(matching_string)) then
            if (isdefined("return_array")) then
                tmp_array = array_append_record(return_array, (/array_string(idx)/), 0)
                delete(return_array)
                return_array = tmp_array
                delete(tmp_array)
            else
                return_array = (/array_string(idx)/)
            end if
        end if
    end do
    if (.not. isdefined("return_array")) then
        print("fatal: Couldn't find attribute '" + str_prefix + "' in SAMonsoon.ncl (prefix_filter())")
        status_exit(1)
    end if
    info_output(">>>>>>>> Leaving prefix_filter", verbosity, 8)
    return return_array
end

undef("get_stored_attribute_names")
function get_stored_attribute_names(vault[1] : logical,
                                    att_prefix[1] : string)
local idx, varatts, requested_atts
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering get_stored_attribute_names", verbosity, 8)
    varatts = getvaratts(vault)
    requested_atts = prefix_filter(att_prefix, varatts)
    info_output(">>>>>>>> Leaving get_stored_attribute_names", verbosity, 8)
    return requested_atts
end

undef("find_dim")
function find_dim(variable, req_dim)
local vardims
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering find_dim", verbosity, 8)
    vardims = getvardims(variable)
    info_output(">>>>>>>> Leaving find_dim", verbosity, 8)
    return ind(vardims .eq. req_dim)
end

undef("flatten_to_string")
function flatten_to_string(array_string[*] : string)
local idx, return_string
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering flatten_to_string", verbosity, 8)
    return_string = ""
    do idx = 0, dimsizes(array_string) - 1
        return_string = return_string + idx + "=" + array_string(idx) + "___"
    end do
    info_output(">>>>>>>> Leaving flatten_to_string", verbosity, 8)
    return return_string
end

undef("string_to_array")
function string_to_array(model_names[1] : string)
local models, strings, idx, model
begin
    strings = str_split(model_names, "___")
    models = new(dimsizes(strings), string)
    do idx = 0, dimsizes(strings) - 1
        model = str_split(strings(idx), "=")
        models(idx) = model(1)
    end do
    return models
end

undef("retrieve_from_vault")
function retrieve_from_vault(storage_vault[1] : logical,
                             storage_record[1] : string)
;; Arguments:
;; storage_vault  -- "pointer" to storage vault
;; storage_record -- name of the variable to fetch from vault
;; Description:
;;     Fetches the full entry 'storage_record' from the
;;     'storage_vault'-"pointer"
local varatts, idx, vardims, variable_dim_name
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering retrieve_from_vault", verbosity, 6)
    ;; Get variable values
    return_var = storage_vault@$storage_record$

    ;; Get variable dimensions
    vardims_string = storage_record + "_dimensions"
    vardims = retrieve_storage_record(storage_vault,\
                                      vardims_string)



    do idx = 0, dimsizes(vardims) - 1
        dim_name = vardims(idx)
        vardim_string = storage_record + "_" + dim_name
        curr_dim = storage_vault@$vardim_string$

        ;; Get variable dimensions attributes
        vardimatts_string = storage_record + "_" + dim_name + "_attributes"
        vardimatts = retrieve_storage_record(storage_vault,\
                                             vardimatts_string)
        do idx_att = 0, dimsizes(vardimatts) - 1
            vardimatt_string = storage_record + "_"\
                             + dim_name + "_"\
                             + vardimatts(idx_att)
            if (vardimatts(idx_att) .eq. "orig_FillValue") then
                vardimatts(idx_att) = "_FillValue"
            end if
            curr_dim@$vardimatts(idx_att)$ =  storage_vault@$vardimatt_string$
        end do
        return_var!idx = dim_name
        return_var&$dim_name$ = curr_dim
        delete(curr_dim)
        delete(vardimatts)
    end do

    varatts_string = storage_record + "_attributes"
    varatts = retrieve_storage_record(storage_vault,\
                                      varatts_string)
    do idx = 0, dimsizes(varatts) - 1
        varatt_string = storage_record + "_" + varatts(idx)
        if (varatts(idx) .eq. "orig_FillValue") then
            varatts(idx) = "_FillValue"
        end if
        return_var@$varatts(idx)$ = storage_vault@$varatt_string$
    end do
    info_output(">>>>>>>> Leaving retrieve_from_vault", verbosity, 6)
    return return_var
end

undef("vault_to_netcdf")
procedure vault_to_netcdf(storage_vault[1] : logical,\
                          diag_script_base [1] : string,\
                          var [1] : string,\
                          field [1] : string)
local filehandle,     filename,     idx,\
      idx_mod,        model_dim,    model_name,\
      models_string,  netcdf_dir,   no_entries,\
      stored_models,  stored_name,  stored_variable_names,\
      tmp,            varname,      verbosity,\
      ref,            no_models
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering vault_to_netcdf", verbosity, 6)

    stored_variable_names = get_stored_attribute_names(storage_vault,\
                                                       "stored_variables")

    netcdf_dir = plot_dir + "/netcdf/"
    system("mkdir -p " + netcdf_dir)

    do idx = 0, dimsizes(stored_variable_names) - 1
        stored_name = stored_variable_names(idx)
        varname = storage_vault@$stored_name$

        tmp = retrieve_from_vault(storage_vault, varname)

        if (isatt(tmp, "idx")) then
            idx_mod = tmp@idx
            model_name = tmp@model_name
            if (model_name .ne. models@name(idx_mod)) then
                print("fatal: (vault_to_netcdf) -- inconsistent model name in ncl.interface and netcdf vault")
                status_exit(1)
            end if
        else
            idx_mod = -1
        end if

        filename = interface_get_figure_filename(diag_script_base,\
                                                 var,\
                                                 field,\
                                                 varname,\
                                                 idx_mod)
        filename = netcdf_dir + "/" + filename + ".nc"
        if (isfilepresent_esmval(filename)) then
            system("rm -f " + filename)
        end if
        filehandle = addfile(filename, "c")

        ;; Work-around:
        ;; 'model' is a variable string dimension and doesn't fit
        ;; in the netCDF data model
        if (isdim(tmp, "model")) then
            model_dim = find_dim(tmp, "model")
            no_entries = dimsizes(tmp&model)
            models_string = flatten_to_string(tmp&model)
            delete(tmp&model)
            tmp@model_names = models_string
        end if

        filehandle->$varname$ = tmp
        delete(tmp)
    end do
    info_output(">>>>>>>> Leaving vault_to_netcdf", verbosity, 6)
end

undef("store_variable")
procedure store_variable(variable : numeric,
                         storage_vault[1] : logical,
                         variable_name[1] : string)
;; Arguments:
;; variable         -- the actual variable to store
;; storage_vaulti   -- logical "pointer" to acces stored variables
;; variable_name    -- string name of the variable to add
;; Description:
;;     Maps 'variable' onto  various 'storage_vault' attributes. These
;;     attributes are prefixed 'variable_name', e.g., global attributes as
;;        storage_vault@variable_name_attribute_name_1 = attribute1
;;        storage_vault@variable_name_attribute_name_2 = attribute2
;;        ...
;;        storage_vault@variable_name_dimensions_1 = dimension1
;;        storage_vault@variable_name@dimensions_1_attribute_1 = dimension1_attribute1
local varatts, idx, vardims, variable_dim_name
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering store_variable", verbosity, 6)

    ;; Add an entry to the list of variables stored in 'storage_vault'
    add_storage_record(storage_vault, "stored_variables", variable_name)

    ;; Copy and store global attributes
    storage_vault@$variable_name$ = variable
    varatts = list_attributes(variable)
    storage_vault = attach_attributes(storage_vault,\
                                      variable_name,\
                                      varatts)

    ;; Add entries to the list with varible attributes
    add_storage_record(storage_vault, variable_name + "_attributes", getvaratts(varatts))

    ;; ----------------------------------------
    ;; Copy dimensions and dimension attributes
    ;; ----------------------------------------
    vardims = getvardims(variable)

    ;; Add entries to the list with varible dimensions
    add_storage_record(storage_vault, variable_name + "_dimensions", vardims)

    do idx = 0, dimsizes(vardims) - 1
        ;; Add dimensions
        variable_dim_name = variable_name + "_" + vardims(idx)
        storage_vault@$variable_dim_name$ = variable&$vardims(idx)$
        delete(varatts)

        ;; Add dimension attributes
        varatts = list_attributes(variable&$vardims(idx)$)
        storage_vault = attach_attributes(storage_vault,\
                                          variable_dim_name,\
                                          varatts)
        add_storage_record(storage_vault,\
                           variable_dim_name + "_attributes",\
                           getvaratts(varatts))
    end do
    info_output(">>>>>>>> Leaving store_variable", verbosity, 6)
end

undef("netcdf_to_vault")
function netcdf_to_vault(storage_vault[1] : logical,
                         varname[1] : string,
                         diag_script_base[1] : string,
                         var[1] : string,
                         field[1] : string,
                         idx_mod[1] : integer)
local fhandle,       filename,  model_dim,\
      models_array,  nc_var,    netcdf_dir,\
      verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering netcdf_to_vault", verbosity, 6)
    netcdf_dir = plot_dir + "/netcdf/"

    filename = interface_get_figure_filename(diag_script_base,\
                                             var,\
                                             field,\
                                             varname,\
                                             idx_mod)

    filename = netcdf_dir + "/" + filename + ".nc"
    fhandle = addfile(filename, "r")
    nc_var = fhandle->$varname$

    ;; Work-around:
    ;; 'model' is a variable string dimension and doesn't fit
    ;; in the netCDF data model
    if (isdim(nc_var, "model")) then
        model_dim = find_dim(nc_var, "model")
        models_array = string_to_array(nc_var@model_names)
        nc_var&model = models_array
    end if

    store_variable(nc_var, storage_vault, varname)

    info_output(">>>>>>>> Leaving netcdf_to_vault", verbosity, 6)
    return storage_vault
end


undef("band_pass_time_lat")
function band_pass_time_lat(x[*][*][*]:numeric,\
                            spd[1]:numeric,\
                            bpf[3]:numeric,\
                            opt:logical)

local nMsg, dimx, ntim, nlat, mlon, dNam, n, saveFillV,\
    bpfStrt, bpfLast, bpfNwgt, fca, fcb, ihp, sigma,\
    bpfWgt, wgtx, work, WORK, cf, WCF, wx
begin
    tmp_lon = x&lon
    dims_lon = dimsizes(tmp_lon)
    wx = new(dims_lon, float)
    wx(:) = 1.0

    nMsg = num(ismissing(x))            ; error check
    if (nMsg .gt. 0) then
        print("fatal: band_pass_time_lat: currently, missing data not allowed: nMsg=" + nMsg)
        status_exit(1)
    end if

    dimx = dimsizes(x)
    dimwx = dimsizes(wx)                ; size wx
    if (dimwx .gt. 1 .and. dimx(2) .ne. dimwx(0)) then
        print("fatal: band_pass_time_lat: sizes of y/lat dimension do not match")
        print("fatal:                       dimsizes(wx)=" + dimwx)
        print("fatal:                       dimx(2)     =" + dimx(2))
        status_exit(1)
    end if

    if (dimwx .eq. 1) then                  ; scalar
        wgtx = new(dimwx, typeof(wx), getFillValue(wx))
    end if
    wgtx = wx

    ntim = dimx(0)
    nlat = dimx(1)
    mlon = dimx(2)

    bpfStrt = bpf(0)                     ; days
    bpfLast = bpf(1)
    bpfNwgt = bpf(2)                     ; effective # weights

    bpfNwgt = (bpfNwgt / 2) * 2 + 1          ; make sure it is odd

    if (bpfStrt .gt. bpfLast) then         ; safety check
        dumy = bpfLast
        bpfLast = bpfStrt
        bpfStrt = dumy
    end if

    fca = 1.0 / (spd * bpfLast)             ; frq corresponding to
    fcb = 1.0 / (spd * bpfStrt)             ; bpfStrt and bpfLast

    dNam = getvardims(x)                 ; get dimension names
    do n = 0, 2                             ; only used if detrending
        if (ismissing(dNam(n))) then
            x!n = "dim" + n
            dNam = x!n
        end if
    end do

    ;;; Weighted average at each (lat, time)
    work = dim_avg_wgt_n_Wrap(x, wgtx, 0, 2)

    work = work - avg(work)              ; Remove overall mean

    ;; By default ... no detrend in time
    if (opt .and. (isatt(opt, "detrend") .and. opt@detrend)) then
        if (isatt(work, "_FillValue")) then
            saveFillV = work@_FillValue  ; avoid annoying warning msg
            delete(work@_FillValue)
        end if

        work = dtrend(work, False)     ; detrend in time

        if (isvar("saveFillV")) then
            work@_FillValue = saveFillV  ; reassign
        end if
        work@detrend = "data detrended in time"
    end if

    ihp = 2                                 ; bpf=>band pass filter
    sigma = 1.0                             ; Lanczos sigma
    bpfwgt = filwgts_lanczos(bpfNwgt, ihp, fca, fcb, sigma)

    work = wgt_runave_n_Wrap(work, bpfwgt, 0, 0)  ; apply filter to time

    work@band_pass_start = bpfStrt
    work@band_pass_last = bpfLast
    work@band_pass_Nwgts = bpfNwgt

    return(work)  ; (time, lon)
end

undef("set_missing_years")
procedure set_missing_years(data_array[*][*]: numeric,
                            data_years[2] : numeric,
                            all_years[*] : numeric,
                            idx[1] : integer)
local missing_first_years, missing_last_years
begin
    missing_first_years = ind(data_years(0) .gt. all_years)
    missing_last_years = ind(data_years(1) .lt. all_years)
    if (.not. any(ismissing(missing_first_years))) then
        data_array(idx, missing_first_years) = default_fillvalue(typeof(data_array))
    end if
    if (.not. any(ismissing(missing_last_years))) then
        data_array(idx, missing_last_years) = default_fillvalue(typeof(data_array))
    end if
end

undef("time_xcor_lag")
function time_xcor_lag(xio[*]:numeric, y[*][*]:numeric,\
                       dates[*]:integer, mxlag[1]:integer,\
                       opt[1]:logical)

; driver to calculate mean seasonal lags
local dimy,       iLast,      iSea,\
    iStrt,      last_sept,  mmdd,\
    mmddStart,  my,         nDay,\
    nLag,       ns,         nSea,\
    ntim,       ntmy,       oct_idx,\
    rLag,       rtyp,       rz,\
    yyyy
begin
    ntim = dimsizes(xio)
    dimy = dimsizes(y)
    ntmy = dimy(0)
    my   = dimy(1)

    if (ntim .ne. ntmy) then
        print("fatal: mjo_xcor_lag: time dimensions do not match")
        print("fatal:               ntim=" + ntim + "  ntmy=" + ntmy)
        status_exit(1)
    end if

    ;; JJAS hardcoded, handles non-365 calendars
    mmddStart = 601
    yyyy = dates / 10000
    mmdd = dates - yyyy * 10000
    oct_idx = ind(mmdd.eq. 1001)
    last_sept = oct_idx - 1

    nSea = num(mmdd .eq. mmddStart)             ; number of 'seaNam' seasons
    iSea = ind(mmdd .eq. mmddStart)
    nDay = last_sept - iSea + 1

    rtyp = "float"
    if (typeof(xio) .eq. "double" .or. typeof(y) .eq. "double") then
        rtyp = "double"
    end if

    nLag = 2 * mxlag + 1
    rz   = new((/nLag, my/), rtyp, getFillValue(y))
    rz   = 0.0

    do ns = 0, nSea - 1                           ; seasonal averaging
        iStrt = iSea(ns)
        iLast = iSea(ns) + nDay(ns) - 1
        if (iLast .gt. (ntim - 1)) then
            break
        end if

        rLag = mjo_xcor_lag_season(xio(iStrt:iLast), y(iStrt:iLast, :), mxlag, opt)
        rLag = where(rLag .gt. 0.999, 0.999, rLag)   ; avoid rounding issue
                                                     ; so (1 - rLag) .ne. 0.0

        rz    = rz + 0.5 * log((1 + rLag) / (1 - rLag))  ; Fischer z-transform
    end do
    rz = rz / nSea                           ; mean Fischer-z
    rz = (exp(2 * rz) - 1) / (exp(2 * rz) + 1)       ; transform back

    copy_VarMeta(rLag, rz)
    return(rz)
end

undef("mjo_xcor_lag_ovly_local")
function mjo_xcor_lag_ovly_local(wks [1]:graphic, ccr_a[*][*]:numeric,\
                            opt[1]:logical)
local pltPath, pltTypeLocal, res1, res2, CCR1, CCR2,\
    plt1, plt2, resP, plot
begin
    if (opt .and. isatt(opt, "colorTable") .or. isatt(opt, "pltColorTable")) then
        if (isatt(opt, "colorTable")) then
            gsn_define_colormap(wks, opt@colorTable)
        else
            gsn_define_colormap(wks, opt@pltColorTable)
        end if
    else
        gsn_define_colormap(wks, "BlWhRe")
    end if

    opt@cnLinesOn            = True
    opt@cnLevelSelectionMode = "ExplicitLevels"  ; set manual contour levels

    opt@cnLineLabelsOn       = True
    opt@gsnContourZeroLineThicknessF = 0.     ; Eliminate 0 line
    opt@gsnContourNegLineDashPattern = 1      ; negative contours dash pattern

    opt@lbLabelBarOn         = False           ; turn off individual cb's
    opt@tiYAxisString        = ""
    if (isatt(opt, "tmYLLabelsOn")) then
        if (opt@tmYLLabelsOn) then
            opt@tiYAxisString        = "lag (days)"
        end if
    end if

    plot = new(1, "graphic")

    CCR1 = ccr_a             ; possible smooth and delete of attribute

    if (opt .and. isatt(opt, "smth9") .and. abs(opt@smth9) .eq. 0.25) then
        CCR1 = smth9(CCR1, 0.50, opt@smth9, False)
    end if
    delete(CCR1@long_name)
    plt1 = gsn_csm_contour(wks, CCR1, opt)       ; contour the variable

    return plt1
end

; =====================================
undef("clmDayTLL_calendar")
function clmDayTLL_calendar(x[*][*][*]:numeric, yyyyddd:integer)
; Re-coded version of the $NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl
; version to handle 360_day-calendars
;
;
; calculate the mean Annual Cycle from daily data.
; The return array will gave the raw climatology at each grid point
;
;              x(time, lat, lon)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "lat"
;              x!2 = "lon"
;
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905
;
;    Leap
;              1908001 =>  Jan  1, 1908]
;              1908032 =>  Feb  1, 1908]
;              1908059 =>  Feb 28, 1908]
;              1908060 =>  Feb 29, 1908]
;              1908061 =>  Mar  1, 1908]
;              1908366 =>  Dec 31, 1908]
;
; Usage:  xClmDay = clmDAY_TLL (x, yyyyddd)
; -------

local dimx, ntim, nlat, mlon, ndys, days, clmDay, ndy, indx, year_day, nFill
begin
    dimx = dimsizes(x)

    ntim = dimx(0)
    nlat = dimx(1)
    mlon = dimx(2)
    if (isatt(yyyyddd, "calendar")) then
        if (yyyyddd@calendar .eq. "360" .or.yyyyddd@calendar .eq. "360_day") then
            ndys = 360
        else
        ndys = 366                                     ; allow for leap year
        end if
    else
        ndys = 366                                     ; allow for leap year
    end if

    days = yyyyddd - (yyyyddd / 1000) * 1000           ; strip year info[yyyy]

    ;; Daily climatology
    clmDay = new((/ndys, nlat, mlon/), typeof(x), getFillValue(x))
;
; Compute averages for each sequential day of the year.
; This uses dimension swapping.
;
    do ndy = 0, ndys - 2                   ; ndy=0->364 ==> day of year 1->365
        indx = ind(days .eq. (ndy + 1))     ; indx:indx Mar 2012, handle 1 element
        if (.not. ismissing(indx(0))) then
            nindx = dimsizes(indx)
            if (nindx .eq. 1) then                       ; force 3rd dimension
            ;; clmDay(ndy,:,:) = dim_avg(x(lat|:, lon|:, time|indx:indx))
            ;; 5 Nov 2009 => dim_avg_n
                clmDay(ndy, :, :) = dim_avg_n(x(indx:indx, :, :), 0)
            else
                ;clmDay(ndy,:,:) = dim_avg(x(lat|:, lon|:, time|indx))
                clmDay(ndy, :, :) = dim_avg_n(x(indx, :, :), 0)  ; 5 Nov 2009 => dim_avg_n
            end if
        end if
        delete(indx)
    end do
                                                     ; nominal day 366
                                                     ; ave(31 Dec + 1 Jan)=leap
    clmDay(ndys - 1, :, :) = (clmDay(0, :, :) + clmDay(ndys - 2, :, :)) * 0.5

    nFill = num(ismissing(clmDay))
    if (nFill .eq. 0) then
        delete(clmDay@_FillValue)
    end if

    clmDay@long_name   = "Daily Climatology"
    if (isatt(x, "long_name")) then
        clmDay@long_name = clmDay@long_name + ": " + x@long_name
    end if
    if (isatt(x, "units")) then
        clmDay@units     = x@units
    end if
    clmDay@information = "Raw daily averages across all years"
    clmDay@smoothing   = "None"

    year_day           = ispan(1, ndys, 1)
    year_day@long_name = "day of year"
    year_day@units     = "ddd"

    clmDay!0           = "year_day"
    clmDay&year_day    =  year_day

    copy_VarCoords(x(0, :, :), clmDay(0, :, :))   ; trick
    delete(clmDay@year_day)                   ; clean up
    return (clmDay)
end

; =====================================
undef("smthClmDayTLL_calendar")
function smthClmDayTLL_calendar(clmDay[*][*][*]:numeric, nHarm:integer)
;
; Re-coded version of the $NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl
; version to handle 360_day-calendars
;
local nFill, dn, z, cf, clmDaySmth
begin

    if (isatt(clmDay, "_FillValue")) then
        nFill = num(ismissing(clmDay))
        if (nFill .gt. 0) then
            print("fatal: smthClmDay_TLL: No missing values allowed: ezfftf does not allow")
            print("fatal: smthClmDay_TLL: nFill=" + nFill)
            status_exit(1)
        end if
    end if

    dn = getvardims(clmDay)                   ; get dimension names
    if (dn(0) .ne. "year_day") then
        print("warning: smthClmDayTLL: Warning: Usually expect year_day to be the dimension name")
    end if

    ; reorder make time fastest varying dimension
    z = clmDay($dn(1)$|:, $dn(2)$|:, $dn(0)$|:)
    cf = ezfftf(z)                   ; [2] x[nlat] x[mlon] x[183]

                                       ; remember NCL is 0-based
                                       ; cf(:, 0:nHarm-1) are retained unaltered
    cf(:, :, :, nHarm) = 0.5 * cf(:, :, :, nHarm)  ; mini-taper
    cf(:, :, :, nHarm + 1:) = 0.0                  ; set all higher coef to 0.0

    z = ezfftb(cf, cf@xbar)                 ; reconstructed series

    clmDaySmth = z($dn(0)$|:, $dn(1)$|:, $dn(2)$|:)
    clmDaySmth@information = "Smoothed daily climatological averages"
    clmDaySmth@smoothing   = "FFT: " + nHarm + " harmonics were retained."

    return(clmDaySmth)
end

; =====================================
undef("calcDayAnomTLL_calendar")
function calcDayAnomTLL_calendar(x[*][*][*]:numeric,\
                        yyyyddd:integer,\
                        clmDay[*][*][*]:numeric)

; Re-coded version of the $NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl
; version to handle 360_day-calendars
;
; Remove the Annual Cycle from daily data.
; On return x will consist of deviations from each day's long term mean.
;
;              x(time, lat, lon)  <==== input dimension order
;              x!0 = "time"     <==== time is in days
;              x!1 = "lat"
;              x!2 = "lon"
;
;
;    non-Leap  yyyyddd
;              1905001 =>  Jan  1, 1905
;              1905032 =>  Feb  1, 1905
;              1905059 =>  Feb 28, 1905
;              1905060 =>  Mar  1, 1905
;              1905365 =>  Dec 31, 1905
;
;    Leap
;              1908001 =>  Jan  1, 1908]
;              1908032 =>  Feb  1, 1908]
;              1908059 =>  Feb 28, 1908]
;              1908060 =>  Feb 29, 1908]
;              1908061 =>  Mar  1, 1908]
;              1908366 =>  Dec 31, 1908]

; Usage:  xAnom = calcDayAnomTLL (x, yyyyddd, clmDay)

local dimx, ntim, nlat, mlon, ndys, days, xAnom, nt
begin

    dimx = dimsizes(x)

    ntim = dimx(0)
    nlat = dimx(1)
    mlon = dimx(2)
    if (isatt(yyyyddd, "calendar")) then
        if (yyyyddd@calendar .eq. "360" .or.yyyyddd@calendar .eq. "360_day") then
            ndys = 360
        else
        ndys = 366                                     ; allow for leap year
        end if
    else
        ndys = 366                                     ; allow for leap year
    end if

    days = yyyyddd - (yyyyddd / 1000) * 1000      ; strip year info[yyyy]

; quick check[debug]
  ;if ( .not. all(days(0:ndys-2) .eq. clmDay&year_day(0:ndys-2))) then
  ;    print("calcDayAnomTLL: day mismatch")
  ;end if

; loop on every day and compute difference.

    xAnom = (/x/)                                 ; create xAnom
    do nt = 0, ntim - 1
        ;; -1 for 0-based subscript
        xAnom(nt, :, :) = x(nt, :, :) - clmDay(days(nt) - 1, :, :)
    end do

    if (isatt(x, "long_name")) then
        xAnom@long_name = "Anomalies: " + x@long_name
    else
        xAnom@long_name = "Anomalies from Daily Climatology"
    end if
    if (isatt(x, "units")) then
        xAnom@units = x@units
    end if

    copy_VarCoords(x, xAnom)

    return(xAnom)
end

undef("write_tuple_with_dash")
function write_tuple_with_dash(tuple [2]:numeric)
begin
    return tostring(tuple(0)) + "-" + tostring(tuple(1))
end

undef("no_unique_vault_entries")
function no_unique_vault_entries(storage_vault[1]:logical,
                                 ref_name [1]:string,
                                 ref_entry [1]:integer)
;                   returns no_entries [1]:integer
; Parameters:
;    storage_vault - all processed data fields attached as attributes
;    ref_name      - fetch this attribute name from the storage vault
;    ref_entry     - for the given ref_name, use this entry
; Description:
;    uses the first part of a storage vault string, e.g. the "mean" of
;    "mean-0:model_name" to count the number of model entries in the
;    storage vault
local idx,         long_names,  no_entries,\
    no_var_fields, ref_var,     short_names,\
    split_name,    upper_limit, stored_models
begin
    storage_names = get_stored_attribute_names(storage_vault, "model")

    long_names = str_get_field(storage_names, 1, ":")  ; array with (/"mean-0:model_name", mean-1...
    no_var_fields = str_fields_count(long_names, "-")
    short_names = new(dimsizes(long_names), string)  ; soon: (/"mean", "mean", ...
    do idx = 0, dimsizes(long_names) - 1
        split_name = str_split(long_names(idx), "-")
        upper_limit = no_var_fields(idx) - 2
        short_names(idx) = str_join(split_name(0:upper_limit), "-")
        delete(split_name)
    end do
    ref_var = short_names(ref_entry)  ; Pick one of the shortened entries as reference
    no_entries = dimsizes(ind(ref_var .eq. short_names))
    return no_entries
end

undef("use_speed_zonal_or_meridional")
function use_speed_zonal_or_meridional(switch[1]:string)
local use_for_contour
begin
    if (switch .eq. "speed") then
        use_for_contour = "mean"
    else if (switch .eq. "zonal") then
        use_for_contour = "mean-ua"
    else if (switch .eq. "meridional") then
        use_for_contour = "mean-va"
    else
        error_msg("fatal",\
                  "SAMonsoon.ncl",\
                  "use_speed_zonal_or_meridional()",\
                  "'switch' not properly defined (=" + switch + ")")
        status_exit(1)
    end if
    end if
    end if
    return use_for_contour
end

undef("empty_string_if_missing_attribute")
function empty_string_if_missing_attribute(res [1]: logical,
                                         prefix[1]: string)
local ret_str
begin
    if (isatt(res, prefix)) then
        ret_str = res@$prefix$
    else
        ret_str = ""
    end if
    return ret_str
end

undef("keep_only_n_cols")
function keep_only_n_cols(res [1]: logical,
                       long_str[1]: string)
;                returns shorter_string [1]:string
; Parameters:
;    res   - resource with attribute 'var_base_lenth', the number of cols
;            in string to keep
; long_str - string to be shortened
local shortened_string
begin
    if (isatt(res, "var_base_length")) then
        shortened_string = str_get_cols(long_str, 0, res@var_base_length - 1)
    else
        shortened_string = long_str
    end if
    return shortened_string
end
