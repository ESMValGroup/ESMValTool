;;#############################################################################
;; GENERAL ROUTINES FOR MODEL ENSEMBLE OPERATIONS
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function get_start_year
;;    function get_end_year
;;    function multi_model_stats
;;
;;#############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/logging.ncl"

load "./diag_scripts/shared/regridding.ncl"
load "./diag_scripts/shared/statistics.ncl"

;;#############################################################################
undef("get_start_year")
function get_start_year(model_index[1]:numeric)
;;
;; Arguments:
;;    model_index: index of array "models"; this function is supposed to be
;;                 called inside a model loop
;;
;; Return value
;;    An integer containing the first year
;;
;; Description
;;    Chooses first year of model input with respect to
;;    diag_script_info@range_option
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140128-A_senf_da: written.
;;
local funcname, scriptname, model_index
begin

    funcname = "get_start_year"
    scriptname = "diag_scripts/shared/ensemble.ncl"
    enter_msg(scriptname, funcname)

    if (.not.isatt(diag_script_info, "range_option")) then
        error_msg("f", "ensemble.ncl", funcname, "diag_script_info must " + \
                  "contain the range_option attribute")
    end if
    if (diag_script_info@range_option.eq.0) then
        start_year = toint(model_info@start_year(model_index))
    else if (diag_script_info@range_option.eq.1) then
        start_year = max(toint(model_info@start_year))
    else
        error_msg("f", "ensemble.ncl", funcname, "range_option not " + \
                  "properly set, must be either 0 or 1")
    end if
    end if

    leave_msg(scriptname, funcname)
    return(start_year)

end

;;#############################################################################
undef("get_end_year")
function get_end_year(model_index:numeric)
;;
;; Arguments
;;    model_index: index of array "models"; this function is supposed to be
;;                 called inside a model loop
;;
;; Return value
;;    An integer containing the last year
;;
;; Description
;;    Chooses last year of model input with respect to
;;    diag_script_info@range_option
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;    20140128-A_senf_da: written.
;;
local funcname, scriptname, model_index
begin

    funcname = "get_end_year"
    scriptname = "diag_scripts/shared/ensemble.ncl"
    enter_msg(scriptname, funcname)

    if (.not.isatt(diag_script_info, "range_option")) then
        error_msg("f", "ensemble.ncl", funcname, "diag_script_info must " + \
                  "contain the range_option attribute")
    end if
    if (diag_script_info@range_option.eq.0) then
        end_year = toint(model_info@end_year(model_index))
    else if (diag_script_info@range_option.eq.1) then
        end_year = min(toint(model_info@end_year))
    else
        error_msg("f", "ensemble.ncl", funcname, "range_option not " + \
                  "properly set, must be either 0 or 1")
    end if
    end if

    leave_msg(scriptname, funcname)
    return(end_year)

end


;;#############################################################################
undef("multi_model_stats")
function multi_model_stats(idx[*]:integer,
                           var:string,
                           field:string,
                           grd[1]:string,
                           statistics[*]:string,
                           opt[1])
;;
;; Arguments
;;    idx: a 1-D array of model indexes to be considered.
;;    var: variable name.
;;    field: field type.
;;    grd: type of grid to be selected:
;;           "coarsest": returns the lowest resolution grid.
;;           "finest": returns the highest resolution grid.
;;    statistics: a 1-D array of strings for selecting the statistical
;;                fields to return. Currently implemented:
;;                "model_mean" -> multi model mean
;;                "stddev" -> standard deviation over "models" dimension
;;                "model_grd" -> models regridded to a common grid
;;    opt: not used yet, set to 0
;;
;; Return value
;;    A field with target grid spatial dimensions (lat, lon) or
;;    (plev, lat, lon), time coordinate, models coordinate.
;;    The "models" coordinate contains the strings from model_info@model and
;;    statistics.
;;    Covers the overlapping time period of all models only.
;;
;; Description
;;    Determines statistical fields of all input models, after regridding them
;;    to a common grid. Statistics is over the "models" coordinate only.
;;    See function find_destination_grid (regridding.ncl) for details of
;;    target grid selection.
;;
;; Caveats
;;    Works for rectilinear grids only, but selection of regridding routines
;;    might be extended (e.g. select automatically, depending on input).
;;    Might be extended to return multi year monthly mean values.
;;    Wrap to add result to the structure returned by function read_data?
;;    Probably not very efficient: try to re-use intermediate data etc.
;;
;; References
;;
;; Modification history:
;;    20140515-A_gott_kl: written.
;;
local funcname, scriptname, dest_grd, monotonicity, field_type, rank_grd, \
    start_year, end_year, idx_mod, var_tmp, var_olap, var_grd, data_mmm, \
    plev_grd, ireturn, tmp, model, Fill, stats, nstats, tmp, iMod, \
    interpol_flag, nvalid, var_plev, result
begin

    funcname = "multi_model_stats"
    scriptname = "diag_scripts/shared/ensemble.ncl"
    enter_msg(scriptname, funcname)

    ;; determine destination grid (regridding.ncl)
    field_type = str_get_cols(field, 0, 1)
    if(field_type.eq."T2" .or. field_type.eq."T3") then
        dest_grd = find_destination_grid(idx, var, field, grd)
        rank_grd = dimsizes(dimsizes(dest_grd))
        ;; check monotonicity of pressure coordinate
        if(rank_grd.eq.3) then
            plev_grd = profile1d&plev
            monotonicity = isMonotonic(plev_grd)
            if (monotonicity.eq.0) then  ; crash
                error_msg("f", "ensemble.ncl", funcname, \
                          "plev_dat is not monotonic")
            end if
        end if
    else
        error_msg("f", "ensemble.ncl", funcname, \
                  "accepts T2* & T3*, but found " + field_type)
    end if

    ;; determine overlapping time period
    start_year = max(toint(model_info@start_year))
    end_year = min(toint(model_info@end_year))
    log_info("overlapping time period: " + tostring(start_year) + \
                " - " + tostring(end_year))

    ;; create data structure for mmm models
    do i = 0, dimsizes(idx) - 1
        idx_mod = idx(i)
        data_mmm = read_data(idx_mod, var, field)
    end do

    ;; loop over mmm models
    do i = 0, dimsizes(idx) - 1
        model = model_info@model(idx(i))

        ;; extract full data
        var_tmp = extract_data(i, data_mmm, -1, 0, 0)

        ;; extract overlapping time period
        var_olap = time_operations(var_tmp, start_year, end_year, \
                                   "extract", "", 0)
        delete(var_tmp)

        ;; initialize target field for regridded data *************************
        if(i.eq.0) then
            if(.not.isatt(var_olap, "_FillValue")) then
                Fill = default_fillvalue(typeof(var_olap))
            else
                Fill = var_olap@_FillValue
            end if
            if(rank_grd.eq.3) then
                var_grd = new((/dimsizes(idx) + 2, dimsizes(var_olap&time), \
                              dimsizes(dest_grd&plev), dimsizes(dest_grd&lat),\
                              dimsizes(dest_grd&lon)/), typeof(var_olap), Fill)
                var_grd!2 = "plev"
                var_grd!3 = "lat"
                var_grd!4 = "lon"
                var_grd&plev = dest_grd&plev
            else if(rank_grd.eq.2) then
                var_grd = new((/dimsizes(idx) + 2, dimsizes(var_olap&time), \
                              dimsizes(dest_grd&lat), \
                              dimsizes(dest_grd&lon)/), typeof(var_olap), Fill)
                var_grd!2 = "lat"
                var_grd!3 = "lon"
            end if
            end if
            var_grd!0 = "models"
            var_grd!1 = "time"
            var_grd&time = var_olap&time
            var_grd&lat = dest_grd&lat
            var_grd&lon = dest_grd&lon

            ;; implemented statistics, hardcoded positions
            stats = (/"model_mean", "stddev"/)
            nstats = dimsizes(stats)
            var_grd&models = array_append_record(stats, \
                                                 model_info@model(idx), 0)

            ;; select return fields
            tmp = new(dimsizes(var_grd&models), integer)
            if(any("model_mean" .eq. statistics)) then
                tmp(0) = 1
            end if
            if(any("stddev" .eq. statistics)) then
                tmp(1) = 1
            end if
            if(any("model_grd" .eq.  statistics)) then
                tmp(nstats:) = 1
            end if
            ireturn = ind(tmp.eq.1)  ; selected indices of models coordinate
            delete(tmp)
        end if  ; target field initialization *********************************
        iMod = nstats + i

        ;; ensure consistent missing values
        var_olap = where(ismissing(var_olap), Fill, var_olap)
        var_olap@_FillValue = Fill

        ;; check that all data are on the same time coordinate
        if(dimsizes(var_olap&time).ne.dimsizes(var_grd&time)) then
            error_msg("f", "ensemble.ncl", funcname, \
                      "inconsistent number of time steps for " + \
                      var_grd&models(iMod))
        end if
        if(any(var_olap&time.ne.var_grd&time)) then
            error_msg("f", "ensemble.ncl", funcname, \
                      "inconsistent time coordinate for " + \
                      var_grd&models(iMod))
        end if

        ;; pressure interpolation *********************************************
; NOT TESTED YET
; make this a function: used in Emmons.ncl, perfmetrics?
; more flexible dimensions, e.g. time missing?
        if(rank_grd.eq.3) then
            ;; check consistency of pressure levels
            if(isMonotonic(var_olap&plev).eq.0) then  ; crash
                error_msg("f", "ensemble.ncl", funcname, \
                          "plev is not monotonic for " + var_grd&models(iMod))
            end if
            if(isMonotonic(var_olap&plev).ne.monotonicity) then
                var_olap = var_olap(:, ::-1, :, :)  ; reverse
            end if

            ;; interpolate to  plev_grd
            if(isatt(diag_script_info, "interpolation")) then
                interpol_flag = diag_script_info@interpolation
                log_debug("Interpolation set by diag_script_info")
            else  ; default
                valid_plev = plev_grd(ind(.not.ismissing(plev_grd)))
                nvalid = dimsizes(valid_plev)
                if(max(var_olap&plev).gt.valid_plev(1) .and.\
                   min(var_olap&plev).lt.valid_plev(nvalid - 2)) then
                    interpol_flag = -1
                    log_debug("Linear interpolation with extrapolation")
                else
                    interpol_flag = 1
                    log_debug("Linear interpolation, no extrapolation")
                end if
            end if
            var_plev = int2p_n_Wrap(var_olap&plev, var_olap,\
                                    plev_grd, interpol_flag, 0)
            delete(var_olap)
            var_olap = var_plev
            copy_VarMeta(var_plev, var_olap)
            delete([/var_plev, valid_plev, nvalid/])
        end if  ; pressure interpolation **************************************

        ;; lat-lon interpolation, mmm, stdev **********************************
        if(grd.eq."coarsest") then
            if(rank_grd.eq.3) then
                var_grd(iMod, :, :, :, :) = \
                    area_hi2lores_Wrap(var_olap&lon, var_olap&lat, \
                                       var_olap, True, 1, dest_grd&lon, \
                                       dest_grd&lat, False)
            else if(rank_grd.eq.2) then
                var_grd(iMod, :, :, :) = \
                    area_hi2lores_Wrap(var_olap&lon, var_olap&lat, \
                                       var_olap, True, 1, dest_grd&lon, \
                                       dest_grd&lat, False)
            else
                error_msg("f", "ensemble.ncl", funcname, \
                          "rank_grd = " + tostring(rank_grd) + \
                          "not implemented")
            end if
            end if
        else
            error_msg("f", "ensemble.ncl", funcname, \
                      "interpolation option " + grd + " not yet implemented")
        end if  ; lat-lon interpolation ***************************************

        ;; clean up
        delete([/var_olap/])

    end do  ; i (loop over mmm models)

    ;; statistics
    if(rank_grd.eq.3) then
        var_grd(0, :, :, :, :) = \
            dim_avg_n_Wrap(var_grd(nstats:, :, :, :, :), 0)
        var_grd(1, :, :, :, :) = \
            dim_stddev_n_Wrap(var_grd(nstats:, :, :, :, :), 0)
        result = var_grd(ireturn, :, :, :, :)
    else if(rank_grd.eq.2) then
        var_grd(0, :, :, :) = dim_avg_n_Wrap(var_grd(nstats:, :, :, :), 0)
        var_grd(1, :, :, :) = dim_stddev_n_Wrap(var_grd(nstats:, :, :, :), 0)
        result = var_grd(ireturn, :, :, :)
    else
        error_msg("f", "ensemble.ncl", funcname, \
                  "multi_model_mean, rank_grd = " + \
                  tostring(rank_grd) + "not implemented")
    end if
    end if

;; DEBUG+++
;;     var_grd@diag_script = "void"
;;     var_grd@var = var_grd@standard_name
;;     printVarSummary(var_grd)
;;     outfile = ncdf_write(var_grd(0,:,:,:), "test0.nc")
;;     outfile = ncdf_write(var_grd(1,:,:,:), "test1.nc")
;;     outfile = ncdf_write(var_grd(2,:,:,:), "test2.nc")
;;     outfile = ncdf_write(var_grd(3,:,:,:), "test3.nc")
;;     outfile = ncdf_write(var_grd(4,:,:,:), "test4.nc")
;;     print("wrote "+outfile)
;; DEBUG---

    leave_msg(scriptname, funcname)
    return(result)

end
