;;#############################################################################
;; GENERAL ROUTINES FOR TRANSFORMING GO DATA
;;#############################################################################
;; Please consider using or extending existing routines before adding new ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;    function GO_month_to_season_extended
;;    function GO_apply_transforms
;;
;;#############################################################################

load "./diag_scripts/lib/ncl/statistics.ncl"
undef("GO_month_to_season_extended")
function GO_month_to_season_extended(indata[*][*][*]:float,
                                     season[1]:string)
;;
;; Arguments
;;    indata: a [lat][lon][time] array.
;;    season: compute the average for this season; include average where
;;            some points missing
;;
;; Return value
;;    An array with the seasonal average for each year.
;;
;; Description
;;    For each year in the input data, averages indata over the given season.
;;
;; Caveats
;;
;; References
;;
;; Modification history
;;
local funcname, scriptname, verbosity, season_indices, dim_season_indices, \
    runaveragedata, start_index, averagedata
begin
    funcname = "month_to_season_extended"
    scriptname = "diag_script/lib/ncl/statistics.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 8)

    if (dimsizes(indata&time) % 12 .ne. 0) then
        error_msg("f", scriptname, funcname, "time dimension must " + \
                  "be divisible by 12" + indata&time % 12)
    end if

    season_indices = calc_season_index(season)
    dim_season_indices = dimsizes(season_indices)

    nyear = dimsizes(indata&time) / 12
    runaveragedata = indata
    runaveragedata = indata@_FillValue

    ;; calculate dim_season_indices-mo running average
    do qq = 1, dimsizes(indata&time) - (dim_season_indices - 1)
        ;; missing data   allowed
        runaveragedata(qq, :, :) = dim_avg_n_Wrap(indata(qq - 1:qq + dim_season_indices - 2, :, :), 0)
    end do

    ;; By picking the correct 'start_index' in the running average array we
    ;; will retrieve the average over the indicated season
    ;; (see runave-documentation for details)
    if (dim_season_indices .eq. dimsizes(indata&time)) then
        start_index = 1
    else if (dim_season_indices % 2 .eq. 0) then
        start_index = season_indices(0) + (dim_season_indices - 2) / 2
    else
        start_index = season_indices(0) + (dim_season_indices - 1) / 2
    end if
    end if

    ;; Extract seasonal average for every year
    averagedata = (/runaveragedata(time|start_index::12, lat|:, lon|:)/)

    averagedata!0 = "time"
    averagedata&time = runaveragedata&time(start_index::12)

    averagedata!1 = "lat"
    averagedata&lat = runaveragedata&lat

    averagedata!2 = "lon"
    averagedata&lon = runaveragedata&lon

    copy_VarAtts(indata, averagedata)

    leave_msg(scriptname, funcname, 8)
    return(averagedata)

end
undef("GO_apply_transforms")
function GO_apply_transforms(transform[1] : logical,
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
;;    Additonal to transforms available in SAMonsoon;
;;    GO_season suitable for data with missing values
local data_regrid,  data_tmp,   data_trans,  dims_3d,\
    idx_trans,    newdims,      ntimes,      Opt,\
    t,            test_regrid,  verbosity,   SrcMask2D,\
    weight,       source,       destination
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering GO_apply_transforms", verbosity, 6)
    data_trans = data
    dim_trans = dimsizes(transform@apply_these)

    do idx_trans = 0, dim_trans - 1

        ;; Restrict data set to a certain season
        if(transform@apply_these(idx_trans) .eq. "GO_season") then
            data_tmp = GO_month_to_season_extended(data, transform@GO_season)
            delete(data_trans)
            data_trans = data_tmp
            delete(data_tmp)
        end if

         ;; Carry out linear detrending
        if(transform@apply_these(idx_trans) .eq. "GO_detrending") then
            data_tmp = dtrend_msg_n(data_trans&time, data_trans, True, False, 0)
            data_tmp!0 = "time"
            data_tmp!1 = "lat"
            data_tmp!2 = "lon"
            data_tmp&time = data_trans&time
            data_tmp&lat = data_trans&lat
            data_tmp&lon = data_trans&lon
            delete(data_trans)
            data_trans = data_tmp
            delete(data_tmp)
        end if
    end do
    info_output(">>>>>>>> Leaving GO_apply_transforms", verbosity, 6)
    return data_trans
end
