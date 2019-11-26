; #############################################################################
; GENERAL ROUTINES FOR STATISTICS
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function dim_stddev_wgt_Wrap
;    function time_operations
;    function calc_season_index
;    function extract_season
;    function month_to_season_extended
;    function coswgt_areaave
;    function coswgt_arearmse
;    function coswgt_pattern_cor
;    function interannual_variability
;    function calculate_metric
;    function normalize_metric
;    function distrib_stats
;    function lognormal_dist
;    function filter121
;    function get_average
;
; #############################################################################

load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/regridding.ncl"

; #############################################################################
undef("dim_stddev_wgt_Wrap")
function dim_stddev_wgt_Wrap(field[*]:numeric,
                             ww[*]:numeric,
                             opt[1]:integer)
;
; Arguments
;    field: a one-dimensional numeric array.
;    ww: a one-dimensional numeric array of the same size of field.
;    opt: a scalar, it has the same meaning as in the corresponding NCL
;    function dim_avg_wgt_Wrap
;
; Return value
;    A float or a double depending on the type of input
;
; Description
;    Calculates the (unbiased) weighted standard deviation consistently with
;    the NCL function dim_std_dev (i.e. it divides by N-1 instead of N). For
;    the weighted case this means applying a correction factor:
;      sum(w_i)/[sum(w_i)^2 - sum(w_i^2)]
;    Missing values are ignored.
;
; Caveats
;
; References
;    en.wikipedia.org/wiki/Weighted_arithmetic_mean#Weighted_sample_variance
;
; Modification history
;    20150511-lauer_axel: modified routine "calculate_metric":
;                         added "no weigth" (nowgt) option
;    20141215-righi_mattia: written.
;
local funcname, scriptname, wavg, wgt, v1, v2, d2, arg
begin

  funcname = "dim_stddev_wgt_Wrap"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  ; Copy to local
  wgt = ww
  wgt@_FillValue = default_fillvalue(typeof(wgt))

  if (dimsizes(field).lt.2) then
    error_msg("f", scriptname, funcname, "input field contains " + \
              "only 1 element, cannot calculate standard deviation")
  end if

  ; Calculate weighted mean
  wavg = dim_avg_wgt_Wrap(field, ww, opt)

  ; Filter missing values
  if (.not.all(ismissing(field))) then
    v1 = sum(where(ismissing(field), wgt@_FillValue, wgt))
    v2 = sum(where(ismissing(field), wgt@_FillValue, wgt) ^ 2)
  else
    out = 1.  ; initialize
    out@_FillValue = field@_FillValue
    out = field@_FillValue
    return(out)
  end if

  ; Calculate weighted standard deviation
  d2 = (field - wavg) ^ 2
  arg = dim_sum_wgt_Wrap(d2, ww, opt)
  out = sqrt(v1 / (v1 ^ 2 - v2) * arg)

  leave_msg(scriptname, funcname)
  return(out)

end

; #############################################################################
undef("time_operations")
function time_operations(field:numeric,
                         y1[1]:integer,
                         y2[1]:integer,
                         oper[1]:string,
                         opt[1]:string,
                         l_wgt[1]:logical)
;
; Arguments
;    field: a numeric array of rank 1 to 4, first dimension must be time.
;    y1: start year of the time period to be averaged (-1 for full range).
;    y2: end year of the time period to be averaged (-1 for full range).
;    oper: type of operations:
;            "extract": no average, just extract selected period.
;            "average": average.
;            "stddev": (unbiased) standard deviation.
;    opt: operation options (has no effect is oper = extract):
;           "annualclim": annual climatology.
;           "seasonalclim": seasonal climatology for the standard seasons DJF,
;                           MAM, JJA, SON.
;           "monthlyclim": monthly climatology jan-dec.
;           "yearly": time average over every year in [y1:y2].
;           [month strings]: climatology of selected (consecutive) months
;                            (e.g., "MAM", "SONDJ").
;           [1, 12]: climatology of the selected month ("1"=Jan, "2"=Feb, ...,
;                    "12"=Dec).
;    l_wgt: if True, calculate weighted average, with days-per-month as
;           weights (has no effect is opt = "extract").
;
; Return value
;    An array of the same rank as field or of rank-1, depending on oper/opt.
;
; Description
;    Performs differnt types of time average, standard deviation or extraction
;    of a selected time period. Weighted average (with days-per-month as
;    weights) can be optionally applied.
;
; Caveats
;    The weighted standard deviation is not yet implmented for all cases
;    The weighted standard deviation is calculated using the unbiased estimator
;    This should take into account missing values and exclude the w_i for
;    which the field contains only missing values. This feature is not
;    implemented yet.
;
; References
;
; Modification history
;    20190503-righi_mattia: removed obsolete option "mymm" (used only in
;                           reformat_obs, now outdated).
;    20140703-gottschaldt_klaus-dirk: added option "mymm".
;    20140312-righi_mattia: extended with standard deviation.
;    20140109-righi_mattia: written.
;
local funcname, scriptname, monthstr, date, year, month, idx1, idx2, loc_y1, \
  loc_y2, rank, subfield, weights, idx, idx_win, idx_spr, idx_sum, idx_aut, \
  mm, idx_1st, idx_arr, p1, p2, d2, arg, v1, v2, ym_in, years, nyear, ym, \
  dims, timec, FillValue, index
begin

  funcname = "time_operations"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  ; Check arguments
  if (all(oper.ne.(/"extract", "average", "stddev"/))) then
    error_msg("f", scriptname, funcname, "unrecognized operation " + oper)
  end if

  ; Check for time dimension
  if (field!0.ne."time") then
    error_msg("f", scriptname, funcname, "the first dimension " + \
              "of input is not time")
  end if

  ; Check for calendar attribute
  if (.not.isatt(field&time, "calendar")) then
    error_msg("f", scriptname, funcname, "time dimension of " + \
              "input must have a calendar attribute")
  end if

  ; Define months string
  monthstr = "JFMAMJJASOND"
  monthstr = monthstr + monthstr

  ; Define flags
  l_ext = oper.eq."extract"
  l_avg = oper.eq."average"
  l_std = oper.eq."stddev"

  ; Calculate date from time coordinate
  date := cd_calendar(field&time, 0)
  year := date(:, 0)
  month := date(:, 1)

  ; Determine indexes for the requested time range
  if (y1.eq.-1) then
    idx1 = 0
    loc_y1 = toint(min(date(:, 0)))
  else
    idx1 = min(ind(year.eq.y1))
    loc_y1 = y1
  end if
  if (y2.eq.-1) then
    idx2 = dimsizes(field&time) - 1
    loc_y2 = toint(max(date(:, 0)))
  else
    idx2 = max(ind(year.eq.y2))
    loc_y2 = y2
  end if
  if (ismissing(idx1).or.ismissing(idx2)) then
    error_msg("f", scriptname, funcname, "the selected time " + \
              "period is out of range")
  end if
  delete(date)
  delete(year)
  delete(month)

  ; Extract requested time range
  rank = dimsizes(dimsizes(field))
  if (rank.eq.4) then
    subfield = field(idx1:idx2, :, :, :)
  end if
  if (rank.eq.3) then
    subfield = field(idx1:idx2, :, :)
  end if
  if (rank.eq.2) then
    subfield = field(idx1:idx2, :)
  end if
  if (rank.eq.1) then
    subfield = field(idx1:idx2)
  end if

  ; Re-calculate date for subfield
  date := cd_calendar(subfield&time, 0)
  year := date(:, 0)
  month := date(:, 1)
  rank := dimsizes(dimsizes(subfield))

  ; Define weights as days-per-month
  if (l_wgt) then
    weights = days_in_month(toint(year), toint(month))
  else
    weights = tofloat(subfield&time)
    weights = 1.
  end if

  ; Extract only
  if (l_ext .and. opt.eq."") then
    leave_msg(scriptname, funcname)
    return(subfield)
  end if

  ; Calculate time average/standard deviation according to the opt argument

  ; Multi-year average
  if (opt.eq."annualclim") then
    if (l_avg) then
      out = dim_avg_wgt_n_Wrap(subfield, weights, 1, 0)
    end if
    if (l_std) then
      error_msg("f", scriptname, funcname, "feature not yet implemented")
    end if
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Year average
  if (opt.eq."yearly") then
    ny = loc_y2 - loc_y1 + 1
    if (rank.eq.4) then
      out = subfield(0:ny - 1, :, :, :)  ; Copy metadata
      do yy = loc_y1, loc_y2
        idx = ind(year.eq.yy)
        if (l_avg) then
          out(yy - loc_y1, :, :, :) = \
            dim_avg_wgt_n_Wrap(subfield(idx, :, :, :), weights(idx), 1, 0)
        end if
        if (l_std) then
          error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
        delete(idx)
      end do
    end if
    if (rank.eq.3) then
      out = subfield(0:ny - 1, :, :)  ; Copy metadata
      do yy = loc_y1, loc_y2
        idx = ind(year.eq.yy)
        if (l_avg) then
          out(yy - loc_y1, :, :) = \
            dim_avg_wgt_n_Wrap(subfield(idx, :, :), weights(idx), 1, 0)
        end if
        if (l_std) then
          error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
        delete(idx)
      end do
    end if
    if (rank.eq.2) then
      out = subfield(0:ny - 1, :)  ; Copy metadata
      do yy = loc_y1, loc_y2
        idx = ind(year.eq.yy)
        if (l_avg) then
          out(yy - loc_y1, :) = \
            dim_avg_wgt_n_Wrap(subfield(idx, :), weights(idx), 1, 0)
        end if
        if (l_std) then
          error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
        delete(idx)
      end do
    end if
    if (rank.eq.1) then
      out = subfield(0:ny - 1)  ; Copy metadata
      do yy = loc_y1, loc_y2
        idx = ind(year.eq.yy)
        if (l_avg) then
          out(yy - loc_y1) = \
            dim_avg_wgt_Wrap(subfield(idx), weights(idx), 1)
        end if
        if (l_std) then
          error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
        delete(idx)
      end do
    end if
    out!0 = "year"
    delete(out&year)
    out&year = ispan(loc_y1, loc_y2, 1)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Season average
  if (opt.eq."seasonalclim") then
    idx_win = ind(month.eq.1.or.month.eq.2.or.month.eq.12)
    idx_spr = ind(month.eq.3.or.month.eq.4.or.month.eq.5)
    idx_sum = ind(month.eq.6.or.month.eq.7.or.month.eq.8)
    idx_aut = ind(month.eq.9.or.month.eq.10.or.month.eq.11)
    if (rank.eq.4) then
      out = subfield(0:3, :, :, :)
      if (l_avg) then
        out(0, :, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_win, :, :, :), \
                             weights(idx_win), 1, 0)
        out(1, :, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_spr, :, :, :), \
                             weights(idx_spr), 1, 0)
        out(2, :, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_sum, :, :, :), \
                             weights(idx_sum), 1, 0)
        out(3, :, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_aut, :, :, :), \
                             weights(idx_aut), 1, 0)
      end if
      if (l_std) then
        error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.3) then
      out = subfield(0:3, :, :)
      if (l_avg) then
        out(0, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_win, :, :), weights(idx_win), 1, 0)
        out(1, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_spr, :, :), weights(idx_spr), 1, 0)
        out(2, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_sum, :, :), weights(idx_sum), 1, 0)
        out(3, :, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_aut, :, :), weights(idx_aut), 1, 0)
      end if
      if (l_std) then
        error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.2) then
      out = subfield(0:3, :)
      if (l_avg) then
        out(0, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_win, :), weights(idx_win), 1, 0)
        out(1, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_spr, :), weights(idx_spr), 1, 0)
        out(2, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_sum, :), weights(idx_sum), 1, 0)
        out(3, :) = \
          dim_avg_wgt_n_Wrap(subfield(idx_aut, :), weights(idx_aut), 1, 0)
      end if
      if (l_std) then
        error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.1) then
      out = subfield(0:3)
      if (l_avg)  then
        out(0) = dim_avg_wgt_Wrap(subfield(idx_win), weights(idx_win), 1)
        out(1) = dim_avg_wgt_Wrap(subfield(idx_spr), weights(idx_spr), 1)
        out(2) = dim_avg_wgt_Wrap(subfield(idx_sum), weights(idx_sum), 1)
        out(3) = dim_avg_wgt_Wrap(subfield(idx_aut), weights(idx_aut), 1)
      end if
      if (l_std) then
        error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    delete(out&time)
    out!0 = "season"
    out&season = (/"DJF", "MAM", "JJA", "SON"/)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Annual cycle
  if (opt.eq."monthlyclim") then
    if (rank.eq.4) then
      out = subfield(0:11, :, :, :)  ; Copy metadata
      do mm = 0, 11
        if (l_avg)  then
          out(mm, :, :, :) = \
            dim_avg_wgt_n_Wrap(subfield(mm::12, :, :, :), \
                               weights(mm::12), 1, 0)
        end if
        if (l_std) then
          error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
      end do
    end if
    if (rank.eq.3) then
      out = subfield(0:11, :, :)  ; Copy metadata
      do mm = 0, 11
        if (l_avg)  then
          out(mm, :, :) = \
            dim_avg_wgt_n_Wrap(subfield(mm::12, :, :), weights(mm::12), 1, 0)
        end if
        if (l_std) then
          error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
      end do
    end if
    if (rank.eq.2) then
      out = subfield(0:11, :)  ; Copy metadata
      do mm = 0, 11
        if (l_avg) then
          out(mm, :) = \
            dim_avg_wgt_n_Wrap(subfield(mm::12, :), weights(mm::12), 1, 0)
        end if
        if (l_std) then
          error_msg("f", scriptname, funcname, "feature not yet implemented")
        end if
      end do
    end if
    if (rank.eq.1) then
      out = subfield(0:11)  ; Copy metadata
      do mm = 0, 11
        if (l_avg) then
          out(mm) = dim_avg_wgt_Wrap(subfield(mm::12), weights(mm::12), 1)
        end if
        if (l_std) then
          out(mm) = dim_stddev_wgt_Wrap(subfield(mm::12), weights(mm::12), 1)
        end if
      end do
    end if
    out!0 = "month"
    delete(out&month)
    out&month = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"/)

    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Months string (at least 2 consecutive months): define indexes
  if (.not.ismissing(str_match_ind_ic(monthstr, opt)).and. \
      strlen(opt).ge.2.and.strlen(opt).le.12) then
    idx_1st =  str_index_of_substr(monthstr, str_upper(opt), 1)
    idx_arr = new(strlen(opt), integer)
    do ii = 0, strlen(opt) - 1
      idx_arr(ii) = idx_1st + ii
    end do
    idx_arr = where(idx_arr.ge.12, idx_arr - 12, idx_arr)  ; Periodicity
    idx_arr = idx_arr + 1  ; From 0-based to month number
    do ii = 0, dimsizes(idx_arr) - 1
      if (.not.isdefined("idx"))  then
        idx = ind(month.eq.idx_arr(ii))
      else
        tmp = array_append_record(idx, ind(month.eq.idx_arr(ii)), 0)
        delete(idx)
        idx = tmp
        delete(tmp)
      end if
    end do
    delete(idx_1st)
    delete(idx_arr)
  end if

  ; Specific-month average: define indexes
  if (any(opt.eq.tostring(ispan(1, 12, 1)))) then
    idx = ind(month.eq.toint(opt))
  end if

  ; Extract or average over the above indexes
  if (isdefined("idx")) then
    if (rank.eq.4) then
      if (l_ext) then
        out = subfield(idx, :, :, :)
      end if
      if (l_avg) then
        if (dimsizes(idx) .eq. 1) then
          out = \
            dim_avg_wgt_n_Wrap(subfield(idx:idx, :, :, :), weights(idx), 1, 0)
        else
          out = dim_avg_wgt_n_Wrap(subfield(idx, :, :, :), weights(idx), 1, 0)
        end if
      end if
      if (l_std) then
        error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.3) then
      if (l_ext) then
        out = subfield(idx, :, :)
      end if
      if (l_avg) then
        if (dimsizes(idx) .eq. 1) then
          out = dim_avg_wgt_n_Wrap(subfield(idx:idx, :, :), weights(idx), 1, 0)
        else
          out = dim_avg_wgt_n_Wrap(subfield(idx, :, :), weights(idx), 1, 0)
        end if
      end if
      if (l_std) then
        error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.2) then
      if (l_ext) then
        out = subfield(idx, :)
      end if
      if (l_avg) then
        if (dimsizes(idx) .eq. 1) then
          out = dim_avg_wgt_n_Wrap(subfield(idx:idx, :), weights(idx), 1, 0)
        else
          out = dim_avg_wgt_n_Wrap(subfield(idx, :), weights(idx), 1, 0)
        end if
      end if
      if (l_std) then
        error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    if (rank.eq.1) then
      if (l_ext) then
        out = subfield(idx)
      end if
      if (l_avg) then
        if (dimsizes(idx) .eq. 1) then
          out = dim_avg_wgt_n_Wrap(subfield(idx:idx), weights(idx), 1, 0)
        else
          out = dim_avg_wgt_n_Wrap(subfield(idx), weights(idx), 1, 0)
        end if
      end if
      if (l_std) then
        error_msg("f", scriptname, funcname, "feature not yet implemented")
      end if
    end if
    leave_msg(scriptname, funcname)
    return(out)
  end if

  error_msg("f", scriptname, funcname, "unrecognized option " + opt)

end

; #############################################################################
undef("calc_season_index")
function calc_season_index(season[1]:string)
;
; Arguments
;    season: the season in upper case.
;
; Return value
;    The indices to the months in season, e.g. "JFM" returns (/0, 1, 2/).
;
; Description
;    Given the "season", i.e., any substring from "JFMAMJJASOND", retrieves
;    the corresponding indices. Crashes if given substring is not unique or
;    does not exist.
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, current_search_set, DEBUG, i, indices, months, \
  start_index, strIndex, stringmonths, string_search_set, subStringLength
begin

  funcname = "calc_season_index"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  DEBUG = False

  ; The months with a wrap around of six months
  ; (so we can find "DJF", etc..)
  stringmonths = "JFMAMJJASONDJFMAMJ"
  months = stringtochar(stringmonths)

  ; Loop using twelve months at a time from 'stringmonths', i.e, first
  ; iteration uses Jan   -> Dec, next iteration uses April -> March, ...
  do start_index = 0, 6, 3
    current_search_set = months(start_index:start_index + 11)
    string_search_set = charactertostring(current_search_set)
    if (DEBUG) then
      error_msg("f", scriptname, funcname, "string_search_set = " + \
                string_search_set)
    end if
    strIndex = str_index_of_substr(string_search_set, season, 0)
    if (DEBUG) then
      error_msg("f", scriptname, funcname, "strIndex = " + strIndex)
    end if

    ; Exit immediately if our substring 'season' is not unique
    ; (e.g., season="J")
    if (dimsizes(strIndex) .gt. 1) then
      error_msg("f", scriptname, funcname, "multiple occurences of substring")
    end if

    ; If there is no match, grab a new set of twelve months
    ; from the wrap-around 'stringmonths'
    if (ismissing(strIndex)) then
      continue
    end if
    break
  end do

  ; Exit if there is no match
  if (ismissing(strIndex)) then
    error_msg("f", scriptname, funcname, "could not find substring")
  end if

  ; Compute the indices for the requested season
  subStringLength = sizeof(stringtochar(season)) - 1
  if (DEBUG) then
    error_msg("f", scriptname, funcname, "subStringLength = " + \
              subStringLength)
  end if
  indices = ispan(strIndex, strIndex + subStringLength - 1, 1)
  indices = indices + start_index
  if (DEBUG) then
    error_msg("f", scriptname, funcname, "indices=" + indices)
  end if

  ; Subtract indices larger than 11 (i.e., we have a wrap around case)
  do i = 0, dimsizes(indices) - 1
    if (indices(i) .gt. 11) then
      indices(i) = indices(i) - 12
    end if
  end do

  leave_msg(scriptname, funcname)
  return(indices)

end

; #############################################################################
undef("extract_season")
function extract_season(data:numeric,
                        season[1]:string)
;
; Arguments
;    data: a numeric field with time dimension.
;    season:  the season in upper case.
;
; Return value
;    The temporal subset of indata defined by the 'season' string.
;
; Description
;    Given the "season", i.e., any substring from "JFMAMJJASOND", retrieves
;    the corresponding months from data.
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, idx, subset_selection, months, newtime, \
  seasonal_selection_ind, season_monthly_indices, sizes, start_of_year, \
  tmp, local_indata, ndims
begin

  funcname = "extract_season"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  season_monthly_indices = calc_season_index(season)
  months = tointeger(cd_convert(data&time, "months since 1850-01-01 00:00"))

  local_data = data
  ; Ensure first month is January
  if (months(0) % 12 .ne. 0) then
    start_of_year = ind(tointeger(months) % 12 .eq. 0)
    sizes = dimsizes(local_data&time)
    tmp = local_data(time|start_of_year(0):sizes - 1, lat|:, lon|:)
    delete(local_data)
    local_data = tmp
    delete(tmp)
    delete(months)
    months = tointeger( \
      cd_convert(local_data&time, "months since 1800-01-01 00:00"))
  end if
  months = months - min(months)

  ; Extract the given months
  seasonal_selection_ind = ind(months % 12 .eq. season_monthly_indices(0))
  do idx = 1, dimsizes(season_monthly_indices) - 1
    tmp = array_append_record(seasonal_selection_ind, \
                              ind(months % 12.eq. \
                                  season_monthly_indices(idx)), 0)
    delete(seasonal_selection_ind)
    seasonal_selection_ind = tmp
    delete(tmp)
  end do
  qsort(seasonal_selection_ind)
  ndims = dimsizes(dimsizes(local_data))
  if (ndims .eq. 4) then  ; Assume time, plev, lat, lon
    subset_selection = \
      local_data(time|seasonal_selection_ind, plev|:, lat|:, lon|:)
  else if (ndims .eq. 3) then  ; Assume time, lat, lon
    subset_selection = local_data(time|seasonal_selection_ind, lat|:, lon|:)
  else if (ndims .eq. 2) then  ; Assume time, lat
    subset_selection = local_data(time|seasonal_selection_ind, lat|:)
  else if (ndims .eq. 1) then  ; Assume time is only dimension
    subset_selection = local_data(time|seasonal_selection_ind)
  else
    printVarSummary(local_data)
    error_msg("fatal", scriptname, funcname,\
              "ndims (=" + ndims + ") is assumed 1 <= ndims <= 4")
    status_exit(1)
  end if
  end if
  end if
  end if

  leave_msg(scriptname, funcname)
  return(subset_selection)

end

; #############################################################################
undef("month_to_season_extended")
function month_to_season_extended(indata:float,
                                  season[1]:string)
;
; Arguments
;    indata: a [lat][lon][time] or.
;            a [lat][lon][plev|[time] array
;    season: compute the average for this season.
;
; Return value
;    An array with the seasonal average for each year.
;
; Description
;    For each year in the input data, averages indata over the given season.
;
; Caveats
;
; References
;
; Modification history
;
local funcname, scriptname, season_indices, dim_season_indices, \
  runaveragedata, start_index, averagedata, dim
begin

  funcname = "month_to_season_extended"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  indata_size = dimsizes(indata)
  if (dimsizes(indata&time) % 12 .ne. 0) then
    error_msg("f", scriptname, funcname, "time dimension must " + \
              "be divisible by 12" + indata&time % 12)
  end if

  season_indices = calc_season_index(season)
  dim_season_indices = dimsizes(season_indices)

  ; Compute average over 'dim_season_indices' number of months with the help
  ; of a running average
  if (dimsizes(indata_size) .eq. 3) then
    runaveragedata = \
      runave_Wrap(indata({lat|:}, {lon|:}, time|:), dim_season_indices, 0)
  else if(dimsizes(indata_size) .eq. 4) then
    runaveragedata = \
      runave_Wrap(indata({plev|:}, {lat|:}, {lon|:}, time|:), \
                  dim_season_indices, 0)
  else
    error_msg("f", scriptname, funcname, "wrong number of dimensions: " + \
              dimsizes(indata_size) + ", only arrays with 3 or 4 " + \
              "dimensions are supported")
  end if
  end if

  ; By picking the correct 'start_index' in the running average array we
  ; will retrieve the average over the indicated season
  ; (see runave-documentation for details)
  if (dim_season_indices % 2 .eq. 0) then
    start_index = season_indices(0) + (dim_season_indices - 2) / 2
  else
    start_index = season_indices(0) + (dim_season_indices - 1) / 2
  end if

  ; Extract seasonal average for every year
  if (dimsizes(indata_size) .eq. 3) then
    averagedata = (/runaveragedata(time|start_index::12, lat|:, lon|:)/)
  else if(dimsizes(indata_size) .eq. 4) then
    averagedata = \
      (/runaveragedata(time|start_index::12, plev|:, lat|:, lon|:)/)
  else
    error_msg("f", scriptname, funcname, "wrong number of dimensions: " + \
              dimsizes(indata_size) + ", only arrays with 3 or 4 " + \
              "dimensions are supported")
  end if
  end if

  dim = 0
  averagedata!dim = "time"
  averagedata&time = runaveragedata&time(start_index::12)
  dim = dim + 1

  if(dimsizes(indata_size) .eq. 4) then
    averagedata!dim = "plev"
    averagedata&plev = runaveragedata&plev
    dim = dim + 1
  end if

  averagedata!dim = "lat"
  averagedata&lat = runaveragedata&lat
  dim = dim + 1

  averagedata!dim = "lon"
  averagedata&lon = runaveragedata&lon

  copy_VarAtts(indata, averagedata)

  leave_msg(scriptname, funcname)
  return(averagedata)

end

; #############################################################################
undef("coswgt_areaave")
function coswgt_areaave(field:numeric)
;
; Arguments
;    field: numeric field.
;
; Return value
;    The area average using cosine lat weights.
;
; Description
;    Computes the area average using cosine lat weights and lon weights=1.
;
; Caveats
;
; References
;
; Modification history
;    20131209-evaldsson_martin: written.
;
local funcname, scriptname, lat, wgt_lat, lon, lon_size, wgt_lon, ave
begin

  funcname = "coswgt_areaave"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  lat = field&lat
  wgt_lat = tofloat(NormCosWgtGlobe(lat))

  lon = field&lon
  lon_size = dimsizes(lon)
  wgt_lon = new((/lon_size(0)/), float)
  wgt_lon = 1.0

  ave = wgt_areaave_Wrap(field, wgt_lat, wgt_lon, 0)

  leave_msg(scriptname, funcname)
  return(ave)

end

; #############################################################################
undef("coswgt_arearmse")
function coswgt_arearmse(field1:numeric,
                         field2:numeric)
;
; Arguments
;    field1: numeric field
;    field2: numeric field
;
; Return value
;    Area rmse average using cosine lat weights.
;
; Description
;    Computes area rmse areage using cosine lat weights and lon weights=1.
;
; Caveats
;
; References
;
; Modification history
;    20131209-evaldsson_martin: written.
;
local funcname, scriptname, lat, wgt_lat, lon, lon_size, wgt_lon, rmse, \
  local_field1, local_field2
begin

  funcname = "coswgt_arearmse"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  field1_grid_size = guestimate_average_grid_area(field1)
  field2_grid_size = guestimate_average_grid_area(field2)

  if (field1_grid_size .gt. field2_grid_size) then
    local_field2 = rect2rect_interp(field2, field1)
    local_field1 = field1
  else
    local_field1 = rect2rect_interp(field1, field2)
    local_field2 = field2
  end if

  lat = local_field1&lat
  wgt_lat = tofloat(NormCosWgtGlobe(lat))

  lon = local_field1&lon
  lon_size = dimsizes(lon)
  wgt_lon = new((/lon_size(0)/), float)
  wgt_lon = 1.0

  rmse = wgt_arearmse(local_field1, local_field2, wgt_lat, wgt_lon, 0)

  leave_msg(scriptname, funcname)
  return(rmse)

end

; #############################################################################
undef("coswgt_pattern_cor")
function coswgt_pattern_cor(field1:numeric,
                            field2:numeric)
;
; Arguments
;    field1: numeric field.
;    field2: numeric field.
;
; Return value
;    Pattern correlation cosine lat weights.
;
; Description
;
; Caveats
;
; References
;
; Modification history
;    20140115-evaldsson_martin: written.
;
local funcname, scriptname, lat, wgt_lat, lon, lon_size, wgt_lon, \
  pattern_correlation, local_field1, local_field2
begin

  funcname = "coswgt_pattern_cor"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  field1_grid_size = guestimate_average_grid_area(field1)
  field2_grid_size = guestimate_average_grid_area(field2)

  if (field1_grid_size .gt. field2_grid_size) then
    local_field2 = rect2rect_interp(field2, field1)
    local_field1 = field1
  else
    local_field1 = rect2rect_interp(field1, field2)
    local_field2 = field2
  end if

  lat = local_field1&lat
  wgt_lat = tofloat(NormCosWgtGlobe(lat))

  lon = local_field1&lon
  lon_size = dimsizes(lon)
  wgt_lon = new((/lon_size(0)/), float)
  wgt_lon = 1.0

  pattern_correlation = pattern_cor(local_field1, local_field2, wgt_lat, 0)

  leave_msg(scriptname, funcname)
  return(pattern_correlation)

end

; #############################################################################
undef("interannual_variability")
function interannual_variability(field: numeric,
                                 y1[1]: integer,
                                 y2[1]: integer,
                                 opt[1]: string,
                                 dtr[1]: string)
;
; Arguments
;    field: a numeric array of rank 1 to 4, first dimension must be time.
;    y1: start year of the time period to be averaged (-1 for full range).
;    y2: end year of the time period to be averaged (-1 for full range).
;    opt: operation options (same as time_operations):
;           "annualclim": annual climatology.
;           "seasonalclim": seasonal climatology for the standard seasons DJF,
;                           MAM, JJA, SON.
;           "monthlyclim": monthly climatology jan-dec.
;           [month strings]: climatology of selected (consecutive) months
;                            (e.g. "MAM", "SONDJ").
;           [1, 12]: climatology of the selected month ("1"=Jan, "2"=Feb, ...,
;                    "12"=Dec).
;    dtr: detrending option:
;           "None": no detrending before standard deviation is calculated
;           "linear": linear detrending using dtrend
;           "quadratic": quadratic detrending using dtrend_quadratic
;
; Return value
;    An array of the same rank as field or of rank-1, depending on opt.
;
; Description
;    Calculates the standard deviation with respect to interannual
;    variability, to be used as input for statistical tests.
;
; Caveats
;    The standard deviation is not weighted, being w.r.t. interannual
;    variability for which all years have the same weight.
;
; Reference
;
; Modification history
;    20181022-lorenz_ruth: added option dtr for possible detrending of data v2
;    20140314-righi_mattia: written.
;
local funcname, scriptname, monthstr, rank, field_avg, field_djf, field_mam, \
  field_jja, field_son, field_avg_djf, field_avg_mam, field_avg_jja, \
  field_avg_son
begin

  funcname = "interannual_variability"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  ; Check for time dimension
  if (field!0.ne."time") then
    error_msg("f", scriptname, funcname, "the first dimension " + \
              "of input is not time")
  end if

  ; Check for calendar attribute
  if (.not.isatt(field&time, "calendar")) then
    error_msg("f", scriptname, funcname, "time dimension of " + \
              "input must have a calendar attribute")
  end if

  ; Define months string
  monthstr = "JFMAMJJASOND"
  monthstr = monthstr + monthstr

  ; Define rank
  rank = dimsizes(dimsizes(field))

  ; Annual climatology
  if (opt.eq."annualclim") then
    field_avg = time_operations(field, y1, y2, "average", "yearly", True)
    if (dtr.eq."None") then
      out = dim_stddev_n_Wrap(field_avg, 0)
    else if (dtr.eq."linear") then
      field_dtr = dtrend_n(field_avg, False, 0)
      tmp = dim_avg_n_Wrap(field_avg, 0)
      tmp_conf = conform(field_dtr, tmp, (/1, 2/))
      field_dtr = field_dtr + tmp_conf
      copy_VarCoords(field_avg, field_dtr)
      out = dim_stddev_n_Wrap(field_dtr, 0)
    else if (dtr.eq."quadratic") then
      field_dtr = dtrend_quadratic_msg_n(field_avg, False, False, 0)
      copy_VarCoords(field_avg, field_dtr)
      out = dim_stddev_n_Wrap(field_dtr, 0)
    end if
    end if
    end if
    out = dim_stddev_n_Wrap(field_avg, 0)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Seasonal climatology
  if (opt.eq."seasonalclim") then
    field_djf = time_operations(field, y1, y2, "extract", "DJF", True)
    field_mam = time_operations(field, y1, y2, "extract", "MAM", True)
    field_jja = time_operations(field, y1, y2, "extract", "JJA", True)
    field_son = time_operations(field, y1, y2, "extract", "SON", True)

    field_avg_djf = time_operations(field_djf, y1, y2, "average", \
                                    "yearly", True)
    field_avg_mam = time_operations(field_mam, y1, y2, "average", \
                                    "yearly", True)
    field_avg_jja = time_operations(field_jja, y1, y2, "average", \
                                    "yearly", True)
    field_avg_son = time_operations(field_son, y1, y2, "average", \
                                    "yearly", True)
    if (dtr.eq."None") then
      field_avg_djf_dtr = field_avg_djf
      field_avg_mam_dtr = field_avg_mam
      field_avg_jja_dtr = field_avg_jja
      field_avg_son_dtr = field_avg_son
    else if (dtr.eq."linear") then
      field_avg_djf_dtr = dtrend_n(field_avg_djf, False, 0)
      tmp = dim_avg_n_Wrap(field_avg_djf, 0)
      tmp_conf = conform(field_avg_djf_dtr, tmp, (/1, 2/))
      field_avg_djf_dtr = field_avg_djf_dtr + tmp_conf
      delete([/tmp, tmp_conf/])
      field_avg_mam_dtr = dtrend_n(field_avg_mam, False, 0)
      tmp = dim_avg_n_Wrap(field_avg_mam, 0)
      tmp_conf = conform(field_avg_mam_dtr, tmp, (/1, 2/))
      field_avg_mam_dtr = field_avg_mam_dtr + tmp_conf
      delete([/tmp, tmp_conf/])
      field_avg_jja_dtr = dtrend_n(field_avg_jja, False, 0)
      tmp = dim_avg_n_Wrap(field_avg_jja, 0)
      tmp_conf = conform(field_avg_jja_dtr, tmp, (/1, 2/))
      field_avg_jja_dtr = field_avg_jja_dtr + tmp_conf
      delete([/tmp, tmp_conf/])
      field_avg_son_dtr = dtrend_n(field_avg_son, False, 0)
      tmp = dim_avg_n_Wrap(field_avg_son, 0)
      tmp_conf = conform(field_avg_son_dtr, tmp, (/1, 2/))
      field_avg_son_dtr = field_avg_son_dtr + tmp_conf
      delete([/tmp, tmp_conf/])
    else if (dtr.eq."quadratic") then
      field_avg_djf_dtr = dtrend_quadratic_msg_n(field_avg_djf, False, \
                                                 False, 0)
      field_avg_mam_dtr = dtrend_quadratic_msg_n(field_avg_mam, False, \
                                                 False, 0)
      field_avg_jja_dtr = dtrend_quadratic_msg_n(field_avg_jja, False, \
                                                 False, 0)
      field_avg_son_dtr = dtrend_quadratic_msg_n(field_avg_son, False, \
                                                 False, 0)
    end if
    end if
    end if
    if (rank.eq.1) then
      out = field_djf(0:3)  ; save metadata
      out(0) = dim_stddev_Wrap(field_avg_djf_dtr)
      out(1) = dim_stddev_Wrap(field_avg_mam_dtr)
      out(2) = dim_stddev_Wrap(field_avg_jja_dtr)
      out(3) = dim_stddev_Wrap(field_avg_son_dtr)
    end if
    if (rank.eq.2) then
      out = field_djf(0:3, :)  ; save metadata
      out(0, :) = dim_stddev_n_Wrap(field_avg_djf_dtr, 0)
      out(1, :) = dim_stddev_n_Wrap(field_avg_mam_dtr, 0)
      out(2, :) = dim_stddev_n_Wrap(field_avg_jja_dtr, 0)
      out(3, :) = dim_stddev_n_Wrap(field_avg_son_dtr, 0)
    end if
    if (rank.eq.3) then
      out = field_djf(0:3, :, :)  ; save metadata
      out(0, :, :) = dim_stddev_n_Wrap(field_avg_djf_dtr, 0)
      out(1, :, :) = dim_stddev_n_Wrap(field_avg_mam_dtr, 0)
      out(2, :, :) = dim_stddev_n_Wrap(field_avg_jja_dtr, 0)
      out(3, :, :) = dim_stddev_n_Wrap(field_avg_son_dtr, 0)
    end if
    if (rank.eq.4) then
      out = field_djf(0:3, :, :, :)  ; save metadata
      out(0, :, :, :) = dim_stddev_n_Wrap(field_avg_djf_dtr, 0)
      out(1, :, :, :) = dim_stddev_n_Wrap(field_avg_mam_dtr, 0)
      out(2, :, :, :) = dim_stddev_n_Wrap(field_avg_jja_dtr, 0)
      out(3, :, :, :) = dim_stddev_n_Wrap(field_avg_son_dtr, 0)
    end if
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Monthly climatology
  if (opt.eq."monthlyclim") then
    out = time_operations(field, y1, y2, "stddev", opt, True)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Month string
  if (.not.ismissing(str_match_ind_ic(monthstr, opt)).and. \
      strlen(opt).ge.2.and.strlen(opt).le.12) then
    field_ext = time_operations(field, y1, y2, "extract", opt, True)
    field_avg = time_operations(field_ext, y1, y2, "average", "yearly", True)
    if (dtr.eq."None") then
      field_avg_dtr = field_avg
    else if (dtr.eq."linear") then
      field_avg_dtr = dtrend_n(field_avg, False, 0)
      tmp = dim_avg_n_Wrap(field_avg, 0)
      tmp_conf = conform(field_avg_dtr, tmp, (/1, 2/))
      field_avg_dtr = field_avg_dtr + tmp_conf
      delete([/tmp, tmp_conf/])
    else if (dtr.eq."quadratic") then
      field_avg_dtr = dtrend_quadratic_msg_n(field_avg, False, False, 0)
    end if
    end if
    end if
    if (rank.eq.1) then
      out = dim_stddev_Wrap(field_avg_dtr)
    else
      out = dim_stddev_n_Wrap(field_avg_dtr, 0)
    end if
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Specific-month
  if (any(opt.eq.tostring(ispan(1, 12, 1)))) then
    out = time_operations(field, y1, y2, "stddev", opt, True)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  error_msg("f", scriptname, funcname, "unrecognized option " + opt)

end

; #############################################################################
undef("calculate_metric")
function calculate_metric(var:numeric,
                          ref:numeric,
                          metric:string)
;
; Arguments
;    var: a 1-D or 2-D numerical array.
;    ref: a numerical array of the same dimensionality of var.
;    metric: a string with the metric to calculate:
;             "RMSD": root-mean square difference.
;             "RMSDxy": root-mean square difference for each grid cell.
;             "BIAS": mean bias.
;             "stddev_ratio": ratio of standard deviations of var and ref
;                             (to be used in Taylor diagram).
;             "correlation": pattern correlation for var and ref
;                            (to be used in Taylor diagram).
;
; Return value
;    A scalar float representing the calculated grading metric.
;
; Description
;    Calculate a grading metrics given two input variables of the same
;    dimensionality.
;
; Modification history
;    20140313-righi_mattia: implemented weights calculation within the
;                           function, depending on dimensionality.
;    20140120-frank_franziska: written.
;
local funcname, scriptname, dims_var, dims_ref, ii, dim_names, mdays, sdays, \
  weights, var1d, ref1d, wgt1d, avg_var, avg_ref, p1, p2, p3, var3d, wgt3d, \
  tmean_m, tmean_o, tmp_m, tmp_o
begin

  funcname = "calculate_metric"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  ; Check for dimensions consistency
  dims_var = dimsizes(var)
  dims_ref = dimsizes(ref)
  if (dimsizes(dims_var).ne.dimsizes(dims_ref)) then
    error_msg("f", scriptname, funcname, "input variables must " + \
              "have the same dimensionality")
  end if
  do ii = 0, dimsizes(dims_var) - 1
    if (dims_var(ii).ne.dims_ref(ii)) then
      error_msg("f", scriptname, funcname, "inconsistent " + \
                "dimension size in input variables for dimension " + ii)
    end if
  end do

  dim_names = getVarDimNames(var)

  mdays = (/31, 28.25, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
  sdays = (/92, 92, 91, 90.25/)

  ; Calculate weights for time case
  if (dimsizes(dim_names).eq.1) then

    ; Monthly climatology
    if (dim_names.eq."month" .and. dims_var.eq.12) then
      weights = mdays
    end if

    ; Seasonal climatology
    if (dim_names.eq."season" .and. dims_var.eq.4) then
      weights = sdays
    end if

  end if

  ; Calculate weights for lat-lon case
  if (dimsizes(dim_names).eq.2) then
    if (dim_names(0).eq."lat" .and. dim_names(1).eq."lon") then
      weights = map_area(var&lat, var&lon)
    end if
  end if

  ; Calculate weights for time-lat-lon case
  if (dimsizes(dim_names).eq.3) then

    ; Monthly climatology
    if (dim_names(0).eq."month" .and. dims_var(0).eq.12) then
      time_weights = mdays
    end if

    ; Seasonal climatology
    if (dim_names(0).eq."season" .and. dims_var(0).eq.4) then
      time_weights = sdays
    end if

    ; Annual-mean time-series
    if (dim_names(0).eq."year") then
      time_weights = new(dims_var(0), float)
      time_weights = 1.
    end if

    if (dim_names(1).eq."lat" .and. dim_names(2).eq."lon") then
      area_weights = map_area(var&lat, var&lon)
    end if

    if (dim_names(1).eq."plev" .and. dim_names(2).eq."lat") then
      areas = map_area(ref&lat, (/1.0, 2.0/))
      nlev = dimsizes(ref&plev)
      ptop = ref&plev(nlev - 1) - \
        0.5 * (ref&plev(nlev - 2) - ref&plev(nlev - 1))
      delta_p = dpres_plevel(ref&plev, 101325., ptop, 0)
      area_weights = new((/nlev, dimsizes(ref&lat)/), float)
      wdims = dimsizes(area_weights)
      area_weights = conform_dims(wdims, delta_p, 0) * \
        conform_dims(wdims, areas(:, 0), 1)
    end if

    if (isdefined("time_weights").and.isdefined("area_weights")) then
      weights = new(dimsizes(var), float)
      do ii = 0, dimsizes(time_weights) - 1
        weights(ii, :, :) = time_weights(ii) * area_weights
      end do
    end if
  end if

  ; Other cases not implemented yet
  if (.not.isdefined("weights")) then
    error_msg("f", scriptname, funcname, "weighting for this " + \
              "variable dimensionality not yet implemented")
  end if

  ; Convert to 1-D arrays
  var1d = ndtooned(var)
  ref1d = ndtooned(ref)
  wgt1d = ndtooned(weights)

  ; optional: no weights --> reset weights
  if (isStrSubset(metric, "nowgt")) then
    wgt1d = 1.0
  end if

  ; Calculate weighted averages
  avg_var = dim_avg_wgt_Wrap(var1d, wgt1d, 1)
  avg_ref = dim_avg_wgt_Wrap(ref1d, wgt1d, 1)

  ; RMSD
  if (metric.eq."RMSD") then
    out = sqrt(dim_avg_wgt_Wrap((var1d - ref1d) ^ 2, wgt1d, 1))
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; BIAS
  if (metric.eq."BIAS") then
    out = dim_avg_wgt_Wrap((var1d - ref1d), wgt1d, 1)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Ratio of standard deviations
  if (isStrSubset(metric, "stddev_ratio")) then
    out = sqrt(dim_sum_wgt_Wrap((var1d - avg_var) ^ 2, wgt1d, 1)) / \
      sqrt(dim_sum_wgt_Wrap((ref1d - avg_ref) ^ 2, wgt1d, 1))
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Pattern correlation
  if (isStrSubset(metric, "correlation")) then
    p1 = dim_sum_wgt_Wrap((var1d - avg_var) * (ref1d - avg_ref), wgt1d, 1)
    p2 = sqrt(dim_sum_wgt_Wrap((var1d - avg_var) ^ 2, wgt1d, 1))
    p3 = sqrt(dim_sum_wgt_Wrap((ref1d - avg_ref) ^ 2, wgt1d, 1))
    out = p1 / (p2 * p3)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Single Model Performance Index
  if (metric.eq."SMPI") then
    nyears = dimsizes(var&year)
    out = new(diag_script_info@smpi_n_bootstrap + 1, float)
    do ibootstrap = 0, diag_script_info@smpi_n_bootstrap
      if (ibootstrap.eq.0) then
        bootvect = ispan(0, nyears - 1, 1)
      else
        icnt = 0
        do while (icnt .le. 10)
          bootvect = generate_sample_indices(nyears, 1)
          icnt = icnt + 1
          if (.not.all(bootvect(:).eq.bootvect(0))) then
            break
          end if
        end do
        if (all(bootvect(:).eq.bootvect(0))) then
          error_msg("f", scriptname, funcname, \
                    "Number of years too small for bootstrapping. Abort.")
        end if
      end if
      obs = ref(bootvect, :, :)
      mod1D = ndtooned(dim_avg_n(var, 0))
      ref1D = ndtooned(dim_avg_n(obs, 0))
      sig1D = ndtooned(dim_stddev_n_Wrap(obs, 0))
      sig1D@_FillValue = default_fillvalue(typeof(sig1D))
      sig1D = where(sig1D.eq.0, sig1D@_FillValue, sig1D)

      delete(weights)
      delete(wgt1d)
      if (isdim(obs, "lon").and.isdim(obs, "lat")) then
        weights = map_area(obs&lat, obs&lon)
      elseif (isdim(obs, "plev").and.isdim(obs, "lat")) then
        areas = map_area(obs&lat, (/1.0, 2.0/))
        nlev = dimsizes(obs&plev)
        ptop = \
          obs&plev(nlev - 1) - 0.5 * (obs&plev(nlev - 2) - obs&plev(nlev - 1))
        delta_p = dpres_plevel(obs&plev, 101325., ptop, 0)
        weights = new((/dimsizes(obs&plev), dimsizes(obs&lat)/), float)
        wdims = dimsizes(weights)
        weights = \
          conform_dims(wdims, delta_p, 0) * conform_dims(wdims, areas(:, 0), 1)
      else
        error_msg("f", diag_script, "", "Unknown dimensions in variable obs.")
      end if

      wgt1d = ndtooned(weights)
      out(ibootstrap) = \
        dim_avg_wgt_Wrap((mod1D - ref1D) ^ 2 / sig1D ^ 2, wgt1d, 1)

    end do
    leave_msg(scriptname, funcname)
    return(out)
  end if

  error_msg("f", scriptname, funcname, "metric " + metric + " not available")

end

; #############################################################################
undef("normalize_metric")
function normalize_metric(var:numeric,
                          opt:string)
;
; Arguments
;    var: numerical array.
;    opt: option determining the used normalization:
;           "max": normalization with max error.
;           "mean": normalization with mean.
;           "median": normalization with median.
;           "stddev_mean": normalization with substracting the mean and
;                          dividing by the standard deviation.
;           "centered_median": substracting and dividing by the median.
;
;
; Return value
;    A numerical array of the same dimensionality as var.
;
; Description
;    Normalizes an array of metrics according to opt.
;
; Caveats
;    Treatment of missing values not explicitely specified (yet).
;
; Reference
;
; Modification history
;    20140609-righi_mattia: absolute value added to "mean" normalization.
;    20140120-frank_franziska: written.
;
local funcname, scriptname, val_var, norm, stdv, p1
begin

  funcname = "normalize_metric"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  val_var = ndtooned(var)

  if (opt.eq."mean") then
    norm = var / dim_avg(abs(val_var))
  end if

  if (opt.eq."median") then
    norm = var / dim_median(val_var)
  end if

  if (opt.eq."stddev_mean") then
    stdv = dim_stddev(val_var)
    if (stdv.gt.0) then
      p1 = ispan(0, dimsizes(dimsizes(var)) - 1, 1)
      norm = dim_rmvmean_n_Wrap(var, p1) / stdv
    else
      error_msg("f", scriptname, funcname, "zero standard deviation")
    end if
  end if

  if (opt.eq."centered_median") then
    norm = (var / dim_median(val_var)) - 1
  end if

  if (opt.eq."maximum") then
    norm = 1. - var / max(var)
  end if

  if (.not.isdefined("norm")) then
    error_msg("f", scriptname, funcname, "no valid normalization defined")

  end if

  copy_VarMeta(var, norm)
  leave_msg(scriptname, funcname)
  return(norm)

end

; #############################################################################
undef("distrib_stats")
function distrib_stats(var[*]:numeric,
                       opt:string)
;
; Arguments
;    var: a one-dimensional input array.
;    opt: type of statistic:
;           "N": number of elements.
;           "mean": mean.
;           "median": median.
;           "min": minimum.
;           "max": maximum.
;           "stddev": standard deviation.
;           [value]: percentile (a value between 0 and 100).
;
; Return value
;    A scalar value.
;
; Description
;    Calculates the relevant statistics for an input one-dimensional
;    distribution. Missing values are ignored.
;
; Caveats
;
; Reference
;
; Modification history
;    20140526-righi_mattia: written.
;
local funcname, scriptname, lvar, nsize, idx
begin

  funcname = "distrib_stats"
  scriptname = "diag_scripts/shared/statistics.ncl"

  enter_msg(scriptname, funcname)

  ; Check input
  if (all(ismissing(var))) then
    error_msg("f", scriptname, funcname, "the input array " + \
              "contains only missing values")
  end if

  ; Exclude missing values
  lvar = var(ind(.not.ismissing(var)))

  ; Number of elements
  if (opt.eq."N") then
    out = dimsizes(lvar)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Minimum
  if (opt.eq."min") then
    out = min(lvar)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Maximum
  if (opt.eq."max") then
    out = max(lvar)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Mean
  if (opt.eq."mean") then
    out = dim_avg(lvar)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Standard deviation
  if (opt.eq."stddev") then
    out = dim_stddev(lvar)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Median
  if (opt.eq."median") then
    out = dim_median(lvar)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  ; Percentiles
  if (toint(opt).gt.0 .and. toint(opt).lt.100) then
    qsort(lvar)  ; sort
    nsize = dimsizes(lvar) / 100.
    idx = round(toint(opt) * nsize, 3) - 1
    idx = where(idx.lt.0, 0, idx)
    idx = where(idx.gt.dimsizes(lvar) - 1, dimsizes(lvar) - 1, idx)
    out = lvar(idx)
    leave_msg(scriptname, funcname)
    return(out)
  end if

  error_msg("f", scriptname, funcname, "unrecognized option " + opt)

end

; #############################################################################
undef("lognormal_dist")
function lognormal_dist(nn:numeric,
                        dg:numeric,
                        sig[1]:numeric,
                        darr[*]:numeric)
;
; Arguments
;    nn: particle number concentration, can be a scalar or 1-D array.
;    dg: median diameter, same dimensionality of nn
;    sig: geometric standard deviation, a scalar
;    darr: array of diameters.
;
; Return value
;    An array of type float, with the same dimensionality of nn, plus the darr
;    dimension on the right, and with the same units of nn.
;
; Description
;    Calculate a lognormal distribution given the three paramters and an array
;    of diameters.
;
; Caveats
;    dg and darr must have the same units.
;
; Reference
;    Seinfeld and Pandis, Atmospheric chemistry and physics, JohnWiley & Sons,
;    New York, US, 1998.
;
; Modification history
;    20130528-righi_mattia: written.
;
local funcname, scriptname, pi, sqrt2pi, dd, ee
begin

  funcname = "lognormal_dist"
  scriptname = "diag_scripts/shared/statistics.ncl"
  enter_msg(scriptname, funcname)

  if (dimsizes(dimsizes(nn)).ne.1 .or. dimsizes(dimsizes(dg)).ne.1) then
    error_msg("f", scriptname, funcname, "input arguments " + \
              "must be scalars or 1-D arrays")
    status_exit(1)
  end if

  rank = dimsizes(nn)

  ; Define output array
  if (rank.eq.1) then
    out = new(dimsizes(darr), double)
  else
    out = new(array_append_record(dimsizes(nn), dimsizes(darr), 0), double)
  end if

  ; Loop over diameter array and calculate size distribution
  pi = 2. * acos(0)
  sqrt2pi = sqrt(2. * pi)
  do dd = 0, dimsizes(darr) - 1
    if (rank.eq.1) then
      ee = (log(darr(dd)) - log(dg)) ^ 2 / (2. * log(sig) ^ 2)
      out(dd) = nn  ; save metadata
      out(dd) = nn / (sqrt2pi * log(sig)) * exp(-ee)
    else
      ee = (log(darr(dd)) - log(dg)) ^ 2 / (2. * log(sig) ^ 2)
      out(:, dd) = nn  ; save metadata
      out(:, dd) = nn / (sqrt2pi * log(sig)) * exp(-ee)
    end if
  end do
  out = where(out.lt.1.d-36, 0., out)
  fout = tofloat(out)
  rank = dimsizes(dimsizes(fout)) - 1
  fout!rank = "diam"
  fout&diam = darr

  leave_msg(scriptname, funcname)
  return(fout)

end

; #############################################################################
undef("filter121")
function filter121(var:numeric, \
                   iter:numeric)
;
; Arguments
;     var: a 1-D array.
;     iter: number of iterations.
;
; Return value
;     Array of size and type of var.
;
; Description
;     Smoothes a time series by 1-2-1 filter in iter time steps.
;
; Modification history
;     20180807-schlund_manuel: ported to v2.0
;     20140721-wenzel_sabrina: written.
;
local yy_sm, y_sm0, n_max
begin

  funcname = "filter121"
  scriptname = "diag_scripts/shared/statistics.ncl"
  enter_msg(scriptname, funcname)

  rank = dimsizes(dimsizes(var))
  n_max = dimsizes(var)

  if (rank .eq. 1) then
    yy_sm = var
    yy_sm0 = yy_sm

    ; Iteration
    do ia = 0, iter
      yy_sm0 = yy_sm
      do it = 1, n_max - 2
        yy_sm(it) = (yy_sm0(it - 1) + 2.0 * yy_sm0(it) + yy_sm0(it + 1)) / 4.0
      end do
    end do
    out = yy_sm
    delete([/yy_sm, yy_sm0/])
  else
    error_msg("f", scriptname, funcname, "smooting is currently not " + \
              "implemented for more than 1 dimension")
  end if

  leave_msg(scriptname, funcname)
  return(out)
end

; #############################################################################
undef("get_average")
function get_average(items: list)

;
; Arguments
;    items: list of input_file_info metadata
;
; Description
;    Calculates average over M different datasets given by items. All given
;    datasets must be of the same shape (X1, X2, X3, ..., Xn). This function
;    creates an intermediate array of shape (M, X1, X2, X3, ..., Xn) with
;    'dataset' as first dimension. Then, the unweighted mean over this first
;    dimension is returned.
;
; Caveats
;
; References
;
; Modification history
;    20181217-schlund_manuel: written
;
local funcname, scriptname, dim_dat, dim_array, rank, data_array, data
begin

  funcname = "get_average"
  scriptname = "diag_scripts/shared/statistics.ncl"
  enter_msg(scriptname, funcname)

  ; Check input list
  dim_dat = ListCount(items)
  if (dim_dat .lt. 1) then
    error_msg("f", scriptname, funcname, "expected at least one dataset, " + \
              "got empty list")
  end if

  ; Collect data
  do idat = 0, dim_dat - 1
    data = read_data(items[idat])

    ; Create array at first iteration
    if (.not. isvar("data_array")) then
      dim_array = dimsizes(data)
      rank = dimsizes(dim_array)
      data_array = new(array_append_record(dim_dat, dim_array, 0), \
                       typeof(data))
    end if
    if (rank .eq. 4) then
      copy_VarCoords(data, data_array(0, :, :, :, :))
      data_array(idat, :, :, :, :) = (/data/)
    elseif (rank .eq. 3) then
      copy_VarCoords(data, data_array(0, :, :, :))
      data_array(idat, :, :, :) = (/data/)
    elseif (rank .eq. 2) then
      copy_VarCoords(data, data_array(0, :, :))
      data_array(idat, :, :) = (/data/)
    elseif (rank .eq. 1) then
      copy_VarCoords(data, data_array(0, :))
      data_array(idat, :) = (/data/)
    else
      error_msg("f", scriptname, funcname, "unsupported rank " + rank)
    end if
  end do

  ; Average
  data_array!0 = "dataset"
  data_array := dim_avg_n_Wrap(data_array, 0)
  copy_VarMeta(data, data_array)
  log_info("Averaged over " + dim_dat + " dataset(s)")

  leave_msg(scriptname, funcname)
  return(data_array)

end
