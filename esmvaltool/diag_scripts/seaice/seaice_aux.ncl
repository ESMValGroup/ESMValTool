; #############################################################################
; # EXTERNAL SEA ICE FUNCTIONS
; #############################################################################
;  Please consider using of extending existing routines before adding new ones.
;  Check the header of each routine for documentation.
;
;  Contents:
;     function calc_yod
;     function sea_ice_area
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/latlon.ncl"

; #############################################################################
undef("sea_ice_area")
function sea_ice_area(sic[*][*][*]:numeric,
                      area[*][*]: numeric,
                      arealat: numeric,
                      l_irreg[1]: logical,
                      region[1]:string,
                      start_year[1]:integer,
                      end_year[1]:integer,
                      month[1]:string,
                      operation[1]:string,
                      fill_pole_hole[1]:logical)

;  Arguments:
;     sic            : sea ice concentration variable
;     area           : area variable
;     arealat        : array with the latitude coordinate of the area
;     l_irreg        : logical for irregular grids
;     region         : "Arctic" or "Antarctic"
;     start_year     : first year averaged time period
;     end_year       : last year of averaged time period
;     month          : [1, 12] = select month (e.g., "3" = Mar, "9" = Sep)
;                      "A" = annual mean
;                      "all" = no average, take all input
;     operation      : "area" or "extent"
;     fill_pole_hole : True = fill northern polar hole with sic = 1
;
;  Return value:
;     1-d array (time) containing sea ice area or sea ice extent
;
;  Description:
;     Calculates sea ice area or extent from sea ice concentration
;
;  Caveats
;
;  Modification history:
;     20180226-senftleben_daniel: written.
;
local funcname, scriptname, verbosity, l_irreg, aux_vars, area_temp, date, \
  month_ind, sic_month, check_month, smonths, month_string, latstep, max_lat, \
  lat1D, actual_lat, data1D, i1D, i2D, dim_i2D, npts, ring, SIA_tmp, SIA_tmp2
begin

  funcname = "sea_ice_area"
  scriptname = "diag_scripts/seaice/seaice_aux.ncl"
  enter_msg(scriptname, funcname)

  ;  Check region and units
  if (region.ne."Arctic" .and. region.ne."Antarctic") then
    error_msg("f", scriptname, funcname, \
              "region can be either 'Arctic' or 'Antarctic'")
  end if
  if (sic@units.eq."%") then
    sic = sic / 100.
    sic@units = "1"
  end if

  ;  Fill pole hole
  if (fill_pole_hole .and. region.eq."Arctic") then
    if (.not.l_irreg) then
      lat1D = ndtooned(conform(sic(0, :, :), arealat, 0))
    else
      lat1D = ndtooned(arealat)  ; already conformal by definition
    end if
    latstep = 0.1
    max_lat = 90.

    do itime = 0, dimsizes(sic&time) - 1  ; hole depends on time
      actual_lat = max_lat
      sic1D = ndtooned(sic(itime, :, :))
      if (all(ismissing(sic1D))) then  ; time step missing
        continue
      end if
      do while (True)  ; will break when hole is filled
        idx1D = ind(lat1D.ge.actual_lat .and. lat1D.lt.actual_lat + latstep)
        if (.not. all(ismissing(idx1D))) then
          idx2D = ind_resolve(ind(lat1D.ge.actual_lat), dimsizes(sic(0, :, :)))
          dim_idx2D = dimsizes(idx2D)
          npts = dim_idx2D(0)
          if (any(ismissing(sic1D(idx1D)))) then
            do n = 0, npts - 1
              ring = sic(itime, idx2D(n, 0), idx2D(n, 1))
              sic(itime, idx2D(n, 0), idx2D(n, 1)) = \
                where(ismissing(ring), 1., ring)
              delete(ring)
            end do
          else  ; out of the hole!
            delete([/idx1D, idx2D, dim_idx2D, npts/])
            break
          end if
          delete([/idx2D, dim_idx2D, npts/])
        end if
        delete([/idx1D/])
        actual_lat = actual_lat - latstep
      end do  ; while
      delete(sic1D)
    end do  ; itime
  end if

  ;  Calculate time coordinate
  date = cd_calendar(sic&time, 0)  ; Create calendar array
  if (month.eq."A") then
    sic_month = time_operations(sic, start_year, end_year, \
                                "average", "yearly", True)
    sic_month!0 = "time"
    month_string = "Annual mean "
  elseif (month.eq."all") then
    sic_month = sic
    month_string = ""
  else  ; "month" must be a number between 1 and 12
    check_month = tostring(ispan(1, 12, 1))
    if (.not.any(check_month.eq.month)) then
      error_msg("f", scriptname, funcname, "month must be between 1 and 12")
    end if
    smonths = (/"January", "February", "March", "April", "May", \
                "June", "July", "August", "September", "October", \
                "November", "December"/)
    month_string = smonths(ind(check_month.eq.month)) + " "
    month_ind = ind(date(:, 1).eq.toint(month))
    sic_month = sic(month_ind, :, :)
  end if

  ; Calculate sea ice area or extent
  if (operation.eq."extent") then
    sic_month = where(sic_month.ge.0.15, 1., 0.)  ; make sic binary for SIE
  end if

  if (l_irreg) then
    if (region.eq."Arctic") then
      area_tmp = where(arealat.gt.0., area * 1.e-12, 0.)
    else
      area_tmp = where(arealat.lt.0., area * 1.e-12, 0.)
    end if
    SIA_tmp  = conform(sic_month, area_tmp, (/1, 2/)) * sic_month
    copy_VarMeta(sic_month, SIA_tmp)
    SIA_tmp2 = dim_sum_Wrap(SIA_tmp)
    SIA      = dim_sum_Wrap(SIA_tmp2)
  else
    area_tmp = area / 1.e12  ; Million square km
    SIA_tmp  = conform(sic_month, area_tmp, (/1, 2/)) * sic_month
    copy_VarMeta(sic_month, SIA_tmp)
    if (region.eq."Arctic") then
      SIA = area_operations(SIA_tmp, 0, 90, 0, 360, "sum", False)
    else
      SIA = area_operations(SIA_tmp, -90, 0, 0, 360, "sum", False)
    end if
  end if

  SIA@units = "10~S~6~N~ km~S~2~N~"
  if (operation.eq."extent") then
    SIA@long_name = month_string + region + " Sea Ice Extent"
    SIA@standard_name = "sea_ice_extent"
  elseif (operation.eq."area") then
    SIA@long_name = month_string + region + " Sea Ice Area"
    SIA@standard_name = "sea_ice_area"
  end if

  return(SIA)

end

; #############################################################################
undef("calc_yod")
function calc_yod(SIE:float,
                  threshold[1]:float,
                  cons_years[1]:integer)
;
; Arguments:
;   SIE:        sea ice extent (one value per year), must be of type float
;   threshold:  threshold value of sea ice extent in million km^2 (usually 1.)
;   cons_years: number of consecutive years in which SIE must be < threshold
;
; Return value:
;   first year of cons_years years in which SIE < threshold, type of SIE
;
; Description:
;   Calculates year of near disappearance of (Arctic) sea ice (YOD). YOD is
;     defined as the first year of five consecutive years in which September
;     Arctic sea ice extent is below 1 million square kilometres. Returns
;     missing value if no YOD occurs in the given timeseries.
;     This function allows the user to specify a different threshold and
;     number of consecutive years.
;
; Caveats
;   (September Arctic) sea ice extent as input expected, must be calculated
;     from sic beforehand. Input must be of one-dimensional: 'time' or 'year'.
;   If function returns a missing value, it is set to be of type float. Make
;     sure that your script recognizes the output correctly as missing.
;
; Modification history:
;    20180827-senftleben_daniel: written.
;
local funcname, scriptname, verbosity, SIE_smoo, a, itime, iyod

begin

  funcname = "calc_yod"
  scriptname = "diag_scripts/seaice/seaice_aux.ncl"
  enter_msg(scriptname, funcname)

  ; check input
  if (SIE!0.ne."year") then
    if (SIE!0.eq."time") then  ; need 'year'
      utc_date = cd_calendar(SIE&time, 0)
      SIE!0    = "year"
      SIE&year = utc_date(:, 0)
      delete(utc_date)
    else
      error_msg("f", scriptname, "", \
                "Input array must have time coordinate" + \
                " 'time' or 'year'")
    end if
  end if

  ; calculate 5-year running mean of input data
  SIE_smoo = runave(SIE, 5, 0)
  copy_VarMeta(SIE, SIE_smoo)
  copy_VarCoords(SIE, SIE_smoo)

  ; calculate YOD
  a = 0
  do itime = 0, dimsizes(SIE_smoo&year) - 1
    if (ismissing(SIE_smoo(itime))) then
      a = 0
      continue
    end if
    if (SIE_smoo(itime).lt.threshold) then
      a = a + 1
      if (a.eq.cons_years) then
        iyod = itime - cons_years + 1
        break
      else
        continue
      end if
    else
      a = 0
    end if
  end do

  ; return YOD
  if (isdefined("iyod")) then
    return_value = tofloat(SIE_smoo&year(iyod))
  else
    return_value = SIE@_FillValue
  end if

  return(return_value)
end
