; #############################################################################
; GENERAL ROUTINES FOR OPERATIONS ON A (RECTILINEAR) GRID
; #############################################################################
; Please consider using of extending existing routines before adding new ones.
; Check the header of each routine for documentation.
;
; Contents:
;    function roi
;    function extract_area
;    function gridcell_area
;    function map_area
;    function area_operations
;    function select_region
;    function make_latlon2D
;    function cdo_remapdis
;
; #############################################################################


load "$diag_scripts/../interface_scripts/auxiliary.ncl"
load "$diag_scripts/../interface_scripts/constants.ncl"
load "$diag_scripts/../interface_scripts/data_handling.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"
; load "./reformat_scripts/obs/reformat_obs_func.ncl"  ; write_nc ; FIX-ME

; #############################################################################
undef("roi")
function roi(latlon_dat[4]:float,
             latlon_roi[4]:float)
;
; Arguments
;    latlon_dat: (/latmin,latmax,lonmin,lonmax/) of the (data) region to check
;    latlon_roi: (/latmin,latmax,lonmin,lonmax/) of the region of interest.
;                Ranges: -90 < lat < 90, 0 < lon < 360
;
; Return value
;    A string array containing all true statements from the following list:
;    contained / center / overlap / encloses / outside (e.g. "contained"
;    implies "center" and "overlap")
;
; Description
;    Checks if the area described by latlon_dat is contained in / has its
;    center in / overlaps with / encloses / is outside the area described by
;    latlon_roi, for use with Emmons.ncl
;
; Caveats
;    Not fully tested for areas containing Greenwich or the date line
;
; References
;
; Modification history
;    20140129-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, latlon_cases, coord_checks, datmin, \
  datmax, roimin, roimax, latlon_dat, latlon_roi, latmin, latmax, lonmin, \
  lonmax, result, result_all, a1, a3
begin

  funcname = "roi"
  scriptname = "diag_scripts/shared/latlon.ncl"
  enter_msg(scriptname, funcname)

  log_debug("  roi: latmin=" + latlon_roi(0) + \
            " ; latmax=" + latlon_roi(1) + \
            " ; lonmin=" + latlon_roi(2) + \
            " ; lonmax=" + latlon_roi(3))
  log_debug("  dat: latmin=" + latlon_dat(0) + \
            " ; latmax=" + latlon_dat(1) + \
            " ; lonmin=" + latlon_dat(2) + \
            " ; lonmax=" + latlon_dat(3))

  ; Avoid errors for areas that contain the Greenwich meridian
  ; -> switch everything to -180 < lon < 180 for that case
  if (latlon_dat(2).gt.latlon_dat(3) .or. latlon_roi(2).gt.latlon_roi(3)) then
    if (latlon_dat(2).gt.180.) then
      latlon_dat(2) = 180. - latlon_dat(2)
    end if
    if (latlon_dat(3).gt.180.) then
      latlon_dat(3) = 180. - latlon_dat(3)
    end if
    if (latlon_roi(2).gt.180.) then
      latlon_roi(2) = 180. - latlon_roi(2)
    end if
    if (latlon_roi(3).gt.180.) then
      latlon_roi(3) = 180. - latlon_roi(3)
    end if
  end if

  ; Check each lat & lon individually for the strongest statement
  ; from this list contained / overlap / encloses / outside
  latlon_cases = new(2, string)  ; overall result for lat & lon
  coord_checks = new(4, logical)  ; results of the 4 checks for 1 coordinate
  do i = 0, 2, 2
    datmin = latlon_dat(i)
    datmax = latlon_dat(i + 1)
    roimin = latlon_roi(i)
    roimax = latlon_roi(i + 1)
    coord_checks = (/datmin.le.roimin, datmin.le.roimax, \
                    datmax.le.roimin, datmax.le.roimax/)
    if (all(coord_checks.eq.(/False, True, False, True/))) then
        latlon_cases(i/2) = "contained"
    else if (all(coord_checks.eq.(/True, True, False, False/))) then
      latlon_cases(i/2) = "enclosing"
    else if (all(coord_checks.eq.True).or.all(coord_checks.eq.False)) then
      latlon_cases(i/2) = "outside"
    else
      latlon_cases(i/2) = "overlap"
    end if
    end if
    end if
  end do

  ; Combine both coordinate directions into one statement
  if (any(latlon_cases.eq."outside")) then
    result = "outside"
  else if (latlon_cases(0).eq.latlon_cases(1)) then
    result = latlon_cases(0)
  else
    result = "overlap"
  end if
  end if

  ; Check "overlap" if the majority of latlon_dat is within latlon_roi
  if (result.eq."overlap") then
    ; area of latlon_dat
    a1 = gridcell_area(latlon_dat(3) - latlon_dat(2), \
                       latlon_dat(0), latlon_dat(1))
    ; area of overlap with latlon_roi
    latmin = max((/latlon_dat(0), latlon_roi(0)/))
    latmax = min((/latlon_dat(1), latlon_roi(1)/))
    lonmin = max((/latlon_dat(2), latlon_roi(2)/))
    lonmax = min((/latlon_dat(3), latlon_roi(3)/))
    a3 = gridcell_area(lonmax - lonmin, latmin, latmax)
    ; compare
    if (2. * a3.gt.a1) then
      result = "center"
    end if
    ; log_debug("warning: check for degree of overlap " + \
    ;             "not implemented yet")
    ; log_debug("warning: -> returning 'center' just for case " + \
    ;             "'contained'")
  end if
  ; add other true, but less strong statements to the output vector
  if (result.eq."contained") then
    result_all = (/"contained", "center", "overlap"/)
  else if (result.eq."center") then
    result_all = (/"center", "overlap"/)
  else if (result.eq."encloses") then
    result_all = (/"encloses", "center", "overlap"/)
  else if (result.eq."overlap") then
    result_all = (/"overlap"/)
  else
    result_all = (/"outside"/)
  end if
  end if
  end if
  end if

  leave_msg(scriptname, funcname)
  return(result_all)

end

; #############################################################################
undef("extract_area")
function extract_area(index[1]:integer,
                      data_pointer[1]:logical,
                      requested_param[1]:string,
                      parent_var[1]:string)
;
; Arguments
;    index: see interface_scripts/data_handling.ncl.
;    data_pointer: see interface_scripts/data_handling.ncl.
;    requested_param: the parameter to fetch for parent_var.
;    parent_var: name of the variable for which the requested params belong.
;
; Return value
;    The field corresponding to the variable var.
;
; Description
;    Extracts additional data from the preproc-file.
;    The parameters for irregular grids (e.g. area,lat,lon) are written to
;    preproc-files by the reformat routines.
;    The functionality is similar to "extract_data", but reduced (e.g.
;    fetching of var sub-sets is not implemented here).
;    This new function was written to avoid modifying "extract_data", which is
;    designed to extract the target variable only.
;
; Caveats
;
; References
;
; Modification history
;    2013????-gottschaldt_klaus-dirk: written.
;
local funcname, scriptname, cn, cnt, datasetfile, fm, val
begin

  funcname = "extract_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  enter_msg(scriptname, funcname)

  ; 'cn' is the basic building block for the keys used in
  ; the 'data_pointer' dictionary
  ; 'dictkeys' is loaded from 'interface_data/ncl.interface'
  cn = interface_get_dictkeys(parent_var, index)
  cnt = cn + "_filepath"
  datasetfile = data_pointer@$cnt$
  print("datasetfile = " + datasetfile)
  fm = addfile(datasetfile, "r")
  val = fm->$requested_param$   ; Fetch entire variable, including metadata

  leave_msg(scriptname, funcname)
  return(val)

end

; #############################################################################
undef("gridcell_area")
function gridcell_area(deltax[1]: numeric,
                       lat_lo[1]: numeric,
                       lat_hi[1]: numeric)
;
; Arguments
;    deltax: longitude resolution [deg].
;    lat_lo: lower limit of the box [deg].
;    lat_hi: upper limit of the box [deg].
;
; Return value
;    The area of the element in units of [m^2].
;
; Description
;    Calculates the area of a grid cell on the sphere.
;
; Modification history
;    20121211-righi_mattia: written.
;
local funcname, scriptname, deg2rad
begin

  funcname = "gridcell_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  enter_msg(scriptname, funcname)

  ; Convert to radiants
  deg2rad = acos(-1.0) / 180.

  deltay = sin(lat_hi * deg2rad) - sin(lat_lo * deg2rad)
  out = abs(EARTHradius ^ 2 * deltay * deltax * deg2rad)
  out@units = EARTHarea@units

  leave_msg(scriptname, funcname)
  return(out)

end

; #############################################################################
undef("map_area")
function map_area(lat[*]:numeric,
                  lon[*]:numeric)
;
; Arguments
;    lat: the latitude coordinate of the map [deg].
;    lon: the longitude coordinate of the map [deg].
;
; Return value
;    A 2D lat-lon area with the area of each gridbox in units of [m^2].
;
; Description
;    Calculates the area of each grid cell on a global map.
;
; Caveats
;    Assumes a constant resolution in longitude.
;
; Modification history
;    20140819-righi_mattia: modified to support non-global input.
;    20121211-righi_mattia: written.
;
local funcname, scriptname, loclat, loclon, nlat, lati, xresol
begin

  funcname = "map_area"
  scriptname = "diag_scripts/shared/latlon.ncl"
  enter_msg(scriptname, funcname)

  ; Type conversion
  loclat = tofloat(lat)
  loclon = tofloat(lon)
  nlat = dimsizes(loclat)

  ; Check ordering (must be S-N)
  if (isMonotonic(loclat).ne.1) then
    error_msg("f", scriptname, funcname, "input latitude is not " + \
              "monotonically increasing")
  end if

  ; Calculate longitude resolution (assumed constant)
  xresol = abs(loclon(1) - loclon(0))
  dsize = dimsizes(loclon)
  if (dsize(0).gt.2) then
    if ((xresol - abs(loclon(2) - loclon(1))).gt.1.0e-4) then
      error_msg("f", scriptname, funcname, "input longitude " + \
                "resolution is not constant")
    end if
  end if

  ; Calculate latitude interfaces
  lati = new(dimsizes(loclat) + 1, float)
  lati(0) = max((/(3 * loclat(0) - loclat(1)) / 2., -90./))
  do ii = 1, dimsizes(lati) - 2
    lati(ii) = 0.5 * (loclat(ii - 1) + loclat(ii))
  end do
  lati(dimsizes(lati) - 1) = \
    min((/(3 * loclat(nlat - 1) - loclat(nlat - 2)) / 2., 90./))

  ; Calculate area
  area = new((/dimsizes(loclat), dimsizes(loclon)/), float)
  do ii = 0, dimsizes(loclat) - 1
    area(ii, :) = gridcell_area(xresol, lati(ii), lati(ii + 1))
  end do

  ; Assign coordinates
  area!0   = "lat"
  area&lat = lat
  area!1   = "lon"
  area&lon = lon

  leave_msg(scriptname, funcname)
  return(area)

end

; #############################################################################
undef("area_operations")
function area_operations(field:numeric,
                         latmin[1]:numeric,
                         latmax[1]:numeric,
                         lonmin[1]:numeric,
                         lonmax[1]:numeric,
                         opt[1]:string,
                         l_wgt[1]:logical)
;
; Arguments
;    field: a numeric array of rank at least 2; second-to-last and last.
;          dimension must be lat and lon, respectively.
;    latmin: minimum latitude boundary of the region to be selected.
;    latmax: maximum latitude boundary of the region to be selected.
;    lonmin: minimum longitude boundary of the region to be selected.
;    lonmax: maximum longitude boundary of the region to be selected.
;    opt: type of operation:
;           "extract": extracts selected region.
;           "average": averages over the selected region.
;           "sum": integrate over the selected region.
;    l_wgt: if True, calculates area-weighted average/sum (has no effect for
;           opt = "extract").
;
; Return value
;    An array of the same rank as field, of rank-1 or of rank-2, depending on
;    opt and on the region boundaries.
;
; Description
;    Extracts a selected region or point on a global map.
;    Performs the (weighted) average over a selected region on a global map.
;    Performs the (weighted) sum over a selected region on a global map.
;
; Caveats
;    This function assumes that the input field is a global map.
;    Mind the order for lonmin and lonmax (e.g., 60,120 is different from.
;    120,60: the order is meant eastwards).
;    To consider the global domain, use -90,90,0,360 as arguments.
;    If lonmin is outside the field boundaries, it is assigned to 0.
;    If lonmax is outside the field boundaries, it is assigned to max(lon).
;    For latmin=latmax and lonmin=lonmax the single gridbox is extracted, no
;    average/sum is possible in this case.
;
; Modification history
;    20140116-righi_mattia: written.
;
local funcname, scriptname, rank, dimnames, wgt, conwgt, f_latmin, f_latmax, \
  f_lonmin, f_lonmax, delta, idx_latmin, idx_latmax, latidx, idx_lonmin, \
  idx_lonmax, lonidx, subfield, subwgt, dnames, didx_r, didx_l, norm
begin

  funcname = "area_operations"
  scriptname = "diag_scripts/shared/latlon.ncl"
  enter_msg(scriptname, funcname)

  ; Check rank
  if (dimsizes(dimsizes(field)) .lt. 2) then
    error_msg("f", scriptname, funcname, "rank of input field <2")
  end if

  ; Check dimensions
  rank = dimsizes(dimsizes(field))
  dimnames = getvardims(field)
  if (dimnames(rank - 1) .ne. "lon") then
    error_msg("f", scriptname, funcname, "the rightmost dimension " + \
              "of the input field must be lon")
  end if
  if (dimnames(rank - 2) .ne. "lat") then
    error_msg("f", scriptname, funcname, "the second rightmost " + \
              "dimension of the input field must be lat")
  end if
  delete(dimnames)

  ; Check that grid is regular
  if (dimsizes(dimsizes(field&lat)) .ne. 1 .or. \
      dimsizes(dimsizes(field&lon)) .ne. 1) then
    error_msg("f", scriptname, funcname, "irregular grids not supported")
  end if

  ; Calculate weights
  wgt = map_area(field&lat, field&lon)
  if (.not.l_wgt) then
    wgt = 1.
  end if

  ; Conform weights and filter for missing values in the input field
  if (rank.eq.4) then
    conwgt = conform(field, wgt, (/2, 3/))
  end if
  if (rank.eq.3) then
    conwgt = conform(field, wgt, (/1, 2/))
  end if
  if (rank.eq.2) then
    conwgt = wgt
  end if
  conwgt = where(ismissing(field), conwgt@_FillValue, conwgt)
  copy_VarCoords(field, conwgt)
  delete(wgt)

  ; Copy arguments to local (swap lat limits if necessary)
  f_latmin = todouble(min((/latmin, latmax/)))
  f_latmax = todouble(max((/latmin, latmax/)))
  f_lonmin = todouble(where(lonmin .lt. 0., lonmin + 360., lonmin))
  f_lonmax = todouble(where(lonmax .lt. 0., lonmax + 360., lonmax))

  ; If single gridbox lon = 360, replace with 0 to account for periodicity
  if (f_lonmin .eq. 360 .and. f_lonmax .eq. 360) then
    f_lonmin = 0.d0
    f_lonmax = 0.d0
  end if

  ; Find nearest points to latitude boundaries
  ; In the (unlikely) case that the selected boundary point is at the
  ; midpoint between two gridboxes, ind will return two points.
  ; The minimum index (= minumum value) is chosen in such case,
  ; consistently with NCL coordinate subscripting
  ; http://www.ncl.ucar.edu/Document/Language/subscript_coord.shtml
  delta = abs(field&lat - f_latmin)
  idx_latmin = min(ind(delta .eq. min(delta)))
  delta = abs(field&lat - f_latmax)
  idx_latmax = min(ind(delta .eq. min(delta)))
  delete(delta)
  if (ismissing(idx_latmin) .or. ismissing(idx_latmax)) then
    error_msg("f", scriptname, funcname, "no values found in the " + \
              "given latitude range")
  end if
  if (latmin .eq. latmax .and. idx_latmin .ne. idx_latmax) then
    error_msg("f", scriptname, funcname, "multiple indexes found " + \
              "for single gridbox latitude selection")
  end if

  ; Define latitude indexes
  latidx = ispan(idx_latmin, idx_latmax, 1)  ; array of indexes
  delete(idx_latmin)
  delete(idx_latmax)

  ; Find neartest points to longitude boundaries
  ; In the (unlikely) case that the selected boundary point is at the
  ; midpoint between two gridboxes, ind will return two points.
  ; The minimum index (= minumum value) is chosen in such case,
  ; consistently with NCL coordinate subscripting
  ; http://www.ncl.ucar.edu/Document/Language/subscript_coord.shtml
  if (f_lonmin.gt.max(field&lon)) then
    idx_lonmin = 0
  else
    delta = abs(field&lon - f_lonmin)
    idx_lonmin = min(ind(delta .eq. min(delta)))
    delete(delta)
  end if
  if (f_lonmax.gt.max(field&lon)) then
    idx_lonmax = dimsizes(field&lon) - 1
  else
    delta = abs(field&lon - f_lonmax)
    idx_lonmax = min(ind(delta .eq. min(delta)))
    delete(delta)
  end if
  if (ismissing(idx_lonmin) .or. ismissing(idx_lonmax)) then
    error_msg("f", scriptname, funcname, "no values found in the " + \
              "given longitude range")
  end if

  ; Special case: single gridbox selection beteen max(field&lon) and 0
  if (lonmin .eq. lonmax .and. lonmin.gt.max(field&lon)) then
    idx_lonmin = where(abs(360. - lonmin) .lt. abs(max(field&lon) - lonmin), \
                       0, dimsizes(field&lon) - 1)
    idx_lonmax = idx_lonmin
  end if
  if (lonmin .eq. lonmax .and. idx_lonmin .ne. idx_lonmax) then
    error_msg("f", scriptname, funcname, "multiple indexes found " + \
              "for single gridbox longitude selection")
  end if

  ; Define longitude indexes
  if (idx_lonmin .le. idx_lonmax) then  ; Greenwich outside
    lonidx = ispan(idx_lonmin, idx_lonmax, 1)
  else  ; Greenwich inside
    lonidx = \
      array_append_record(ispan(idx_lonmin, dimsizes(field&lon) - 1, 1), \
                          ispan(0, idx_lonmax, 1), 0)
  end if
  delete(idx_lonmin)
  delete(idx_lonmax)

  ; Extract region
  if (rank .eq. 4) then
    subfield = field(:, :, latidx, lonidx)
    subwgt = conwgt(:, :, latidx, lonidx)
  end if
  if (rank .eq. 3) then
    subfield = field(:, latidx, lonidx)
    subwgt = conwgt(:, latidx, lonidx)
  end if
  if (rank .eq. 2) then
    subfield = field(latidx, lonidx)
    subwgt = conwgt(latidx, lonidx)
  end if
  delete(conwgt)

  if (opt .eq. "extract") then
    leave_msg(scriptname, funcname)
    return(subfield)
  end if

  ; Find dimensions to be averaged or summed
  dnames = getVarDimNames(subfield)
  didx_r = ind(dnames.eq."lat" .or. dnames.eq."lon")
  didx_l = ind(dnames.ne."lat" .and. dnames.ne."lon")
  if (all(ismissing(didx_r))) then
    error_msg("f", scriptname, funcname, "cannot perform " + \
              "average/sum on single gridcell, use the extract option")
  end if

  ; Normalize weights
  if (opt .eq. "average") then
    if (dimsizes(didx_r).eq.dimsizes(dnames)) then
      subwgt = subwgt / sum(subwgt)
    else
      norm = dim_sum_n_Wrap(subwgt, didx_r)
      subwgt = subwgt / conform(subwgt, norm, didx_l)
    end if
  end if

  ; Calculate average/sum. The same code applies, as the only difference
  ; is in the normalization of the weights, which is calculated above
  ; only for opt="average"
  if (opt .eq. "average" .or. opt .eq. "sum") then
    subfield = subfield * subwgt
    if (dimsizes(didx_r).eq.dimsizes(dnames)) then
      out = dim_sum_Wrap(dim_sum_Wrap(subfield))
    else
      out = dim_sum_n_Wrap(subfield, didx_r)
    end if
    leave_msg(scriptname, funcname)
    return(out)
  end if

  error_msg("f", scriptname, funcname, "unrecognized option " + opt)

end

; #############################################################################
undef("select_region")
function select_region(region:string)
;
; Arguments
;    region: a string specifying the region to be selected.
;
; Return value
;    An array with the region boundary as (latmin, latmax, lonmin, lonmax)
;    with the name of the region as a string attribute @name.
;
; Description
;    Translates a region specification into lat/lon boundaries and a region
;    name as an attribute.
;
; Modification history
;    20190405-righi_mattia: change region names to match provenance.
;    20141205-gottschaldt_klaus-dirk: adjusted names to Righi et al. (2015).
;    20140410-frank_franziska: extended to midlat, equatorial and polar.
;    20140129-frank_franziska: written.
;
local funcname, scriptname, region
begin

  funcname = "select_region"
  scriptname = "diag_scripts/shared/latlon.ncl"
  enter_msg(scriptname, funcname)

  if (region.eq."global") then
    region_array = (/-90., 90., 0., 360./)
    region_array@name = "Global"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."trop") then
    region_array = (/-20., 20., 0., 360./)
    region_array@name = "Tropics"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."nhext") then
    region_array = (/20., 90., 0., 360./)
    region_array@name = "Northern extratropics"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."shext") then
    region_array = (/-90., -20., 0., 360./)
    region_array@name = "Southern extratropics"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."nhtrop") then
    region_array = (/0., 30., 0., 360./)
    region_array@name = "Northern tropics"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."shtrop") then
    region_array = (/-30., 0., 0., 360./)
    region_array@name = "Southern tropics"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."nh") then
    region_array = (/20., 90., 0., 360./)
    region_array@name = "Northern hemisphere"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."sh") then
    region_array = (/-90., -20., 0., 360./)
    region_array@name = "Southern hemisphere"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."nhmidlat") then
    region_array = (/35., 60., 0., 360./)
    region_array@name = "Northern midlatitudes"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."shmidlat") then
    region_array = (/-60., -35., 0., 360./)
    region_array@name = "Southern midlatitudes"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."nhpolar") then
    region_array = (/60., 90., 0., 360./)
    region_array@name = "Northern polar"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."shpolar") then
    region_array = (/-90., -60., 0., 360./)
    region_array@name = "Southern polar"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  if (region.eq."eq") then
    region_array = (/-10., 10., 0., 360./)
    region_array@name = "Equatorial"
    leave_msg(scriptname, funcname)
    return(region_array)
  end if

  error_msg("f", scriptname, funcname, "unrecognized region " + region)

end

; #############################################################################
undef("make_latlon2D")
function make_latlon2D(var[1]:string,
                       index[1]:integer,
                       data_pointer[1]:logical,
                       field:numeric)
;
; Arguments:
;    var             = name of the variable
;    index           = see interface_scripts/data_handling.ncl
;    data_pointer    = see interface_scripts/data_handling.ncl
;    field           = a numeric array of rank 3 (time, lon, lat) or
;                      rank 4 (time, [p]lev, lon, lat)
;
; Return value:
;    the same field, but with 2D latitude and longitude
;
; Description:
;    Get 2d lat & lon coordinates for variables (see e.g.
;    http://www.ncl.ucar.edu/Applications/Scripts/ice_3.ncl)
;
; Modification history
;    20141023-vanulft_bert: written based on code in SeaIce_polcon_diff
;
local funcname, scriptname, area_lat, area_lon, area_temp, out_field, rank
begin

  funcname = "make_latlon2D"
  scriptname = "diag_scripts/shared/latlon.ncl"
  enter_msg(scriptname, funcname)

  if (all(iscoord(field, (/"lat", "lon"/)))) then

    log_debug(var + " is on a lat-lon grid")

    ; Check rank
    rank = dimsizes(dimsizes(field))
    if (rank .lt. 3 .or. rank .gt. 4) then
      error_msg("f", "latlon.ncl", funcname, \
                "rank of input variable must be 3 or 4")
    end if

    ; Order coordinates
    if (rank .eq. 3) then
      out_field = field(time|:, lon|:, lat|:)
    end if
    if (rank .eq. 4) then
      if (any(field .eq. "plev")) then
        out_field = field(time|:, plev|:, lon|:, lat|:)
      else
        out_field = field(time|:, lev|:, lon|:, lat|:)
      end if
    end if

    ; Create 2d lat & lon
    area_temp = map_area(out_field&lat, out_field&lon)

    ; Make 2D to be compatible with irregular grid handling
    area_lat = conform_dims( \
      (/dimsizes(out_field&lon), dimsizes(out_field&lat)/), out_field&lat, 1)
    area_lon = conform_dims( \
      (/dimsizes(out_field&lon), dimsizes(out_field&lat)/), out_field&lon, 0)
    area_lat!0 = "i"
    area_lat&i = out_field&lon
    area_lat!1 = "j"
    area_lat&j = out_field&lat
    out_field!(rank - 1) = "i"
    out_field!(rank - 2) = "j"
    copy_VarMeta(area_lat, area_lon)
    delete(area_temp)
  else

    ; irregular grid
    log_debug(var + " is not on a lat-lon grid")
    out_field = field

    ; Get lat-lon from another file
    ; (dataset specific: ~/reformat_scripts folder)
    area_lon = extract_area(index, data_pointer, "lon", var)
    area_lat = extract_area(index, data_pointer, "lat", var)

  end if

  out_field@lat2d = area_lat
  out_field@lon2d = area_lon
  delete(area_lat)
  delete(area_lon)

  leave_msg(scriptname, funcname)
  return(out_field)

end

; #############################################################################
undef("cdo_remapdis")
function cdo_remapdis(var[1]:string,
                      field:numeric,
                      index[1]:integer,
                      data_pointer[1]:logical,
                      dst_grid[1]:string,
                      opt[1]:string)
;
; Arguments:
;    var          = name of the variable
;    field        = a numeric array of rank 2 (lon, lat) or
;                   rank 3 (time, lon, lat), higher ranks not tested yet
;    index        = index of current dataset in data_pointer
;    data_pointer = see interface_scripts/data_handling.ncl
;    dst_grid     = path to file containing destination grid information
;    opt          = unused
;
; Return value:
;    The field regridded to the destination grid specified in dst_grid
;
; Description:
;    Regrids the field to the destination field using the
;      'cdo remapdis' distance-weighted interpolation method
;
; Caveats
;    Climate Data Operators (CDO) package must be installed on the system
;
; Modification history
;    20170120-senftleben_daniel: written.
;
local regfile_in_tmp, regfile_out_tmp, lat, lon, gAtt, field_out, quote
begin

  ; Necessary files
  if (.not.isfilepresent(work_dir + "/CDO/")) then
    system("mkdir -p " + work_dir + "/CDO/")
  end if
  regfile_in_tmp  = work_dir + "/CDO/regrid_in.nc"   ; temporary -
  regfile_out_tmp = work_dir + "/CDO/regrid_out.nc"  ; regrid-files
  if (isfilepresent(regfile_in_tmp)) then
    system("rm " + regfile_in_tmp)
  end if
  if (isfilepresent(regfile_out_tmp)) then
    system("rm " + regfile_out_tmp)
  end if

  ; Attach coordinate information to array and file
  lat          = extract_area(index, data_pointer, "lat", var)
  lon          = extract_area(index, data_pointer, "lon", var)
  field@lat    = lat
  field@lon    = lon
  gAtt         = True
  gAtt@purpose = "Temporary file used for regridding " + \
    "with cdo, generated by cdo_remapdis (latlon.ncl)"
;    write_nc(regfile_in_tmp, var, field, gAtt)  ; FIX-ME
  field_out      = addfile(regfile_in_tmp, "w")
  field_out->lat = lat
  field_out->lon = lon

  ; Call 'CDO remapdis' with system call
  quote = inttochar(34)
  system("ncatted -a coordinates," + var + ",c,c," + quote + \
         "lon lat" + quote + " " + regfile_in_tmp)
  system("cdo -s remapdis," + dst_grid + " -selname," + \
         var + " " + regfile_in_tmp + " " + regfile_out_tmp)

  ; Read regridded field and return it
  var_reg = ncdf_read(regfile_out_tmp, var)
  return(var_reg)

end
