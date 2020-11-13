; #############################################################################
; FUNCTIONS FOR THE reformat_obs_*.ncl SCRIPTS
; #############################################################################
; General purpose functions called by the reformat_obs_*.ncl scripts.
;
; Contents
;    function create_timec
;    procedure format_time
;    procedure format_plev
;    procedure format_lev
;    procedure format_alt40
;    procedure format_lat
;    procedure format_lon
;    procedure format_coords
;    function read_cmor
;    function format_variable
;    function guess_bounds_time
;    function guess_bounds_lev
;    function guess_bounds_alt40
;    function guess_bounds_lat
;    function guess_bounds_lon
;    function guess_coord_bounds
;    function set_global_atts
;    procedure write_nc
;    procedure write_nc_profile
;    function set_size_array
;    function process_EBAS_data
;
; #############################################################################

; Time units
TUNITS = "days since 1950-01-01 00:00:00"

; CMOR FillValue
FILL = 1.e+20

; #############################################################################
undef("create_timec")
function create_timec(y1:integer,
                      y2:integer)
;
; Arguments
;    y1: start year of the time range.
;    y2: end year of the time range.
;
; Return value
;    A one-dimensional array of size 12*(y2-y1+1).
;
; Description
;    Create a monthly time coordinate for the given time range.
;
; Modification history
;    20140124-righi_mattia: written.
;
local funcname, scriptname, yy, mm, out
begin

  funcname = "create_timec"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  out = new(12 * (y2 - y1 + 1), double)
  do yy = y1, y2
    do mm = 1, 12
      out(12 * (yy - y1) + mm - 1) = \
        cd_inv_calendar(yy, mm, 15, 0, 0, 0, TUNITS, 0)
    end do
  end do

  return(out)

end

; #############################################################################
undef("format_time")
procedure format_time(var:numeric,
                      start_date:string,
                      end_date:string,
                      frequency:string)
;
; Arguments
;    var: input variable
;    start_date: start date as YYYYMMDD
;    end_date: end date as YYYYMMDD
;    frequency: time frequency ("3hr", "6hr", "day", "mon", "yr")
;
; Description
;    Check the time range of the time coordinate, set the values depending on
;    the given frequency and set the standard CMOR attributes.
;
; Caveats
;
; References
;
; Modification history
;    20200518-righi_mattia: improve timesteps consistency check.
;    20190216-righi_mattia: written.
;
local funcname, scriptname, ctime, ntime, year1, month1, day1, year2, month2, \
  day2, calendar, date, exp_ntime, yy, mm, m1, m2, dd, opt, newtime
begin

  funcname = "format_time"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Check supported frequency
  if (all(frequency.ne.(/"3hr", "6hr", "day", "mon", "yr"/))) then
    error_msg("f", scriptname, funcname, "unsupported frequency " + frequency)
  end if

  ; Read coordinate
  ctime = var&time
  ntime = dimsizes(ctime)

  ; Check monotonicity
  if (ntime.ne.1 .and. isMonotonic(ctime).ne.1) then
    error_msg("f", scriptname, funcname, \
              "non-monotonically-increasing time coordinate")
  end if

  ; Extract expected start/end date
  year1  = toint(str_get_cols(start_date, 0, 3))
  month1 = toint(str_get_cols(start_date, 4, 5))
  day1   = toint(str_get_cols(start_date, 6, 7))
  year2  = toint(str_get_cols(end_date, 0, 3))
  month2 = toint(str_get_cols(end_date, 4, 5))
  day2   = toint(str_get_cols(end_date, 6, 7))

  ; Set array of zeros with the same dimensionality for cd_inv_calendar
  zero = year1
  zero = 0

  ; Set calendar
  if (isatt(ctime, "calendar")) then
    calendar = ctime@calendar

    ; Special treatment for proleptic_gregorian calendars
    ; (not supported by 'days_in_month(...)'
    if (calendar .eq. "proleptic_gregorian") then
      if (year1 .lt. 1582) then
        error_msg("w", scriptname, funcname, \
                  "changing calendar from proleptic_gregorian " + \
                  "to gregorian prior to 1582")
      end if
      calendar = "gregorian"
    end if

  else
    calendar = "standard"
  end if

  ; Actual date, to be compared with the expected date
  date = cd_calendar(ctime, 0)
  delete(ctime)

  ; Yearly frequency
  if (frequency.eq."yr") then

    ; Check size
    exp_ntime = year2 - year1 + 1
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + \
                exp_ntime + " expected, " + ntime + " found")
    end if

    ; Reset date (1st of July)
    do yy = 0, ntime - 1
      date(yy, 0) = year1 + yy
    end do
    date(:, 1) = 7
    date(:, 2) = 1
    date(:, 3:5) = 0

  end if

  ; Monthly frequency
  if (frequency.eq."mon") then

    ; Check size
    exp_ntime = 12 * (year2 - year1 + 1)
    exp_ntime = exp_ntime - (month1 - 1)
    exp_ntime = exp_ntime - (12 - month2)
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + \
                exp_ntime + " expected, " + ntime + " found")
    end if

    ; Reset date (middle of the month)
    tt = 0
    do yy = year1, year2
      date(tt, 0) = yy
      yy@calendar = calendar
      m1 = where(yy.eq.year1, month1, 1)
      m2 = where(yy.eq.year2, month2, 12)
      do mm = m1, m2
        date(tt, 1) = mm
        dm = days_in_month(yy, mm) / 2. + 1
        date(tt, 2) = toint(dm)
        date(tt, 3) = toint(24 * (dm - toint(dm)))
        tt = tt + 1
        delete(dm)
      end do
    end do
    date(:, 4:5) = 0

  end if

  ; Daily frequency
  if (frequency.eq."day") then

    ; Check size
    opt = 0
    opt@calendar = calendar
    exp_ntime = \
      cd_inv_calendar(year2, month2, day2, zero, zero, zero, TUNITS, opt) - \
      cd_inv_calendar(year1, month1, day1, zero, zero, zero, TUNITS, opt) + 1
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + \
                exp_ntime + " expected, " + ntime + " found")
    end if
    delete(opt)

    ; Reset date (middle of the day)
    tt = 0
    do yy = year1, year2
      date(tt, 0) = yy
      m1 = where(yy.eq.year1, month1, 1)
      m2 = where(yy.eq.year2, month2, 12)
      do mm = m1, m2
        date(tt, 1) = mm
        d1 = where(yy.eq.year1 .and. mm.eq.m1, day1, 1)
        d2 = where(yy.eq.year2 .and. mm.eq.m2, day2, days_in_month(yy, mm))
        do dd = d1, d2
          date(tt, 2) = dd
          tt = tt + 1
        end do
      end do
    end do
    date(:, 3) = 12
    date(:, 4:5) = 0

  end if

  ; 6-hourly frequency (check size only)
  if (frequency.eq."6hr") then

    ; Check size
    opt = 0
    opt@calendar = calendar
    exp_ntime = \
      cd_inv_calendar(year2, month2, day2, zero, zero, zero, TUNITS, opt) - \
      cd_inv_calendar(year1, month1, day1, zero, zero, zero, TUNITS, opt) + 1
    exp_ntime = 4 * exp_ntime
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + \
                exp_ntime + " expected, " + ntime + " found")
    end if
    delete(opt)

  end if

  ; 3-hourly frequency (check size only)
  if (frequency.eq."3hr") then

    ; Check size
    opt = 0
    opt@calendar = calendar
    exp_ntime = \
      cd_inv_calendar(year2, month2, day2, zero, zero, zero, TUNITS, opt) - \
      cd_inv_calendar(year1, month1, day1, zero, zero, zero, TUNITS, opt) + 1
    exp_ntime = 8 * exp_ntime
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + \
                exp_ntime + " expected, " + ntime + " found")
    end if
    delete(opt)

  end if

  ; Set calendar
  opt = 0
  opt@calendar = calendar

  ; Redefine time coordinate
  coord = cd_inv_calendar(toint(date(:, 0)), toint(date(:, 1)), \
                          toint(date(:, 2)), toint(date(:, 3)), \
                          toint(date(:, 4)), toint(date(:, 5)), \
                          TUNITS, opt)

  ; Set standard attributes
  newtime = todouble(coord)  ; this also removes attributes
  copy_VarCoords(coord, newtime)
  newtime@bounds = "time_bnds"
  newtime@calendar = calendar
  newtime@long_name = "time"
  newtime@axis = "T"
  newtime@units = TUNITS
  newtime@standard_name = "time"
  if (isatt(newtime, "_FillValue")) then
    delete(newtime@_FillValue)
  end if

  ; Reset time coordinate
  delete(var&time)
  var&time = newtime

end

; #############################################################################
undef("format_plev")
procedure format_plev(var:numeric)
;
; Arguments
;    var: input variable
;
; Description
;    Check the monotonicity of the plev (pressure) coordinate and set the
;     standard CMOR attributes.
;
; Caveats
;
; References
;
; Modification history
;    20190216-righi_mattia: written.
;
local funcname, scriptname, rank, cplev, newplev
begin

  funcname = "format_plev"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Set rank
  rank = dimsizes(dimsizes(var))

  ; Check monotonicity
  if (isMonotonic(var&plev) .eq. 0) then
    error_msg("f", scriptname, funcname, "non-monotonic vertical coordinate")
  end if
  if (isMonotonic(var&plev).eq.1) then  ; must be monotonically decreasing
    if (rank.eq.4) then
      var = var(:, ::-1, :, :)
    elseif (rank.eq.3) then
      var = var(:, ::-1, :)
    elseif (rank.eq.2) then
      var = var(:, ::-1)
    end if
  end if

  ; Read coordinate
  cplev = var&plev

  ; Set standard attributes
  newplev = todouble(cplev)  ; this also removes attributes
  copy_VarCoords(cplev, newplev)
  newplev@positive = "down"
  newplev@long_name = "pressure"
  newplev@axis = "Z"
  newplev@units = "Pa"
  newplev@standard_name = "air_pressure"
  if (isatt(newplev, "_FillValue")) then
    delete(newplev@_FillValue)
  end if

  ; Reset plev coordinate
  delete(var&plev)
  var&plev = newplev

end

; #############################################################################
undef("format_lev")
procedure format_lev(var:numeric)
;
; Arguments
;    var: input variable
;
; Description
;    Check the monotonicity of the lev (depth) coordinate and set the
;    standard CMOR attributes.
;
; Caveats
;
; References
;
; Modification history
;    20190216-righi_mattia: written.
;
local funcname, scriptname, rank, clev, newlev
begin

  funcname = "format_lev"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Set rank
  rank = dimsizes(dimsizes(var))

  ; Check monotonicity
  if (isMonotonic(var&lev) .eq. 0) then
    error_msg("f", scriptname, funcname, "non-monotonic vertical coordinate")
  end if
  if (isMonotonic(var&lev).eq.-1) then  ; must be monotonically increasing
    if (rank.eq.4) then
      var = var(:, ::-1, :, :)
    elseif (rank.eq.3) then
      var = var(:, ::-1, :)
    elseif (rank.eq.2) then
      var = var(:, ::-1)
    end if
  end if

  ; Read coordinate
  clev = var&lev

  ; Set standard attributes
  newlev = todouble(clev)  ; this also removes attributes
  copy_VarCoords(clev, newlev)
  newlev@bounds = "lev_bnds"
  newlev@positive = "down"
  newlev@long_name = "ocean depth coordinate"
  newlev@axis = "Z"
  newlev@units = "m"
  newlev@standard_name = "depth"
  if (isatt(newlev, "_FillValue")) then
    delete(newlev@_FillValue)
  end if

  ; Reset lev coordinate
  delete(var&lev)
  var&lev = newlev

end

; #############################################################################
undef("format_alt40")
procedure format_alt40(var:numeric)
;
; Arguments
;    var: input variable
;
; Description
;    Check the monotonicity of the alt40 (altitude) coordinate and set the
;    standard CMOR attributes.
;
; Caveats
;
; References
;
; Modification history
;    20200204-lauer_axel: written.
;
local funcname, scriptname, rank, clev, newlev
begin

  funcname = "format_alt40"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Set rank
  rank = dimsizes(dimsizes(var))

  ; Check monotonicity
  if (isMonotonic(var&alt40) .eq. 0) then
    error_msg("f", scriptname, funcname, "non-monotonic vertical coordinate")
  end if
  if (isMonotonic(var&alt40).eq.-1) then  ; must be monotonically increasing
    if (rank.eq.4) then
      var = var(:, ::-1, :, :)
    elseif (rank.eq.3) then
      var = var(:, ::-1, :)
    elseif (rank.eq.2) then
      var = var(:, ::-1)
    end if
  end if

  ; Read coordinate
  calt = var&alt40

  ; Set standard attributes
  newalt = todouble(calt)  ; this also removes attributes
  copy_VarCoords(calt, newalt)
  newalt@bounds = "alt40_bnds"
  newalt@positive = "up"
  newalt@long_name = "altitude"
  newalt@axis = "Z"
  newalt@units = "m"
  newalt@standard_name = "altitude"
  if (isatt(newalt, "_FillValue")) then
    delete(newalt@_FillValue)
  end if

  ; Reset lev coordinate
  delete(var&alt40)
  var&alt40 = newalt

end

; #############################################################################
undef("format_lat")
procedure format_lat(var:numeric)
;
; Arguments
;    var: input variable
;
; Description
;    Check the monotonicity of the latitude coordinat (S->N) and set the
;    standard CMOR attributes.
;
; Caveats
;
; References
;
; Modification history
;    20190216-righi_mattia: written.
;
local funcname, scriptname, rank, dims, dpos, lcheck, clat, newlat
begin

  funcname = "format_lat"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Set rank
  rank = dimsizes(dimsizes(var))
  dims = getvardims(var)
  dpos = ind(dims.eq."lat")

  ; Check monotonicity but not for point variables
  if (isMonotonic(var&lat) .eq. 0 .and. \
      dimsizes(var&lat) .gt.  1) then
    error_msg("f", scriptname, funcname, "non-monotonic latitude coordinate")
  end if
  if (isMonotonic(var&lat) .eq. -1) then  ; must be S->N
    lcheck = False
    if (rank.eq.4) then
      if (dpos.eq.2) then
        var = var(:, :, ::-1, :)
        lcheck = True
      end if
    end if
    if (rank.eq.3) then
      if (dpos.eq.1) then
        var = var(:, ::-1, :)
        lcheck = True
      end if
      if (dpos.eq.2) then
        var = var(:, :, ::-1)
        lcheck = True
      end if
    end if
    if (rank.eq.2) then
      if (dpos.eq.0) then
        var = var(::-1, :)
        lcheck = True
      end if
    end if
    if (.not.lcheck) then
      error_msg("f", scriptname, funcname, "cannot locate latitude position")
    end if
  end if

  ; Read coordinate
  clat = var&lat

  ; Set standard attributes
  newlat = todouble(clat)  ; this also removes attributes
  copy_VarCoords(clat, newlat)
  newlat@bounds = "lat_bnds"
  newlat@long_name = "latitude"
  newlat@axis = "Y"
  newlat@units = "degrees_north"
  newlat@standard_name = "latitude"
  if (isatt(newlat, "_FillValue")) then
    delete(newlat@_FillValue)
  end if

  ; Reset lat coordinate
  delete(var&lat)
  var&lat = newlat

end

; #############################################################################
undef("format_lon")
procedure format_lon(var:numeric)
;
; Arguments
;    var: input variable
;
; Description
;    Check the lon coordinate (0:360) and set the standard CMOR attributes.
;
; Caveats
;
; References
;
; Modification history
;    20190216-righi_mattia: written.
;
local funcname, scriptname, clon, newlon
begin

  funcname = "format_lon"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Check monotonicity but not for point variables
  if (isMonotonic(var&lon) .eq. 0 .and. \
      dimsizes(var&lon) .gt.  1) then
    error_msg("f", scriptname, funcname, "non-monotonic longitude coordinate")
  end if

  ; Check that lon is 0:360
  if (any(var&lon.lt.0.)) then
    var = lonFlip(var)
  end if

  ; Read coordinate
  clon = var&lon

  ; Set standard attributes
  newlon = todouble(clon)  ; this also removes attributes
  copy_VarCoords(clon, newlon)
  newlon@bounds = "lon_bnds"
  newlon@long_name = "longitude"
  newlon@axis = "X"
  newlon@units = "degrees_east"
  newlon@standard_name = "longitude"
  if (isatt(newlon, "_FillValue")) then
    delete(newlon@_FillValue)
  end if

  ; Reset lon coordinate
  delete(var&lon)
  var&lon = newlon

end

; #############################################################################
undef("format_coords")
procedure format_coords(var:numeric,
                        date1:string,
                        date2:string,
                        frequency:string)
;
; Arguments
;    var: input variable with named dimensions and assiociated coordinates.
;    date1: start date as YYYYMMDD
;    date2: end date as YYYYMMDD
;    frequency: time frequency ("3hr", "6hr", "day", "mon", "yr")
;
; Description
;    Format the coordinate according to the CF/CMOR standard.
;
; Caveats
;
; References
;
; Modification history
;    20190216-righi_mattia: written.
;
begin

  funcname = "format_coords"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Get variable dimensions
  dnames = getvardims(var)

  ; Loop over dimensions, call formatting procedures for each coordinate
  do dd = 0, dimsizes(dnames) - 1

    found = False

    if (dnames(dd).eq."time") then
      format_time(var, date1, date2, frequency)
      found = True
    end if

    if (dnames(dd).eq."plev") then
      format_plev(var)
      found = True
    end if

    if (dnames(dd).eq."lev") then
      format_lev(var)
      found = True
    end if

    if (dnames(dd).eq."alt40") then
      format_alt40(var)
      found = True
    end if

    if (dnames(dd).eq."lat") then
      format_lat(var)
      found = True
    end if

    if (dnames(dd).eq."lon") then
      format_lon(var)
      found = True
    end if

    if (.not.found) then
      error_msg("f", scriptname, funcname, "cannot format coordinate " + \
                dnames(dd))
    end if

  end do

end

; #############################################################################
undef("read_cmor")
function read_cmor(name:string,
                   table:string)
;
; Arguments
;    name: standard variable name.
;    string: full path to the CMOR table of the variable.
;
; Return value
;    A logical variable with the CMOR table attached as attributes.
;
; Description
;    Read variable attributes from the CMOR tables (cmor/<name>.cmor).
;
; Caveats
;
; References
;
; Modification history
;    20190107-righi_mattia: modify to read standard CMIP5 tables
;    20130528-righi_mattia: written.
;
local funcname, scriptname, data, idxu, idxd, attn, attv, out
begin

  funcname = "read_cmor"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Read attributes from cmor table
  if (.not.fileexists(table)) then
    error_msg("f", scriptname, funcname, \
              "cannot find CMOR table " + table)
  end if

  if (isStrSubset(table, ".json")) then  ; CMIP6 tables

    error_msg("f", scriptname, funcname, \
              "use of CMIP6 CMOR tables not supported")

  else  ; CMIP5 and custom tables

    data = readAsciiTable(table, 1, "string", 0)

    ; Extract variable block
    idxu = ind(data(:, 0).eq."variable_entry:    " + name)
    if (any(ismissing(idxu))) then
      error_msg("f", scriptname, funcname, \
                "table for variable " + name + " not found in table " + \
                table)
    end if
    tmp = ind(str_get_field(data(:, 0), 1, ":").eq."variable_entry")
    if (dimsizes(tmp).gt.1) then
      next = min(ind(tmp.gt.idxu))
      if (.not.ismissing(next))
        idxd = tmp(min(ind(tmp.gt.idxu))) - 2
      else
        idxd = dimsizes(data(:, 0)) - 1
      end if
      data := data(idxu:idxd, 0)
      delete(idxd)
    else
      data := data(:, 0)  ; just 1 variable in this table
    end if
    delete(idxu)
    delete(tmp)

    ; Extract attributes
    idxu = ind(str_get_field(data, 1, ":").eq."! Variable attributes") + 2
    idxd = ind(str_get_field(data, 1, ":").eq. \
               "! Additional variable information") - 2
    attn = str_squeeze(str_get_field(data(idxu:idxd), 1, ":"))
    attv = str_squeeze(str_get_field(data(idxu:idxd), 2, ":"))

  end if

  out = True
  do ii = 0, dimsizes(attn) - 1
    out@$attn(ii)$ = attv(ii)
  end do

  return(out)

end

; #############################################################################
undef("format_variable")
function format_variable(var:numeric,
                         name:string,
                         table:string)
;
; Arguments
;    var: input variable.
;    name: standard name of the input variable.
;    string: full path to the CMOR table containing the variable.
;
; Return value
;    An array of the same dimensionality of var.
;
; Description
;    Set standard variable attributes according to the given CMOR table.
;
; Caveats
;
; References
;
; Modification history
;    20190107-righi_mattia: add extra argument for CMOR table
;    20161202-lauer_axel: preserve attribute "coordinates" if present
;    20130528-righi_mattia: written.
;
local funcname, scriptname, coordattr, out, tmp, att, ii
begin

  funcname = "var_attrib"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Set fill value first
  if(isatt(var, "_FillValue")) then
    var = where(var.eq.var@_FillValue, FILL, var)
  end if
  var@_FillValue = FILL

  if (isatt(var, "coordinates")) then
    coordattr = var@coordinates
  end if

  ; Remove attributes
  delete_VarAtts(var, -1)

  ; Convert to float
  if (typeof(var).ne."float") then
    out = tofloat(var)
    copy_VarCoords(var, out)
  else
    out = var
  end if

  ; Append attributes
  out@_FillValue = FILL
  tmp = read_cmor(name, table)
  att = getvaratts(tmp)
  do ii = 0, dimsizes(att) - 1
    out@$att(dimsizes(att) - 1 - ii)$ = tmp@$att(dimsizes(att) - 1 - ii)$
  end do

  if (isvar("coordattr")) then
    out@coordinates = coordattr
  end if

  return(out)

end

; #############################################################################
undef("guess_bounds_time")
function guess_bounds_time(coord[*]:double,
                           frequency:string)
;
; Arguments
;    coord: input time coordinate.
;    frequency: time frequency ("3hr", "6hr", "day", "mon", "yr")
;
; Return value
;    A two dimensional array, with the first dimension of the same size of the
;    input coordinate and the second dimension of size 2.
;
; Description
;    Calculate the boundaries of the time coordinate given the frequency.
;
; Caveats
;
; References
;
; Modification history
;    20190217-righi_mattia: written.
;
local funcname, scriptname, date, year, month, day, opt, units, tyear, tmonth
begin

  funcname = "guess_bounds_time"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; This function assumes that units of days are used
  if (.not.isStrSubset(coord@units, "days since")) then
    error_msg("f", scriptname, funcname, "unsupported units " + \
              coord@units + ", expected " + TUNITS)
  end if

  ; Define boundaries
  bounds = new((/dimsizes(coord), 2/), double)
  bounds!0 = "time"
  bounds&time = coord
  bounds!1 = "bnds"
  delete(bounds@_FillValue)

  ; Get date
  date = cd_calendar(coord, 0)
  year = date(:, 0)
  month = date(:, 1)
  day = date(:, 2)

  ; Set array of constants with the same dimensionality for cd_inv_calendar
  zero = year
  zero = 0
  one = year
  one = 1

  ; Settings for calendar
  units = coord@units
  opt = 0
  copy_VarAtts(coord, opt)

  if (frequency.eq."yr") then
    ; 1st day of the year
    bounds(:, 0) = \
      (/cd_inv_calendar(year, one, one, zero, zero, zero, units, opt)/)
    ; 1st day of the next year
    bounds(:, 1) = \
      (/cd_inv_calendar(year + 1, one, one, zero, zero, zero, units, opt)/)
    return(bounds)
  end if

  if (frequency.eq."mon") then
    ; 1st day of the month
    bounds(:, 0) = \
      (/cd_inv_calendar(year, month, one, zero, zero, zero, units, opt)/)
    ; 1st day of the next month (special case December!)
    tmonth = where(month.eq.12, 1, month + 1)
    tyear = where(month.eq.12, year + 1, year)
    bounds(:, 1) = \
      (/cd_inv_calendar(tyear, tmonth, one, zero, zero, zero, units, opt)/)
    return(bounds)
    delete([/tyear, tmonth/])
  end if

  if (frequency.eq."day") then
    ; Shift of half a day
    bounds(:, 0) = coord - 0.5
    bounds(:, 1) = coord + 0.5
    return(bounds)
  end if

  if (frequency.eq."6hr") then
    ; Shift of 3 hours (day/8)
    bounds(:, 0) = coord - 0.125
    bounds(:, 1) = coord + 0.125
    return(bounds)
  end if

  if (frequency.eq."3hr") then
    ; Shift of 1.5 hours (day/16)
    bounds(:, 0) = coord - 0.0625
    bounds(:, 1) = coord + 0.0625
    return(bounds)
  end if

  error_msg("f", scriptname, funcname, "unsupported frequency " + frequency)

end

; #############################################################################
undef("guess_bounds_lev")
function guess_bounds_lev(coord[*]:double)
;
; Arguments
;    coord: input level coordinate.
;
; Return value
;    A two dimensional array, with the first dimension of the same size of the
;    input coordinate and the second dimension of size 2.
;
; Description
;    Calculate the boundaries of the level coordinate as midpoints between
;    the input levels. The first (top) boundary is set to a maximum of zero.
;
; Caveats
;
; References
;
; Modification history
;    20190217-righi_mattia: written.
;
local funcname, scriptname, size, top
begin

  funcname = "guess_bounds_lev"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  bounds = new((/dimsizes(coord), 2/), double)
  bounds!0 = "lev"
  bounds&lev = coord
  bounds!1 = "bnds"
  delete(bounds@_FillValue)

  size = dimsizes(coord)

  bounds(1:size - 1, 0) = 0.5 * (coord(0:size - 2) + coord(1:size - 1))
  bounds(0:size - 2, 1) = bounds(1:size - 1, 0)

  ; Set top and bottom separately
  top = coord(0) - 0.5 * (coord(1) - coord(0))
  bounds(0, 0) = where(top.ge.0., top, 0.)
  bounds(size - 1, 1) = \
    coord(size - 1) + 0.5 * (coord(size - 1) - coord(size - 2))

  return(bounds)

end

; #############################################################################
undef("guess_bounds_alt40")
function guess_bounds_alt40(coord[*]:double)
;
; Arguments
;    coord: input level coordinate.
;
; Return value
;    A two dimensional array, with the first dimension of the same size of the
;    input coordinate and the second dimension of size 2.
;
; Description
;    Calculate the boundaries of the altitude40 coordinate as midpoints between
;    the input levels. The first (top) boundary is set to a maximum of zero.
;
; Caveats
;
; References
;
; Modification history
;    20200204-lauer_axel: written.
;
local funcname, scriptname, size, top
begin

  funcname = "guess_bounds_alt40"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  bounds = new((/dimsizes(coord), 2/), double)
  bounds!0 = "alt40"
  bounds&lev = coord
  bounds!1 = "bnds"
  delete(bounds@_FillValue)

  size = dimsizes(coord)

  bounds(1:size - 1, 0) = 0.5 * (coord(0:size - 2) + coord(1:size - 1))
  bounds(0:size - 2, 1) = bounds(1:size - 1, 0)

  ; Set top and bottom separately
  top = coord(0) - 0.5 * (coord(1) - coord(0))
  bounds(0, 0) = where(top.ge.0., top, 0.)
  bounds(size - 1, 1) = \
    coord(size - 1) + 0.5 * (coord(size - 1) - coord(size - 2))

  return(bounds)

end

; #############################################################################
undef("guess_bounds_lat")
function guess_bounds_lat(coord[*]:double)
;
; Arguments
;    coord: input latitude coordinate.
;
; Return value
;    A two dimensional array, with the first dimension of the same size of the
;    input coordinate and the second dimension of size 2.
;
; Description
;    Calculate the boundaries of the latitude coordinate as midpoints between
;    the input values. The first and last boundary is set to 90 S and 90 N,
;    respectively.
;
; Caveats
;    This function works only with regular grids.
;
; References
;
; Modification history
;    20190217-righi_mattia: written.
;
local funcname, scriptname
begin

  funcname = "guess_bounds_lat"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  bounds = new((/dimsizes(coord), 2/), double)
  bounds!0 = "lat"
  bounds&lat = coord
  bounds!1 = "bnds"
  delete(bounds@_FillValue)

  size = dimsizes(coord)

  bounds(1:size - 1, 0) = 0.5 * (coord(0:size - 2) + coord(1:size - 1))
  bounds(0:size - 2, 1) = bounds(1:size - 1, 0)
  bounds(0, 0) = -90.
  bounds(size - 1, 1) = 90.

  return(bounds)

end

; #############################################################################
undef("guess_bounds_lon")
function guess_bounds_lon(coord[*]:double)
;
; Arguments
;    coord: input longitude coordinate.
;
; Return value
;    A two dimensional array, with the first dimension of the same size of the
;    input coordinate and the second dimension of size 2.
;
; Description
;    Calculate the boundaries of the longitude coordinate as midpoints between
;    the input values.
;
; Caveats
;    This function works only with regular grids.
;
; References
;
; Modification history
;    20190217-righi_mattia: written.
;
local funcname, scriptname
begin

  funcname = "guess_bounds_lon"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  bounds = new((/dimsizes(coord), 2/), double)
  bounds!0 = "lon"
  bounds&lon = coord
  bounds!1 = "bnds"
  delete(bounds@_FillValue)

  size = dimsizes(coord)

  bounds(1:size - 1, 0) = 0.5 * (coord(0:size - 2) + coord(1:size - 1))
  bounds(0:size - 2, 1) = bounds(1:size - 1, 0)
  bounds(0, 0) = coord(0) - 0.5 * (coord(1) - coord(0))
  bounds(size - 1, 1) = \
    coord(size - 1) + 0.5 * (coord(size - 1) - coord(size - 2))

  return(bounds)

end

; #############################################################################
undef("guess_coord_bounds")
function guess_coord_bounds(var:numeric,
                            frequency:string)
local funcname, scriptname, time_bnds, plev_bnds, lev_bnds, lat_bnds, lon_bnds
begin

  funcname = "guess_coord_bounds"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Check supported frequency
  if (all(frequency.ne.(/"3hr", "6hr", "day", "mon", "yr", "fx"/))) then
    error_msg("f", scriptname, funcname, "unsupported frequency " + frequency)
  end if

  ; Get variable dimensions
  dnames = getvardims(var)

  ; Initialize list
  bounds_list = NewList("fifo")

  ; Loop over dimensions, call formatting procedures for each coordinate
  do dd = 0, dimsizes(dnames) - 1

    if (dnames(dd).eq."time") then
      time_bnds = guess_bounds_time(var&time, frequency)
      ListPush(bounds_list, time_bnds)
      continue
    end if

    if (dnames(dd).eq."plev") then
      continue  ; plev bounds are not required
    end if

    if (dnames(dd).eq."lev") then
      lev_bnds = guess_bounds_lev(var&lev)
      ListPush(bounds_list, lev_bnds)
      continue
    end if

    if (dnames(dd).eq."alt40") then
      alt40_bnds = guess_bounds_lev(var&alt40)
      ListPush(bounds_list, alt40_bnds)
      continue
    end if

    if (dnames(dd).eq."lat") then
      lat_bnds = guess_bounds_lat(var&lat)
      ListPush(bounds_list, lat_bnds)
      continue
    end if

    if (dnames(dd).eq."lon") then
      lon_bnds = guess_bounds_lon(var&lon)
      ListPush(bounds_list, lon_bnds)
      continue
    end if

    if (.not.found) then
      error_msg("f", scriptname, funcname, "cannot guess bounds for " + \
                "coordinate " + dnames(dd))
    end if

  end do

  return(bounds_list)

end

; #############################################################################
undef("set_global_atts")
function set_global_atts(obsname:string,
                         tier:integer,
                         source:string,
                         reference:string,
                         comment:string)
;
; Argument
;    obsname: name of the observational dataset.
;    reference: reference for the dataset, or leave empty if not available.
;    source: link to the data source.
;    tier: tier number (2 or 3).
;    comment: additional information if required, or leave empty.
;
; Return value
;    A logical containing the arguments as attributes.
;
; Description
;    Generate the global attribute for the output file by combining user
;    provided information with default ones (author, host, date, etc.).
;
; Modification history
;    20190202-righi_mattia: written.
;
local funcname, scriptname,
    dim_unlim, ii
begin

  funcname = "set_global_atts"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  global = True
  global@title = obsname + " data reformatted for the ESMValTool v2.0"
  global@tier = tier
  global@source = source
  if (strlen(str_squeeze(reference)).ne.0) then
    global@reference = reference
  else
    global@reference = "not available"
  end if
  if (strlen(str_squeeze(comment)).ne.0) then
    global@comment = comment
  end if
  global@user = systemfunc("echo $USER")
  global@host = systemfunc("hostname -f")
  global@history = "Created on " + systemfunc("date")
  global@conventions = "CF/CMOR"

  return(global)

end

; #############################################################################
undef("write_nc")
procedure write_nc(outfile:string,
                   name:string,
                   var:numeric,
                   bounds:list,
                   gAtt:logical)
;
; Arguments
;    outfile: the name of the file to be written, including its path.
;    name: the variable standard name.
;    var: the variable array.
;    bounds: a list containing the bounds of the variable coordinates.
;    gAtt: a logical variable, whose attributes are appended as file
;          attributes.
;
; Description
;    Write the given variable to the given NetCDF file, together with its
;    coordinates and boundaries, and append the provided global attributes.
;
; Modification history
;    20190218_righi_mattia: extend with coordinate bounds.
;    20140123-righi_mattia: written.
;
local funcname, scriptname, w, gAtt, dim_names, ndims, dim_sizes, dim_types, \
    dim_unlim, ii
begin

  funcname = "write_nc"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Open file
  if (fileexists(outfile)) then
    system("rm -f " + outfile)
  end if
  w = addfile(outfile, "c")
  setfileoption(w, "DefineMode", True)

  ; Get coordinates
  dim_names = getvardims(var)
  ndims = dimsizes(dim_names)
  dim_sizes = new(ndims, integer)
  dim_types = new(ndims, string)
  dim_unlim = new(ndims, logical)
  do ii = 0, ndims - 1
    dim_sizes(ii) = dimsizes(var&$dim_names(ii)$)
    dim_types(ii) = typeof(var&$dim_names(ii)$)
    dim_unlim(ii) = False
  end do

  ; Time coordinate must be unlimited
  if (any(dim_names.eq."time")) then
    dim_sizes(ind(dim_names.eq."time")) = -1
    dim_unlim(ind(dim_names.eq."time")) = True
  end if

  ; Define dimensions
  filedimdef(w, dim_names, dim_sizes, dim_unlim)
  filedimdef(w, "bnds", 2, False)
  do ii = 0, ndims - 1
    ; Coordinates
    filevardef(w, dim_names(ii), dim_types(ii), dim_names(ii))
    filevarattdef(w, dim_names(ii), var&$dim_names(ii)$)
    ; Bounds
    bnds_name = dim_names(ii) + "_bnds"
    listidx = ListIndexFromName(bounds, bnds_name)
    if (listidx.ne.-1) then
      filevardef(w, bnds_name, dim_types(ii), (/dim_names(ii), "bnds"/))
    end if
    ; No attributes for the bounds
  end do

  ; Define variable
  filevardef(w, name, "float", dim_names)
  filevarattdef(w, name, var)

  ; Append global attributes
  fileattdef(w, gAtt)

  ; Write coordinate and corresponding bounds (if available)
  setfileoption(w, "DefineMode", False)
  do ii = 0, ndims - 1
    w->$dim_names(ii)$ = (/var&$dim_names(ii)$/)
    bnds_name = dim_names(ii) + "_bnds"
    listidx = ListIndexFromName(bounds, bnds_name)
    if (listidx.ne.-1) then
      bound = bounds[listidx]
      w->$bnds_name$ = (/bound/)
      delete(bound)
    end if
  end do

  ; Write variable
  w->$name$ = (/var/)

end

; #############################################################################
undef("write_nc_profile")
procedure write_nc_profile(outfile:string,
                           name:string,
                           var:numeric,
                           gAtt:logical)
;
; Arguments
;    outfile: the name of the file to be written, including its path.
;    name: the variable name.
;    var: the variable field.
;    gAtt: a logical variable, whose attributes are appended as file
;          attributes.
;
; Description
;    Write the given variable to the given NetCDF file, appending also the
;    provided global attributes.
;    Designed to write multiple variables for the vertical profiles data.
;
; Modification history
;    20140422-righi_mattia: written.
;
local funcname, scriptname, w, coords, cc, jj, locname, locvar, cname
begin

  funcname = "write_nc_profile"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Open file
  if (fileexists(outfile)) then
    system("rm -f " + outfile)
  end if
  w = addfile(outfile, "c")
  setfileoption(w, "DefineMode", True)

  ; Attach global attributes
  fileattdef(w, gAtt)

  ; Write dimensions
  coords = getvardims(var)
  do cc = 0, dimsizes(coords) - 2  ; skip column
    cname = coords(cc)
    filedimdef(w, cname, dimsizes(var&$cname$), False)
    filevardef(w, cname, typeof(var&$cname$), cname)
  end do

  ; Write variable
  do jj = 0, dimsizes(var&column) - 1

    ; Extract given column
    locname = str_sub_str(name + "_" + var&column(jj), "%", "")
    if (isdim(var, "case")) then
      locvar = var(:, :, jj)
    else
      locvar = var(:, jj)
    end if
    if (var&column(jj).eq."N") then
      locvar@units = "1"
    end if

    ; Define dimensions
    filevardef(w, locname, "float", coords(0: dimsizes(coords) - 2))
    do cc = 0, dimsizes(coords) - 2
      cname = coords(cc)
      filevarattdef(w, cname, locvar&$cname$)
    end do
    filevarattdef(w, locname, locvar)

    ; Write
    setfileoption(w, "DefineMode", False)
    do cc = 0, dimsizes(coords) - 2
      cname = coords(cc)
      w->$cname$ = (/locvar&$cname$/)
    end do
    w->$locname$ = (/locvar/)
    delete(locvar)
    delete(locname)

  end do

end

; #############################################################################
undef("set_size_array")
function set_size_array()
;
; Arguments
;
; Return value
;    An array of type double.
;
; Description
;    Set a logarithmic array of sizes to be used for particle size
;    distribution calculations.
;
; Caveats
;
; References
;
; Modification history
;    20130528-righi_mattia: written.
;
local funcname, scriptname, minsize, maxsize, nbins, bin, out
begin

  funcname = "set_size_array"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; Size range (0.5 nm - 10 um)
  minsize = 0.5e-9
  maxsize = 10.e-6
  nbins = 100

  ; Generate array
  out = new(nbins, double)
  bin = 10. ^ (log10(maxsize / minsize) / (nbins - 1))
  out(0) = minsize
  do ii = 1, nbins - 1
    out(ii) = out(ii - 1) * bin
  end do

  return(out)

end

; #############################################################################
undef("process_EBAS_data")
function process_EBAS_data(in_vars[*]:string,
                           in_units[*]: string,
                           in_matrix[*]:string,
                           in_compon[*]:string,
                           in_column[*]:string,
                           indir[1]:string,
                           st_code:string,
                           y1:integer,
                           y2:integer)
;
; Arguments
;    in_vars: variables standard name.
;    in_units: variables units in the raw data.
;    in_matrix: variables matrix in the raw data.
;    in_compon: variables name in the raw data.
;    in_column: variables name in the header.
;    indir: the input directory for raw data.
;    stcode: the code of the station to be processed (used for
;            cross-checking)
;    y1: start year of the considered time interval.
;    y2: end year of the considered time interval.
;
; Return value
;    A two-dimensional array (time, variable) with the monthly mean time
;    series of for each of the processed variables.
;
; Description
;    Process the data from the EBAS database (e.g., EANET, EMEP).
;
; Caveats
;    For the time coordinate in the input data, only units of days are
;    currently accepted.
;
; Modification history
;    20150413-righi_mattia: improved time selection.
;    20140124-righi_mattia: written.
;
local timec, datec, vID, fID, bn, en, head, hh, cline, syear, smonth, sday, \
  scode, comp, matr, unit, scale, fills, lline, cols, data_col, flag_col, \
  start_col, end_col, data, value, flag, start_time, end_time, mm, sidx, \
  monthind, stday, enday, nd, pt1, pt2, data_arr
begin

  funcname = "process_EBAS_data"
  scriptname = "esmvaltool/cmorizers/obs/utilities.ncl"

  ; EBAS flags for valid measurements
  ; (see http://www.nilu.no/projects/ccc/flags/index.html)
  validflags = (/798, 797, 782, 781, 780, 771, 770, 741, 740, 680, 679, \
                678, 676, 675, 668, 665, 662, 660, 657, 656, 655, 654, \
                653, 652, 651, 650, 649, 648, 644, 640, 559, 558, 557, \
                556, 555, 532, 531, 521, 499, 498, 476, 475, 470, 458, \
                457, 450, 440, 420, 410, 394, 392, 390, 382, 380, 370, \
                299, 298, 276, 275, 258, 257, 250, 249, 248, 247, 220, \
                211, 210, 191, 190, 189, 188, 187, 186, 185, 147, 120, \
                111, 110, 103, 102, 101, 100, 000/)

  ; Create time coordinate
  timec = create_timec(y1, y2)
  datec = cd_calendar(timec, -1)

  ; Create output array
  data_arr = new((/dimsizes(timec), dimsizes(in_vars)/), float)
  data_arr!0 = "time"
  data_arr&time = timec
  data_arr@_FillValue = FILL

  ; Create a temporary arrays for time averages and weights
  temp_arr = new(dimsizes(timec), float)
  temp_arr!0 = "time"
  temp_arr&time = timec
  temp_arr@_FillValue = FILL
  ndays_arr = new(dimsizes(timec), float)
  ndays_arr!0 = "time"
  ndays_arr&time = timec
  ndays_arr@_FillValue = FILL

  ; Loop over variables
  do vID = 0, dimsizes(in_vars) - 1

    log_info("  Processing variable " + in_compon(vID))

    ; Initialize
    temp_arr = 0.
    ndays_arr = 0.

    ; Read file list
    cstr = "find " + indir + " -type f -name '" + \
      st_code + ".*." + in_compon(vID) + "." + in_matrix(vID) + "*.nas'"
    in_files = systemfunc(cstr)
    if (all(ismissing(in_files))) then
      delete(in_files)
      continue
    end if
    in_bnames = systemfunc(cstr + " -exec basename {} " + inttochar(92) + ";")
    sy = str_get_cols(str_get_field(in_bnames, 2, "."), 0, 5)
    delete(cstr)

    sqsort(in_files)
    sqsort(in_bnames)

    ; Check for duplicates
    if (dimsizes(UNIQ(sy)).ne.dimsizes(sy)) then
      log_info("Duplicated data in input files")
      do fID = 0, dimsizes(in_files) - 1
        log_info("  " + in_files(fID))
      end do
      log_info("Remove duplicated files considering the following criteria")
      log_info("  most recent revision date")
      log_info("  most complete time coverage")
      log_info("  same instrument in different years")
      error_msg("f", scriptname, funcname, \
                "rerun this station after removing duplicates")
    end if
    delete(sy)

    ; Loop over input files
    do fID = 0, dimsizes(in_files) - 1

      log_info("   Reading file " + in_bnames(fID))

      ; Read header
      head = readAsciiHead(in_files(fID), "starttime")

      ; Extract and check starting date
      syear = toint(str_get_field(head(6), 1, " "))
      smonth = toint(str_get_field(head(6), 2, " "))
      sday = toint(str_get_field(head(6), 3, " "))

      ; Get time units
      utime = str_sub_str(head(8), "file reference point", "")
      if (.not.isStrSubset(utime, "days")) then
        error_msg("f", scriptname, funcname, "unexpected time units")
      end if
      utime = utime + syear + "-" + smonth + "-" + sday
      delete(syear)
      delete(smonth)
      delete(sday)

      ; Use first file units as reference
      if (fID.eq.0) then
        ref_utime = utime
      end if

      ; Check units consistency
      do hh = 0, dimsizes(head) - 1
        if (isStrSubset(head(hh), "Unit:")) then
          unit = str_squeeze(str_get_field(head(hh), 2, ":"))
          if (unit .ne. in_units(vID) .and. unit.ne."ug/m3") then
            error_msg("f", scriptname, funcname, \
                      "units in the file not as expected " + \
                      "(" + unit + " vs. " + in_units(vID) + ")")
          end if
          delete(unit)
        end if
      end do

      ; Get specific fill values and scale factors
      scale = tofloat(str_get_field(head(10), 2, " "))
      fills = tofloat(str_get_field(head(11), 2, " "))

      ; Get column names
      lline = head(dimsizes(head) - 1)
      ncols = str_fields_count(lline, " ")
      cols = new(ncols, string)
      do cc = 0, ncols - 1
        cols(cc) = str_get_field(lline, cc + 1, " ")
      end do
      data_col = min(ind(cols.eq.in_column(vID))) + 1
      flag_col = (ind(cols.eq."flag_" + in_column(vID).or.cols.eq."flag")) + 1
      start_col = ind(cols.eq."starttime") + 1
      end_col = ind(cols.eq."endtime") + 1
      delete(cols)
      delete(ncols)

      ; Read data
      data = readAsciiTable(in_files(fID), 1, "string", dimsizes(head))
      delete(head)

      ; Read data (for the given month)
      value = tofloat(str_get_field(data(:, 0), data_col, " "))
      value@_FillValue = -999.

      ; Read flag
      flag = toint(1000 * tofloat(str_get_field(data(:, 0), flag_col, " ")))

      ; Filter for valid values
      value = where(value.eq.fills, value@_FillValue, value)
      value = where(value.lt.0, value@_FillValue, value)
      do jj = 0, dimsizes(value) - 1
        if (all(flag(jj).ne.validflags)) then
          value(jj) = value@_FillValue
        end if
      end do
      delete(flag)
      delete(fills)

      ; Apply scaling
      if (scale.ne.1) then
        value = value * scale
      end if
      delete(scale)

      ; Get start and end time
      stt = todouble(str_get_field(data(:, 0), start_col, " "))
      stt@units = utime
      ent = todouble(str_get_field(data(:, 0), end_col, " "))
      ent@units = utime
      delete(data)

      ; Convert to reference time units
      stt = cd_convert(stt, ref_utime)
      ent = cd_convert(ent, ref_utime)

      ; Create time series
      if (fID.eq.0) then
        start_time = stt
        end_time = ent
        var = value
      else
        tmp = array_append_record(start_time, stt, 0)
        delete(start_time)
        start_time = tmp
        delete(tmp)
        tmp = array_append_record(end_time, ent, 0)
        delete(end_time)
        end_time = tmp
        delete(tmp)
        tmp = array_append_record(var, value, 0)
        delete(var)
        var = tmp
        delete(tmp)
      end if
      delete(stt)
      delete(ent)
      delete(value)

    end do
    delete(in_files)
    delete(in_bnames)

    ; Check monotonicity
    if (isMonotonic(start_time).ne.1) then
      error_msg("f", scriptname, funcname, \
                "non-monotonically increasing time-series, possible " + \
                "duplicated data in input")
    end if

    ; Calculate monthly mean from the time series
    do vv = 0, dimsizes(var) - 1

      if (ismissing(var(vv))) then
        continue
      end if

      pstart = start_time(vv)
      pend = -1000.d0

      do while (pend.lt.end_time(vv))

        wdate = cd_calendar(pstart, -5)
        wdatec = cd_calendar(pstart, -1)

        ; Find beginning of next month
        if (wdate(0, 1).eq.12) then
          wyear = wdate(0, 0) + 1
          wmonth = 1
        else
          wyear = wdate(0, 0)
          wmonth = wdate(0, 1) + 1
        end if
        pend = cd_inv_calendar(wyear, wmonth, 1, 0, 0, 0, ref_utime, 0)

        if (pend.gt.end_time(vv)) then
          pend = (/end_time(vv)/)
        end if

        didx = ind(wdatec.eq.datec)
        if (wdate(0, 0).lt.y1 .or. wdate(0, 0).gt.y2) then
          pstart = pend
          continue
        end if
        nd = tofloat(pend - pstart)
        temp_arr(didx) = temp_arr(didx) + var(vv) * nd
        ndays_arr(didx) = ndays_arr(didx) + nd

        ; DEBUG+++
        ; print(start_time(vv) +"  "+end_time(vv) + "  " + "(" + \
        ;       cd_calendar(start_time(vv), -2) + "-" + \
        ;       cd_calendar(end_time(vv), -2)+") " + datec(didx) + \
        ;       " nd="+nd)
        ; DEBUG---

        pstart = pend

      end do
      delete(pstart)
      delete(pend)

    end do

    delete(var)
    delete(start_time)
    delete(end_time)

    ; Calculate monthly mean
    temp_arr = where(temp_arr.eq.0, temp_arr@_FillValue, temp_arr)
    ndays_arr = where(ndays_arr.eq.0, ndays_arr@_FillValue, ndays_arr)
    temp_arr = temp_arr / ndays_arr

    ; Assign to global data array
    idx_nm = ind(.not.ismissing(temp_arr))
    if (all(ismissing(idx_nm))) then
      delete(idx_nm)
      continue
    end if
    do ii = 0, dimsizes(idx_nm) - 1
      data_arr(idx_nm(ii), vID) = temp_arr(idx_nm(ii))
    end do
    delete(idx_nm)

  end do  ; variables

  return(data_arr)

end
