; #############################################################################
; FUNCTIONS FOR THE reformat_obs_*.ncl SCRIPTS
; #############################################################################
; General purpose functions called by the reformat_obs_*.ncl scripts.
;
; Contents
;    procedure format_time
;    procedure format_plev
;    procedure format_lev
;    procedure format_lat
;    procedure format_lon
;    procedure format_coords

;    function plev_attrib
;    function lev_attrib
;    function lat_attrib
;    function lon_attrib
;    function read_cmor
;    function var_attrib
;    function set_size_array
;    function create_timec
;    function set_global_atts
;    procedure write_nc
;    procedure write_nc_profile
;    function process_EBAS_data
;
; #############################################################################

; Time units
TUNITS = "days since 1950-01-01 00:00:00"

; CMOR FillValue
FILL = 1.e+20

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
;    20190216-A_righ_ma: written.
;
local funcname, scriptname, ctime, ntime, year1, month1, day1, year2, month2, \
  day2, calendar, date, exp_ntime, yy, mm, m1, m2, dd, opt, newtime
begin

  funcname = "format_time"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

  if (all(frequency.ne.(/"3hr", "6hr", "day", "mon", "yr"/))) then
    error_msg("f", scriptname, funcname, "unsupported frequency " + frequency)
  end if

  ; Read coordinate
  ctime = var&time
  ntime = dimsizes(ctime)

  ; Check monotonicity
  if (isMonotonic(ctime).ne.1) then
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
                "incorrect number of timesteps in input data: " + exp_ntime + \
                "expected, " + ntime + " found")
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
    if (year1.eq.year2) then
      exp_ntime = month2 - month1 + 1
    else
      exp_ntime = (12 - month1 + 1) + month2 + \
        12 * where((year2 - year1 - 1).gt.0, year2 - year1 - 1, 0)
    end if
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + exp_ntime + \
                "expected, " + ntime + " found")
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
        dm = days_in_month(yy, mm) / 2.
        date(tt, 2) = toint(dm)
        date(tt, 3) = toint(24 * (dm - toint(dm)))
        tt = tt + 1
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
      cd_inv_calendar(year2, month2, day2, 0, 0, 0, TUNITS, opt) - \
      cd_inv_calendar(year1, month1, day1, 0, 0, 0, TUNITS, opt) + 1
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + exp_ntime + \
                "expected, " + ntime + " found")
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
        do dd = day1, day2
          date(tt, 2) = dd
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
      cd_inv_calendar(year2, month2, day2, 0, 0, 0, TUNITS, opt) - \
      cd_inv_calendar(year1, month1, day1, 0, 0, 0, TUNITS, opt) + 1
    exp_ntime = 4 * exp_ntime
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + exp_ntime + \
                "expected, " + ntime + " found")
    end if
    delete(opt)

  end if

  ; 3-hourly frequency (check size only)
  if (frequency.eq."3hr") then

    ; Check size
    opt = 0
    opt@calendar = calendar
    exp_ntime = \
      cd_inv_calendar(year2, month2, day2, 0, 0, 0, TUNITS, opt) - \
      cd_inv_calendar(year1, month1, day1, 0, 0, 0, TUNITS, opt) + 1
    exp_ntime = 8 * exp_ntime
    if (ntime.ne.exp_ntime) then
      error_msg("f", scriptname, funcname, \
                "incorrect number of timesteps in input data: " + exp_ntime + \
                "expected, " + ntime + " found")
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

  leave_msg(scriptname, funcname)

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
;    20190216-A_righ_ma: written.
;
local funcname, scriptname, cplev, newplev
begin

  funcname = "format_plev"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

  ; Check monotonicity
  if (isMonotonic(var&plev) .eq. 0) then
    error_msg("f", DIAG_SCRIPT, "", "non-monotonic vertical coordinate")
  end if
  if (isMonotonic(var&plev).eq.1) then
    var = var(:, ::-1, :, :)
  end if

  ; Read coordinate
  cplev = var&plev

  ; Set standard attributes
  newplev = todouble(cplev)  ; this also removes attributes
  copy_VarCoords(cplev, newplev)
  newplev@bnds = "plev_bnds"
  newplev@positive = "down"
  newplev@axis = "Z"
  newplev@long_name = "pressure"
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
;    20190216-A_righ_ma: written.
;
begin

  funcname = "format_coords"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

  dnames = getvardimnames(var)

  ; Loop over dimensions, call formatting procedures
  do dd = 0, dimsizes(dnames) - 1

    found = False

    if (dnames(dd).eq."time") then
      format_time(var, date1, date2, frequency)
      found = True
    end if

    if (dnames(dd).eq."plev") then
;      format_plev(var)
      found = True
    end if

    if (dnames(dd).eq."lev") then
;      format_lev(var)
      found = True
    end if

    if (dnames(dd).eq."lat") then
;      format_lat(var)
      found = True
    end if

    if (dnames(dd).eq."lon") then
;      format_lon(var)
      found = True
    end if

    if (.not.found) then
      error_msg("f", scriptname, funcname, "cannot format coordinate " + \
                dnames(dd))
    end if

  end do

  leave_msg(scriptname, funcname)

end


; #############################################################################
undef("lev_attrib")
function lev_attrib(coord:numeric)
;
; Arguments
;    coord: input ocean level coordinate.
;
; Return value
;    An array of the same dimensionality of coord.
;
; Description
;    Assigns attributes to the ocean level coordinate.
;
; Caveats
;
; References
;
; Modification history
;    20150415-A_righ_ma: written.
;
local funcname, scriptname, out_coord
begin

  funcname = "lev_attrib"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

  out_coord = todouble(coord)  ; this also removes attributes
  copy_VarCoords(coord, out_coord)
  out_coord@long_name = "ocean depth coordinate"
  out_coord@positive = "down"
  out_coord@axis = "Z"
  out_coord@units = "m "
  out_coord@standard_name = "depth"
  if (isatt(out_coord, "_FillValue")) then
    delete(out_coord@_FillValue)
  end if

  leave_msg(scriptname, funcname)
  return(out_coord)

end

; #############################################################################
undef("lat_attrib")
function lat_attrib(coord)
;
; Arguments
;    coord: input latitude coordinate.
;
; Return value
;    An array of the same dimensionality of coord.
;
; Description
;    Assigns attributes to the latitude coordinate.
;
; Caveats
;
; References
;
; Modification history
;    20130528-A_righ_ma: written.
;
local funcname, scriptname, out_coord
begin

  funcname = "lat_attrib"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

  out_coord = todouble(coord)  ; this also removes attributes
  copy_VarCoords(coord, out_coord)
  out_coord@long_name = "latitude"
  out_coord@axis = "Y"
  out_coord@units = "degrees_north"
  out_coord@standard_name = "latitude"
  if (isatt(out_coord, "_FillValue")) then
    delete(out_coord@_FillValue)
  end if

  leave_msg(scriptname, funcname)
  return(out_coord)

end

; #############################################################################
undef("lon_attrib")
function lon_attrib(coord)
;
; Arguments
;    coord: input longitude coordinate.
;
; Return value
;    An array of the same dimensionality of coord.
;
; Description
;    Assigns attributes to the longitude coordinate.
;
; Caveats
;
; References
;
; Modification history
;    20130528-A_righ_ma: written.
;
local funcname, scriptname, out_coord
begin

  funcname = "lon_attrib"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

  out_coord = todouble(coord)  ; this also removes attributes
  copy_VarCoords(coord, out_coord)
  out_coord@long_name = "longitude"
  out_coord@axis = "X"
  out_coord@units = "degrees_east"
  out_coord@standard_name = "longitude"
  if (isatt(out_coord, "_FillValue")) then
    delete(out_coord@_FillValue)
  end if

  leave_msg(scriptname, funcname)
  return(out_coord)

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
;    20190107-A_righ_ma: modify to read standard CMIP5 tables
;    20130528-A_righ_ma: written.
;
local funcname, scriptname, data, idxu, idxd, attn, attv, out
begin

  funcname = "read_cmor"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

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

  leave_msg(scriptname, funcname)
  return(out)

end

; #############################################################################
undef("var_attrib")
function var_attrib(var:numeric,
                    name:string,
                    table:string)
;
; Arguments
;    var: the variable's field.
;    name: the variable's standard name.
;    string: full path to the CMOR table of the variable.
;
; Return value
;    An array of the same dimensionality of var.
;
; Description
;    Sets standard variable attributes.
;
; Caveats
;
; References
;
; Modification history
;    20190107-A_righ_ma: add extra argument for CMOR table
;    20161202-A_laue_ax: preserve attribute "coordinates" if present
;    20130528-A_righ_ma: written.
;
local funcname, scriptname, coordattr, out, tmp, att, ii
begin

  funcname = "var_attrib"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

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

  leave_msg(scriptname, funcname)
  return(out)

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
;    Sets a logarithmic array of sizes to be used for particle size
;    distribution calculations.
;
; Caveats
;
; References
;
; Modification history
;    20130528-A_righ_ma: written.
;
local funcname, scriptname, minsize, maxsize, nbins, bin, out
begin

  funcname = "set_size_array"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

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

  leave_msg(scriptname, funcname)
  return(out)

end

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
;    20140124-A_righ_ma: written.
;
local funcname, scriptname, yy, mm, out
begin

  funcname = "create_timec"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

  out = new(12 * (y2 - y1 + 1), double)
  do yy = y1, y2
    do mm = 1, 12
      out(12 * (yy - y1) + mm - 1) = \
        cd_inv_calendar(yy, mm, 15, 0, 0, 0, TUNITS, 0)
    end do
  end do

  leave_msg(scriptname, funcname)
  return(out)

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
;    20190202-A_righ_ma: written.
;
local funcname, scriptname,
    dim_unlim, ii
begin

  funcname = "set_global_atts"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

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

  leave_msg(scriptname, funcname)

end


; #############################################################################
undef("write_nc")
procedure write_nc(outfile:string,
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
;    Writes the given variable to the given NetCDF file, appending also the
;    provided global attributes.
;
; Modification history
;    20140123-A_righ_ma: written.
;
local funcname, scriptname, w, gAtt, dim_names, ndims, dim_sizes, dim_types, \
    dim_unlim, ii
begin

  funcname = "write_nc"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

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
  do ii = 0, ndims - 1
    filevardef(w, dim_names(ii), dim_types(ii), dim_names(ii))
  end do
  filevardef(w, name, "float", dim_names)
  do ii = 0, ndims - 1
    filevarattdef(w, dim_names(ii), var&$dim_names(ii)$)
  end do
  filevarattdef(w, name, var)

  ; Append global attributes
  fileattdef(w, gAtt)

  ; Write
  setfileoption(w, "DefineMode", False)
  do ii = 0, ndims - 1
    w->$dim_names(ii)$ = (/var&$dim_names(ii)$/)
  end do
  w->$name$ = (/var/)

  leave_msg(scriptname, funcname)

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
;    Writes the given variable to the given NetCDF file, appending also the
;    provided global attributes.
;    Designed to write multiple variables for the vertical profiles data.
;
; Modification history
;    20140422-A_righ_ma: written.
;
local funcname, scriptname, w, coords, cc, jj, locname, locvar, cname
begin

  funcname = "write_nc_profile"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

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

  leave_msg(scriptname, funcname)

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
;    This function is used to process the data from the EBAS database (e.g.,
;    EANET, EMEP).
;
; Caveats
;    For the time coordinate in the input data, only units of days are
;    currently accepted.
;
; Modification history
;    20150413-A_righ_ma: improved time selection.
;    20140124-A_righ_ma: written.
;
local timec, datec, vID, fID, bn, en, head, hh, cline, syear, smonth, sday, \
  scode, comp, matr, unit, scale, fills, lline, cols, data_col, flag_col, \
  start_col, end_col, data, value, flag, start_time, end_time, mm, sidx, \
  monthind, stday, enday, nd, pt1, pt2, data_arr
begin

  funcname = "process_EBAS_data"
  scriptname = "esmvaltool/utils/cmorizers/obs/utilities.ncl"
  enter_msg(scriptname, funcname)

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

  leave_msg(scriptname, funcname)
  return(data_arr)

end
