; #############################################################################
; ESMValTool CMORizer for ESACCI-AEROSOL data
; #############################################################################
;
; Tier
;   Tier 2: other freely-available dataset.
;
; Source
;   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/fire/data/
;
; Last access
;   20190124
;
; Download and processing instructions
;   Download the data from:
;     burned_area/MERIS/grid/v4.1/
;   Put all files in input_dir_path (no subdirectories with years).
;
; Modification history
;   20190124-A_righ_ma: written based on a python script by muel_bn.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ESACCI-FIRE.ncl"

  ; Source name
  OBSNAME = "ESACCI-FIRE"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 2005
  YEAR2 = 2011

  ; Field
  FIELD = "T2Ms"

  ; Selected variable (standard name)
  VARS = "burntArea"

  ; Name in the raw data
  NAME = "burned_area"

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/cmip5/Tables/CMIP5_Lmon"

end

begin

  time = create_timec(YEAR1, YEAR2)
  date = cd_calendar(time, 1)

  do yy = YEAR1, YEAR2
    do mm = 1, 12

      ldate = yy + sprinti("%0.2i", mm)

      files = systemfunc("ls " + input_dir_path + "/" + \
                         ldate + "??-ESACCI-L4_FIRE-BA-MERIS-fv4.1.nc")
      f = addfiles(files, "r")

      xx = f[:]->$NAME$

      ; Calculate area
      if (.not.isdefined("area")) then
        deg2rad = acos(-1.0) / 180.
        lat = f[0]->lat
        lon = f[0]->lon
        nlat = dimsizes(lat)
        deltax = abs(lon(1) - lon(0))
        lati = new(dimsizes(lat) + 1, float)
        lati(0) = max((/(3 * lat(0) - lat(1)) / 2., -90./))
        do ii = 1, dimsizes(lati) - 2
          lati(ii) = 0.5 * (lat(ii - 1) + lat(ii))
        end do
        lati(dimsizes(lati) - 1) = \
          min((/(3 * lat(nlat - 1) - lat(nlat - 2)) / 2., 90./))
        area = new((/dimsizes(lat), dimsizes(lon)/), float)
        do ii = 0, dimsizes(lat) - 1
          deltay = sin(lati(ii + 1) * deg2rad) - sin(lati(ii) * deg2rad)
          area(ii, :) = abs(6371000. ^ 2 * deltay * deltax * deg2rad)
        end do
        delete([/lat, lon, nlat, deltax, lati, deltay/])
      end if

      ; Calculate fraction
      xx = xx / (/conform(xx, area, (/1, 2/))/)

      ; Assign to global array
      if (.not.isdefined("output")) then
        dims = array_append_record(dimsizes(time), dimsizes(xx(0, :, :)), 0)
        output = new(dims, float)
        output!0 = "time"
        output&time = time
        output!1 = "lat"
        output&lat = f[0]->lat
        output!2 = "lon"
        output&lon = f[0]->lon
      end if
      output(ind(toint(ldate).eq.date), :, :) = dim_sum_n(xx, 0)

      delete(files)
      delete(f)

    end do
  end do

  ; Convert units [1] --> [%]
  output = output * 100.

  ; Format time coordinate
  ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
  delete(output&time)
  output&time = ctime
  delete(ctime)

  ; Format latitude coordinate
  output!1 = "lat"
  if (isMonotonic(output&lat) .eq. 0) then
    error_msg("f", DIAG_SCRIPT, "", "non-monotonic latitude coordinate")
  end if
  if (isMonotonic(output&lat) .eq. -1) then
    output = output(:, ::-1, :)
  end if
  clat = lat_attrib(output&lat)
  delete(output&lat)
  output&lat = clat
  delete(clat)

  ; Format longitude coordinate
  output!2 = "lon"
  if (isMonotonic(output&lon) .eq. 0) then
    error_msg("f", DIAG_SCRIPT, "", "non-monotonic longitude coordinate")
  end if
  if (any(output&lon.lt.0.)) then
    output = lonFlip(output)
  end if
  clon = lon_attrib(output&lon)
  delete(output&lon)
  output&lon = clon
  delete(clon)

  ; Set variable attributes
  tmp = var_attrib(output, VARS, CMOR_TABLE)
  delete(output)
  output = tmp
  delete(tmp)

    ; Set global attributes
    gAtt = True
    gAtt@history = "Created on " + systemfunc("date")
    gAtt@host    = systemfunc("echo $HOST")
    gAtt@user    = systemfunc("echo $USER")
    gAtt@period  = YEAR1 + "-" + YEAR2
    gAtt@field   = FIELD
    gAtt@tier    = TIER
    gAtt@source  = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/fire/data/"
    gAtt@reference = ""
    gAtt@title = OBSNAME + " data reformatted for the ESMValTool"
    gAtt@conventions = "CF/CMOR"

    ; Outfile
    fout = output_dir_path + "/OBS_" + OBSNAME + "_sat_L4-BA-MERIS-fv4.1_" + \
      FIELD + "_" + VARS + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ; Write variable
    write_nc(fout, VARS, output, gAtt)
    delete(gAtt)
    delete(output)

end
