; #############################################################################
; ESMValTool CMORizer for ESACCI-FIRE data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/fire/data/
;
; Last access
;    20190124
;
; Download and processing instructions
;    Download the data from:
;      burned_area/MERIS/grid/v4.1/
;    Put all files in input_dir_path (no subdirectories with years).
;
; Modification history
;    20190124-righi_mattia: written based on a python script by Ben Mueller.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "esacci_fire.ncl"

  ; Source name
  OBSNAME = "ESACCI-FIRE"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 2005
  YEAR2 = 2011

  ; Selected variable (standard name)
  VAR = "burntArea"

  ; Name in the raw data
  NAME = "burned_area"

  ; MIP
  MIP = "Lmon"

  ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_Lmon"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L4-BA-MERIS-fv4.1"

  ; Global attributes
  SOURCE = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/fire/data/"
  REF = "Chuvieco et al.ESA Fire Climate Change Initiative (Fire_cci): " + \
    "Burned Area Grid Product Version 4.1. " + \
    "Centre for Environmental Data Analysis, " + \
    "doi:10.5285/D80636D4-7DAF-407E-912D-F5BB61C142FA, 2016."
  COMMENT = ""

end

begin

  time = create_timec(YEAR1, YEAR2)
  date = cd_calendar(time, 1)

  do yy = YEAR1, YEAR2
    do mm = 1, 12

      ldate = yy + sprinti("%0.2i", mm)

      files = systemfunc("ls " + input_dir_path + ldate + \
                         "??-ESACCI-L4_FIRE-BA-MERIS-fv4.1.nc")
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

  ; Format coordinates
  output!0 = "time"
  output!1 = "lat"
  output!2 = "lon"
  format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ)

  ; Set variable attributes
  tmp = format_variable(output, VAR, CMOR_TABLE)
  delete(output)
  output = tmp
  delete(tmp)

  ; Calculate coordinate bounds
  bounds = guess_coord_bounds(output, FREQ)

  ; Set global attributes
  gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

  ; Output file
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE, VERSION, \
               MIP, VAR, DATESTR/), "_") + ".nc"

  ; Write variable
  write_nc(fout, VAR, output, bounds, gAtt)
  delete(gAtt)
  delete(output)
  delete(bounds)

end
