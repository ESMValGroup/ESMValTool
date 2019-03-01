; #############################################################################
; ESMValTool CMORizer for HadCRUT data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://www.esrl.noaa.gov/psd/data/gridded/data.ghcngridded.html
;
; References
;    Jones, P. D., and A. Moberg, 2003: Hemispheric and Largescale Surface Air
;      Temperature Variations: An extensive Revision and an Update to 2001,
;      J. Climate, 16, 206-223.
;    Peterson, Thomas C. and Russell S. Vose, 1997: An Overview of the Global
;      Historical Climatology Network Temperature Data Base, Bulletin of the
;      American Meteorological Society, 78, 2837-2849.
;
; Last access
;    20190227
;
; Download and processing instructions
;    Download the dataset "precip.mon.total.nc" (precipitation, total, surface,
;    1900-2015 on a 5x5 grid).
;
; Caveats
;
; Modification history
;    20190227-A_bock_ls: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Source name
  OBSNAME = "GHCN"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1900
  YEAR2 = 2014

  ; Selected variable (standard name)
  VAR = "pr"

  ; MIPS
  MIP = "Amon"

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/cmip5/Tables/CMIP5_Amon"

end

begin

  ; Read file
  fname = input_dir_path + "precip.mon.total.nc"
  f = addfile(fname, "r")
  setfileoption("nc", "MissingToFillValue", False)

  ; Read absolute precipitation without last incomplete year
  output = f->precip(time|0:1379, lat|:, lon|:)

  ; Calculate days per month
  date = cd_calendar(output&time, 0)
  dpm = days_in_month(toint(date(:, 0)), toint(date(:, 1)))
  dpmc = conform(output, dpm, 0)

  ; Check time range
  if (dimsizes(date(:, 0)).ne.12 * (YEAR2 - YEAR1 + 1)) then
    error_msg("f", cmorization_script, "", "incorrect number of timesteps")
  end if

  ; Convert units [mm/month] --> [kg/m2/s]
  output = output / (24 * 3600 * dpmc)

  log_info("  Climatology range: " + min(output) + \
           " kg/m2/s to " + max(output) + " kg/m2/s")

  dims = dimsizes(output)

  ; Format time coordinate
  ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
  delete(output&time)
  output&time = ctime
  delete(ctime)

  ; Format latitude coordinate
  output!1 = "lat"
  if (isMonotonic(output&lat) .eq. 0) then
    error_msg("f", cmorization_script, "", "non-monotonic latitude coordinate")
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
    error_msg("f", cmorization_script, "", \
              "non-monotonic longitude coordinate")
  end if
  if (any(output&lon.lt.0.)) then
    output1 = lonFlip(output)
  end if
  clon = lon_attrib(output&lon)
  delete(output&lon)
  output&lon = clon
  delete(clon)

  ; Set variable attributes
  tmp = var_attrib(output, VAR, CMOR_TABLE)
  output = tmp
  delete(tmp)

  ; Write temperature time-series
  gAtt = set_global_atts( \
    OBSNAME, TIER, \
    "https://www.esrl.noaa.gov/psd/data/gridded/data.ghcngridded.html", \
    "Jones, P. D., and A. Moberg, J. Climate, 2003", \
    "Peterson, Thomas C. and Russell S. Vose, Bull. Amer. Meteor. Soc., 1997")
  fout = output_dir_path + "OBS_GHCN_ground_1_" + MIP + "_" + VAR + "_" + \
    YEAR1 + "01-" + YEAR2 + "12.nc"
  write_nc(fout, VAR, output, gAtt)
  delete(gAtt)
  delete(output)

end
