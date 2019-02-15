; #############################################################################
; ESMValTool CMORizer for HadCRUT data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://crudata.uea.ac.uk/cru/data/temperature/
;
; Last access
;    20190208
;
; Download and processing instructions
;    Download the dataset "HadCRUT4" (median temperature anomalies) and
;    the dataset "Absolute" (absolute temperatures for the base period
;    1961-90 on a 5x5 grid).
;
; Caveats
;    In contrast to the HadCRUT3 reformat script which produces temperature
;    anomalies (relative to the 1961-90 climatology), this script calculates
;    absolute tas by adding the climatology ("absolute.nc") to the anomalies
;    ("HadCRUT.4.6.0.0.median.nc"). It creates two output, one with the
;    temperature time-series and one with the climatology of the temperature
;    anomaly.
;
; Modification history
;    20190208-A_righ_ma: added output for climatology and adapted to v2.
;    20180222-A_laue_ax: bug fix (added swapping of latitudes if needed).
;    20160203-A_laue_ax: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Source name
  OBSNAME = "HadCRUT4"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1850
  YEAR2 = 2018

  ; Selected variable (standard name)
  VAR = "tas"

  ; MIPS
  MIP = "Amon"

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/cmip5/Tables/CMIP5_Amon"

end

begin

  ; Read file
  fname1 = input_dir_path + "HadCRUT.4.6.0.0.median.nc"
  fname2 = input_dir_path + "absolute.nc"
  f1 = addfile(fname1, "r")
  setfileoption("nc", "MissingToFillValue", False)
  f2 = addfile(fname2, "r")

  ; Read anomaly
  anomaly = f1->temperature_anomaly

  ; Read absolute temperature
  tmp = f2->tem
  clim = tofloat(tmp * tmp@scale_factor) + 273.15
  copy_VarCoords(tmp, clim)
  delete(tmp)

  ; Swap latitudes
  if (isMonotonic(anomaly&latitude).eq.-1) then
    anomaly = anomaly(:, ::-1, :)
  end if

  if (isMonotonic(clim&lat).eq.-1) then
    clim = clim(:, ::-1, :)
  end if

  log_info("  Climatology range: " + min(clim) + \
           " K to " + max(clim) + " K")
  log_info("  Anomaly range: " + min(anomaly) + \
           " K to " + max(anomaly) + " K")

  output1 = anomaly
  output2 = clim
  dims = dimsizes(output1)

  ; Add absolute temperature to anomaly
  do yr = 0, dims(0) / 12 - 1
    m1 = yr * 12
    m2 = m1 + 11
    output1(m1:m2, :, :) = where(.not.ismissing(clim), \
                                 anomaly(m1:m2, :, :) + clim, \
                                 tofloat(anomaly@_FillValue))
  end do

  ; Format time coordinate
  ctime = time_attrib(output1&time, YEAR1, YEAR2, "M")
  delete(output1&time)
  output1&time = ctime
  delete(ctime)

  ; Set date to center year of the climatology (1961-1990 --> 1976)
  delete(output2&time)
  output2&time = create_timec(1976, 1976)
  ctime = time_attrib(output2&time, 1976, 1976, "M")
  delete(output2&time)
  output2&time = ctime
  delete(ctime)

  ; Format latitude coordinate
  output1!1 = "lat"
  if (isMonotonic(output1&lat) .eq. 0) then
    error_msg("f", DIAG_SCRIPT, "", "non-monotonic latitude coordinate")
  end if
  if (isMonotonic(output1&lat) .eq. -1) then
    output1 = output1(:, ::-1, :)
  end if
  clat = lat_attrib(output1&lat)
  delete(output1&lat)
  output1&lat = clat
  delete(clat)

  output2!1 = "lat"
  if (isMonotonic(output2&lat) .eq. 0) then
    error_msg("f", DIAG_SCRIPT, "", "non-monotonic latitude coordinate")
  end if
  if (isMonotonic(output2&lat) .eq. -1) then
    output2 = output2(:, ::-1, :)
  end if
  clat = lat_attrib(output2&lat)
  delete(output2&lat)
  output2&lat = clat
  delete(clat)

  ; Format longitude coordinate
  output1!2 = "lon"
  if (isMonotonic(output1&lon) .eq. 0) then
    error_msg("f", DIAG_SCRIPT, "", "non-monotonic longitude coordinate")
  end if
  if (any(output1&lon.lt.0.)) then
    output1 = lonFlip(output1)
  end if
  clon = lon_attrib(output1&lon)
  delete(output1&lon)
  output1&lon = clon
  delete(clon)

  output2!2 = "lon"
  if (isMonotonic(output2&lon) .eq. 0) then
    error_msg("f", DIAG_SCRIPT, "", "non-monotonic longitude coordinate")
  end if
  if (any(output2&lon.lt.0.)) then
    output2 = lonFlip(output2)
  end if
  clon = lon_attrib(output2&lon)
  delete(output2&lon)
  output2&lon = clon
  delete(clon)

  ; Set variable attributes
  tmp = var_attrib(output1, VAR, CMOR_TABLE)
  output1 = tmp
  delete(tmp)

  tmp = var_attrib(output2, VAR, CMOR_TABLE)
  output2 = tmp
  delete(tmp)

  ; Add height coordinate
  output1@coordinates = "height"
  output2@coordinates = "height"
  height = 2.d
  height!0 = "ncl_scalar"
  height@units = "m"
  height@axis = "Z"
  height@positive = "up"
  height@long_name = "height"
  height@standard_name = "height"

  ; Write temperature time-series
  gAtt = set_global_atts( \
    OBSNAME, TIER, \
    "https://crudata.uea.ac.uk/cru/data/temperature/", \
    "Morice et al., J. Geophys. Res., doi:10.1029/2011JD017187, 2012", \
    "Temperature time-series calculated from the anomaly time-series by " + \
    "adding the temperature climatology for 1961-1990")
  fout = output_dir_path + "OBS_HadCRUT4_ground_1_" + MIP + "_" + VAR + "_" + \
    YEAR1 + "01-" + YEAR2 + "12.nc"
  write_nc(fout, VAR, output1, gAtt)
  w = addfile(fout, "w")
  w->height = height
  delete(w)
  delete(gAtt)
  delete(output1)

  ; Write temperature climatology
  gAtt = set_global_atts( \
    OBSNAME, TIER, \
    "https://crudata.uea.ac.uk/cru/data/temperature/", \
    "Jones et al., Rev. Geophys., doi:10.1029/1999RG900002, 1999", \
    "Climatology 1961-1990")
  fout = output_dir_path + "OBS_HadCRUT4_clim_1_" + MIP + "_" + VAR + "_" + \
    "197601-197612.nc"
  write_nc(fout, VAR, output2, gAtt)
  w = addfile(fout, "w")
  w->height = height
  delete(w)
  delete(gAtt)
  delete(output2)

end
