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
;    20201125
;
; Download and processing instructions
;    Download the dataset "HadCRUT4" (median temperature anomalies) and
;    the dataset "Absolute" (absolute temperatures for the base period
;    1961-90 on a 5x5 grid).
;    For the 5% to 95% confidence interval of the combined effects of all the
;    uncertainties described in the HadCRUT4 error model (measurement and
;    sampling, bias, and coverage uncertainties) download:
;    HadCRUT.4.6.0.0.annual_ns_avg.txt from
;    https://www.metoffice.gov.uk/hadobs/hadcrut4/data/current/download.html
;
; Caveats
;    In contrast to the HadCRUT3 reformat script which produces temperature
;    anomalies (relative to the 1961-90 climatology), this script calculates
;    absolute tas by adding the climatology ("absolute.nc") to the anomalies
;    ("HadCRUT.4.6.0.0.median.nc"). It creates 2 output, one with the
;    temperature time-series and one with the anomaly time-series
;
; Modification history
;    20201125-bock_lisa: add tas uncertainty
;    20190916-righi_mattia: remove no-longer used climatology output.
;    20190229-righi_mattia: added output for anomaly (tasa).
;    20190208-righi_mattia: added output for climatology and adapted to v2.
;    20180222-lauer_axel: bug fix (added swapping of latitudes if needed).
;    20160203-lauer_axel: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_hadcrut4.ncl"

  ; Source name
  OBSNAME = "HadCRUT4"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1850
  YEAR2 = 2018

  ; Selected variable (standard name)
  VAR = (/"tas", "tasa", "tasConf5", "tasConf95"/)

  ; MIP
  MIP = "Amon"

  ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE1 = getenv("cmor_tables") + \
    "/cmip5/Tables/CMIP5_Amon"
  CMOR_TABLE2 = getenv("cmor_tables") + \
    "/custom/CMOR_tasa.dat"
  CMOR_TABLE3 = getenv("cmor_tables") + \
    "/custom/CMOR_tasConf5.dat"
  CMOR_TABLE4 = getenv("cmor_tables") + \
    "/custom/CMOR_tasConf95.dat"

  ; Version
  VERSION = "1"

  ; Type
  TYPE1 = "ground"
  TYPE2 = "ground"

  ; Global attributes
  SOURCE = "https://crudata.uea.ac.uk/cru/data/temperature/"
  REF1 =  "Morice et al., J. Geophys. Res., doi:10.1029/2011JD017187, 2012"
  REF2 =  "Morice et al., J. Geophys. Res., doi:10.1029/2011JD017187, 2012"
  COMMENT1 = "Temperature time-series calculated from the anomaly " + \
    "time-series by adding the temperature climatology for 1961-1990"
  COMMENT2 = "Temperature anomaly with respect to the period 1961-1990"
  COMMENT3 = "Yearly uncertainty of temperature time-series"

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
  output2 = anomaly
  dims = dimsizes(output1)

  ; Add absolute temperature to anomaly
  do yr = 0, dims(0) / 12 - 1
    m1 = yr * 12
    m2 = m1 + 11
    output1(m1:m2, :, :) = where(.not.ismissing(clim), \
                                 anomaly(m1:m2, :, :) + clim, \
                                 tofloat(anomaly@_FillValue))
  end do

  ; Format coordinates
  output1!0 = "time"
  output1!1 = "lat"
  output1!2 = "lon"
  format_coords(output1, YEAR1 + "0101", YEAR2 + "1231", FREQ)
  output2!0 = "time"
  output2!1 = "lat"
  output2!2 = "lon"
  format_coords(output2, YEAR1 + "0101", YEAR2 + "1231", FREQ)

  ; Calculate coordinate bounds
  bounds1 = guess_coord_bounds(output1, FREQ)
  bounds2 = guess_coord_bounds(output2, FREQ)

  ; Set variable attributes
  tmp = format_variable(output1, VAR(0), CMOR_TABLE1)
  delete(output1)
  output1 = tmp
  delete(tmp)
  tmp = format_variable(output2, VAR(1), CMOR_TABLE2)
  delete(output2)
  output2 = tmp
  delete(tmp)

  ; Add height coordinate
  output1@coordinates = "height"
  height = 2.d
  height!0 = "ncl_scalar"
  height@units = "m"
  height@axis = "Z"
  height@positive = "up"
  height@long_name = "height"
  height@standard_name = "height"

  ; Set global attributes
  gAtt1 = set_global_atts(OBSNAME, TIER, SOURCE, REF1, COMMENT1)
  gAtt2 = set_global_atts(OBSNAME, TIER, SOURCE, REF2, COMMENT2)

  ; Write temperature time-series
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE1, VERSION, \
               MIP, VAR(0), DATESTR/), "_") + ".nc"
  write_nc(fout, VAR(0), output1, bounds1, gAtt1)
  w = addfile(fout, "w")
  w->height = height
  delete(w)
  delete(gAtt1)
  delete(bounds1)
  delete(output1)

  ; Write temperature anomaly time-series
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE2, VERSION, \
               MIP, VAR(1), DATESTR/), "_") + ".nc"
  write_nc(fout, VAR(1), output2, bounds2, gAtt2)
  w = addfile(fout, "w")
  delete(w)
  delete(gAtt2)
  delete(bounds2)
  delete(output2)

  ; -------------------------------------
  ; Uncertainties
  ; -------------------------------------

  ; Read uncertainty file
  fname = input_dir_path + "HadCRUT.4.6.0.0.annual_ns_avg.txt"

  ntime = YEAR2 - YEAR1 + 1
  tmp = asciiread(fname, (/ntime, 12/), "float")

  ntime2 = 12 * (YEAR2 - YEAR1 + 1)
  data1 = new(ntime2, float)
  data2 = new(ntime2, float)

  i = 0
  do yy = 0, ntime - 1
    data1((yy * 12) : (yy * 12 + 11)) = tmp(yy, 1) - tmp(yy, 10)
    data2((yy * 12) : (yy * 12 + 11)) = tmp(yy, 11) - tmp(yy, 1)
  end do

  delete(tmp)

  data1!0 = "time"
  data1&time =  create_timec(YEAR1, YEAR2)
  format_coords(data1, YEAR1 + "0101", YEAR2 + "1231", FREQ)
  data2!0 = "time"
  data2&time = data1&time

  ; Set variable attributes
  tmp = format_variable(data1, VAR(2), CMOR_TABLE3)
  delete(data1)
  data1 = tmp
  delete(tmp)
  tmp = format_variable(data2, VAR(3), CMOR_TABLE4)
  delete(data2)
  data2 = tmp
  delete(tmp)

  ; Calculate coordinate bounds
  bounds = guess_coord_bounds(data1, FREQ)
  bounds = guess_coord_bounds(data2, FREQ)

  ; Set global attributes
  gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF1, COMMENT3)

  ; Write temperature time-series uncertainties
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE1, VERSION, \
               MIP, VAR(2), DATESTR/), "_") + ".nc"
  write_nc(fout, VAR(2), data1, bounds, gAtt)
  w = addfile(fout, "w")
  delete(w)
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE1, VERSION, \
               MIP, VAR(3), DATESTR/), "_") + ".nc"
  write_nc(fout, VAR(3), data2, bounds, gAtt)
  w = addfile(fout, "w")
  delete(w)
  delete(gAtt)
  delete(bounds)
  delete(data1)
  delete(data2)

end
