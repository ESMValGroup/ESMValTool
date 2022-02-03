; #############################################################################
; ESMValTool CMORizer for NOAAGlobalTemp data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    "https://www.ncdc.noaa.gov/noaa-merged-land-ocean-global-surface-"+
;      "temperature-analysis-noaaglobaltemp-v5"
;
; Last access
;    20191101
;
; Download and processing instructions
;    Download
;    "https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v5/"+
;      "access/gridded/NOAAGlobalTemp_v5.0.0_gridded_s188001_e202112_"+
;      "c20220108T133303.nc"
;    and
;    "https://www.ncei.noaa.gov/data/noaa-global-surface-temperature/v5/"+
;       "access/timeseries/aravg.mon.land_ocean.90S.90N.v5.0.0.202112.asc"
;
; Caveats
;    Produces temperature anomalies (relative to the 1971-2000 climatology)
;
; Modification history
;    20210203-A_bock_lisa: prepared for public branch
;    20191101-A_kosaka_yu: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_noaaglobaltemp.ncl"

  ; Source name
  OBSNAME = "NOAAGlobalTemp"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1880
  YEAR2 = 2021

  ; Selected variable (standard name)
  VAR1 = "tasa"
  VAR2 = "tasaga"

  ; MIPS
  MIP = "Amon"

    ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE1 = getenv("cmor_tables") + \
    "/custom/CMOR_tasa.dat"
  CMOR_TABLE2 = getenv("cmor_tables") + \
    "/custom/CMOR_tasaga.dat"

  ; Version
  VERSION = "v5"

  ; Type
  TYPE = "ground"

  ; Global attributes
  SOURCE = "https://www.ncdc.noaa.gov/noaa-merged-land-ocean-global-surface" \
           + "-temperature-analysis-noaaglobaltemp-v5"
  REF1 =  "Zhang et al., Eos, in press, 2019"
  REF2 =  REF1
  COMMENT1 = "Temperature anomaly with respect to the period 1971-2000"
  COMMENT2 = "Global mean temperature anomaly with respect to the period " \
             + "1971-2000"

end

begin

  ; Read file
  fname = input_dir_path + "NOAAGlobalTemp_v5.0.0_gridded_s188001_e202112_" + \
                           "c20220108T133303.nc"
  f = addfile(fname, "r")

  ; Read anomaly
  time = f->time
  yyyymm = cd_calendar(time, -1)
  ts = ind(yyyymm.eq.YEAR1*100+1)
  te = ind(yyyymm.eq.YEAR2*100+12)
  anomaly = f->anom(ts:te, 0, :, :)

  ; Swap latitudes
  if (isMonotonic(anomaly&lat).eq.-1) then
    anomaly = anomaly(:, ::-1, :)
  end if

  log_info("  Anomaly range: " + min(anomaly) + \
           " K to " + max(anomaly) + " K")

  output1 = anomaly
  dims = dimsizes(output1)

  ; Format coordinates
  output1!0 = "time"
  output1!1 = "lat"
  output1!2 = "lon"
  format_coords(output1, YEAR1 + "0101", YEAR2 + "1231", FREQ)

  ; Format time coordinate
  time_save = output1&time

  ; Calculate coordinate bounds
  bounds = guess_coord_bounds(output1, FREQ)

  ; Set variable attributes
  tmp = format_variable(output1, VAR1, CMOR_TABLE1)
  delete(output1)
  output1 = tmp
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

  ; Write temperature time-series
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE, VERSION, \
               MIP, VAR1, DATESTR/), "_") + ".nc"
  write_nc(fout, VAR1, output1, bounds, gAtt1)
  w = addfile(fout, "w")
  w->height = height
  delete(w)
  delete(gAtt1)
  delete(output1)
  delete(bounds)

end

begin

  ; Read file
  len = (YEAR2 - YEAR1 + 1)*12
  fname = input_dir_path + "aravg.mon.land_ocean.90S.90N.v5.0.0.202112.asc"
  data = asciiread(fname, (/len, 10/), "float")
  output2 = data(:, 2)
  output2!0 = "time"
  output2&time = time_save
  delete(data)

  log_info("  Anomaly range: " + min(output2) + \
           " K to " + max(output2) + " K")

  output2!0 = "time"
  output2&time = time_save

  ; Format coordinates
  format_coords(output2, YEAR1 + "0101", YEAR2 + "1231", FREQ)

  ; Calculate coordinate bounds
  bounds = guess_coord_bounds(output2, FREQ)

  ; Set variable attributes
  tmp = format_variable(output2, VAR2, CMOR_TABLE2)
  delete(output2)
  output2 = tmp
  delete(tmp)

  ; Add height coordinate
  output2@coordinates = "height"
  height = 2.d
  height!0 = "ncl_scalar"
  height@units = "m"
  height@axis = "Z"
  height@positive = "up"
  height@long_name = "height"
  height@standard_name = "height"

  ; Set global attributes
  gAtt2 = set_global_atts(OBSNAME, TIER, SOURCE, REF2, COMMENT2)

  ; Write temperature anomaly time-series
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, TYPE, VERSION, \
               MIP, VAR2, DATESTR/), "_") + ".nc"
  write_nc(fout, VAR2, output2, bounds, gAtt2)
  w = addfile(fout, "w")
  w->height = height
  delete(w)
  delete(gAtt2)
  delete(bounds)
  delete(output2)

end
