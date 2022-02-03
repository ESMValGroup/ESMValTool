; #############################################################################
; ESMValTool CMORizer for Kadow data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    "https://github.com/FREVA-CLINT/climatereconstructionAI/tree/master/"+
;      "reconstructions"
;
; Last access
;    20210212
;
; Download and processing instructions
;    Download "http://users.met.fu-berlin.de/~ChristopherKadow/"+
;      "HadCRUT.5.0.1.0.anomalies.Kadow_et_al_2020_20crAI-infilled."+
;      "ensemble_mean_185001-202012.nc"
;
; Modification history
;    20220203-bock_lisa: fixed for the public branch
;    20210212-kosaka_yu: updated to version 5.0.1.0
;    20201105-kosaka_yu: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_kadow.ncl"

  ; Source name
  OBSNAME = "Kadow"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1850
  YEAR2 = 2020

  ; Selected variable (standard name)
  VAR1 = "tasa"

  ; MIP
  MIP = "Amon"

  ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE1 = getenv("cmor_tables") + \
    "/custom/CMOR_tasa.dat"

  ; Version
  VERSION = "5.0.1.0"

  ; Type
  TYPE1 = "ground"

  ; Global attributes
  SOURCE = "http://users.met.fu-berlin.de/~ChristopherKadow/HadCRUT." + \
           "5.0.1.0.anomalies.Kadow_et_al_2020_20crAI-infilled.ensemble" + \
           "_mean_185001-202012.nc"
  REF1 =  "Kadow et al., Nature Geosci, doi:10.1038/s41561-020-0582-5, 2020"
  COMMENT1 = "Temperature anomaly with respect to the period 1961-1990, " + \
             "ensemble mean, infilled by AI"

end

begin

  ; Read file
  fname = input_dir_path + "HadCRUT.5.0.1.0.anomalies.Kadow_et_al_2020_" \
                         + "20crAI-infilled.ensemble_mean_185001-202012.nc"
  f = addfile(fname, "r")

  ; Read anomaly
  len = (YEAR2 - YEAR1 + 1)*12
  anomaly = f->tas(0:len-1, :, :)

  ; Swap latitudes
  if (isMonotonic(anomaly&lat).eq.-1) then
    anomaly = anomaly(:, ::-1, :)
  end if

  log_info("  Anomaly range: " + min(anomaly) + \
           " K to " + max(anomaly) + " K")

  output = anomaly
  dims = dimsizes(output)
  delete(anomaly)

  ; Format coordinates
  output!0 = "time"
  output!1 = "lat"
  output!2 = "lon"
  format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ)
  time_save = output&time

  ; Calculate coordinate bounds
  bounds = guess_coord_bounds(output, FREQ)

  ; Set variable attributes
  tmp = format_variable(output, VAR1, CMOR_TABLE1)
  delete(output)
  output = tmp
  delete(tmp)

  ; Add height coordinate
  output@coordinates = "height"
  height = 2.d
  height!0 = "ncl_scalar"
  height@units = "m"
  height@axis = "Z"
  height@positive = "up"
  height@long_name = "height"
  height@standard_name = "height"

  ; Set global attributes
  gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF1, COMMENT1)

  ; Write temperature anomaly time-series
  DATESTR = YEAR1 + "01-" + YEAR2 + "12"
  fout = output_dir_path + \
    str_join((/"OBS", OBSNAME, "mean", VERSION, \
               MIP, VAR1, DATESTR/), "_") + ".nc"
  write_nc(fout, VAR1, output, bounds, gAtt)
  w = addfile(fout, "w")
  w->height = height
  delete(w)
  delete(gAtt)
  delete(bounds)
  delete(output)

end
