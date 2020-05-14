; #############################################################################
; ESMValTool CMORizer for GHCN data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://www.esrl.noaa.gov/psd/data/gridded/data.ghcngridded.html
;
; Last access
;    20190308
;
; Download and processing instructions
;    Download the dataset "precip.mon.total.nc" (precipitation, total, surface,
;    1900-2015 on a 5x5 grid).
;
; Modification history
;    20190308-righi_mattia: minor changes to include coordinate boundaries.
;    20190227-bock_lisa: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "ghcn.ncl"

  ; Source name
  OBSNAME = "GHCN"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1900
  YEAR2 = 2014

  ; Selected variable (standard name)
  VAR = "pr"

  ; Name in the raw data
  NAME = "precip"

  ; MIP
  MIP = "Amon"

  ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_Amon"

  ; Type
  TYPE = "ground"

  ; Version
  VERSION = "1"

  ; Global attributes
  SOURCE = "https://www.esrl.noaa.gov/psd/data/gridded/data.ghcngridded.html"
  REF = "Jones and Moberg, J. Clim., " + \
    "doi:10.1175/1520-0442(2003)016<0206:HALSSA>2.0.CO;2, 2003."
  COMMENT = ""

end

begin

  ; Read file
  fname = input_dir_path + "precip.mon.total.nc"
  f = addfile(fname, "r")
  setfileoption("nc", "MissingToFillValue", False)

  ; Read absolute precipitation without last incomplete year
  output = f->$NAME$(time|0:1379, lat|:, lon|:)

  ; Calculate days per month
  date = cd_calendar(output&time, 0)
  dpm = days_in_month(toint(date(:, 0)), toint(date(:, 1)))
  dpmc = conform(output, dpm, 0)

  ; Check time range
  if (dimsizes(date(:, 0)).ne.12 * (YEAR2 - YEAR1 + 1)) then
    error_msg("f", DIAG_SCRIPT, "", "incorrect number of timesteps")
  end if

  ; Convert units [mm/month] --> [kg/m2/s]
  output = output / (24 * 3600 * dpmc)

  log_info("  Climatology range: " + min(output) + \
           " kg/m2/s to " + max(output) + " kg/m2/s")

  ; Format coordinates
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

  ; Write temperature time-series
  write_nc(fout, VAR, output, bounds, gAtt)
  delete(gAtt)
  delete(output)
  delete(bounds)

end
