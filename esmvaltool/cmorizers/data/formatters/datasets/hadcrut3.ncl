; #############################################################################
; ESMValTool CMORizer for HadCRUT3 data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    http://www.metoffice.gov.uk/hadobs/hadcrut3/data/download.html
;
; Last access
;    20190221
;
; Download and processing instructions
;    Download the HadCRUT3v.nc file (variance adjusted dataset).
;
; Caveats
;    The HadCRUT3v variance-adjusted dataset for tas is actually the anomaly
;    with respect to the period 1958-2001.
;
; Modification history
;    20190221-righi_mattia: adapted to v2 and renamed to HadCRUT3.
;    20150330-righi_mattia: updated paths and global attributes.
;    20140311-senftleben_daniel: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "hadcrut3.ncl"

  ; Source name
  OBSNAME = "HadCRUT3"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = get_year(start_year, 1850)
  YEAR2 = get_year(end_year, 2013)

  ; Selected variable (standard name)
  VAR = "tasa"

  ; Name in the raw data
  NAME = "temp"

  ; MIP
  MIP = "Amon"

  ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/custom/CMOR_tasa.dat"

  ; Type
  TYPE = "ground"

  ; Version
  VERSION = "1"

  ; Global attributes
  SOURCE = "http://www.metoffice.gov.uk/hadobs/hadcrut3/data/download.html"
  REF = "Brohan et al., J. Geophys. Res., doi:10.1029/2005JD006548, 2006"
  COMMENT = "Temperature anomaly with respect to the period 1958-2001"

end

begin

  ; Read file
  fname = input_dir_path + "HadCRUT3v.nc"
  f = addfile(fname, "r")

  ; Read variable
  output = f->temp

  ; Delete level coordinate (dimension 1)
  tmp = rm_single_dims(output)
  delete(output)
  output = tmp
  delete(tmp)

  ; Extract time period
  output!0 = "time"
  date = cd_calendar(output&time, 0)
  idx = ind(date(:, 0).ge.YEAR1 .and. date(:, 0).le.YEAR2)
  output := output(idx, :, :)

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
  w = addfile(fout, "w")
  delete(w)
  delete(gAtt)
  delete(output)
  delete(bounds)

end
