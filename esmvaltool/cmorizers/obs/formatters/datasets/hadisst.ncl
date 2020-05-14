; #############################################################################
; ESMValTool CMORizer for HadISST data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    http://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
;
; Last access
;    20190208
;
; Download and processing instructions
;    Download and unzip "HadISST_ice.nc.gz" and "HadISST_sst.nc.gz".
;
; Modification history
;    20190208-hassler_birgit: adapted to v2.
;    20180530-righi_mattia: fixed coordinates and metadata.
;    20170217-senftleben_daniel: modified to write SST fields as 'tos'.
;    20150422-lauer_axel: written.
;
; ############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "hadisst.ncl"

  ; Source name
  OBSNAME = "HadISST"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1870
  YEAR2 = 2017

  ; Selected variable (standard name)
  VAR = (/"ts", "tos", "sic"/)

  ; Name in the raw data
  NAME = (/"sst", "sst", "ice"/)

  ; MIP
  MIP = (/"Amon", "Omon", "OImon"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_" + MIP

  ; Type
  TYPE = "reanaly"

  ; Version
  VERSION = "1"

  ; Global attributes
  SOURCE = "http://www.metoffice.gov.uk/hadobs/hadisst/data/download.html"
  REF = "Rayner et al., J. Geophys. Res., doi:10.1029/2002JD002670, 2013"
  COMMENT = ""

end

begin

  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    ; Read variables
    fname = input_dir_path + "HadISST_" + NAME(vv) + ".nc"
    f = addfile(fname, "r")
    if (NAME(vv).eq."sst") then
      output = f->sst
      output@_FillValue = -1000.
    end if
    if (NAME(vv).eq."ice") then
      output = f->sic
    end if

    ; Convert units
    if (isStrSubset(NAME(vv), "sst")) then
      output = output + 273.15  ; [degC] --> [K]
      output@units = "K"
    end if
    if (isStrSubset(NAME(vv), "ice")) then
      output = output * 100.  ; [1] --> [%]
      output@units = "%"
    end if

    ; Extract time period
    date = cd_calendar(output&time, 0)
    idx = ind(date(:, 0).ge.YEAR1 .and. date(:, 0).le.YEAR2)
    output := output(idx, :, :)

    ; Format coordinates
    output!0 = "time"
    output!1 = "lat"
    output!2 = "lon"
    format_coords(output, YEAR1 + "0101", YEAR2 + "1231", FREQ(vv))

    ; Set variable attributes
    tmp = format_variable(output, VAR(vv), CMOR_TABLE(vv))
    delete(output)
    output = tmp
    delete(tmp)

    ; Calculate coordinate bounds
    bounds = guess_coord_bounds(output, FREQ(vv))

    ; Set global attributes
    gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

    ; Output file
    DATESTR = YEAR1 + "01-" + YEAR2 + "12"
    fout = output_dir_path + \
      str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                 MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

    ; Write variable
    write_nc(fout, VAR(vv), output, bounds, gAtt)
    delete(gAtt)
    delete(output)
    delete(bounds)

  end do

end
