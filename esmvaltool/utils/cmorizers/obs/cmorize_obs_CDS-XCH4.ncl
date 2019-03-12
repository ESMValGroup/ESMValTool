; #############################################################################
; ESMValTool CMORizer for HadISST data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-methane?tab=overview
;
; Last access
;    20190311
;
; Download and processing instructions
;    Download data file directly as nc file from the CDS
;
; Modification history
;    20190311-A_hass_bg: written.
;
; ############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_CDS-XCH4.ncl"

  ; Source name
  OBSNAME = "CDS-XCH4"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 2003
  YEAR2 = 2016

  ; Selected variable (standard name)
  VAR = (/"xch4"/)

  ; Name in the raw data
  NAME = (/"xch4"/)

  ; MIP
  MIP = (/"Amon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/custom/CMOR_xch4.dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3"

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-methane?tab=overview"
  REF = ""
  COMMENT = ""

end

begin

  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    ; Read variables
    fname = input_dir_path + "OBS_" + OBSNAME + "_sat_L3_" + NAME(vv) + "_200301-201612.nc"
    setfileoption("nc", "MissingToFillValue", False)
    f = addfile(fname, "r")
    

    output = f->xch4

    ; Extract time period
    date = cd_calendar(output&time, 0)
    idx = ind(date(:, 0).ge.YEAR1 .and. date(:, 0).le.YEAR2)
    output := output(idx, :, :)

    ;; fix the longitude
    ;if (any(OUTPUT&lon.lt.0.)) then
    ;  OUTPUT = lonFlip(OUTPUT)
    ;end if
    ;clon = lon_attrib(OUTPUT&lon)
    ;delete(OUTPUT&lon)
    ;OUTPUT&lon = clon
    ;delete(clon)
    
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
