; #############################################################################
; ESMValTool CMORizer for CDS-XCH4 data
; #############################################################################
;
; Tier
;    Tier 3 other freely-available dataset.
;
; Source
;    https://cds.climate.copernicus.eu/cdsapp#!/dataset/
;    satellite-methane?tab=form
;
; Last access
;    20190311
;
; Download and processing instructions
;    A registration is required to download the data.
;
; Modification history
;    20190318-A_crez_ba: written.
;
; ############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_CDS-SATELLITE-ALBEDO.ncl"

  ; Source name
  OBSNAME = "CDS-SATELLITE-ALBEDO-025deg"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 2003
  YEAR2 = 2003

  ; Selected variable (standard name)
  VAR = (/"bhalb"/)

  ; Name in the raw data
  NAME = (/"AL_BH_BB"/)
  ; Name in the raw data filename
  NAME_FILENAME = (/"albh"/)


  ; MIP
  MIP = (/"Lmon"/)

  ; Frequency ; BAS: actually, 10-daily, let's see how to deal with that. 
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/custom/CMOR_bhalb.dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3"

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/dataset/" + \
    "satellite-albedo?tab=form"
  REF = ""
  COMMENT = ""

end

begin

  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    ; Read variables
    fname = input_dir_path + "OBS_" + OBSNAME + "_sat_L3_10day_" + NAME_FILENAME(vv) + \
      "_"+YEAR1+"01-"+YEAR2+"12.nc"
    setfileoption("nc", "MissingToFillValue", False)
    f = addfile(fname, "r")
    output = f->AL_BH_BB

    ; Format coordinates
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
