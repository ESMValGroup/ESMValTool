; #############################################################################
; ESMValTool CMORizer for CDS-XCH4 data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    https://cds.climate.copernicus.eu/cdsapp#!/dataset/
;    satellite-methane?tab=form
;
; Last access
;    20190311
;
; Download and processing instructions
;    Select Processing level "Level 3", variable "Column-average dry-air mole
;    fraction of atmospheric methane (XCH4) and related variables", Sensor and
;    algorithm "MERGED and OBS4MIPS".
;    A registration is required to download the data.
;
; Modification history
;    20190311-hassler_birgit: written.
;
; ############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cds_xch4.ncl"

  ; Source name
  OBSNAME = "CDS-XCH4"

  ; Tier
  TIER = 3

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
  CMOR_TABLE = getenv("cmor_tables") + "/custom/CMOR_xch4.dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3"

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/dataset/" + \
    "satellite-methane?tab=form"
  REF = "Buchwitz et al., Adv. Astronaut. Sci. Technol., " + \
    "doi:10.1007/s42423-018-0004-6, 2018."
  COMMENT = ""

end

begin

  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    ; Read variables
    fname = input_dir_path + "OBS_" + OBSNAME + "_sat_L3_" + NAME(vv) + \
      "_200301-201612.nc"
    setfileoption("nc", "MissingToFillValue", False)
    f = addfile(fname, "r")
    output = f->xch4

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
