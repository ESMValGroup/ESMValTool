; #############################################################################
; ESMValTool CMORizer for CDS-SICONC data
; #############################################################################
;
; Tier
;    Tier 2 - other freely-available dataset.
;
; Source
;    https://cds.climate.copernicus.eu/cdsapp#!/dataset/
;    satellite-sea-ice?tab=form
;
; Last access
;    20190311
;
; Download and processing instructions
;    Select Processing level "Level 2", variable "sea_ice_area_fraction", Sensor and
;    algorithm - please refer to documentation.
;    A registration is required to download the data.
;
; Modification history
;    20190320-rana_ar: written.
;
; ############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_CDS-SICONC.ncl"

  ; Source name
  OBSNAME = "CDS-SICONC"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1979
  YEAR2 = 2015

  ; Selected variable (standard name)
  VAR = (/"siconc"/)

  ; Name in the raw data
  NAME = (/"siconc"/)

  ; MIP
  MIP = (/"SImon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/custom/CMOR_siconc.dat"

  ; Type
  TYPE = "sat"

  ; Hemisphere -  this needs to be changed for southern hemisphere
  HEMI = NH

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/dataset/" + \
    "satellite-sea-ice?tab=form"
  REF = ""
  COMMENT = ""

end

begin

  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    ; Read variables
    fname = input_dir_path + "OBS_" + OBSNAME + "_sat_${HEMI}_" + NAME(vv) +\
      "_${YEAR1}01-${YEAR2}12.nc"
    setfileoption("nc", "MissingToFillValue", False)
    f = addfile(fname, "r")
    output = f->siconc

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
