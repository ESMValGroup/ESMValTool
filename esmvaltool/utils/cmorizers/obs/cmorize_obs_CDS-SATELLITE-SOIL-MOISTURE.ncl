; #############################################################################
; ESMValTool CMORizer for CDS-SATELLITE-SOIL-MOISTURE data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-soil-moisture?tab=form
;
; Last access
;    20190314
;
; Download and processing instructions
;    - Download the data from source as follows:
;          - Variable: Volumetric surface soil moisture
;          - Time aggregation: Day average
;          - Year: as you wish
;          - Month: as you wish
;          - Day: as you wish
;          - Format: 'tar.gz' (recommended)
;          - Type of sensor: Combined passive and active
;          - Type of record: CDR
;          - Version: v201812.0.0
;    - Decompress the files and put them in a single directory (no subdirectories with years)
;
; Modification history
;    20190314-A_crez_ba: adapted from file cmorize_obs_ESACCI-SOILMOISTURE.ncl
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_CDS-SATELLITE-SOIL-MOISTURE.ncl"

  ; Source name
  OBSNAME = "CDS-SATELLITE-SOIL-MOISTURE" ;

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1979
  YEAR2 = 2018

  ; Selected variable (standard name)
  VAR = (/"sm", "smStderr"/)

  ; Name in the raw data
  NAME = (/"sm", "sm_uncertainty"/)

  ; MIP
  MIP = (/"Lmon", "Lmon"/)

  ; Frequency
  FREQ = (/"day", "day"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/custom/CMOR_" + VAR + ".dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3S-SSMV-COMBINED-v4.2"; BAS should be changed

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/dataset/satellite-soil-moisture"
  REF = ""
  COMMENT = ""

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    do yy = YEAR1, YEAR2

      ; Set list of files
      files = systemfunc("ls " + input_dir_path + \
                         "C3S-SOILMOISTURE-L3S-SSMV-" + \
                         "COMBINED-DAILY-" + yy + "????000000-TCDR-v201812.0.0.nc")
      f = addfiles(files, "r")
      delete(files)

      ; Read data
      xx = f[:]->$NAME(vv)$
      if (isatt(xx, "scale_factor")) then
        tmp = tofloat(xx * xx@scale_factor)
        copy_VarAtts(xx, tmp)
        copy_VarCoords(xx, tmp)
        delete(xx)
        xx = tmp
        delete(tmp)
      end if
      delete(f)

      ; Append to time-series
      if (.not.isdefined("output")) then
        output = xx
      else
        output := array_append_record(output, xx, 0)
      end if
      delete(xx)

    end do

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
