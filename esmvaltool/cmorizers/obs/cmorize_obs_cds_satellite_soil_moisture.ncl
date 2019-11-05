; #############################################################################
; ESMValTool CMORizer for CDS-SATELLITE-SOIL-MOISTURE data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    https://cds.climate.copernicus.eu/cdsapp#!/dataset/
;    satellite-soil-moisture?tab=form
;
; Last access
;    20190314
;
; Download and processing instructions
;    Download the data from source as follows:
;      Variable: Volumetric surface soil moisture
;      Time aggregation: Day average
;      All Years, Month and Days
;      Format: 'tar.gz' (recommended)
;      Type of sensor: Combined passive and active
;      Type of record: CDR
;      Version: v201812.0.0
;    Decompress the files and put them in a single directory (no subdirectories
;    with years)
;
; Modification history
;    20190314-crezee_bas: written based on cmorize_obs_esacci_soilmoisture.ncl.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_cds_satellite_soil_moisture.ncl"

  ; Source name
  OBSNAME = "CDS-SATELLITE-SOIL-MOISTURE"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1979
  YEAR2 = 2018

  ; Selected variable (standard name)
  VAR = (/"sm", "smStderr"/)

  ; Name in the raw data
  NAME = (/"sm", "sm_uncertainty"/)

  ; MIP
  MIP = (/"day", "day"/)

  ; Frequency
  FREQ = (/"day", "day"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/custom/CMOR_" + VAR + ".dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "COMBINED-TCDR-v201812.0.0"

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/" + \
    "dataset/satellite-soil-moisture"
  REF = "Gruber et al., Earth Syst. Sci. Data, " + \
    "doi: 10.5194/essd-11-717-2019, 2019."
  COMMENT = ""

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    do yy = YEAR1, YEAR2

      ; Set list of files
      files = systemfunc("ls " + input_dir_path + \
                         "C3S-SOILMOISTURE-L3S-SSMV-" + \
                         "COMBINED-DAILY-" + yy + "????000000-" + \
                         "TCDR-v201812.0.0.nc")
      f = addfiles(files, "r")
      delete(files)

      ; Read data
      output = f[:]->$NAME(vv)$
      delete(f)

      ; Format coordinates
      output!0 = "time"
      output!1 = "lat"
      output!2 = "lon"
      format_coords(output, yy + "0101", yy + "1231", FREQ(vv))

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
      DATESTR = yy + "01-" + yy + "12"
      fout = output_dir_path + \
        str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                   MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

      ; Write variable
      write_nc(fout, VAR(vv), output, bounds, gAtt)
      delete(gAtt)
      delete(output)
      delete(bounds)

    end do

  end do

end
