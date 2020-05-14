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
;    Use the download script to retrieve the data.
;
; Modification history
;    20190314-crezee_bas: written based on esacci_soilmoisture.ncl.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cds_satellite_soil_moisture.ncl"

  ; Source name
  OBSNAME = "CDS-SATELLITE-SOIL-MOISTURE"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1978
  YEAR2 = 2020

  ; Selected variables (standard name)
  VAR = (/"sm", "sm"/)

  ; Name in the raw data
  NAME = (/"sm", "sm"/)

  ; MIP
  MIP = (/"day", "Lmon"/)

  ; Frequency
  FREQ = (/"day", "mon"/)

  ; Name of frequency in filename
  FREQNAME = (/"-DAILY-", "-MONTHLY-"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/custom/CMOR_" + VAR + ".dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "CUSTOM-TCDR-ICDR-20200602"

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
                         "C3S-SOILMOISTURE-L3S-SSMV-COMBINED" + \
                         FREQNAME(vv) + yy + \
                         "????000000-?CDR-v201?12.0.0.nc")

      ; Retrieve start and end date (YYYYMM)
      start_date = systemfunc("basename " + files(0))
      start_date := str_get_cols(str_get_field(start_date, 7, "-"), 0, 7)
      end_date = systemfunc("basename " + files(dimsizes(files) - 1))
      end_date := str_get_cols(str_get_field(end_date, 7, "-"), 0, 7)

      ; Read data
      f = addfiles(files, "r")
      delete(files)
      output = f[:]->$NAME(vv)$
      delete(f)

      ; Format coordinates
      output!0 = "time"
      output!1 = "lat"
      output!2 = "lon"
      format_coords(output, start_date, end_date, FREQ(vv))

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
      DATESTR = start_date + "-" + end_date
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
