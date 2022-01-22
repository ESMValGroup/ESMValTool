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
;    20201207
;
; Download and processing instructions
;    - Use the download script to retrieve the data.
;    - Use the configuration file to choose what data to process.
;    - Available years: 1992-2019 (ACTIVE) or 1979-2019 (others).
;    - Versions: v201812.0.0 and v201912.0.0
;
; Modification history
;    20201207-kazeroni-remi: fix ACTIVE version
;    20200722-crezee_bas: added PASSIVE/ACTIVE versions
;    20200602-crezee_bas: include ICDR
;    20200525-righi_mattia: fix time checker.
;    20190522-crezee_bas: time period extended.
;    20190314-crezee_bas: written based on cmorize_obs_esacci_soilmoisture.ncl.
;
; Caveats
;    - Daily data has to consist of complete months
;    - Version ACTIVE can only deal with complete years for both
;      monthly and daily data
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")
; Read in variable JOBS which defines what combination of
; [time resolution] x [sensor] will be processed.
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/cmor_config/" + \
           "CDS-SATELLITE-SOIL-MOISTURE.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_cds_satellite_soil_moisture.ncl"

  ; Source name
  OBSNAME = "CDS-SATELLITE-SOIL-MOISTURE"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1979  ; 1992
  YEAR2 = 2020  ; 2019

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

  ; Global attributes
  SOURCE = "https://cds.climate.copernicus.eu/cdsapp#!/" + \
    "dataset/satellite-soil-moisture"
  REF = "Gruber et al., Earth Syst. Sci. Data, " + \
    "doi: 10.5194/essd-11-717-2019, 2019."
  COMMENT = ""

  SENSORS_FILEBASENAMES = (/"C3S-SOILMOISTURE-L3S-SSMV-COMBINED", \
                            "C3S-SOILMOISTURE-L3S-SSMV-PASSIVE", \
                            "C3S-SOILMOISTURE-L3S-SSMS-ACTIVE"/)
  ; Version
  VERSIONS = (/"COMBINED", "PASSIVE", "ACTIVE"/)
  PATTERNS = (/"????000000-?CDR-v201?12.0.0.nc", \
               "????000000-?CDR-v201?12.0.0.nc", \
               "????000000-?CDR-v201?12.0.0.nc"/)

end

begin

  do xx = 0, dimsizes(VERSIONS) - 1

    do vv = 0, dimsizes(FREQ) - 1

      if JOBS(xx, vv) then
        log_info("Processing " + VAR(vv) + " " + VERSIONS(xx) + \
                 " (" + MIP(vv) + ") as requested.")

        do yy = YEAR1, YEAR2

          ; Set list of files
          files = systemfunc("ls " + input_dir_path + \
                             SENSORS_FILEBASENAMES(xx) + \
                             FREQNAME(vv) + yy + \
                             PATTERNS(xx))

          if all(ismissing(files)) then
            log_info("Skipping year " + yy + " since no data.")
            continue
          end if

          ; Retrieve start and end date (YYYYMM(DD))
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
            str_join((/"OBS", OBSNAME, TYPE, VERSIONS(xx), \
                       MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

          ; Write variable
          write_nc(fout, VAR(vv), output, bounds, gAtt)
          delete(gAtt)
          delete(output)
          delete(bounds)

        end do

      else
        log_info("Skipping " + VAR(vv) + " " + VERSIONS(xx) + \
                 " (" + MIP(vv) + ") as requested.")
      end if

    end do

  end do

end
