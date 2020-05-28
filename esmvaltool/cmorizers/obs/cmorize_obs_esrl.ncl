; ######################################################################
;  ESMValTool CMORizer for ESRL observational data
; ######################################################################
; 
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    http://www.esrl.noaa.gov/gmd/dv/data/index.php
;
; Last access
;    20140115
;
; Download and processing instructions
;    Select the following filters:
;        Category: Greenhouse Gases.
;        Parameter name: Carbon Dioxide.
;        Type: Insitu.
;        Frequency: Monthly Averages.
;        Site: BRW, MLO, SMO, SPO.
;    Convert the data to unix format using dos2unix.
;    Get the list of stations at:
;        http://www.esrl.noaa.gov/gmd/obop/
;    This list must be manually formatted in a comma-separated table, 
;    given as
;   FILE_STAT below.
;    Check the period covered by the data and give it below as YEAR1 and
;    YEAR2.
;
; Modification history
;    20200508-zechlau_sabrina: adapted to v2.
;    20140122-righi_mattia: written.
;
; ######################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_esrl.ncl"

  ; Source name
  OBSNAME = "ESRL"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1973
  YEAR2 = 2012

  ; Selected variable (standard name)
  VAR = "co2"

  ; Field
  MIP = "Amon"

  ; Frequency
  FREQ = "yr"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_" + MIP

  ; Type
  TYPE = "ground"

  ; Stations list
  FILE_STAT = input_dir_path + OBSNAME + "_stations.csv"

  ; Global attributes
  SOURCE = "http://www.esrl.noaa.gov/gmd/dv/data/index.php"
  REF = "Dr. Pieter Tans, NOAA/ESRL " + \
        "(www.esrl.noaa.gov/gmd/ccgg/trends/) and " + \
        "Dr. Ralph Keeling, Scripps Institution of Oceanography " + \
        "(scrippsco2.ucsd.edu/)"
  COMMENT = ""

end

begin

  ; Create time coordinate
  timec = create_timec(YEAR1, YEAR2)
  datec = ut_calendar(timec, 1)

  ; Read station info
  stat = readAsciiTable(FILE_STAT, 1, "string", 0)
  st_code = str_squeeze(str_get_field(stat(:, 0), 1, ","))
  st_name = str_squeeze(str_get_field(stat(:, 0), 2, ","))
  st_lat  = tofloat(str_get_field(stat(:, 0), 3, ","))
  st_lon  = tofloat(str_get_field(stat(:, 0), 4, ","))
  st_alt  = toint(str_get_field(stat(:, 0), 5, ","))
  nstats  = dimsizes(st_code)
  delete(stat)

  ; Loop over stations
  do stID = 0, nstats - 1

    infile = input_dir_path + "/" + OBSNAME + "/" + \
             VAR + "_" + str_lower(st_code(stID)) + \
             "_surface-insitu_1_ccgg_month.txt"

    if (.not.isfilepresent(infile)) then

      head = readAsciiHead(infile, "# data_field")
      data = readAsciiTable(infile, 1, "string", dimsizes(head))
      delete(head)

      ; Read site name
      tmp = str_get_field(data(:, 0), 1, " ")
      if (.not.all(tmp(0).eq.tmp)) then
          log_info("Inconsistent site name in " + files(ii))
      else
        site = tmp(0)
      end if
      delete(tmp)

      ; Read time and extract full years only
      year = str_squeeze(str_get_field(data(:, 0), 2, " "))
      month = str_squeeze(str_get_field(data(:, 0), 3, " "))
      date = toint(year + month)
      var = tofloat(str_get_field(data(:, 0), 4, " "))
      delete(data)

      ; Check time range
      if (min(toint(year)).lt.YEAR1 .or. max(toint(year)).gt.YEAR2) then
        error_msg("f", diag_script, "", "given time range (" + \
                  YEAR1 + "-" + YEAR2 + ")" + \
                  " not consistent with data (" + \
                  min(toint(year)) + "-" + max(toint(year)) + ")")
      end if
      delete(year)
      delete(month)

      ; Create time coordinate
      timec = create_timec(YEAR1, YEAR2)
      datec = cd_calendar(timec, 1)

      ; Create output variable
      output = new(dimsizes(timec), float)
      output!0 = "time"
      output&time = timec
      output@_FillValue = FILL

      do ii = 0, dimsizes(var) - 1
        output(ind(datec.eq.date(ii))) = (/var(ii)/)
      end do
      delete(var)
      delete(date)
      output = where(output.lt.0., output@_FillValue, output)

      ; Set variable attributes (CMOR)
      tmp = format_variable(output, VAR, CMOR_TABLE)
      delete(output)
      output = tmp
      delete(tmp)

      ; Calculate coordinate bounds
      bounds = guess_coord_bounds(output, FREQ)

      ; Reformat time coordinate
      ctime = time_attrib(output&time)
      delete(output&time)
      output&time = ctime
      delete(ctime)

      ; Set global attributes
      gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

      ; Outfile
      VERSION = st_code(stID)
      DATESTR = YEAR1 + "01-" + YEAR2 + "12"

      fout = output_dir_path + \
             str_join((/"OBS", OBSNAME, TYPE, VERSION, \
             MIP, VAR, DATESTR/), "_") + ".nc"
      
      ; Write variable
      write_nc(fout, VAR, output, bounds, gAtt)
      delete(gAtt)
      delete(output)

    end if

  end do

end
