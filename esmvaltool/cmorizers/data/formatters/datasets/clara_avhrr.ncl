; #############################################################################
; ESMValTool CMORizer for CM SAF CLARA-AHRR v2 data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    https://wui.cmsaf.eu/
;
; Last access
;    2021-03-22
;
; Download and processing instructions
;    1) Create ("register") an user account at
;       https://wui.cmsaf.eu/safira/action/viewLogin?menuName=NUTZER_HOME
;    2) login (same URL as above)
;    3) Search data using search form at
;    https://wui.cmsaf.eu/safira/action/viewProduktHome?menuName=PRODUKT_HOME
;
;      - Product group: Climate Data Records
;      - Product family: CLARA-A ed. 2.1
;      - Product name: CFC - Factional cloud cover
;                      IWP - Ice water path
;                      LWP - Liquid water path
;      - Area: Global
;      - Temporal resolution: Monthly
;
;    4) Select "CLARA-A ed. 2.1 AVHRR on polar orbiting satellites" from
;       list of results.
;    5) Click on "Add to order cart"
;    6) Follow download instructions in automatic email received when data
;       are ready for download.
;    7) Untar all .tar files into a single directory.
;
; Modification history
;    20210506-lauer_axel: output of lwp instead of clwvi
;    20210323-lauer_axel: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/data/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_clara_avhrr.ncl"

  ; Source name
  OBSNAME = "CLARA-AVHRR"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1982
  YEAR2 = 2018

  ; Selected variable (standard name)
  VAR = (/"clt", "clivi", "lwp"/)

  ; Name in the raw data
  NAME = (/"cfc", "iwp_allsky", "lwp_allsky"/)

  ; Filename base
  FNBASE = (/"CFCmm", "IWPmm", "LWPmm"/)

  ; Conversion factor
  ; Remark: total cloud cover (CFC) is reported as "1" but is actually "%"
  ; IWP and LWP use scale_factor to convert to kg/m2
  ; CONV = (/1., 1., 1./)

  ; MIP
  MIP = (/"Amon", "Amon", "Amon"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/custom/CMOR_lwp.dat"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "V002_01"

  ; Global attributes
  SOURCE = "https://wui.cmsaf.eu/"
  REF = "https://doi.org/10.5676/EUM_SAF_CM/CLARA_AVHRR/V002_01"
  COMMENT = "The CM SAF data are owned by EUMETSAT and are available to " \
    + "all users free of charge and with no conditions to use. If you wish " \
    + "to use these products, EUMETSAT's copyright credit must be shown by " \
    + "displaying the words 'Copyright (c) (2020) EUMETSAT' under/in each " \
    + "of these SAF Products used in a project or shown in a publication " \
    + "or website. Please follow the citation guidelines given at " \
    + "https://doi.org/10.5676/EUM_SAF_CM/CLARA_AVHRR/V002_01 and also " \
    + "register as a user at http://cm-saf.eumetsat.int/ to receive latest " \
    + "information on CM SAF services and to get access to the CM SAF User " \
    + "Help Desk."

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    time = create_timec(YEAR1, YEAR2)
    date = cd_calendar(time, 1)

    ; Create timeseries
    do yy = YEAR1, YEAR2

      syear = sprinti("%i", yy)
      do mm = 1, 12

        smonth = sprinti("%0.2i", mm)

        ; Read file
        fname = systemfunc("ls " + input_dir_path + FNBASE(vv) + \
                           syear + smonth + "01*.nc")

        ; No files found
        if (ismissing(fname)) then
          log_info("Warning: no input data found for variable " + VAR(vv) + \
                   " (" + syear + smonth + ")")
          continue
        end if

        ; Extract data
        f = addfile(fname, "r")
        val = f->$NAME(vv)$
        if (isatt(val, "scale_factor")) then
          scalefac = tofloat(val@scale_factor)
        else
          scalefac = 1.0
        end if
        if (isatt(val, "add_offset")) then
          offset = tofloat(val@add_offset)
        else
          offset = 0.0
        end if
        xx = tofloat(val) * scalefac + offset
        delete(val)

        ; Assign to global array
        if (.not.isdefined("output")) then
          dims = dimsizes(xx)
          dims(0) = dimsizes(time)
          output = new(dims, float)
          output!0 = "time"
          output&time = time
          output!1 = "lat"
          output&lat = f->lat
          output!2 = "lon"
          output&lon = f->lon
        end if
        output(ind(toint(yy * 100 + mm).eq.date), :, :) = (/xx/)

        delete(fname)
        delete(f)

      end do
    end do

    ; Set fill value
    output = where(output.eq.xx@_FillValue, output@_FillValue, output)

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
