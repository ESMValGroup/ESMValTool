; #############################################################################
; ESMValTool CMORizer for ESACCI-WV data (WATER VAPOR)
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    https://wui.cmsaf.eu/safira/action/viewDoiDetails?acronym=HOAPS_V002
;
; Last access
;    20190201
;
; Download and processing instructions
;    Download the data (see source) and put all files into the same directory.
;
; Modification history
;    20190508-A_laue_ax: written (base on reformat_obs_ESACCI-CLOUD.ncl).
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ESACCI-WV.ncl"

  ; Source name
  OBSNAME = "ESACCI-WV"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1988
  YEAR2 = 2014

  ; Selected variable (standard name)
  VAR = (/"prw"/)

  ; Name in the raw data
  NAME = (/"wvpa"/)

  ; Conversion factor
  CONV = (/1./)

  ; MIP
  MIP = (/"Amon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + "/cmor/tables/" + \
    (/"cmip5/Tables/CMIP5_Amon"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "HOAPS-v4.0"

  ; Global attributes
  SOURCE = "https://wui.cmsaf.eu/safira/action/viewDoiDetails?" + \
           "acronym=HOAPS_V002"
  REF = "Andersson et al., DOI:10.5676/EUM_SAF_CM/HOAPS/V002, https://" + \
        "doi.org/10.5676/EUM_SAF_CM/HOAPS/V002, 2017"
  COMMENT = ""

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
        fname = systemfunc("ls " + input_dir_path + "HTWmm" + syear + \
                           smonth + "01000000213SCPOS01GL.nc")

        ; No files found
        if (all(ismissing(fname))) then
          continue
        end if

        ; Extract data
        f = addfile(fname, "r")
        xx = f->$NAME(vv)$

        ; Convert units
        xx = xx * CONV(vv)

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
    output = where(output.eq.-999, output@_FillValue, output)

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
