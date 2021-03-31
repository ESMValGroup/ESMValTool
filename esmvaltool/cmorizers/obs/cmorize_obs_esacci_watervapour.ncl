; #############################################################################
; ESMValTool CMORizer for (preliminary) ESACCI-WATERVAPOUR data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    ftp://ftp.brockmann-consult.de
;
; Last access
;    20210331
;
; Download and processing instructions
;    Contact Marc Schröder (Marc.Schroeder@dwd.de) to get access to the data.
;    Go to
;    ftp://ftp.brockmann-consult.de/data/tcwv/dataset3_1/CDR-2/monthlies/yyyy
;    and download 0.5x0.5 deg files
;    ESACCI-WATERVAPOUR-L3S-TCWV-*-05deg-yyyymm-fv3.1.nc into the single
;    directory ${RAWOBS}/Tier2/ESACCI-WATERVAPOUR (no subdirectories with
;    years)

; Modification history
;    20210331-laue_axel: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_esacci_watervapour.ncl"

  ; Source name
  OBSNAME = "ESACCI-WATERVAPOUR"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 2003
  YEAR2 = 2017

  ; Selected variable (standard name)
  VAR = (/"prw"/)

  ; Name in the raw data
  NAME = (/"tcwv"/)

  ; Conversion factor
  CONV = (/1./)

  ; MIP
  MIP = (/"Amon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_Amon"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L3S-COMBI-TCWV-05deg-fv3.1"

  ; Global attributes
  SOURCE = "ftp://ftp.brockmann-consult.de/" \
           + "data/tcwv/dataset3_1/CDR-2/monthlies/"
  REF = "DWD"
  COMMENT = "prelimiary data"

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
        fname = systemfunc("ls " + input_dir_path + \
                           "ESACCI-WATERVAPOUR-L3S-TCWV-*-05deg-" + \
                           syear + smonth + "-fv3.1.nc")

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
;    output = where(output.eq.-999, output@_FillValue, output)

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
