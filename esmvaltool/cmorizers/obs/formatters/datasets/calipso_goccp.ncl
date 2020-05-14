; #############################################################################
; ESMValTool CMORizer for CALIPSO-GOCCP data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    ftp://ftp.climserv.ipsl.polytechnique.fr/
;
; Last access
;    20200127
;
; Download and processing instructions
;    Download the data from:
;      ftp://ftp.climserv.ipsl.polytechnique.fr/
;      path: /cfmip/GOCCP_v3/3D_CloudFraction/grid_2x2xL40/yyyy/avg/
;    Put all files under a single directory (no subdirectories with years).
;      3D_CloudFraction330m_yyyymm_avg_CFMIP2_sat_3.1.2.nc
;
; Modification history
;    20200127-laue_axel: written
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "calipso-goccp.ncl"

  ; Source name
  OBSNAME = "CALIPSO-GOCCP"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 2007
  YEAR2 = 2015

  ; Selected variable (standard name)
  VAR = (/"clcalipso"/)

  ; Names in the raw data
  NAME = (/"clcalipso"/)

  ; Conversion factor
  CONV = (/100./)

  ; MIP
  MIP = (/"cfMon"/)

  ; Frequency
  FREQ = (/"mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_cfMon"/)

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "3.1.2"

  ; Global attributes
  SOURCE = "ftp://ftp.climserv.ipsl.polytechnique.fr/" + \
    "/cfmip/GOCCP_v3/3D_CloudFraction/grid_2x2xL40/yyyy/avg/"
  REF = "Chepfer et al., J. Geophys. Res., doi:10.1029/2009JD012251, 2010"
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
        fname = input_dir_path + "3D_CloudFraction330m_" + syear + smonth + \
          "_avg_CFMIP2_" + TYPE + "_" + VERSION + ".nc"

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
          output!1 = "alt40"
          ; km --> m and round to 1 m
          output&alt40 = round(f->alt_mid * 1000.0, 0)
          output!2 = "lat"
          output&lat = f->latitude
          output!3 = "lon"
          output&lon = f->longitude
        end if
        output(ind(toint(yy * 100 + mm).eq.date), :, :, :) = (/xx/)
        delete(fname)
        delete(f)
      end do
    end do

    ; Set fill value
    output = where(output.eq.-999, output@_FillValue, output)

    ; Format coordinates
    output!0 = "time"
    output!1 = "alt40"
    output!2 = "lat"
    output!3 = "lon"
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
