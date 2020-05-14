; #############################################################################
; ESMValTool CMORizer for UWisc data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    Data provided by Ralf Bennartz.
;
; Last access
;    20150415
;
; Download and processing instructions
;    Contact Ralf Bennartz (Earth and Environmental Sciences, Vanderbilt
;    University, USA).
;
; Modification history
;    20190208-righi_mattia: adapted to v2.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")
begin

  ; Script name (for logger)
  DIAG_SCRIPT = "UWISC.ncl"

  ; Source name
  OBSNAME = "UWisc"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1988
  YEAR2 = 2007

  ; Selected variable (standard name)
  VAR = (/"lwp", "lwpStderr"/)

  ; Name in the raw data
  NAME = (/"LWP", "LWP_ERROR"/)

  ; Conversion factor
  CONVERSION = (/1.e-3, 1.e-3/)

  ; MIP
  MIP = (/"Amon", "Amon"/)

  ; Frequency
  FREQ = (/"mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/custom/CMOR_" + VAR + ".dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "v2"

  ; Global attributes
  SOURCE = "Data provided by Ralf Bennartz (Vanderbilt University, USA)"
  REF = "O'Dell et al., J. Clim., doi:10.1175/2007JCLI1958.1, 2008"
  COMMENT = ""

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    do yr = YEAR1, YEAR2

      fname = input_dir_path + "UWisc_LWPMEAN_" + yr + "_v2.nc"

      f = addfile(fname, "r")
      output = (/f->$NAME(vv)$/)

      lat = (/f->lat/)
      lon = (/f->lon/)

      ; Convert units
      output = output * CONVERSION(vv)  ; for clivi this will be equal 0

      ; Format coordinates
      output!0 = "time"
      output!1 = "lat"
      output!2 = "lon"
      output&time = create_timec(yr, yr)
      output&lat = lat
      output&lon = lon
      format_coords(output, yr + "0101", yr + "1231", FREQ(vv))

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
      DATESTR = yr + "01-" + yr + "12"
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
