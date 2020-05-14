; #############################################################################
; ESMValTool CMORizer for MAC-LWP data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    https://search.earthdata.nasa.gov/search?q=MAC-LWP
;
; Last access
;    20200130
;
; Download and processing instructions
;    Select the MACLWP-mean dataset.
;    Download the script file for both projects.
;    Download the data using http either by selecting each granule
;    individually or by using the option "download all".
;    Data is freely available, but a registration is required.
;
; Modification history
;    20200130-hassler_birgit: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")
begin

  ; Script name (for logger)
  DIAG_SCRIPT = "MAC-LWP.ncl"

  ; Source name
  OBSNAME = "MAC-LWP"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1988
  YEAR2 = 2016

  ; Selected variable (standard name)
  VAR = (/"lwp", "lwpStderr"/)

  ; Name in the raw data
  NAME = (/"cloudlwp", "cloudlwp_error"/)

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
  VERSION = "v1"

  ; Global attributes
  SOURCE = "https://search.earthdata.nasa.gov/search?q=MAC-LWP"
  REF = "Elsaesser et al., J. Clim., doi:10.1175/JCLI-D-16-0902.1, 2017"
  COMMENT = ""

end

begin

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    do yr = YEAR1, YEAR2

      fname = input_dir_path + "maclwp_cloudlwpave_" + yr + "_v1.nc4"

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
