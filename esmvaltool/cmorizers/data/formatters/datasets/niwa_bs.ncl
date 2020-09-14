; #############################################################################
; ESMValTool CMORizer for NIWA-BS data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    http://www.bodekerscientific.com/data/total-column-ozone
;
; Last access
;    20190207
;
; Download and processing instructions
;    To get the access data send an email to datasets@bodekerscientific.com
;    Download all files from
;      ftp://ftp.bodekerscientific.com/CombinedTCOV3.3/Monthly/Patched/NetCDF/
;    Newer versions may become available over time, but make sure to download
;    the patched one. Only complete years should be downloaded.
;
; Modification history
;    20190207-righi_mattia: renamed to NIWA-BS and adapted to v2.
;    20140528-gottschaldt_klaus-dirk: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "niwa_bs.ncl"

  ; Source name
  OBSNAME = "NIWA-BS"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = get_year(start_year, 1979)
  YEAR2 = get_year(end_year, 2016)

  ; Selected variable (standard name)
  VAR = (/"toz", "tozStderr"/)

  ; Name in the raw data
  NAME = (/"tco", "tco_uncert"/)

  ; MIP
  MIP = (/"Amon", "Amon"/)

  ; Frequency
  FREQ = (/"mon", "mon"/)

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + "/custom/CMOR_" + VAR + ".dat"

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "V3.3"

  ; Global attributes
  SOURCE = "http://www.bodekerscientific.com/data/total-column-ozone"
  REF = "Bodeker et al., Atmos. Chem. Phys., doi:10.5194/acp-5-2603-2005, 2005"
  COMMENT = ""

end

begin

  files = systemfunc("ls " + input_dir_path + \
                     "NIWA-BS_CombinedTCO_" + VERSION + \
                     "_????_Monthly_Patched.nc")

  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    f = addfiles(files, "r")
    output = f[:]->$NAME(vv)$

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
      str_join((/"OBS", OBSNAME, TYPE, str_sub_str(VERSION, "V", "v"), \
                 MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

    ; Write variable
    write_nc(fout, VAR(vv), output, bounds, gAtt)
    delete(gAtt)
    delete(output)
    delete(bounds)

  end do

end
