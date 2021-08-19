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
;    20210817
;
; Download and processing instructions
;    Download all files from
;      https://zenodo.org/record/4535247
;    Newer versions may become available over time. 
;    Only complete years should be downloaded.
;
; Modification history
;    20210819-hassler_birgit: update for version 3.5.1 of the dataset
;    20190207-righi_mattia: renamed to NIWA-BS and adapted to v2.
;    20140528-gottschaldt_klaus-dirk: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_niwa_bs.ncl"

  ; Source name
  OBSNAME = "BS"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1979
  YEAR2 = 2019

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
  VERSION = "V3.5.1"

  ; Global attributes
  SOURCE = "https://zenodo.org/record/4535247"
  REF = "Bodeker et al. (2021), BS Filled Total Column Ozone Database V3.5.1 \
         (3.5.1) [Data set], Zenodo, https://doi.org/10.5281/zenodo.4535247; \
         Bodeker et al., Earth Sys. Sci. Data., doi:10.5194/essd-13-3885-2021, \
         2021"
  COMMENT = ""

end

begin

  files = systemfunc("ls " + input_dir_path + \
                     "BSFilledTCO_" + VERSION + \
                     "_????_Monthly.nc")

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
