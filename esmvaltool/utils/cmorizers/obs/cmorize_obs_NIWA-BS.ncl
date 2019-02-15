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
;    20190207-A_righ_ma: renamed to NIWA-BS and adapted to v2.
;    20140528-A_gott_kl: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Source name
  OBSNAME = "NIWA-BS"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1979
  YEAR2 = 2016

  ; Selected variable (standard name)
  VARS = (/"toz", "tozStderr"/)

  ; Name in the raw data
  NAME = (/"tco", "tco_uncert"/)

  ; MIPS
  MIPS = "Amon"

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + "/cmor/tables/custom/CMOR_" + \
    VARS + ".dat"

  ; Version
  VERSION = "3.3"

end

begin

  files = systemfunc("ls " + input_dir_path + \
                     "NIWA-BS_CombinedTCO_V" + VERSION + \
                     "_????_Monthly_Patched.nc")

  do vv = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vv))

    f = addfiles(files, "r")
    output = f[:]->$NAME(vv)$

    ; Format time coordinate
    ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ; Format latitude coordinate
    output!1 = "lat"
    if (isMonotonic(output&lat) .eq. 0) then
      error_msg("f", diag_script, "", "non-monotonic latitude coordinate")
    end if
    if (isMonotonic(output&lat) .eq. -1) then
      output = output(:, :, ::-1, :)
    end if
    clat = lat_attrib(output&lat)
    delete(output&lat)
    output&lat = clat
    delete(clat)

    ; Format longitude coordinate
    output!2 = "lon"
    if (isMonotonic(output&lon) .eq. 0) then
      error_msg("f", diag_script, "", "non-monotonic longitude coordinate")
    end if
    if (any(output&lon.lt.0.)) then
      output = lonFlip(output)
    end if
    clon = lon_attrib(output&lon)
    delete(output&lon)
    output&lon = clon
    delete(clon)

    ; Set variable attributes
    tmp = var_attrib(output, VARS(vv), CMOR_TABLE(vv))
    delete(output)
    output = tmp
    delete(tmp)

    ; Set global attributes
    gAtt = set_global_atts( \
      OBSNAME, TIER, \
      "http://www.bodekerscientific.com/data/total-column-ozone", \
      "Bodeker et al., Atmos. Chem. Phys., " + \
      "doi:10.5194/acp-5-2603-2005, 2005", "")

    ; Outfile
    fout = output_dir_path + "OBS_" + OBSNAME + "_sat_v" + VERSION + "_" + \
      MIPS + "_" + VARS(vv) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ; Write variable
    write_nc(fout, VARS(vv), output, gAtt)
    delete(gAtt)
    delete(output)

  end do

end
