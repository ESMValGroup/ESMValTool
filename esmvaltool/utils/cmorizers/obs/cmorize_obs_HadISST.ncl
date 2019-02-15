; #############################################################################
; REFORMAT SCRIPT FOR HadISST REANALYSIS DATA
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    http://www.metoffice.gov.uk/hadobs/hadisst/data/download.html
;
; Last access
;    02/2019
;
; Download and processing instructions
;    Download and unzip "HadISST_ice.nc.gz" and "HadISST_sst.nc.gz".
;
; Modification history
;    20190208-A_hass_bg: adapted to v2.
;    20180530-A_Righ_ma: fixed coordinates and metadata
;    20170217-A_senf_da: modified so that SST fields are also written as 'tos'
;    20150422-A_laue_ax: written
;
; ############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_HadISST.ncl"

  ; Source name
  OBSNAME = "HadISST"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1870
  YEAR2 = 2017

  ; MIPS
  MIPS = (/"Amon", "Omon", "OImon"/)

  ; Selected variable (standard name)
  VARS = (/"ts", "tos", "sic"/)

  ; Name in the raw data
  NAME = (/"sst", "sst", "ice"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/cmip5/Tables/CMIP5_" + MIPS

end

begin

  ; Loop over variables
  do vID = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vID))

    fname = input_dir_path + "HadISST_" + NAME(vID) + ".nc"
    f = addfile(fname, "r")
    if (NAME(vID).eq."sst") then
      output = f->sst
      output@_FillValue = -1000.
    end if
    if (NAME(vID).eq."ice") then
      output = f->sic
    end if

    if (isStrSubset(NAME(vID), "sst")) then
      output = output + 273.15  ; convert from degrees C to K
      output@units = "K"
    end if
    if (isStrSubset(NAME(vID), "ice")) then
      output = output * 100.  ; convert from 1 to %
      output@units = "%"
    end if

    date = cd_calendar(output&time, 0)
    idx = ind(date(:, 0).ge.YEAR1 .and. date(:, 0).le.YEAR2)
    output := output(idx, :, :)

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
      output = output(:, ::-1, :)
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
    tmp = var_attrib(output, VARS(vID), CMOR_TABLE(vID))
    delete(output)
    output = tmp
    delete(tmp)

    ; Set global attributes
    gAtt = set_global_atts( \
      OBSNAME, TIER, \
      "http://www.metoffice.gov.uk/hadobs/hadisst/data/download.html", \
      "Rayner et al., J. Geophys. Res., doi:10.1029/2002JD002670, 2013", "")

    ; Outfile
    fout = output_dir_path + "OBS_" + OBSNAME + "_reanaly_1_" + MIPS(vID) + \
      "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ; Write variable
    write_nc(fout, VARS(vID), output, gAtt)
    delete(gAtt)
    delete(output)

  end do

end
