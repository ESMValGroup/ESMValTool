; #############################################################################
; ESMValTool CMORizer for ESACCI-OZONE data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/ozone/data/
;
; Last access
;    20190201
;
; Download and processing instructions
;    Download the data from:
;      total_columns/l3/merged/v0100/
;    Put all files under a single directory (no subdirectories with years).
;
; Modification history
;    20190201-A_righ_ma: adapted to v2 and replace NaN/inf with FillValue.
;    20160224-A_wenz_sa: written based on reformat_obs_ESACCI-AEROSOL.ncl.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ESACCI-OZONE.ncl"

  ; Source name
  OBSNAME = "ESACCI-OZONE"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = (/1997, 1997, 2007, 2007/)
  YEAR2 = (/2010, 2010, 2008, 2008/)

  ; MIP
  MIP = "Amon"

  ; Selected variable (standard name)
  VARS = (/"toz", "tozStderr", "tro3prof", "tro3profStderr"/)

  ; Name in the raw data
  NAME = (/"atmosphere_mole_content_of_ozone", \
           "atmosphere_mole_content_of_ozone_standard_error", \
           "merged_ozone_vmr", \
           "uncertainty_of_merged_ozone"/)

  ; File name
  FNAME = (/"ESACCI-OZONE-L3S-TC-MERGED-DLR_1M-_DATE_??-fv0100.nc", \
            "ESACCI-OZONE-L3S-TC-MERGED-DLR_1M-_DATE_??-fv0100.nc", \
            "ESACCI-OZONE-L3-LP-MERGED-MZM-_DATE_-fv0002.nc", \
            "ESACCI-OZONE-L3-LP-MERGED-MZM-_DATE_-fv0002.nc"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/custom/CMOR_" + VARS + ".dat"

end

begin

  do vID = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vID))

    time = create_timec(YEAR1(vID), YEAR2(vID))
    date = cd_calendar(time, 1)

    ; Create timeseries
    do yy = YEAR1(vID), YEAR2(vID)
      do mm = 1, 12

        ldate = yy + sprinti("%0.2i", mm)

        ; File name
        fname = systemfunc("ls " + input_dir_path + \
                           str_sub_str(FNAME(vID), "_DATE_", ldate))

        ; Check
        if (all(ismissing(fname))) then
          error_msg("f", DIAG_SCRIPT, "", "no file found for date " + ldate)
        end if

        ; Extract data
        f = addfile(fname(0), "r")
        xx = f->$NAME(vID)$
        xx@_FillValue = FILL
        xx@missing_value = xx@_FillValue
        xx = where(xx.lt.0., xx@_FillValue, xx)
        xx = where(xx.gt.1e35, xx@_FillValue, xx)  ; get rid of infinity values
        replace_ieeenan(xx, xx@_FillValue, 0)

        ; Assign to global array
        dimnames = getvardimnames(xx)
        if (.not.isdefined("output")) then
          dims = array_append_record(dimsizes(time), dimsizes(xx), 0)
          output = new(dims, typeof(xx))
          output!0 = "time"
          output&time = time
          if (dimnames(0).eq."latitude") then
            output!1 = "lat"
            output&lat = f->latitude
          end if
          if (dimnames(0).eq."air_pressure") then
            output!1 = "plev"
            output&plev = f->air_pressure
          end if
          output!2 = "lon"
          output&lon = f->$dimnames(1)$
        end if
        output(ind(toint(ldate).eq.date), :, :) = (/xx/)
        delete(fname)
        delete(xx)

      end do
    end do

    ; Format time coordinate
    ctime = time_attrib(output&time, YEAR1(vID), YEAR2(vID), "M")
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ; Format pressure level coordinate
    if (output!1.eq."plev") then
      if (isMonotonic(output&plev) .eq. 0) then
        error_msg("f", diag_script, "", "non-monotonic level coordinate")
      end if
      cplev = plev_attrib(output&plev)
      delete(output&plev)
      output&plev = cplev
      delete(cplev)
    end if

    ; Format latitude coordinate
    if (output!1.eq."lat") then
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
    end if

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
      "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/ozone/data/", \
      "Loyola et al., Int. J. Remote Sens. " + \
      "doi:10.1080/01431160902825016, 2009.", "")

    ; Outfile
    fout = output_dir_path + "OBS_" + OBSNAME + "_sat_L3_" + MIP + \
      "_" + VARS(vID) + "_" + YEAR1(vID) + "01-" + YEAR2(vID) + "12.nc"

    ; Write variable
    write_nc(fout, VARS(vID), output, gAtt)
    delete(gAtt)
    delete(output)
    delete(time)
    delete(date)

  end do

end
