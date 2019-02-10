; #############################################################################
; REFORMAT SCRIPT FOR MODIS Level 3 OBSERVATIONAL DATA
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    https://ladsweb.modaps.eosdis.nasa.gov/search/order
;
; Last access
;    20190209
;
; Download and processing instructions
;    In Products: select "MODIS Aqua", "Collection 6.1" and
;                 "L3 Atmosphere Product", click on MYD08_M3.
;    In Time: select from 2000-01-01 to today.
;    In Location: skip, the global domain will be applied.
;    In Files: select all.
;    Submit the order.
;    A registration is required to download the data.
;
; Caveats
;    clwvi and clivi data are in-cloud values whereas CMIP5 models provide
;    grid-box averages --> multiply MODIS clwvi and clivi values with cloud
;    fraction as a first guess
;
; Modification history
;    20180209-A_hass_bg: adapted to v2.
;    20180810-A_righ_ma: fix minor calendar issue
;    20180806-A_righ_ma: code cleaning
;    20170116-A-laue_ax: using cirrus fraction to convert lwp and iwp
;                        from in-cloud values to gridbox averages
;                        (test)
;    20160408-A-laue_ax: added processing of uncertainties
;                        removed suffixes for variables names
;    20151118-A-laue_ax: bugfix: added unit conversion
;                        clivi, clwvi: g/m2 --> kg/m2
;                        clt: frac --> %
;    20150430-eval_ma: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_MODIS.ncl"

  ; Source name
  OBSNAME = "MODIS"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 2003
  YEAR2 = 2018

  ; Selected variable (standard name)
  VARS = (/"clt",       \
           "clwvi",     \
           "clivi",     \
           "lwpStderr", \
           "iwpStderr", \
           "od550aer"/)

  ; Name in the raw data
  NAME = (/"Cloud_Fraction_Mean_Mean",                \
           "Cloud_Water_Path_Liquid_Mean_Mean",        \
           "Cloud_Water_Path_Ice_Mean_Mean",           \
           "Cloud_Water_Path_Liquid_Mean_Uncertainty", \
           "Cloud_Water_Path_Ice_Mean_Uncertainty",    \
           "AOD_550_Dark_Target_Deep_Blue_Combined_Mean_Mean"/)

  ; MIPS
  MIPS = (/"Amon", "Amon", "Amon", "Amon", "Amon", "aero"/)

  ; Version
  VERSION = "MYD08_M3"

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + "/cmor/tables/" + \
    (/"cmip5/Tables/CMIP5_Amon", \
      "cmip5/Tables/CMIP5_Amon", \
      "cmip5/Tables/CMIP5_Amon", \
      "custom/CMOR_lwpStderr.dat", \
      "custom/CMOR_iwpStderr.dat", \
      "cmip5/Tables/CMIP5_aero"/)

end

begin

  ; List of files
  FILES = systemfunc("ls -1 " + input_dir_path + VERSION + ".A*.hdf")

  ; Create time coordinate
  time = create_timec(YEAR1, YEAR2)
  date = cd_calendar(time, 0)

  ; Define output array
  f = addfile(FILES(0), "r")
  if (.not.isdefined("output")) then
    lat = f->YDim
    lon = f->XDim
    output = new((/dimsizes(time), dimsizes(lat), dimsizes(lon)/), float)
    output!0 = "time"
    output&time = time
    output!1 = "lat"
    output&lat = lat
    output!2 = "lon"
    output&lon = lon
  end if

  ; Loop over variable
  do vv = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vv))

    output = 0.

    ; Loop over files
    do ff = 0, dimsizes(FILES) - 1

      ; Get date
      infile = systemfunc("basename " + FILES(ff))
      ymd   = yyyyddd_to_yyyymmdd(toint(str_get_cols(infile, 10, 16)))
      year   = toint(str_get_cols(tostring(ymd), 0, 3))
      month  = toint(str_get_cols(tostring(ymd), 4, 5))
      if (year.lt.YEAR1 .or. year.gt.YEAR2) then
        continue
      end if
      tidx = ind(date(:, 0).eq.year .and. date(:, 1).eq.month)
      f = addfile(FILES(ff), "r")

      tmp = f->$NAME(vv)$
      xx = tofloat(tmp * tmp@scale_factor + tmp@add_offset)
      delete(tmp)

      ; Special case clwvi as the sum lwp + iwp
      if (VARS(vv).eq."clwvi") then

        ; Read cirrus cloud fraction and estimate liquid cloud fraction
        tmp = f->Cirrus_Fraction_SWIR_FMean
        cirrus_frac = tofloat(tmp * tmp@scale_factor + tmp@add_offset)
        liquid_frac = 1. - cirrus_frac
        delete(tmp)

        ; Read ice water path
        tmp = f->Cloud_Water_Path_Ice_Mean_Mean
        iwp = tofloat(tmp * tmp@scale_factor + tmp@add_offset)

        ; Convert in-cloud to gridbox
        iwp = iwp * liquid_frac  ; convert lwp in-cloud value to gridbox avg
        xx = xx * cirrus_frac  ; convert iwp in-cloud value to grid-box avg
        delete(liquid_frac)
        delete(cirrus_frac)

        ; Add lwp and iwp
        xx = xx + iwp
        delete(iwp)

      end if

      ; Input lwp and iwp are in-cloud values, convert to grid-box averages
      ; by multiplying with average cloud fraction (not optimal but best we
      ; can do at the moment)
      if (any(VARS(vv).eq. \
              (/"clwvi", "clivi", "iwpStderr", "lwpStderr"/))) then

        ; Read cirrus cloud fraction
        tmp = f->Cirrus_Fraction_SWIR_FMean
        cirrus_frac = tofloat(tmp * tmp@scale_factor + tmp@add_offset)
        delete(tmp)

        ; Convert to grid-box average
        if (any(VARS(vv).eq.(/"clwvi", "lwpStderr"/))) then
          cirrus_frac = 1.0 - cirrus_frac
        end if
        xx = xx * cirrus_frac
        delete(cirrus_frac)
        xx = 0.001 * xx  ; [g/m2] --> [kg/m2]

      end if

      if (VARS(vv).eq."clt") then
        xx = 100. * xx  ; [1] --> [%]
      end if

      ; Assign to global array
      output(tidx, :, :) = (/xx/)
      delete(xx)

    end do

    ; Format time coordinate
    ctime = time_attrib(output&time, YEAR1, YEAR2, "M")
    delete(output&time)
    output&time = ctime
    delete(ctime)

    ; Format latitude coordinate
    output!1 = "lat"
    if (isMonotonic(output&lat) .eq. 0) then
      error_msg("f", DIAG_SCRIPT, "", "non-monotonic latitude coordinate")
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
      error_msg("f", DIAG_SCRIPT, "", "non-monotonic longitude coordinate")
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
      "https://ladsweb.modaps.eosdis.nasa.gov/search/order", "", "")

    ; Outfile
    fout = output_dir_path + "OBS_" + OBSNAME + "_sat_" + \
      str_sub_str(VERSION, "_", "-") + "_" + MIPS(vv) + "_" + VARS(vv) + \
      "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

    ; Write variable
    write_nc(fout, VARS(vv), output, gAtt)

  end do

end
