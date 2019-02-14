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
;    20180810-A_righ_ma: fix minor calendar issue.
;    20180806-A_righ_ma: code cleaning.
;    20170116-A-laue_ax: using cirrus fraction to convert lwp and iwp
;                        from in-cloud values to gridbox averages
;                        (test).
;    20160408-A-laue_ax: added processing of uncertainties
;                        removed suffixes for variables names.
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

  ; Selected variable (standard name)
  VARS = (/"clwvi",     \
           "clivi",     \
           "clt",       \
           "lwpStderr", \
           "iwpStderr", \
           "od550aer"/)

  ; Name in the raw data
  NAME = (/"Cloud_Water_Path_Liquid_Mean_Mean",        \
           "Cloud_Water_Path_Ice_Mean_Mean",           \
           "Cloud_Fraction_Mean_Mean",                 \
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

  do ff = 0, dimsizes(FILES) - 1

    fin = addfile(FILES(ff), "r")

    ; Get time
    infile = systemfunc("basename " + FILES(ff))
    date   = yyyyddd_to_yyyymmdd(toint(str_get_cols(infile, 10, 16)))
    year   = toint(str_get_cols(tostring(date), 0, 3))
    month  = toint(str_get_cols(tostring(date), 4, 5))
    day    = toint(str_get_cols(tostring(date), 6, 7)) + 14  ; center in month

    time = cd_inv_calendar(year, month, day, 0, 0, 0, TUNITS, 0)
    time@calendar = "standard"
    time@long_name = "time"
    time@axis = "T"
    time@standard_name = "time"
    if (isatt(time, "_FillValue")) then
      delete(time@_FillValue)
    end if

    ; Loop over variables to fetch from input file
    do vv = 0, dimsizes(VARS) - 1

      invar = fin->$NAME(vv)$
      invar_fv = invar@_FillValue
      invar_coords = invar
      invar := tofloat(invar)
      invar := where(invar.eq.tofloat(invar_fv), \
                     default_fillvalue("float"), invar)

      ; Special case clwvi as the sum lwp + iwp
      if (VARS(vv).eq."clwvi") then
        if (NAME(vv).ne."Cloud_Water_Path_Liquid_Mean_Mean") then
          error_msg("f", DIAG_SCRIPT, "", "Cannot calculate clwvi")
        end if

        ; Read cirrus fraction
        ; cfin = fin->Cirrus_Fraction_SWIR_FMean
        cfin = fin->Cirrus_Fraction_Infrared_FMean
        cif = tofloat(cfin * cfin@scale_factor + cfin@add_offset)
        ; liquid fraction is estimated assuming random overlap, i.e.
        ; ctot = 1 - (1 - cif) * (1 - lif)
        ; --> lif = 1 - (1 - ctot) / (1 - cif)
        delete(cfin)
        cfin = fin->Cloud_Fraction_Mean_Mean
        ctot = tofloat(cfin * cfin@scale_factor + cfin@add_offset)
        delete(cfin)
        cif = where(cif.gt.0.999, cif@_FillValue, cif)
        lif = 1.0 - (1.0 - ctot) / (1.0 - cif)
        lif = where(lif.lt.0, 0, lif)
        tmpvar = fin->Cloud_Water_Path_Ice_Mean_Mean  ; read ice water path
        tmpvar_fv = tmpvar@_FillValue
        tmpvar := tofloat(tmpvar)
        tmpvar := where(tmpvar.eq.tofloat(tmpvar_fv), \
                        default_fillvalue("float"), \
                        tmpvar)
        tmpvar = tmpvar * cif  ; convert iwp in-cloud value to gridbox avg
        invar = invar * lif  ; convert lwp in-cloud value to grid-box avg
        invar = invar + tmpvar  ; clwvi = lwp + iwp
        delete(tmpvar)
        delete(lif)
        delete(cif)
        invar = 0.001 * invar  ; [g/m2] --> [kg/m2]
      end if

      ; lwp and iwp are in-cloud values
      ; convert lwp/iwp to grid-box averages by multiplying with
      ; average cloud fraction (not optimum but best we can do at the moment)
      if (any((/"clivi", "iwpStderr", "lwpStderr"/) .eq. VARS(vv))) then

        ; Read cirrus fraction (0-1)
        ; cfin = fin->Cirrus_Fraction_SWIR_FMean
        cfin = fin->Cirrus_Fraction_Infrared_FMean
        cf = tofloat(cfin * cfin@scale_factor + cfin@add_offset)
        delete(cfin)
        if (VARS(vv).eq."lwpStderr") then
          cfin = fin->Cloud_Fraction_Mean_Mean
          ctot = tofloat(cfin * cfin@scale_factor + cfin@add_offset)
          delete(cfin)
          cif = where(cf.gt.0.999, cf@_FillValue, cf)
          cf = 1.0 - (1.0 - ctot) / (1.0 - cif)
          cf = where(cf.lt.0, 0, cf)
          delete(cif)
          delete(ctot)
        end if
        invar = invar * cf  ; ; "grid-box average" lwp/iwp
        delete(cf)
        invar = 0.001 * invar  ; [g/m2] --> [kg/m2]
      end if

      invar@_FillValue = default_fillvalue("float")
      copy_VarCoords(invar_coords, invar)
      if (isatt(invar_coords, "scale_factor")) then
        invar = invar * tofloat(invar_coords@scale_factor)
      end if
      if (isatt(invar_coords, "add_offset")) then
        invar = invar + tofloat(invar_coords@add_offset)
      end if

      if (VARS(vv).eq."clt") then
        invar = 100.0 * invar  ; [1] --> [%]
      end if

      ; Create output variable
      lat = fin->YDim
      lon = fin->XDim
      output = new((/1, dimsizes(lat), dimsizes(lon)/), float)

      ; Set variable attributes
      output(0, :, :) = var_attrib(invar, VARS(vv), CMOR_TABLE(vv))

      ; Set reformatted time dimension
      output!0 = "time"
      output&time = time

      ; Set reformatted latitude dimension
      output!1 = "lat"
      output&lat = lat
      if (isMonotonic(output&lat) .eq. 0) then
        error_msg("f", DIAG_SCRIPT, "", "non-monotonic latitude coordinate")
      end if
      if (isMonotonic(output&lat) .eq. -1) then
        output = output(:, ::-1, :)
      end if
      clat := lat_attrib(output&lat)
      output&lat := clat
      delete(lat)

      ; Set reformatted longitude dimension
      output!2 = "lon"
      output&lon = lon

      if (isMonotonic(output&lon) .eq. 0) then
        error_msg("f", DIAG_SCRIPT, "", "non-monotonic longitude coordinate")
      end if
      if (any(output&lon .lt. 0.)) then
        output = lonFlip(output)
      end if
      clon := lon_attrib(output&lon)
      output&lon := clon
      log_info(min(clon) + "  " + max(clon))
      delete(lon)

      ; Set global attributes
      gAtt = set_global_atts( \
        OBSNAME, TIER, \
        "https://ladsweb.modaps.eosdis.nasa.gov/search/order", "", "")

      ; Define output file
      datestr = sprinti("%0.4i", year) + sprinti("%0.2i", month)
      fout = output_dir_path + "OBS_" + OBSNAME + "_sat_" + \
        str_sub_str(VERSION, "_", "-") + "_" + MIPS(vv) + "_" + VARS(vv) + \
        "_" + datestr + "-" + datestr + ".nc"

      ; Write attribute
      write_nc(fout, VARS(vv), output, gAtt)

      delete(invar)
      delete(invar_coords)
      delete(output)

    end do

  end do

end
