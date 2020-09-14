; #############################################################################
; ESMValTool CMORizer for MODIS data
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
;    20180209-righi_mattia: fixed bug in lwpStderr.
;    20180209-hassler_birgit: adapted to v2.
;    20180810-righi_mattia: fix minor calendar issue.
;    20180806-righi_mattia: code cleaning.
;    20170116-lauer_axel: using cirrus fraction to gridbox averages.
;    20160408-lauer_axel: added processing of uncertainties.
;    20151118-lauer_axel: bugfix: added unit conversion.
;    20150430-evaldsson_martin: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "modis.ncl"

  ; Source name
  OBSNAME = "MODIS"

  ; Tier
  TIER = 3

  ; Selected variable (standard name)
  VAR = (/"clwvi",     \
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

  ; MIP
  MIP = (/"Amon", "Amon", "Amon", "Amon", "Amon", "aero"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon", "mon", "mon", "mon"/)

  ; Version
  VERSION = "MYD08_M3"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/cmip5/Tables/CMIP5_Amon", \
      "/custom/CMOR_lwpStderr.dat", \
      "/custom/CMOR_iwpStderr.dat", \
      "/cmip5/Tables/CMIP5_aero"/)

  ; Type
  TYPE = "sat"

  ; Global attributes
  SOURCE = "https://ladsweb.modaps.eosdis.nasa.gov/search/order"
  REF1 = "Platnick et al., IEEE Trans. Geosci. Remote Sens., " + \
    "doi:10.1109/TGRS.2002.808301, 2003."
  REF2 = "Levy et al., Atmos. Meas. Tech., " + \
    "doi:10.5194/amt-6-2989-2013, 2013."
  COMMENT = ""

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
    dm = days_in_month(year, month)

    ; Loop over variables to fetch from input file
    do vv = 0, dimsizes(VAR) - 1

      invar = fin->$NAME(vv)$
      invar_fv = invar@_FillValue
      invar_coords = invar
      invar := tofloat(invar)
      invar := where(invar.eq.tofloat(invar_fv), \
                     default_fillvalue("float"), invar)

      ; Special case clwvi as the sum lwp + iwp
      if (VAR(vv).eq."clwvi") then
        if (NAME(vv).ne."Cloud_Water_Path_Liquid_Mean_Mean") then
          error_msg("f", DIAG_SCRIPT, "", "cannot calculate clwvi")
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
      if (any((/"clivi", "iwpStderr", "lwpStderr"/) .eq. VAR(vv))) then

        ; Read cirrus fraction (0-1)
        ; cfin = fin->Cirrus_Fraction_SWIR_FMean
        cfin = fin->Cirrus_Fraction_Infrared_FMean
        cf = tofloat(cfin * cfin@scale_factor + cfin@add_offset)
        delete(cfin)
        if (VAR(vv).eq."lwpStderr") then
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

      if (VAR(vv).eq."clt") then
        invar = 100.0 * invar  ; [1] --> [%]
      end if

      ; Create output variable
      lat = fin->YDim
      lon = fin->XDim
      output = new((/1, dimsizes(lat), dimsizes(lon)/), float)
      output!0 = "time"
      output!1 = "lat"
      output!2 = "lon"
      output&time = cd_inv_calendar(year, month, 15, 0, 0, 0, TUNITS, 0)
      output&lat = lat
      output&lon = lon
      output(0, :, :) = (/invar/)
      delete(invar)
      delete(invar_coords)

      ; Format coordinates
      format_coords(output, year + sprinti("%0.2i", month) + "01", \
                    year + sprinti("%0.2i", month) + dm, FREQ(vv))

      ; Set variable attributes
      tmp = format_variable(output, VAR(vv), CMOR_TABLE(vv))
      delete(output)
      output = tmp
      delete(tmp)

      ; Calculate coordinate bounds
      bounds = guess_coord_bounds(output, FREQ(vv))

      ; Set global attributes
      if (VAR(vv).ne."od550aer") then
        gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF1, COMMENT)
      else
        gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF2, COMMENT)
      end if

      ; Output file
      DATESTR = \
        year + sprinti("%0.2i", month) + "-" + year + sprinti("%0.2i", month)
      fout = output_dir_path + \
        str_join((/"OBS", OBSNAME, TYPE, str_sub_str(VERSION, "_", "-"), \
                   MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

      ; Write variable
      write_nc(fout, VAR(vv), output, bounds, gAtt)
      delete(gAtt)
      delete(output)
      delete(bounds)

    end do

  end do

end
