; #############################################################################
; ESMValTool CMORizer for NCEP data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    http://www.esrl.noaa.gov/psd/data/gridded/data.ncep.reanalysis.html
;
; Last access
;    20190204
;
; Download and processing instructions
;    To facilitate the download, the links to the ftp server are provided.
;    Since the filenames are sometimes identical across different
;    save the data in two subdirectories in input_dir_path.
;    Subdirectory pressure/:
;      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/pressure/
;        air.mon.mean.nc
;        hgt.mon.mean.nc
;        rhum.mon.mean.nc
;        shum.mon.mean.nc
;        uwnd.mon.mean.nc
;        vwnd.mon.mean.nc
;        omega.mon.mean.nc
;      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/pressure/
;        uwnd.????.nc
;        vwnd.????.nc
;    Subdirectory surface/:
;      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface/
;        air.mon.mean.nc
;      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.derived/surface_gauss/
;        prate.mon.mean.nc
;      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/surface_gauss/
;        prate.sft.gauss.????.nc
;      ftp://ftp.cdc.noaa.gov/Datasets/ncep.reanalysis.dailyavgs/other_gauss/
;        ulwrf.ntat.gauss.????.nc
;
;    Select the section "Pressure" and "Surface" and download the variables
;    listed below. Since raw data on pressure levels and for surface have the
;    same file and variable name, save the data in two different subdirectories
;    "press" and "surf" in input_dir_path.
;    Specify the time range of the data as YEAR1-YEAR2 below, considering only
;    complete years (Jan to Dec).
;
; Modification history
;    20190204-A_righ_ma: merged with NCEP-daily and adapted to v2.
;    20140128-A_righ_ma: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Source name
  OBSNAME = "NCEP"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1948
  YEAR2 = 2018

  ; Selected variable (standard name)
  VARS = (/"ta", "zg", "hur", "hus", "ua", \
           "va", "wap", "tas", "pr", \
           "ua", "va", "pr", "rlut"/)

  ; Name in the raw data
  NAME = (/"air", "hgt", "rhum", "shum", "uwnd", \
           "vwnd", "omega", "air", "prate", \
           "uwnd", "vwnd", "prate", "ulwrf"/)

  ; Subdirectory
  SUBDIR = (/"pressure/", "pressure/", "pressure/", "pressure/", "pressure/", \
             "pressure/", "pressure/", "surface/", "surface/", \
             "pressure/", "pressure/", "surface/", "surface/"/)

  ; Expected units (according to CMOR standard)
  EXPUNITS = (/"K", "m", "%", "1", "m/s", \
               "m/s", "Pascal/s", "K", "Kg/m^2/s", \
               "m/s", "m/s", "Kg/m^2/s", "W/m^2"/)

  ; MIPS
  MIPS = (/"Amon", "Amon", "Amon", "Amon", "Amon", \
           "Amon", "Amon", "Amon", "Amon", \
           "day", "day", "day", "day"/)

  ; CMOR tables
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/cmip5/Tables/CMIP5_" + MIPS

end

begin

  ; Loop over variables
  do vID = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vID))

    do yy = YEAR1, YEAR2

      if (isStrSubset(MIPS(vID), "mon")) then
        fname = input_dir_path + SUBDIR(vID) + NAME(vID) + ".mon.mean.nc"
      end if

      if (isStrSubset(MIPS(vID), "day")) then
        flist = systemfunc("ls " + input_dir_path + SUBDIR(vID) + NAME(vID) + \
                           ".*" + yy + ".nc")
        fname = flist(0)
        delete(flist)
      end if

      f = addfile(fname, "r")
      tmp = f->$NAME(vID)$
      delete(fname)
      delete(f)
      fill_start = num(ismissing(tmp))

      ; Extract time range
      tmp&time@calendar = "standard"
      date = cd_calendar(tmp&time, 0)
      idx = ind(date(:, 0).eq.yy)
      if (dimsizes(dimsizes(tmp)).eq.4) then
        output = tmp(idx, :, :, :)
      else
        output = tmp(idx, :, :)
      end if
      delete(idx)
      delete(tmp)
      delete(date)

      ; Unpack variable according to metadata information
      if (isatt(output, "scale_factor") .or. \
          isatt(output, "add_offset")) then
        tmp = output * output@scale_factor + output@add_offset
        copy_VarMeta(output, tmp)
        delete(output)
        output = tmp
        delete(tmp)
      end if

      ; Convert units
      if (any(VARS(vID).eq.(/"ta", "tas"/)) .and. \
          output@units.eq."degC") then
        output = output + 273.15
        output@units = "K"
      end if
      if (VARS(vID).eq."hus" .and. output@units.eq."grams/kg") then
        output = output / 1000.
        output@units = "1"
      end if
      if (output@units.ne.EXPUNITS(vID)) then
        error_msg("f", diag_script, "", \
                  "possibly wrong input units for " + VARS(vID))
      end if

      rank = dimsizes(dimsizes(output))

      ; Format time coordinate
      if (isStrSubset(MIPS(vID), "mon")) then
        ctime = time_attrib(output&time, yy, yy, "M")
      end if
      if (isStrSubset(MIPS(vID), "day")) then
        ctime = time_attrib(output&time, yy, yy, "D")
      end if
      delete(output&time)
      output&time = ctime
      delete(ctime)

      ; Format vertical coordinate
      if (rank.eq.4) then
        output!1 = "plev"
        if (isMonotonic(output&plev) .eq. 0) then
          error_msg("f", diag_script, "", "non-monotonic vertical coordinate")
        end if
        if (isMonotonic(output&plev).eq.1) then
          output = output(:, ::-1, :, :)
        end if
        if (output&plev@units.eq."millibar") then
          output&plev = output&plev * 100
        else
          error_msg("f", diag_script, "", \
                    "unexpected units for plev " + output&plev@units)
        end if
        cplev = plev_attrib(output&plev)
        delete(output&plev)
        output&plev = cplev
        delete(cplev)
      end if

      ; Format latitude coordinate
      if (rank.eq.4) then
        output!2 = "lat"
        if (isMonotonic(output&lat) .eq. 0) then
          error_msg("f", diag_script, "", "non-monotonic latitude coordinate")
        end if
        if (isMonotonic(output&lat) .eq. -1) then
          output = output(:, :, ::-1, :)
        end if
      end if
      if (rank.eq.3)
        output!1 = "lat"
        if (isMonotonic(output&lat) .eq. 0) then
          error_msg("f", diag_script, "", "non-monotonic latitude coordinate")
        end if
        if (isMonotonic(output&lat) .eq. -1) then
          output = output(:, ::-1, :)
        end if
      end if
      clat = lat_attrib(output&lat)
      delete(output&lat)
      output&lat = clat
      delete(clat)

      ; Format longitude coordinate
      if (rank.eq.4) then
        output!3 = "lon"
      end if
      if (rank.eq.3) then
        output!2 = "lon"
      end if
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
        "http://www.esrl.noaa.gov/psd/data/gridded/" + \
        "data.ncep.reanalysis.html", \
        "Kalnay et al., B. Am. Meteorol. Soc., " + \
        "doi:10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2, 1996", "")

      ; Check fill values
      fill_end = num(ismissing(output))
      if (fill_start.ne.fill_end) then
        error_msg("f", diag_script, "", \
                  "missing values lost during conversion")
      end if
      delete(fill_start)
      delete(fill_end)

      ; Outfile
      fout = output_dir_path + "OBS_" + OBSNAME + "_reanaly_1_" + MIPS(vID) + \
        "_" + VARS(vID) + "_" + yy + "01-" + yy + "12.nc"

      ; Add height coordinate to tas variable (required by the new backend)
      if (VARS(vID).eq."tas") then
        output@coordinates = "height"
      end if

      ; Write variable
      write_nc(fout, VARS(vID), output, gAtt)
      delete(gAtt)
      delete(output)

      ; Add height coordinate to tas variable (required by the new backend)
      if (VARS(vID).eq."tas") then
        height = 2.d
        height!0 = "ncl_scalar"
        height@units = "m"
        height@axis = "Z"
        height@positive = "up"
        height@long_name = "height"
        height@standard_name = "height"
        w = addfile(fout, "w")
        w->height = height
        delete(w)
      end if

    end do

  end do

end
