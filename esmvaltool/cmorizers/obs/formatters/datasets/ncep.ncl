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
;    20190204-righi_mattia: merged with NCEP-daily and adapted to v2.
;    20140128-righi_mattia: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "ncep.ncl"

  ; Source name
  OBSNAME = "NCEP"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1948
  YEAR2 = 2018

  ; Selected variable (standard name)
  VAR = (/"ta", "zg", "hur", "hus", "ua", \
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

  ; MIP
  MIP = (/"Amon", "Amon", "Amon", "Amon", "Amon", \
          "Amon", "Amon", "Amon", "Amon", \
          "day", "day", "day", "day"/)

  ; Frequency
  FREQ = (/"mon", "mon", "mon", "mon", "mon", \
           "mon", "mon", "mon", "mon", \
           "day", "day", "day", "day"/)

  ; CMOR tables
  CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_" + MIP

  ; Type
  TYPE = "reanaly"

  ; Version
  VERSION = "1"

  ; Global attributes
  SOURCE = "http://www.esrl.noaa.gov/psd/data/gridded/" + \
    "data.ncep.reanalysis.html"
  REF = "Kalnay et al., B. Am. Meteorol. Soc., " + \
    "doi:10.1175/1520-0477(1996)077<0437:TNYRP>2.0.CO;2, 1996"
  COMMENT = ""

end

begin

  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP(vv) + ")")

    do yy = YEAR1, YEAR2

      if (isStrSubset(MIP(vv), "mon")) then
        fname = input_dir_path + SUBDIR(vv) + NAME(vv) + ".mon.mean.nc"
      end if

      if (isStrSubset(MIP(vv), "day")) then
        flist = systemfunc("ls " + input_dir_path + SUBDIR(vv) + NAME(vv) + \
                           ".*" + yy + ".nc")
        fname = flist(0)
        delete(flist)
      end if

      f = addfile(fname, "r")
      tmp = f->$NAME(vv)$
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
      if (any(VAR(vv).eq.(/"ta", "tas"/)) .and. \
          output@units.eq."degC") then
        output = output + 273.15
        output@units = "K"
      end if
      if (VAR(vv).eq."hus" .and. output@units.eq."grams/kg") then
        output = output / 1000.
        output@units = "1"
      end if
      if (output@units.ne.EXPUNITS(vv)) then
        error_msg("f", DIAG_SCRIPT, "", \
                  "possibly wrong input units for " + VAR(vv))
      end if

      rank = dimsizes(dimsizes(output))
      output!0 = "time"
      if (rank.eq.4) then
        output!1 = "plev"
        output!2 = "lat"
        output!3 = "lon"
        output&plev = output&plev * 100.  ; [mb] --> [Pa]
      elseif (rank.eq.3)
        output!1 = "lat"
        output!2 = "lon"
      end if

      ; Format coordinates
      format_coords(output, yy + "0101", yy + "1231", FREQ(vv))

      ; Set variable attributes
      tmp = format_variable(output, VAR(vv), CMOR_TABLE(vv))
      delete(output)
      output = tmp
      delete(tmp)

      ; Calculate coordinate bounds
      bounds = guess_coord_bounds(output, FREQ(vv))

      ; Set global attributes
      gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

      ; Check fill values
      fill_end = num(ismissing(output))
      if (fill_start.ne.fill_end) then
        error_msg("f", DIAG_SCRIPT, "", \
                  "missing values lost during conversion")
      end if
      delete(fill_start)
      delete(fill_end)

      ; Output file
      DATESTR = yy + "01-" + yy + "12"
      fout = output_dir_path + \
        str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                   MIP(vv), VAR(vv), DATESTR/), "_") + ".nc"

      ; Add height coordinate to tas variable (required by the new backend)
      if (VAR(vv).eq."tas") then
        output@coordinates = "height"
      end if

      ; Write variable
      write_nc(fout, VAR(vv), output, bounds, gAtt)
      delete(gAtt)
      delete(output)
      delete(bounds)

      ; Add height coordinate to tas variable (required by the new backend)
      if (VAR(vv).eq."tas") then
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
