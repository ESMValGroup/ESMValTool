; #############################################################################
; ESMValTool CMORizer for ESACCI-SST data
; #############################################################################
;
; Tier
;   Tier 2: other freely-available dataset.
;
; Source
;   ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sst/data/
;
; Last access
;   20190201
;
; Download and processing instructions
;   Download the data from:
;     lt/Analysis/L4/v01.1/
;   Put all files under a single directory (no subdirectories with years).
;
; Modification history
;   20190201-A_righ_ma: adapted to v2.
;   20180529-A_righ_ma: modified to have consistent metadata across multiple
;                       files and to speed up processing time
;   20160818-A_laue_ax: written
;
;#############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ESACCI-SST.ncl"

  ; Source name
  OBSNAME = "ESACCI-SST"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1992
  YEAR2 = 2010

  ; Field
  FIELD = "T2Ms"

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + "/cmor/tables/"

end

begin

  ; Save date for consistent history attribute
  today = systemfunc("date")

  ; Loop over time period
  do yy = YEAR1, YEAR2
    do mm = 1, 12

      mo = sprinti("%0.2i", mm)

      ; Select files for this year/month
      files = \
        systemfunc("ls " + input_dir_path + "/" + yy + mo + "??120000-" + \
                   "ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_LT-v02.0-fv01.1.nc")
      nfiles = dimsizes(files)

      if (nfiles.ne.days_in_month(yy, mm)) then
        error_msg("f", DIAG_SCRIPT, "", "Incomplete data in " + yy + mo + \
                  " (" + nfiles + " != " + days_in_month(yy, mm) + \
                  " files found)")
      end if

      ; Create time-series
      f = addfiles(files, "r")
      ListSetType(f, "cat")

      ; Read mask
      lsmask = f[:]->mask
        
      ; Read SST
      tmp = f[:]->analysed_sst
      sst = tmp * tmp@scale_factor + tmp@add_offset
      copy_VarCoords(tmp, sst)
      delete(tmp)

      ; Read error
      tmp = f[:]->analysis_error
      err = tmp * tmp@scale_factor + tmp@add_offset
      copy_VarCoords(tmp, err)
      delete(tmp)

      ; Mask anything that is not open sea water (i.e. land, ice, lakes)
      sst = where(lsmask.eq.1, sst, sst@_FillValue)
      err = where(lsmask.eq.1, err, err@_FillValue)
      delete(lsmask)

      ; Calculate time averages
      sstavg = sst(0:0, :, :)
      sstavg(0, :, :) = dim_avg_n_Wrap(sst, 0)
      erravg = err(0:0, :, :)
      tmp = err ^ 2
      erravg(0, :, :) = sqrt(dim_avg_n_Wrap(tmp, 0))
      delete(tmp)
      delete(sst)
      delete(err)

      ; Format time coordinate
      sstavg!0 = "time"
      erravg!0 = "time"
      delete(sstavg&time)
      delete(erravg&time)
      time = cd_inv_calendar(yy, mm, 15, 0, 0, 0, TUNITS, 0)
      time@calendar = "standard"
      time@long_name = "time"
      time@axis = "T"
      time@units = TUNITS
      time@standard_name = "time"
      sstavg&time = time
      erravg&time = time

      ; Format latitude coordinate
      sstavg!1 = "lat"
      erravg!1 = "lat"
      if (isMonotonic(sstavg&lat).eq.0) then
        error_msg("f", DIAG_SCRIPT, "", "non-monotonic latitude coordinate")
      end if
      if (isMonotonic(sstavg&lat) .eq. -1) then
        sstavg = sstavg(:, ::-1, :)
        erravg = erravg(:, ::-1, :)
      end if
      clat = lat_attrib(sstavg&lat)
      delete(sstavg&lat)
      delete(erravg&lat)
      sstavg&lat = clat
      erravg&lat = clat

      ; Format longitude coordinate
      sstavg!2 = "lon"
      erravg!2 = "lon"
      if (isMonotonic(sstavg&lon).eq.0) then
        error_msg("f", DIAG_SCRIPT, "", "non-monotonic lonitude coordinate")
      end if
      if (any(sstavg&lon.lt.0.)) then
        sstavg = lonFlip(sstavg)
        erravg = lonFlip(erravg)
      end if
      clon = lon_attrib(sstavg&lon)
      delete(sstavg&lon)
      delete(erravg&lon)
      sstavg&lon = clon
      erravg&lon = clon

      ; Format variable metadata
      table = CMOR_TABLE + "cmip5/Tables/CMIP5_Amon"
      tmp = var_attrib(sstavg, "ts", table)
      sstavg := tmp
      delete(tmp)

      table = CMOR_TABLE + "custom/CMOR_tsStderr.dat"
      tmp = var_attrib(erravg, "tsStderr", table)
      erravg := tmp
      delete(tmp)

      ; Set global attributes
      gAtt = True
      gAtt@history = "Created on " + systemfunc("date")
      gAtt@host    = systemfunc("echo $HOST")
      gAtt@user    = systemfunc("echo $USER")
      gAtt@period  = YEAR1 + "-" + YEAR2
      gAtt@field   = FIELD
      gAtt@tier    = TIER
      gAtt@source  = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sst/data/"
      gAtt@reference = "Merchant et al. " + \
         "Geoscience Data Journal, doi:10.1002/gdj3.20, 2014."
      gAtt@title = OBSNAME +  " data reformatted for the ESMValTool"
      gAtt@conventions = "CF/CMOR"

      ; Write output
      fout = output_dir_path + "/OBS_" + OBSNAME + \
        "_sat_L4-GHRSST-SSTdepth-OSTIA-GLOB_" + \
        FIELD + "_ts_" + yy + mo + "-" + yy + mo + ".nc"
      write_nc(fout, "ts", sstavg, gAtt)
      delete(sstavg)

      fout = output_dir_path + "/OBS_" + OBSNAME + \
        "_sat_L4-GHRSST-SSTdepth-OSTIA-GLOB_" + \
        FIELD + "_tsStderr_" + yy + mo + "-" + yy + mo + ".nc"
      write_nc(fout, "tsStderr", erravg, gAtt)
      delete(erravg)

      delete(gAtt)
      delete(files)
      delete(f)

    end do
  end do

end
