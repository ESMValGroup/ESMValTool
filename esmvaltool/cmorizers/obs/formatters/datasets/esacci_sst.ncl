; #############################################################################
; ESMValTool CMORizer for ESACCI-SST data
; #############################################################################
;
; Tier
;    Tier 2: other freely-available dataset.
;
; Source
;    ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sst/data/
;
; Last access
;    20190201
;
; Download and processing instructions
;    Download the data from:
;      lt/Analysis/L4/v01.1/
;    Put all files under a single directory (no subdirectories with years).
;
; Modification history
;    20190201-righi_mattia: adapted to v2.
;    20180529-righi_mattia: modified to have consistent metadata across
;                           multiple files and to speed up processing time.
;    20160818-lauer_axel: written.
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/obs/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "esacci_sst.ncl"

  ; Source name
  OBSNAME = "ESACCI-SST"

  ; Tier
  TIER = 2

  ; Period
  YEAR1 = 1992
  YEAR2 = 2010

  ; MIP
  MIP = "Amon"

  ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables")

  ; Type
  TYPE = "sat"

  ; Version
  VERSION = "L4-GHRSST-SSTdepth-OSTIA-GLOB"

  ; Global attributes
  SOURCE = "ftp://anon-ftp.ceda.ac.uk/neodc/esacci/sst/data/"
  REF = "Merchant et al., Geosci. Data J., doi:10.1002/gdj3.20, 2014"
  COMMENT = ""

end

begin

  ; Save date for consistent history attribute
  today = systemfunc("date")

  ; Loop over time period
  do yy = YEAR1, YEAR2
    do mm = 1, 12

      mo = sprinti("%0.2i", mm)

      dm = days_in_month(yy, mm)

      ; Select files for this year/month
      files = \
        systemfunc("ls " + input_dir_path + yy + mo + "??120000-" + \
                   "ESACCI-L4_GHRSST-SSTdepth-OSTIA-GLOB_LT-v02.0-fv01.1.nc")
      nfiles = dimsizes(files)

      if (nfiles.ne.dm) then
        error_msg("f", DIAG_SCRIPT, "", "incomplete data in " + yy + mo + \
                  " (" + dm + " != " + nfiles + " files found)")
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
      sstavg!1 = "lat"
      sstavg!2 = "lon"
      format_coords(sstavg, yy + mo + "01", yy + mo + dm, FREQ)
      erravg!0 = "time"
      erravg!1 = "lat"
      erravg!2 = "lon"
      format_coords(erravg, yy + mo + "01", yy + mo + dm, FREQ)

      ; Format variable metadata
      table = CMOR_TABLE + "/cmip5/Tables/CMIP5_Amon"
      tmp = format_variable(sstavg, "ts", table)
      delete(sstavg)
      sstavg = tmp
      delete(tmp)

      table = CMOR_TABLE + "/custom/CMOR_tsStderr.dat"
      tmp = format_variable(erravg, "tsStderr", table)
      delete(erravg)
      erravg = tmp
      delete(tmp)

      ; Calculate coordinate bounds
      bounds_sst = guess_coord_bounds(sstavg, FREQ)
      bounds_err = guess_coord_bounds(erravg, FREQ)

      ; Set global attributes
      gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

      ; Write output
      DATESTR = yy + mo + "-" + yy + mo
      fout = output_dir_path + \
        str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                   MIP, "ts", DATESTR/), "_") + ".nc"
      write_nc(fout, "ts", sstavg, bounds_sst, gAtt)
      delete(sstavg)
      delete(bounds_sst)

      fout = output_dir_path + \
        str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                   MIP, "tsStderr", DATESTR/), "_") + ".nc"
      write_nc(fout, "tsStderr", erravg, bounds_err, gAtt)
      delete(erravg)
      delete(bounds_err)

      delete(gAtt)
      delete(files)
      delete(f)

    end do
  end do

end
