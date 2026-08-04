; #############################################################################
; ESMValTool CMORizer for NASA MERRA reanalysis v5.2.0
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset (registration required).
;
; Source
;    EarthData via https://goldsmr3.gesdisc.eosdis.nasa.gov/data/MERRA_MONTHLY/
;
; Last access
;    20230201
;
; Download and processing instructions
; (requires EarthData login; see https://urs.earthdata.nasa.gov/)
;    Use ESMValTool automatic download:
;       esmvaltool data download MERRA
;
; Modification history
;    20230818-lauer_axel: added output of clwvi (iwp + lwp)
;    20230201-lauer_axel: written
;
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/data/formatters/interface.ncl")

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "merra.ncl"

  ; Source name
  OBSNAME = "MERRA"

  ; Tier
  TIER = 3

  ; Period (complete years only)
  YEAR1 = get_year(start_year, 1979)
  YEAR2 = get_year(end_year, 2015)

  ; Selected variable (standard name)
  VAR = (/"zg", "wap", "cli", "clw", "hus", "hur", "ta", "ua", "va",  \
          "ps", "psl",                                                \
          "clivi", "clwvi", "lwp", "prw",                             \
          "clt", "rlut", "rlutcs", "rsdt", "rsut", "rsutcs", "ts",    \
          "tas", "pr"/)

  ; Name in the raw data
  NAME = (/"H", "OMEGA", "QI", "QL", "QV", "RH", "T", "U", "V",  \  ; 3d asm
           "PS", "SLP",                                          \
           "TQI", "TQI", "TQL", "TQV",                           \  ; 2d int
           "CLDTOT", "LWTUP", "LWTUPCLR", "SWTDN", "SWTNT",      \  ; 2d rad
           "SWTNTCLR", "TS",                                     \
           "T2M", "PRECLSC"/)                                       ; 2d chm

  ; unit conversion factor
  CONVFAC = (/1.0, 1.0, 1.0, 1.0, 1.0, 100.0, 1.0, 1.0, 1.0,  \  ; 3d asm
             1.0, 1.0,                                        \
             1.0, 1.0, 1.0, 1.0,                              \  ; 2d int
             100.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,             \  ; 2d rad
             1.0, 1.0/)                                          ; 2d chm

  ; additional attribute "positive" (for radiative fluxes)
  POSITIVE = (/"", "", "", "", "", "", "", "", "", "", "",    \  ; 3d asm
               "", "", "", "",                                \  ; 2d int
               "", "up", "up", "down", "up", "up", "",        \  ; 2d rad
               "", ""/)                                          ; 2d chm

  ; names source files
  SOURCEFILE = (/"instM_3d_asm_Cp.", "instM_3d_asm_Cp.", \  ; 3d asm
                 "instM_3d_asm_Cp.", "instM_3d_asm_Cp.", \
                 "instM_3d_asm_Cp.", "instM_3d_asm_Cp.", \
                 "instM_3d_asm_Cp.", "instM_3d_asm_Cp.", \
                 "instM_3d_asm_Cp.", "instM_3d_asm_Cp.", \
                 "instM_3d_asm_Cp.",                     \
                 "instM_2d_int_Nx.", "instM_2d_int_Nx.", \  ; 2d int
                 "instM_2d_int_Nx.", "instM_2d_int_Nx.", \
                 "tavgM_2d_rad_Nx.", "tavgM_2d_rad_Nx.", \  ; 2d rad
                 "tavgM_2d_rad_Nx.", "tavgM_2d_rad_Nx.", \
                 "tavgM_2d_rad_Nx.", "tavgM_2d_rad_Nx.", \
                 "tavgM_2d_rad_Nx.",                     \
                 "tavgM_2d_chm_Fx.", "tavgM_2d_chm_Fx."/)   ; 2d chm

  ; dataset doi numbers
  DOI = (/"10.5067/YX0AVASQRTNW", "10.5067/YX0AVASQRTNW",  \  ; 3d asm
          "10.5067/YX0AVASQRTNW", "10.5067/YX0AVASQRTNW",  \
          "10.5067/YX0AVASQRTNW", "10.5067/YX0AVASQRTNW",  \
          "10.5067/YX0AVASQRTNW", "10.5067/YX0AVASQRTNW",  \
          "10.5067/YX0AVASQRTNW", "10.5067/YX0AVASQRTNW",  \
          "10.5067/YX0AVASQRTNW",                          \
          "10.5067/QL0PGBK2CYJS", "10.5067/QL0PGBK2CYJS",  \  ; 2d int
          "10.5067/QL0PGBK2CYJS", "10.5067/QL0PGBK2CYJS",  \
          "10.5067/6UX3EDUNVUFK", "10.5067/6UX3EDUNVUFK",  \  ; 2d rad
          "10.5067/6UX3EDUNVUFK", "10.5067/6UX3EDUNVUFK",  \
          "10.5067/6UX3EDUNVUFK", "10.5067/6UX3EDUNVUFK",  \
          "10.5067/6UX3EDUNVUFK",                          \
          "10.5067/IYDN3LNZ63UE", "10.5067/IYDN3LNZ63UE"/)    ; 2d chm

  ; MIP
  MIP = "Amon"

  ; Frequency
  FREQ = "mon"

  ; CMOR table
  CMOR_TABLE = getenv("cmor_tables") + \
    (/"/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \  ; 3d asm
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \
      "/cmip5/Tables/CMIP5_" + MIP,                                \
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \  ; 2d int
      "/custom/CMOR_lwp.dat", "/cmip5/Tables/CMIP5_" + MIP,        \
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \  ; 2d rad
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP,  \
      "/cmip5/Tables/CMIP5_" + MIP,                                \
      "/cmip5/Tables/CMIP5_" + MIP, "/cmip5/Tables/CMIP5_" + MIP/)    ; 2d chm

  ; Type
  TYPE = "reanaly"

  ; Version
  VERSION = "5.2.0"

  ; Global attributes
  SOURCE = "https://goldsmr3.gesdisc.eosdis.nasa.gov/data/MERRA_MONTHLY/"
  COMMENT = "Goddard Earth Sciences Data and Information Services Center " + \
            "(GES DISC)"

end

begin

  ; Loop over variables
  do vv = 0, dimsizes(VAR) - 1

    log_info("Processing " + VAR(vv) + " (" + MIP + ")")
    if (isvar("output")) then
      delete(output)
    end if

    do yy = YEAR1, YEAR2
      do mm = 1, 12

        datestr = tostring(yy)
        if (mm .lt. 10) then
          datestr = datestr + "0"
        end if
        datestr = datestr + tostring(mm)

        fname = systemfunc("ls " + input_dir_path + "MERRA???.prod.assim." + \
                           SOURCEFILE(vv) + datestr + ".hdf")

        f = addfile(fname, "r")
        tmp = f->$NAME(vv)$

        ; Extract time range
        tmp&TIME_EOSGRID@calendar = "standard"
        date = cd_calendar(tmp&TIME_EOSGRID, 0)
        if ((date(0, 0) .ne. yy) .or. (date(0, 1) .ne. mm)) then
          error_msg("f", DIAG_SCRIPT, "", \
                    "date in input file does not match date in filename: " + \
                    fname)
        end if

        delete(date)

        if (.not.isvar("output")) then
          dims = dimsizes(tmp)
          ; overwrite time dimension
          dims(0) = 12
          output = new(dims, float)
          delete(dims)
          rank = dimsizes(dimsizes(output))

          output!0 = "time"
          if (rank.eq.4) then
            output!1 = "plev"
            output!2 = "lat"
            output!3 = "lon"
            output&plev = tmp&Height_EOSGRID * 100.  ; [hPa] --> [Pa]
          elseif (rank.eq.3)
            output!1 = "lat"
            output!2 = "lon"
          end if

          output&time = fspan(1, 12, 12)
          output&time@calendar = "standard"
          output&time@units = "days since 1950-01-01 00:00:00"

          output&lat = tmp&YDim_EOSGRID
          output&lon = tmp&XDim_EOSGRID

        end if

        ; Unpack variable according to metadata information
        if (isatt(tmp, "scale_factor") .or. isatt(tmp, "add_offset")) then
          tmp = tmp * tmp@scale_factor + tmp@add_offset
        end if

        if (rank.eq.4) then
          output(mm - 1, :, :, :) = (/ tmp(0, :, :, :) /)
        else
          output(mm - 1, :, :) = (/ tmp(0, :, :) /)
        end if

        delete(tmp)

        ; calculation of outgoing fluxes: out = in - net
        if ((VAR(vv) .eq. "rsut") .or. (VAR(vv) .eq. "rsutcs")) then
          tmp = f->SWTDN
          if (isatt(tmp, "scale_factor") .or. isatt(tmp, "add_offset")) then
            tmp = tmp * tmp@scale_factor + tmp@add_offset
          end if
          output(mm - 1, :, :) = (/ tmp(0, :, :) /) - output(mm - 1, :, :)

          delete(tmp)
        end if

        ; calculation of total precipitation flux =
        ; large-scale+convective+anvil
        if (VAR(vv) .eq. "pr") then
          tmp = f->PRECCON  ; surface precipitation flux from convection
          if (isatt(tmp, "scale_factor") .or. isatt(tmp, "add_offset")) then
            tmp = tmp * tmp@scale_factor + tmp@add_offset
          end if
          output(mm - 1, :, :) = output(mm - 1, :, :) + tmp(0, :, :)
          delete(tmp)
          tmp = f->PRECANV  ; surface precipitation flux from anvils
          if (isatt(tmp, "scale_factor") .or. isatt(tmp, "add_offset")) then
            tmp = tmp * tmp@scale_factor + tmp@add_offset
          end if
          output(mm - 1, :, :) = output(mm - 1, :, :) + tmp(0, :, :)
          delete(tmp)
        end if

        ; calculation of clwvi
        if (VAR(vv) .eq. "clwvi") then
          tmp = f->TQL
          if (isatt(tmp, "scale_factor") .or. isatt(tmp, "add_offset")) then
            tmp = tmp * tmp@scale_factor + tmp@add_offset
          end if
          output(mm - 1, :, :) = output(mm - 1, :, :) + tmp(0, :, :)
          delete(tmp)
        end if

        delete(f)

      end do  ; loop over months (mm)

      ; Convert units
      if (CONVFAC(vv) .ne. 1.0) then
        output = output * CONVFAC(vv)
      end if

      ; Format coordinates
      format_coords(output, yy + "0101", yy + "1231", FREQ)

      ; Set variable attributes
      tmp = format_variable(output, VAR(vv), CMOR_TABLE(vv))
      delete(output)
      output = tmp
      delete(tmp)

      ; set additional attribute "positive" for radiative fluxes
      if (POSITIVE(vv) .ne. "") then
        output@positive = POSITIVE(vv)
      end if

      ; Calculate coordinate bounds
      bounds = guess_coord_bounds(output, FREQ)

      ; Set global attributes
      REF = "doi: " + DOI(vv)
      gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

      ; Output file
      DATESTR = yy + "01-" + yy + "12"
      fout = output_dir_path + \
        str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                   MIP, VAR(vv), DATESTR/), "_") + ".nc"

      ; Add height coordinate to tas variable (required by the new backend)
      if (VAR(vv).eq."tas") then
        output@coordinates = "height"
      end if

      ; Write variable
      write_nc(fout, VAR(vv), output, bounds, gAtt)
      delete(gAtt)
      delete(output)
      delete(bounds)

      ; Add height coordinate to tas variable (required by CMOR checker)
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

        ; ---------------------------------------------------------------------

    end do  ; loop over years (yy)
  end do  ; loop over variables (vv)

end
