; #############################################################################
; ESMValTool CMORizer for ERA-Interim data
; #############################################################################
;
; Tier
;    Tier 3: restricted dataset.
;
; Source
;    http://apps.ecmwf.int/datasets/data/interim-full-moda/
;
; Last access
;    20190205
;
; Download and processing instructions
;    Select Type of level ("Pressure levels" or "Surface") and
;    ERA Interim Fields ("Monthly means of daily means" or "Daily").
;    For "Daily" select only the step 0 and all times (00, 06, 12, 18).
;    Download the data on a single variable and single year basis, and save
;    them as ERA-Interim_<var>_<mean>_YYYY.nc,
;    where <var> is the ERA-Interim variable name and <mean> is either monthly
;    or daily. Further download "land-sea mask" from the "Invariant" data and
;    save it in ERA-Interim_lsm.nc.
;    It is also possible to download data in an automated way, see:
;      https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets
;      https://confluence.ecmwf.int/display/WEBAPI/Python+ERA-interim+examples
;    A registration is required for downloading the data.
;
; Caveats
;    The data are updated regularly: recent years are added, but also the past
;    years are sometimes corrected. To have a consistent timeseries, it is
;    therefore recommended to download the full timeseries and not just add
;    new years to a previous version of the data.
;
; Modification history
;    20190204-A_righ_ma: adapted to v2.
;    20171023-A_laue_ax: added variables prw, ts
;    20160913-A_laue_ax: added variable tcc
;    20150820-A_laue_ax: added variables lwp, iwp, sftlf
;    20150327-A_righ_ma: merged with ERA-Interim-surf.ncl and
;                        ERA-Interim-surf_daily.ncl.
;    20140128-A_righ_ma: written.
;
; #############################################################################

begin

  ; Script name (for logger)
  DIAG_SCRIPT = "cmorize_obs_ERA-Interim.ncl"

  ; Source name
  OBSNAME = "ERA-Interim"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1979
  YEAR2 = 2018

  ; Selected variable (standard name)
  VARS = (/"zg", "psl", "ta", "zg", "hus", \
           "ua", "va", "tas", "ps", "tos",  \
           "tauu", "tauv", "clwvi", "clivi", "sftlf", \
           "clt", "prw", "ts", "wap", "psl"/)

  ; Name in the raw data
  NAME = (/"z", "msl", "t", "z", "q", \
           "u", "v", "t2m", "sp", "sst", \
           "iews", "inss", "p56.162", "p57.162", "lsm", \
           "tcc", "tcwv", "skt", "w", "msl"/)

  ; Units in the input  data and corresponding conversion facto to CMOR units
  RAWUNITS = (/"m**2 s**-2", "Pa", "K", "m**2 s**-2", "kg kg**-1", \
               "m s**-1", "m s**-1", "K", "Pa", "K", \
               "N m**-2", "N m**-2", "kg m**-2", "kg m**-2", "(0 - 1)", \
               "(0 - 1)", "kg m**-2", "K", "Pa s**-1", "Pa"/)

  CONVERSION = (/1./9.80665, 1., 1., 1./9.80665, 1., \
                1., 1., 1., 1., 1., \
                1., 1., 1., 1., 100., \
                100., 1., 1., 1., 1./)

  ; MIPs
  MIPS = (/"day", "day", "Amon", "Amon", "Amon", \
           "Amon", "Amon", "Amon", "Amon", "Omon", \
           "Amon", "Amon", "Amon", "Amon", "fx", \
           "Amon", "Amon", "Amon", "Amon", "Amon"/)

  ; CMOR table
  CMOR_TABLE = getenv("esmvaltool_root") + \
    "/cmor/tables/cmip5/Tables/CMIP5_" + MIPS

end
; #############################################################################
loadscript(getenv("esmvaltool_root") + "/utils/cmorizers/obs/interface.ncl")

begin

  ; Loop over variables
  do vID = 0, dimsizes(VARS) - 1

    log_info("Processing " + VARS(vID))

    l_day = False
    l_fx = False
    if (isStrSubset(MIPS(vID), "mon")) then
      freq1 = "monthly"
      freq2 = "M"
    end if
    if (isStrSubset(MIPS(vID), "day")) then
      freq1 = "daily"
      freq2 = "D"
      l_day = True
    end if
    if (isStrSubset(MIPS(vID), "fx")) then
      l_fx = True
    end if

    ; Create timeseries
    do yy = YEAR1, YEAR2

      ; fx variables are time invariant
      if (l_fx) then
        fname = input_dir_path + "ERA-Interim_" + NAME(vID) + ".nc"
        f = addfile(fname, "r")
      else
        fname = input_dir_path + \
          "ERA-Interim_" + NAME(vID) + "_" + freq1 + "_" + yy + ".nc"
        f = addfile(fname, "r")

        ; Check time coordinate
        date = cd_calendar(f->time, 0)
        if (any(date(:, 0).ne.yy)) then
          error_msg("f", DIAG_SCRIPT, "", "time coordinate not " + \
                    "corresponding to selected year")
        end if
      end if

      ; Rank of the input data
      rank = dimsizes(getfilevardims(f, NAME(vID)))

      ; For daily data, split in monthly files, otherwise yearly
      if (l_day) then
        mend = 12
      else
        mend = 1  ; dummy loop, the whole year will be processed
      end if

      do mm = 1, mend

        ; Read variable
        if (l_day) then
          if (rank.eq.4) then
            xx = f->$NAME(vID)$(ind(date(:, 1).eq.mm), :, :, :)
          else if (rank.eq.3) then
            xx = f->$NAME(vID)$(ind(date(:, 1).eq.mm), :, :)
          else
            error_msg("f", DIAG_SCRIPT, "", "invalid rank")
          end if
          end if
        else
          xx = f->$NAME(vID)$
        end if

        ; Unpack variable according to metadata information
        tmp1 = todouble(xx)
        tmp1 = tmp1 * xx@scale_factor + xx@add_offset
        output = tofloat(tmp1)
        copy_VarMeta(xx, output)
        delete(xx)
        delete(tmp1)

        ; Calculate clwvi = lwp + iwp (Amon case only)
        if (VARS(vID).eq."clwvi") then

          if (l_day) then
            error_msg("f", DIAG_SCRIPT, "", "clwvi calculation for daily " + \
                      "values not implemented")
          end if

          ; Read 2nd file containing iwp (variable "p57.162")
          iwpname = "p57.162"
          fname2 = input_dir_path + \
            "ERA-Interim_" + iwpname + "_" + freq1 + "_" + yy + ".nc"
          f2 = addfile(fname2, "r")

          ; Read iwp from 2nd file
          xx2 = f2->$iwpname$

          ; Unpack variable according to metadata information
          tmp1a = todouble(xx2)
          tmp1a = tmp1a * xx2@scale_factor + xx2@add_offset
          tmp2a = tofloat(tmp1a)
          copy_VarMeta(xx2, tmp2a)
          delete(xx2)
          delete(tmp1a)

          ; Add iwp to lwp to calculate cmor variable "clwvi"
          output = output + tmp2a
          delete(tmp2a)

        end if

        ; Calculate daily mean
        if (l_day)
          if (rank.eq.4) then
            dmean = output(::4, :, :, :)
            do ii = 0, dimsizes(dmean&time) - 1
              dmean(ii, :, :, :) = \
                dim_avg_n_Wrap(output(4*ii:4*ii + 3, :, :, :), 0)
            end do
          else if (rank.eq.3) then
            dmean = output(::4, :, :)
            do ii = 0, dimsizes(dmean&time) - 1
              dmean(ii, :, :) = dim_avg_n_Wrap(output(4*ii:4*ii + 3, :, :), 0)
            end do
          else
            error_msg(f, DIAG_SCRIPT, "", "invalid rank")
          end if
          end if
          delete(output)
          output = dmean
          delete(dmean)
        end if

        ; Convert units
        if (output@units.ne.RAWUNITS(vID)) then
          error_msg("f", DIAG_SCRIPT, "", "unexpected units for " + NAME(vID))
        else
          output = output * CONVERSION(vID)
        end if

        ; Format time coordinate
        if (l_day) then
          newdate = cd_calendar(output&time, 0)
          newdate(:, 3) = 12.
          newdate(:, 4) = 0.
          newdate(:, 5) = 0.
          output!0 = "time"
          ctime = \
            cd_inv_calendar(newdate(:, 0), newdate(:, 1), newdate(:, 2), \
                            newdate(:, 3), newdate(:, 4), newdate(:, 5), \
                            TUNITS, 0)
          ctime@calendar = "standard"
          ctime@long_name = "time"
          ctime@axis = "T"
          ctime@standard_name = "time"
          if (isatt(ctime, "_FillValue")) then
            delete(ctime@_FillValue)
          end if
          delete(output&time)
          output&time = ctime
          delete(ctime)
          delete(newdate)
        else if (l_fx) then
          output := output(0, :, :)
          rank = 2
        else
          ctime = time_attrib(output&time, yy, yy, freq2)
          delete(output&time)
          output&time = ctime
          delete(ctime)
        end if
        end if

        ; Format vertical coordinate
        if (rank.eq.4) then
          output!1 = "plev"
          if (isMonotonic(output&plev).eq.1) then
            output = output(:, ::-1, :, :)
          end if
          if (output&plev@units.eq."millibars") then
            output&plev = output&plev * 100
          else
            error_msg("f", DIAG_SCRIPT, "", \
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
          if (isMonotonic(output&lat) .eq. -1) then
            output = output(:, :, ::-1, :)
          end if
        else if (rank.eq.3) then
          output!1 = "lat"
          if (isMonotonic(output&lat) .eq. -1) then
            output = output(:, ::-1, :)
          end if
        else if (rank.eq.2) then
          output!0 = "lat"
          if (isMonotonic(output&lat) .eq. -1) then
            output = output(::-1, :)
          end if
        else
          error_msg("f", DIAG_SCRIPT, "", "invalid rank")
        end if
        end if
        end if

        clat = lat_attrib(output&lat)
        delete(output&lat)
        output&lat = clat
        delete(clat)

        ; Format longitude coordinate
        if (rank.eq.4) then
          output!3 = "lon"
        else if (rank.eq.3) then
          output!2 = "lon"
        else if (rank.eq.2) then
          output!1 = "lon"
        else
          error_msg("f", DIAG_SCRIPT, "", "invalid rank")
        end if
        end if
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
          "http://apps.ecmwf.int/datasets/data/interim_full_moda/", \
          "Dee, D. P. et al., Q. J. Roy. Meteor. Soc., " + \
          "doi:10.1002/qj.828, 2011", "")

        ; Outfile
        if (l_fx) then
          datestr = ""
        else if (l_day) then
          datestr = \
            "_" + yy + sprinti("%0.2i", mm) + "-" + yy + sprinti("%0.2i", mm)
        else
          datestr = "_" + yy + "01-" + yy + "12"
        end if
        end if
        fout = output_dir_path + \
          "OBS_" + OBSNAME + "_reanaly_1_" + MIPS(vID) + \
          "_" + VARS(vID) + datestr + ".nc"

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

      end do  ; month (for daily data, otherwise dummy loop)

      ; No year loop need for fx (time invariant)
      if (l_fx) then
        break
      end if

      delete(date)

    end do  ; year

  end do  ; variable

end
