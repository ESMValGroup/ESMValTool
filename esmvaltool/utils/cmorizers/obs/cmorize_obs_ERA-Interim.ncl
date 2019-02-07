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
  YEAR2 = 1980  ; 2018

  ; Selected variable (standard name)
  VARS = (/"psl", "ta", "zg", "hus", "ua", \
           "va", "tas", "ps", "tos", "tauu", \
           "tauv", "clwvi", "clivi", "sftlf", "clt", \
           "prw", "ts", "wap", "psl"/)

  ; Name in the raw data
  NAME = (/"msl", "t", "z", "q", "u", \
           "v", "t2m", "sp", "sst", "iews", \
           "inss", "p56.162", "p57.162", "lsm", "tcc", \
           "tcwv", "skt", "w", "msl"/)

  ; Units in the input  data and corresponding conversion facto to CMOR units
  RAWUNITS = (/"Pa", "K", "m**2 s**-2", "kg kg**-1", "m s**-1", \
               "m s**-1", "K", "Pa", "K", "N m**-2", \
               "N m**-2", "kg m**-2", "kg m**-2", "(0 - 1)", "(0 - 1)",\
               "kg m**-2", "K", "Pa s**-1", "Pa"/)

  CONVERSION = (/1., 1., 1./9.80665, 1., 1., \
                1., 1., 1., 1., 1., \
                1., 1., 1., 100., 100.,\
                1., 1., 1., 1./)

  ; MIPs
  MIPS = (/"day", "Amon", "Amon", "Amon", "Amon", \
           "Amon", "Amon", "Amon", "Omon", "Amon", \
           "Amon", "Amon", "Amon", "fx", "Amon", \
           "Amon", "Amon", "Amon", "Amon"/)

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

    fill_start = 0

    if (isStrSubset(MIPS(vID), "mon")) then
      freq1 = "monthly"
      freq2 = "M"
    end if
    if (isStrSubset(MIPS(vID), "day")) then
      freq1 = "daily"
      freq2 = "D"
    end if

    ; Create timeseries
    do yy = YEAR1, YEAR2

      ; land sea mask is invariant
      if (MIPS(vID).eq."fx") then
        fname = input_dir_path + "ERA-Interim_" + NAME(vID) + ".nc"
        print("path1: " + fname) 
        f = addfile(fname, "r")
      else
        fname = input_dir_path + \
          "ERA-Interim_" + NAME(vID) + "_" + freq1 + "_" + yy + ".nc"
        print("path2: " + fname) 
        f = addfile(fname, "r")

        ; Check time coordinate
        date = cd_calendar(f->time, 0)
        print("Date: " + date)
        if (any(date(:, 0).ne.yy)) then
          error_msg("f", DIAG_SCRIPT, "", "time coordinate not " + \
                    "corresponding to selected year")
        end if
        delete(date)
      end if

      ; Read variable
      xx = f->$NAME(vID)$
      
      ; Add to missing values count
      fill_start = fill_start + num(ismissing(xx))
 
      ; Unpack variable according to metadata information
      tmp1 = todouble(xx)
      tmp1 = tmp1 * xx@scale_factor + xx@add_offset
      tmp2 = tofloat(tmp1)
      copy_VarMeta(xx, tmp2)
      delete(xx)
      delete(tmp1)
      
      ; Calculate clwvi = lwp + iwp
      if (VARS(vID).eq."clwvi") then

        ; Read 2nd file containing iwp (variable "p57.162")
        iwpname = "p57.162"
        fname2 = input_dir_path + \
          "ERA-Interim_" + iwpname + "_" + freq1 + "_" + yy + ".nc"
        f2 = addfile(fname2, "r")

        ; Check time coordinate of 2nd file
        date = cd_calendar(f2->time, 0)
        if (any(date(:, 0).ne.yy)) then
          error_msg("f", DIAG_SCRIPT, "", "time coordinate not " + \
                    "corresponding to selected year")
        end if
        delete(date)
        
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
        tmp2 = tmp2 + tmp2a
        delete(tmp2a)

      end if

      rank = dimsizes(dimsizes(tmp2))

      ; Calculate daily mean
      if (isStrSubset(MIPS(vID), "day"))
        if (rank.eq.4) then
          dmean = tmp2(::4, :, :, :)
          do ii = 0, dimsizes(dmean&time) - 1
            dmean(ii, :, :, :) = \
              dim_avg_n_Wrap(tmp2(4*ii:4*ii + 3, :, :, :), 0)
          end do
        else
          dmean = tmp2(::4, :, :)
          do ii = 0, dimsizes(dmean&time) - 1
            dmean(ii, :, :) = dim_avg_n_Wrap(tmp2(4*ii:4*ii + 3, :, :), 0)
          end do
        end if
        delete(tmp2)
        tmp2 = dmean
        delete(dmean)
      end if

      ; Add to time series
      if (.not.isdefined("output")) then
        output = tmp2
      else
        tmp3 = array_append_record(output, tmp2, 0)
        delete(output)
        output = tmp3
        delete(tmp3)
      end if
      delete(tmp2)
    end do

    ; Convert units
    if (output@units.ne.RAWUNITS(vID)) then
      error_msg("f", DIAG_SCRIPT, "", "unexpected units for " + NAME(vID))
    else
      output = output * CONVERSION(vID)
    end if

    ; Format time coordinate
    if (MIPS(vID).ne."fx") then
      ctime = time_attrib(output&time, yy, yy, freq2)
      delete(output&time)
      output&time = ctime
      delete(ctime)
    else
      output := output(0, :, :)
      rank = 2
    end if

    ; Format vertical coordinate
    if (rank.eq.4) then
      output!1 = "plev"
      if (isMonotonic(output&plev) .eq. 0) then
        error_msg("f", DIAG_SCRIPT, "", "non-monotonic vertical coordinate")
      end if
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
    end if
    end if
    end if
    if (isMonotonic(output&lat) .eq. 0) then
      error_msg("f", DIAG_SCRIPT, "", "non-monotonic latitude coordinate")
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
    end if
    end if
    end if
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
    tmp = var_attrib(output, VARS(vID), CMOR_TABLE(vID))
    delete(output)
    output = tmp
    delete(tmp)

    ; Check fill values
    fill_end = num(ismissing(output))
    if (fill_start .ne. fill_end) then
      error_msg("f", DIAG_SCRIPT, "", "missing values lost during conversion")
    end if
    delete(fill_start)
    delete(fill_end)

    ; Set global attributes
    gAtt = set_global_atts( \
      OBSNAME, TIER, \
      "http://apps.ecmwf.int/datasets/data/interim_full_moda/", \
      "Dee, D. P. et al., Q. J. Roy. Meteor. Soc., doi:10.1002/qj.828, 2011", \
      "")

    ; Outfile
    fout = output_dir_path + \
      "OBS_" + OBSNAME + "_reanaly_1_" + MIPS(vID) + \
      "_" + VARS(vID) + "_" + YEAR1 + "01-" + YEAR2 + "12.nc"

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

end
