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
;    Select "Era Interim Fields":
;      Daily: for daily values
;      Invariant: for time invariant variables (like land-sea mask)
;      Monthly Means of Daily Means: for monthly values
;      Monthly Means of Daily Forecast Accumulation: for accumulated variables
;      like precipitation or radiation fluxes
;    Select "Type of level" (Surface or Pressure levels)
;    Download the data on a single variable and single year basis, and save
;    them as ERA-Interim_<var>_<mean>_YYYY.nc, where <var> is the ERA-Interim
;    variable name and <mean> is either monthly or daily. Further download
;    "land-sea mask" from the "Invariant" data and save it in
;    ERA-Interim_lsm.nc.
;    It is also possible to download data in an automated way, see:
;      https://confluence.ecmwf.int/display/WEBAPI/Access+ECMWF+Public+Datasets
;      https://confluence.ecmwf.int/display/WEBAPI/Python+ERA-interim+examples
;    A registration is required for downloading the data.
;
; Caveats
;    Make sure to select the right steps for accumulated fluxes, see:
;      https://confluence.ecmwf.int/pages/viewpage.action?pageId=56658233
;      https://confluence.ecmwf.int/display/CKB/ERA-Interim%3A+monthly+means
;    for a detailed explanation.
;    The data are updated regularly: recent years are added, but also the past
;    years are sometimes corrected. To have a consistent timeseries, it is
;    therefore recommended to download the full timeseries and not just add
;    new years to a previous version of the data.
;
; Modification history
;    20190311-A_righ_ma: added surface fluxes.
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
  DIAG_SCRIPT = "cmorize_obs_era_interim.ncl"

  ; Source name
  OBSNAME = "ERA-Interim"

  ; Tier
  TIER = 3

  ; Period
  YEAR1 = 1979
  YEAR2 = 2018

  ; Acceleration of gravity [m s-2]
  G = 9.80665

  ; Variable settings as list [VAR, NAME, CONVERSION, MIP, FREQ]
  ; VAR: selected variable (standard name)
  ; NAME: name in the raw data
  ; RAWUNITS: expected units in the input data
  ; CONVERSION: corresponding conversion factor to CMOR units
  ; MIP: mip
  ; FREQ: frequency
  VLIST = \
    [/[/"sftlf",  "lsm",                            100.,   "fx",   "fx"/], \
      [/"pr",     "tp",                             1.e3,   "day",  "day"/], \
      [/"psl",    "msl",                            1.,     "day",  "day"/], \
      [/"tas",    "t2m",                            1.,     "day",  "day"/], \
      [/"tasmax", "mx2t",                           1.,     "day",  "day"/], \
      [/"tasmin", "mn2t",                           1.,     "day",  "day"/], \
      [/"zg",     "z",                              1. / G, "day",  "day"/], \
      [/"clivi",  "p57.162",                        1.,     "Amon", "mon"/], \
      [/"clt",    "tcc",                            100.,   "Amon", "mon"/], \
      [/"clwvi",  (/"p56.162", "p57.162"/),         1.,     "Amon", "mon"/], \
      [/"hfds",   (/"ssr", "str", "slhf", "sshf"/), 1.,     "Omon", "mon"/], \
      [/"hur",    "r",                              1.,     "Amon", "mon"/], \
      [/"hus",    "q",                              1.,     "Amon", "mon"/], \
      [/"pr",     "tp",                             1.e3,   "Amon", "mon"/], \
      [/"prw",    "tcwv",                           1.,     "Amon", "mon"/], \
      [/"ps",     "sp",                             1.,     "Amon", "mon"/], \
      [/"psl",    "msl",                            1.,     "Amon", "mon"/], \
      [/"ta",     "t",                              1.,     "Amon", "mon"/], \
      [/"tas",    "t2m",                            1.,     "Amon", "mon"/], \
      [/"tauu",   "iews",                           1.,     "Amon", "mon"/], \
      [/"tauv",   "inss",                           1.,     "Amon", "mon"/], \
      [/"tos",    "sst",                            1.,     "Omon", "mon"/], \
      [/"ts",     "skt",                            1.,     "Amon", "mon"/], \
      [/"ua",     "u",                              1.,     "Amon", "mon"/], \
      [/"va",     "v",                              1.,     "Amon", "mon"/], \
      [/"wap",    "w",                              1.,     "Amon", "mon"/], \
      [/"zg",     "z",                              1. / G, "Amon", "mon"/], \
      [""]/]

  ; Type
  TYPE = "reanaly"

  ; Version
  VERSION = "1"

  ; Global attributes
  SOURCE =  "http://apps.ecmwf.int/datasets/data/interim_full_moda/"
  REF = "Dee, D. P. et al., Q. J. Roy. Meteor. Soc., doi:10.1002/qj.828, 2011"
  COMMENT = ""

end
; #############################################################################
loadscript(getenv("esmvaltool_root") + \
           "/esmvaltool/cmorizers/obs/interface.ncl")

begin

  ; Loop over variables
  do vv = 0, ListCount(VLIST) - 2

    ; Setting for the current variable
    clist = VLIST[vv]
    VAR = clist[0]
    NAME = clist[1]
    CONVERSION = clist[2]
    MIP = clist[3]
    FREQ = clist[4]
    CMOR_TABLE = getenv("cmor_tables") + "/cmip5/Tables/CMIP5_" + MIP
    delete(clist)

    log_info("Processing " + VAR + " (" + MIP + ")")

    l_day = False
    l_fx = False
    if (isStrSubset(MIP, "mon")) then
      freqkey = "monthly"
    end if
    if (isStrSubset(MIP, "day")) then
      freqkey = "daily"
      l_day = True
    end if
    if (isStrSubset(MIP, "fx")) then
      l_fx = True
    end if

    ; Create timeseries
    do yy = YEAR1, YEAR2

      ; fx variables are time invariant
      if (l_fx) then
        fname = input_dir_path + "ERA-Interim_" + NAME(0) + ".nc"
        f = addfile(fname, "r")
      else
        fname = input_dir_path + \
          "ERA-Interim_" + NAME(0) + "_" + freqkey + "_" + yy + ".nc"
        f = addfile(fname, "r")
      end if

      ; Rank of the input data
      rank = dimsizes(getfilevardims(f, NAME(0)))

      ; For daily data, split in monthly files, otherwise yearly
      if (l_day) then
        mend = 12
      else
        mend = 1  ; dummy loop, the whole year will be processed
      end if

      do mm = 1, mend

        ; Read variable
        if (l_day) then

          ; Check time-index
          time = f->time
          if (any(VAR.eq.(/"tasmin", "tasmax", "pr"/))) then
            time = time - 1  ; Shift back 1h to include 12:00 + step=12h
          end if
          date = cd_calendar(time, 0)
          idxt = ind(date(:, 1).eq.mm)

          ; Extract current month
          if (rank.eq.4) then
            xx = f->$NAME(0)$(idxt, :, :, :)
          elseif (rank.eq.3) then
            xx = f->$NAME(0)$(idxt, :, :)
          else
            error_msg("f", DIAG_SCRIPT, "", "invalid rank")
          end if
          start_date = yy + sprinti("%0.2i", mm) + "01"
          end_date = yy + sprinti("%0.2i", mm) + \
            sprinti("%0.2i", days_in_month(yy, mm))

          ; Shift back 1h to include 12:00 + step=12h
          if (any(VAR.eq.(/"tasmin", "tasmax", "pr"/))) then
            xx&time = xx&time - 1
          end if
          delete(date)
          delete(time)
          delete(idxt)

        else

          xx = f->$NAME(0)$
          start_date = yy + "0101"
          end_date = yy + "1231"

        end if

        ; Unpack variable according to metadata information
        output = short2flt(xx)
        delete(xx)

        ; Multiple input case
        if (dimsizes(NAME).gt.1) then

          if (l_day) then
            error_msg("f", DIAG_SCRIPT, "", "multiple input for daily " + \
                      "values not implemented")
          end if

          do iv = 1, dimsizes(NAME) - 1
            fname2 = input_dir_path + \
              "ERA-Interim_" + NAME(iv) + "_" + freqkey + "_" + yy + ".nc"
            f2 = addfile(fname2, "r")
            xx2 = f2->$NAME(iv)$
            output = output + short2flt(xx2)
            delete(xx2)
            delete(f2)
            delete(fname2)
          end do

        end if

        ; Calculate daily mean
        if (l_day) then

          if (VAR.eq."tasmin") then
            dmean = calculate_daily_values(output, "min", 0, False)
          elseif (VAR.eq."tasmax") then
            dmean = calculate_daily_values(output, "max", 0, False)
          elseif (VAR.eq."pr") then
            dmean = calculate_daily_values(output, "sum", 0, False)
          else
            dmean = calculate_daily_values(output, "avg", 0, False)
          end if
          delete(output)
          output = dmean
          delete(dmean)

        end if

        ; Time-invariant variables
        if (l_fx) then
          output := output(0, :, :)
          rank = 2
        end if

        ; Convert units
        output = output * CONVERSION

        ; Special case: accumulated fluxes
        if (any(VAR.eq.(/"pr", "hfds"/))) then  ; [X] --> [X s-1]
          if (l_day) then
            output = output / 24. / 3600.
          else
            locy = new(12, integer)
            locy = yy
            dm = conform(output, days_in_month(locy, ispan(1, 12, 1)), 0)
            output = output / 24. / 3600. / dm
            delete(locy)
            delete(dm)
          end if
        end if

        ; Format coordinates
        if (rank.eq.4) then
          output!0 = "time"
          output!1 = "plev"
          output!2 = "lat"
          output!3 = "lon"
          output&plev = output&plev * 100  ; [mb] --> [Pa]
        elseif (rank.eq.3) then
          output!0 = "time"
          output!1 = "lat"
          output!2 = "lon"
        elseif (rank.eq.2) then
          output!0 = "lat"
          output!1 = "lon"
        else
          error_msg("f", DIAG_SCRIPT, "", "invalid rank")
        end if
        format_coords(output, start_date, end_date, FREQ)

        ; Set variable attributes
        tmp = format_variable(output, VAR, CMOR_TABLE)
        delete(output)
        output = tmp
        delete(tmp)

        ; Calculate coordinate bounds
        bounds = guess_coord_bounds(output, FREQ)

        ; Set global attributes
        gAtt = set_global_atts(OBSNAME, TIER, SOURCE, REF, COMMENT)

        ; Output file
        if (l_fx) then
          DATESTR = ""
        elseif (l_day) then
          DATESTR = \
            yy + sprinti("%0.2i", mm) + "-" + yy + sprinti("%0.2i", mm)
        else
          DATESTR = yy + "01-" + yy + "12"
        end if
        if (DATESTR.eq."") then
          fout = output_dir_path + \
            str_join((/"OBS", OBSNAME, TYPE, VERSION, MIP, VAR/), "_") + ".nc"

        else
          fout = output_dir_path + \
            str_join((/"OBS", OBSNAME, TYPE, VERSION, \
                       MIP, VAR, DATESTR/), "_") + ".nc"
        end if

        ; Add height coordinate to tas variable (required by the new backend)
        if (any(VAR.eq.(/"tas", "tasmin", "tasmax"/))) then
          output@coordinates = "height"
        end if

        ; Write variable
        write_nc(fout, VAR, output, bounds, gAtt)
        delete(gAtt)
        delete(output)
        delete(bounds)

        ; Add height coordinate to tas variable (required by the new backend)
        if (any(VAR.eq.(/"tas", "tasmin", "tasmax"/))) then
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

    end do  ; year

    delete(NAME)

  end do  ; variable

end
