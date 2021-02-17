; #############################################################################
; YOD - Calculate year of near disappearance of (Arctic) sea ice
; Author: Senftleben, Daniel (DLR, Germany)
; CRESCENDO project
; #############################################################################
; Description
; Calculates the year of near-disappearance of (Arctic) sea ice (YOD). YOD is
; defined as the first of five consecutive years in which the minimum Arctic
; sea ice extent drops below one million km2. The diagnostic creates a time
; series plot of (September Arctic) sea ice extent for each model given in the
; recipe and adds three multi-model statistics: the mean, the standard
; deviation and the YOD. It optionally accepts a list of pre-determined model
; weights and adds to the time series the three weighted multi-model
; statistics.
;
; Required diag_script_info attributes
;   month: selected month (1, 2, ..., 12) or annual mean ("A")
;   region: region to be analyzed ( "Arctic" or "Antarctic")
;
; Optional diag_script_info attributes
;   fill_pole_hole: fill observational hole at North pole (default: False)
;   wgt_file: netCDF containing pre-determined model weights
;
; Caveats
;   - CMOR units expected: areacello in m2, sic in %
;
; Modification history
;    20191030-lauer_axel: rewritten for v2.0
;    20190211-senftleben_daniel: implemented.
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/../interface_scripts/logging.ncl"

load "$diag_scripts/seaice/seaice_aux.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  ; Get metadata items
  att = True
  att@mip = "OImon"
  info = select_metadata_by_atts(input_file_info, att)  ; variable
  var0 = info[0]@short_name
  datasets = metadata_att_as_array(info, "dataset")
  ensembles = metadata_att_as_array(info, "ensemble")
  exps = metadata_att_as_array(info, "exp")
  infiles = metadata_att_as_array(info, "filename")
  dim_MOD = dimsizes(datasets)
  att@mip = "fx"
  areainfo = select_metadata_by_atts(input_file_info, att)  ; area
  areadatasets = metadata_att_as_array(areainfo, "dataset")
  areafiles = metadata_att_as_array(areainfo, "filename")
  delete(att)

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"month", "region"/))

  ; Check fill_pole_hole option
  if (.not.isatt(diag_script_info, "fill_pole_hole")) then
    diag_script_info@fill_pole_hole = False
  else
    if (diag_script_info@fill_pole_hole) then
      if (diag_script_info@region.eq."Antarctic") then
        diag_script_info@fill_pole_hole = False
        error_msg("f", DIAG_SCRIPT, "", "No pole hole at Antarctic " + \
                  "- setting diag_script_info@fill_pole_hole to False")
      end if
    end if
  end if

  ; Create output plot directory
  if (diag_script_info@write_plots) then
    plot_dir = config_user_info@plot_dir
    system("mkdir -p " + plot_dir)
  end if

  ; Create output directory (work)
  system("mkdir -p " + config_user_info@work_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type

  ; Determine time range
  start_years = metadata_att_as_array(info, "start_year")
  end_years   = metadata_att_as_array(info, "end_year")

  ; Check diagnostic attributes
  region = diag_script_info@region
  if (isatt(diag_script_info, "wgt_file")) then
    if (isfilepresent(diag_script_info@wgt_file)) then
      wgt_flag = True  ; do both, weighted and unweighted
      weights = ncdf_read(diag_script_info@wgt_file, "weight")  ; load weights
    else
      error_msg("f", DIAG_SCRIPT, "", "Weights file not found: " + \
                diag_script_info@wgt_file)
    end if
  else
    wgt_flag = False
  end if

  smonths = (/"January", "February", "March", "April", \
              "May", "June", "July", "August", \
              "September", "October", "November", "December"/)

  if (diag_script_info@month .eq. "A") then
    timestring = "annual"
  else
    timestring = smonths(toint(diag_script_info@month) - 1)
  end if

  all_rcps   = exps(str_match_ind_ic(exps, "rcp"))
  which_rcps = get_unique_values(all_rcps)

  YOD     = new((/dimsizes(which_rcps), 2, 2/), float)
  YOD!0   = "rcp"
  YOD!1   = "wgt"
  YOD!2   = "val"
  YOD&rcp = which_rcps
  YOD&wgt = (/"unweighted", "weighted"/)
  YOD&val = (/"avg", "std"/)

  ; for writing provenance

  if (diag_script_info@region.eq."Arctic") then
    domain = "nh"
  else if (diag_script_info@region.eq."Antarctic") then
    domain = "sh"
  end if
  end if

  do ircp = 0, dimsizes(which_rcps) - 1

    ; initialize arrays
    rcp_ind          = ind(exps.eq.which_rcps(ircp))
    rcp_models       = datasets(rcp_ind)
    YOD_tmp          = new((/dimsizes(rcp_models)/), float)
    YOD_tmp!0        = "dataset"
    YOD_tmp&dataset  = rcp_models
    years            = ispan(min(toint(start_years(rcp_ind))), \
                             max(toint(end_years(rcp_ind))), 1)
    SSIE_all         = new((/dimsizes(YOD_tmp), dimsizes(years)/), float)
    SSIE_all!0       = "dataset"
    SSIE_all!1       = "year"
    SSIE_all&dataset = rcp_models
    SSIE_all&year    = years

    ; check if weights are consistent
    if (wgt_flag) then
      if (dimsizes(weights).ne.dimsizes(YOD_tmp&dataset)) then
        error_msg("f", diag_script, "", "weights do not fit to datasets " + \
                  "in recipe for " + which_rcps(ircp))
      else  ; make sure each weight fits to the correct dataset
        sort_ind = new(dimsizes(YOD_tmp&dataset), integer)
        do imod = 0, dimsizes(weights)-1
          sort_ind(imod) = ind(weights&dataset.eq.YOD_tmp&dataset(imod))
        end do
        if (any(ismissing(sort_ind))) then
          error_msg("f", diag_script, "", "weights do not fit to " + \
                    "datasets in recipe for " + which_rcps(ircp))
        else
          weights = weights(sort_ind)
        end if
      end if
    end if

    ; loop over datasets
    do ii = 0, dimsizes(rcp_ind) - 1
      imod = rcp_ind(ii)
      log_info("Processing dataset " + datasets(imod))

      sic = read_data(info[imod])

      ; check grid type
      l_irreg = True
      if (sic!1.eq."lat" .and. sic!2.eq."lon") then
        if (dimsizes(dimsizes(sic&lat)).eq.1 .and. \
            dimsizes(dimsizes(sic&lon)).eq.1) then
          l_irreg = False
        end if
      end if

      ; look for area variable
      if (any(datasets(imod).eq.areadatasets)) then

        ; read areacello
        log_info("Using variable areacello for the area")
        area = read_data(areainfo[imod])
        locinfo = areainfo[imod]
        locinfo@short_name = "lat"  ; read auxiliary latitude coordinate
        arealat = read_data(locinfo)
        delete(locinfo)

        ; check grid consistency
        if (.not.isconform(sic(0, :, :), area)) then
          if (l_irreg) then
            error_msg("f", DIAG_SCRIPT, "", "areacello not conformal")
          else
            delete(area)
            delete(arealat)
            log_info("areacello has inconsistent dimensions, " + \
                     "calculating area for a regular grid instead")
            area = map_area(sic&lat, sic&lon)
            arealat = area&lat
          end if
        end if
      else
        if (.not.l_irreg) then
          log_info("Calculating area for a regular grid")
          area = map_area(sic&lat, sic&lon)
          arealat = area&lat
        else
          error_msg("f", DIAG_SCRIPT, "", "cannot calculate area")
        end if
      end if

      ; calcluate sea ice extent (SIE)
      SIE  = sea_ice_area(sic, area, arealat, l_irreg, \
                          diag_script_info@region, \
                          start_years(imod), end_years(imod), \
                          diag_script_info@month, "extent", \
                          diag_script_info@fill_pole_hole)

      delete([/sic, area, arealat/])

      SSIE = SIE
      if (SSIE!0.eq."time") then  ; need 'year'
        utc_date = cd_calendar(SSIE&time, 0)
        SSIE!0 = "year"
        SSIE&year = utc_date(:, 0)
        delete(utc_date)
      end if

      ; calculate ensemble mean
      if (imod.ne.max(rcp_ind)) then
        if ((datasets(imod) + exps(imod)) .eq. \
            (datasets(imod + 1) + exps(imod + 1))) then
          if (.not.isvar("store")) then
            store = new((/1, dimsizes(SSIE)/), typeof(SSIE))
            store(0, :) = SSIE
            continue
          else
            tmp = new((/1, dimsizes(SSIE)/), typeof(SSIE))
            tmp(0, :) = SSIE
            store := array_append_record(store, tmp, 0)
            continue
          end if
        else
          if (isvar("store")) then
            SSIE_em = dim_avg_n_Wrap(store, 0)
            delete(store)
          else
            SSIE_em = SSIE
          end if
        end if
      else
        if (isvar("store")) then
          SSIE_em = dim_avg_n_Wrap(store, 0)
        else
          SSIE_em = SSIE
        end if
      end if

      ; calculate first year of disappearance (YOD)
      ;  = 1st year when 5-yr smoothed SSIE drops below 10^6 km^2
      ;    for more than 5 consecutive years
      mod_ind = ind(datasets(imod).eq.YOD_tmp&dataset)
      YOD_tmp(mod_ind) = calc_yod(SSIE_em, 1., 5)
      if (YOD_tmp(mod_ind).lt.min(tofloat(start_years(rcp_ind))) .or. \
          YOD_tmp(mod_ind).gt.max(tofloat(end_years(rcp_ind)))) then
        ; make sure missing values are recognized
        YOD_tmp(mod_ind) = YOD_tmp@_FillValue
      end if

      ; store
      SSIE_all(mod_ind, :) = (/SSIE_em/)

      ; clean up
      delete([/SIE, SSIE, SSIE_em, mod_ind/])
    end do  ; imod

    ; calculate YOD from multi-model mean (MMM)
    MMM_uwt = dim_avg_n_Wrap(SSIE_all, 0)
    STD_uwt = dim_stddev_n_Wrap(SSIE_all, 0)
    STD_wgt = STD_uwt  ; initialize
    MMM_wgt = MMM_uwt
    if (wgt_flag) then
      MMM_wgt = dim_avg_wgt_n_Wrap(SSIE_all, weights, 1, 0)
      do iyear = 0, dimsizes(SSIE_all&year) - 1
        STD_wgt(iyear) = dim_stddev_wgt_Wrap(SSIE_all(:, iyear), weights, 1)
      end do
    else
      STD_wgt = STD_wgt@_FillValue
      MMM_wgt = MMM_wgt@_FillValue
    end if

    YOD_MMM_uwt   = calc_yod(MMM_uwt, 1., 5)
    YOD_MMM_uwt_2 = calc_yod(MMM_uwt, 2., 1)

    if (wgt_flag) then
      YOD_MMM_wgt   = calc_yod(MMM_wgt, 1., 5)
      YOD_MMM_wgt_2 = calc_yod(MMM_wgt, 2., 1)
    else
      YOD_MMM_wgt   = 0.0
      YOD_MMM_wgt_2 = 0.0
    end if

    ; =========================================================================
    ; plotting section
    ; =========================================================================

    ; --------------------------------------------------------------------
    ; timeseries of models (ensemble means), MMM (wgt/uwt), stdd (wgt/uwt)
    ; --------------------------------------------------------------------
    outfile_ts = "timeseries_" + which_rcps(ircp)
    wks_ts = get_wks(0, DIAG_SCRIPT, outfile_ts)
    res_ts = True
    res_ts@gsnDraw  = False
    res_ts@gsnFrame = False
    res_ts@gsnMaximize = True
    res_ts@vpHeightF = 0.4
    res_ts@vpWidthF = 0.8
    res_ts@xyDashPattern = 1
    res_ts@trXMinF = min(SSIE_all&year)
    res_ts@trXMaxF = max(SSIE_all&year)
    res_std_uwt = True
    copy_VarMeta(res_ts, res_std_uwt)
    res_std_wgt = True
    copy_VarMeta(res_ts, res_std_wgt)
    res_mmm = True
    copy_VarMeta(res_ts, res_mmm)
    res_ts@tiMainString = which_rcps(ircp)
    plot_ts = gsn_csm_xy(wks_ts, SSIE_all&year, SSIE_all, res_ts)
    res_std_uwt@gsnXYFillColors = (/0., 0., 0., .07/)  ; transparent gray
    res_std_uwt@xyLineColor = "red"
    res_std_uwt@xyDashPattern = 0
    res_std_uwt@xyLineThicknessF = 2.
    plot_std_uwt = gsn_csm_xy(wks_ts, SSIE_all&year, (/MMM_uwt + STD_uwt, \
                              MMM_uwt - STD_uwt/), res_std_uwt)
    overlay(plot_ts, plot_std_uwt)
    if (wgt_flag) then
      res_std_wgt@gsnXYFillColors = (/0., 0., 0., .07/)  ; transparent gray
      res_std_wgt@xyLineColor = "blue"
      res_std_wgt@xyDashPattern = 0
      res_std_wgt@xyLineThicknessF = 2.
      plot_std_wgt = gsn_csm_xy(wks_ts, SSIE_all&year, (/MMM_wgt + STD_wgt, \
                                MMM_wgt - STD_wgt/), res_std_wgt)
      overlay(plot_ts, plot_std_wgt)
    end if
    res_mmm@xyLineColors = (/"red", "blue"/)
    res_mmm@xyLineThicknessF = 4.
    res_mmm@xyDashPattern = 0
    if (wgt_flag) then
      plot_mmm = gsn_csm_xy(wks_ts, SSIE_all&year, \
                            (/MMM_uwt, MMM_wgt/), res_mmm)
    else
      plot_mmm = gsn_csm_xy(wks_ts, SSIE_all&year, MMM_uwt, res_mmm)
    end if
    overlay(plot_ts, plot_mmm)
    res_yod = True
    res_yod@gsLineThicknessF = 3.
    res_yod@gsLineColor = "red"
    res_yod@gsLineDashPattern = 0
    yod_line_uwt = gsn_add_polyline(wks_ts, plot_ts, (/YOD_MMM_uwt, \
                                    YOD_MMM_uwt/), (/0., 12./), res_yod)
    if (wgt_flag) then
      res_yod@gsLineColor = "blue"
      res_yod@gsLineDashPattern = 0
      yod_line_wgt = gsn_add_polyline(wks_ts, plot_ts, (/YOD_MMM_wgt, \
                                      YOD_MMM_wgt/), (/0., 12./), res_yod)
    end if
    res_yod@xyLineThicknessF  = 0.1
    res_yod@gsLineColor       = "black"  ; 1 million km2 reference line
    res_yod@gsLineDashPattern = 0
    yod_ref_line = gsn_add_polyline(wks_ts, plot_ts, (/min(SSIE_all&year), \
                                    max(SSIE_all&year)/), (/1., 1./), res_yod)
    draw(plot_ts)
    frame(wks_ts)

    ; store
    YOD(ircp, 0, 0) = avg(YOD_tmp)
    YOD(ircp, 0, 1) = stddev(YOD_tmp)
    if (wgt_flag) then
      YOD(ircp, 1, 0) = dim_avg_wgt(YOD_tmp, weights, 1)
      YOD(ircp, 1, 1) = dim_stddev_wgt_Wrap(YOD_tmp, weights, 1)
    end if

    ; =========================================================================
    ; Provenance
    ; =========================================================================

    SSIE_all@var = "SSIE"
    SSIE_all@diag_script = DIAG_SCRIPT
    outf = config_user_info@work_dir + outfile_ts + ".nc"

    MMM_uwt@var = "SSIE_MMM_uwt"
    MMM_uwt@diag_script = DIAG_SCRIPT
    STD_uwt@var = "SSIE_MMM_STD_uwt"
    STD_uwt@diag_script = DIAG_SCRIPT

    nc_ext = ncdf_write(SSIE_all, outf)
    outf@existing = "append"
    nc_ext = ncdf_write(MMM_uwt, outf)
    nc_ext = ncdf_write(STD_uwt, outf)

    if (wgt_flag) then
      MMM_wgt@var = "SSIE_MMM_wgt"
      MMM_wgt@diag_script = DIAG_SCRIPT
      STD_wgt@var = "SSIE_MMM_STD_wgt"
      STD_wgt@diag_script = DIAG_SCRIPT

      nc_ext = ncdf_write(MMM_wgt, outf)
      nc_ext = ncdf_write(STD_wgt, outf)
    end if

    caption = "Time series of " + timestring + " " + diag_script_info@region \
              + " sea ice extent for individual models and multi-model mean," \
              + " and multi-model standard deviation of year of dissaperance" \
              + " for scenario " + which_rcps(ircp) + "."

    log_provenance(nc_ext, wks_ts@fullname, caption, \
                   "other", domain, "times", \
                   (/"senftleben_daniel", "lauer_axel"/), \
                   "", infiles(rcp_ind))

    ; clean up
    delete([/rcp_models, SSIE_all, YOD_tmp, MMM_uwt, MMM_wgt, STD_uwt, \
           STD_wgt, YOD_MMM_uwt, YOD_MMM_wgt, res_yod, outf, nc_ext, wks_ts/])
  end do  ; ircp

  ; ---
  ; YOD
  ; ---

  outfile  = "YOD"
  wks = get_wks(0, DIAG_SCRIPT, outfile)
  if (wgt_flag) then
    factor = 2
  else
    factor = 1
  end if
  ds = dimsizes(which_rcps) * factor
  if (ds .gt. 1) then
    x = fspan(1, ds, dimsizes(which_rcps) * factor)
  else
    x = 1
  end if
  if (factor .gt. 1) then
    x(::2) = x(::2) + 0.5
  end if
  y        = new((/dimsizes(YOD&rcp) * factor, 5/), float)
  XBLabels = new((/dimsizes(YOD&rcp) * factor/), string)
  do ircp = 0, dimsizes(YOD&rcp) - 1
    ; unweighted
    y(ircp * factor,     0)     = YOD(ircp, 0, 0) - YOD(ircp, 0, 1)
    y(ircp * factor,     1:3)   = YOD(ircp, 0, 0)
    y(ircp * factor,     4)     = YOD(ircp, 0, 0) + YOD(ircp, 0, 1)
    XBLabels(ircp * factor)     = YOD&rcp(ircp) + "_uwt"
    if (wgt_flag) then
      ; weighted
      y(ircp * factor + 1, 0)     = YOD(ircp, 1, 0) - YOD(ircp, 1, 1)
      y(ircp * factor + 1, 1:3)   = YOD(ircp, 1, 0)
      y(ircp * factor + 1, 4)     = YOD(ircp, 1, 0) + YOD(ircp, 1, 1)
      XBLabels(ircp * factor + 1) = YOD&rcp(ircp) + "_wgt"
    end if
  end do
  res  = True
  res@tmXBLabels      = XBLabels
  res@tmXBLabelAngleF = 90
  res@tiYAxisString = "Year of Disappearance"
  plot = boxplot(wks, x, y, False, res, False)
  draw(plot)
  frame(wks)

  ; ===========================================================================
  ; Provenance
  ; ===========================================================================

  YOD@var = "yod"
  YOD@diag_script = DIAG_SCRIPT

  outf = config_user_info@work_dir + outfile + ".nc"

  nc_ext2 = ncdf_write(YOD, outf)

  caption = "Year of dissaperance of " + timestring + " " \
            + diag_script_info@region + " sea ice for scenarios " \
            + tostring(which_rcps) + "."

  log_provenance(nc_ext2, wks@fullname, caption, \
                 "other", domain, "bar", \
                 (/"senftleben_daniel", "lauer_axel"/), \
                 "", infiles)

  leave_msg(DIAG_SCRIPT, "")

end
