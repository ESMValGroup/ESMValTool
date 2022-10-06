; #############################################################################
; PLOTS TIME SERIES OF PRECIPITATION ANOMALY
; Authors: Lisa Bock (DLR, Germany)
; #############################################################################
;
; Description
;    This script is for plotting a time series of the global and regional
;    anomaly of precipitation.
;
; Required diag_script_info attributes (diagnostic specific)
;    panels: list of variables plotted in each panel
;    start_year: start of time coordinate
;    end_year: end of time coordinate
;
; Optional diag_script_info attributes (diagnostic specific)
;    anomaly: true if anomaly should be calculated
;    ref_start: start year of reference period for anomalies
;    ref_end: end year of reference period for anomalies
;    ref_mask: if true, model fields will be masked by reference fields
;    region: name of domain
;    plot_units: variable unit for plotting
;    header: if true, region name as header
;    stat: statistics for multi model nc-file (MinMax,5-95,10-90)
;    y_min: set min of y-axis
;    y_max: set max of y-axis
;
; Caveats
;
; Modification history
;    20190911-bock_lisa: written
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"

load "$diag_scripts/ipcc_ar6/percentiles.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"panels", "start_year", \
                                           "end_year"/))

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "anomaly", True)
  set_default_att(diag_script_info, "ref_mask", False)
  set_default_att(diag_script_info, "region", "Global")
  set_default_att(diag_script_info, "header", False)
  set_default_att(diag_script_info, "stat", "MinMax")

  ; Create output directories
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  out_path = config_user_info@work_dir
  system("mkdir -p " + out_path)

  if (diag_script_info@anomaly) then
    anom = "anom"
  end if

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Determine time range
  start_year = diag_script_info@start_year
  end_year = diag_script_info@end_year
  all_years = ispan(start_year, end_year, 1)

  ; Create time coordinate
  ntime = end_year - start_year + 1
  time = new(ntime, integer)
  time = ispan(start_year, end_year, 1)

  ; *************************************
  ; General settings for the plot
  ; *************************************

  ; Define workstation
  outfile = config_user_info@plot_dir + "precip_" + anom + "_" \
                                      + start_year + "-" + end_year
  wks = gsn_open_wks(file_type, outfile)

  ; Set resources
  res = True
  res@gsnDraw   = False
  res@gsnFrame  = False

  res@vpXF      = 0.05
  res@vpYF      = 0.7
  res@vpHeightF = 0.4
  res@vpWidthF  = 0.7
  res@pmLegendDisplayMode = "Never"
  res@tmYLLabelFontHeightF = 0.03
  res@tmXBLabelFontHeightF = 0.03
  res@tiXAxisFontHeightF = 0.018
  res@tiYAxisFontHeightF = 0.018

  xmin = start_year
  xmax = end_year
  res@trXMinF = xmin
  res@trXMaxF = xmax
  res@tmXBMode = "Explicit"
  if (xmax - xmin.gt.20) then
    res@tmXBValues = ispan(xmin, xmax, 10)
    res@tmXBLabels = ispan(xmin, xmax, 10)
    res@tmXBMinorValues =  ispan(xmin, xmax, 5)
  else
    res@tmXBValues = ispan(xmin, xmax, 5)
    res@tmXBLabels = ispan(xmin, xmax, 5)
    res@tmXBMinorValues =  ispan(xmin, xmax, 1)
  end if

  res@tmXBMode         = "Manual"
  res@tmXBTickSpacingF = 10

  if (isatt(diag_script_info, "y_min")) then
    res@trYMinF = diag_script_info@y_min
  end if
  if (isatt(diag_script_info, "y_max")) then
    res@trYMaxF = diag_script_info@y_max
  end if
  res@tiYAxisOn = True
  res@tiXAxisString = ""

  resb = True
  resb@vpWidthF  = 0.2
  resb@trXMinF = -0.8
  resb@trXMaxF = 2.5
  resb@trYMinF = -0.002
  resb@trYMaxF = 0.003
  resb@tmXTOn = False
  resb@tmXBOn = False
  resb@tmYRBorderOn = True
  resb@tmXTBorderOn = True
  resb@tmYROn = True
  resb@tmYUseLeft = False
  resb@tmYRLabelsOn = True
  resb@tmYRLabelFontHeightF = 0.03
  resb@tiYAxisOn = True
  resb@tiYAxisSide   = "Right"
  resb@tiYAxisAngleF = 90.
  resb@tiYAxisFontHeightF = 0.018

  llres = True
  llres@gsLineThicknessF = 2.5

  cmap = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                            + "ipcc-ar6_damip_line.rgb")
  opti = True
  opti@boxColors = (/"peru", "darkgreen"/)

  ; ------------------------------------------------------------------
  ; Panel settings
  ; ------------------------------------------------------------------

  nvar = ListCount(variable_info)

  panels = diag_script_info@panels
  tmp := dimsizes(panels)
  if(dimsizes(tmp).gt.1) then
    npanels = tmp(0)
    nexp = tmp(1)
  else
    npanels = 1
    nexp = tmp(0)
  end if
  delete(tmp)
  plots = new((/npanels/), graphic)

  ; ------------------------------------------------------------------
  ; Loop over panels
  ; ------------------------------------------------------------------

  do ivar = 0, npanels-1

    ; create array for statistics
    model_arr_stat = new((/3, 3, ntime/), double)
    model_arr_stat!0 = "exp"
    model_arr_stat!1 = "stat"
    model_arr_stat!2 = "time"
    model_arr_stat&exp = (/"exp1", "exp2", "ref"/)
    model_arr_stat&time = time
    model_arr_stat&stat = (/"mean", "min", "max"/)

    model_arr_ref = new((/3, ntime/), double)
    model_arr_ref!0 = "ref"
    model_arr_ref!1 = "time"
    model_arr_ref&ref = (/"ref1", "ref2", "ref3"/)

    model_arr_box = new((/2, 5/), double)
    model_arr_box!0 = "exp"
    model_arr_box!1 = "stat"
    model_arr_box&exp = (/"exp1", "exp2"/)

    do iexp = 0, 1

      if(npanels .gt. 1) then
        var = panels(ivar, iexp)
      else
        var = panels(iexp)
      end if

      log_info("processing " + var)

      ; Retrieve metadata item
      atts = True
      atts@variable_group = var
      info_items = select_metadata_by_atts(input_file_info, atts)
      delete(atts)

      var0 = info_items[0]@short_name
      exp0 = info_items[0]@exp
      datasetnames := metadata_att_as_array(info_items, "dataset")
      ensembles := metadata_att_as_array(info_items, "ensemble")
      dim_MOD = ListCount(info_items)

      ; get multi-model mean index
      mmm_ind = ind(datasetnames .eq. "MultiModelMean")

      ; Set index of the reference dataset and read it
      if (isatt(info_items[0], "reference_dataset")) then
        ref_ind = ind(datasetnames.eq.info_items[0]@reference_dataset)
        A_ref = read_data(info_items[ref_ind])
      end if

      ; Anomaly
      if (isdefined(anom)) then
        if (isatt(diag_script_info, "ref_start") .and. \
            isatt(diag_script_info, "ref_end")) then
          ref_start = diag_script_info@ref_start
          ref_end   = diag_script_info@ref_end
          if (ref_start.lt.start_year) then
            ref_start = start_year
          end if
          if ((ref_start.lt.start_year) .or. (ref_end.gt.end_year) .or. \
              (ref_end.lt.ref_start)) then
            error_msg("f", DIAG_SCRIPT, "", \
                      "period for reference years is not properly defined")
          end if
        else
          error_msg("f", DIAG_SCRIPT, "", \
                    "period for reference years is not defined " + \
                    "(needed for anomaly)")
        end if
      end if

      ; Create model array
      model_arr = new((/dim_MOD, ntime/), double)
      model_arr!0 = "model"
      model_arr!1 = "time"
      model_arr&model = datasetnames
      model_arr&time = time
      model_arr@_FillValue = 1e+20
      model_arr = model_arr@_FillValue

      model_arr_trend := new(dim_MOD, double)
      model_arr_trend!0 = "model"
      model_arr_trend&model = datasetnames
      model_arr_trend@_FillValue = 1e+20
      model_arr_trend = model_arr_trend@_FillValue

      ; Loop over models
      do imod = 0, dim_MOD - 1

        log_info("Process dataset: " + datasetnames(imod))

        ; Read data
        A0 = read_data(info_items[imod])
        dnames = getVarDimNames(A0)

        ; Convert units for plotting (if required)
        if (isatt(diag_script_info, "plot_units")) then
          A0 = convert_units(A0, diag_script_info@plot_units)
        end if

        ; masking with reference dataset
        if (isatt(info_items[0], "reference_dataset")) then
          if (diag_script_info@ref_mask .and. imod .ne. ref_ind) then
            ; Determine start/end year
            start_year = info_items[imod]@start_year
            end_year = info_items[imod]@end_year
            A_ref_mask = time_operations(A_ref, start_year, \
                                         end_year, "extract", "", 0)
            A0 = where(A_ref_mask.eq.A_ref_mask@_FillValue, \
                       A_ref_mask@_FillValue, A0)
            delete(A_ref_mask)
          end if
        end if

        ; Anomaly
        if (anom .eq. "anom") then
          if (ref_start.lt.start_year) then
            ref_start = start_year
          end if
          tmp = time_operations(A0, ref_start, ref_end, "extract", "", 0)
          A0_monavg_ref = time_operations(tmp, ref_start, ref_end, "average", \
                                          "monthlyclim", True)
          delete(tmp)

          ; calculate anomaly for each grid point
          do i = 0, dimsizes(A0&time) - 1
            A0(i, :, :) = A0(i, :, :) - A0_monavg_ref(mod(i, 12), :, :)
          end do

          delete(A0_monavg_ref)

        end if

        ; Calculate time average of dataset
        A0_timavg = time_operations(A0, -1, -1, "average", \
                                    "yearly", True)
        delete(A0)

        ; calculate global mean of complete time series
        procmod = area_operations(A0_timavg, -90., 90., 0., 360., \
                                  "average", True)
        delete(A0_timavg)

        ; Match time coordinate
        date = procmod&year
        idx1 = ind(date(0).eq.model_arr&time)
        idx2 = ind(date(dimsizes(date) - 1).eq.model_arr&time)
        model_arr(imod, idx1:idx2) = (/procmod/)
        if (imod.eq.0) then
          copy_VarAtts(procmod, model_arr)
        end if
        delete(procmod)
        delete(date)

        ; calculate trend
        tmp = dtrend_msg_n(time, model_arr(imod, :), True, True, 0)
        model_arr_trend(imod) = tmp@slope
        delete(tmp)

      end do

      ; *************************************
      ; calculate means over ensemble members
      ; *************************************

      ndim := dimsizes(datasetnames)
      ensemble_assign := new(ndim, integer)
      ensemble_size := new(ndim, integer)

      enum = 1
      temp = datasetnames(0)
      do gg = 0, ndim - 1
        ; does the model name match what's in temp?
        if (temp.eq.datasetnames(gg)) then
          ensemble_assign(gg) = enum  ; if so, assign it the same number
        else
          enum = enum + 1  ; if not, assign it the next number
          ensemble_assign(gg) = enum
          temp = datasetnames(gg)
        end if
      end do
      ensemble_assign@model = datasetnames
      delete(temp)
      ensemble_assign@models = str_join(datasetnames+"/"+ensembles, ", ")

      model_arr_em := model_arr(:max(ensemble_assign)-1, :)
      model_arr_em = model_arr_em@_FillValue
      model_arr_trend_em := model_arr_trend(:max(ensemble_assign)-1)
      model_arr_trend_em = model_arr_em@_FillValue
      datasetnames_em := datasetnames(:max(ensemble_assign)-1)

      do gg = 1, max(ensemble_assign)         ; calculate ensemble means
        wind := ind(ensemble_assign.eq.gg)
        ensemble_size(wind) = dimsizes(wind)
        if (dimsizes(wind).eq.1) then
          model_arr_em(gg-1, :) = (/ model_arr(wind, :) /)
          model_arr_trend_em(gg-1) = (/ model_arr_trend(wind) /)
          datasetnames_em(gg-1) = datasetnames(wind)
        else
          model_arr_em(gg-1, :) = (/ dim_avg_n(model_arr(wind, :), 0) /)
          model_arr_trend_em(gg-1) = (/ dim_avg_n(model_arr_trend(wind), 0) /)
          datasetnames_em(gg-1) = datasetnames(wind(0))
        end if
      end do

      ; --------------------------------------------------------------------
      ; calculate multi model statistics
      ; --------------------------------------------------------------------

      if (isatt(info_items[0], "reference_dataset")) then
        ind_mod := ind(datasetnames_em.ne."MultiModelMean" .and. \
                       datasetnames_em.ne.info_items[0]@reference_dataset \
                       .and. datasetnames_em.ne."GHCN" .and. \
                       datasetnames_em.ne."GPCP-SG" .and. \
                       datasetnames_em.ne."CRU")
      else
        ind_mod := ind(datasetnames_em.ne."MultiModelMean" .and. \
                       datasetnames_em.ne."GHCN" .and. \
                       datasetnames_em.ne."GPCP-SG" .and. \
                       datasetnames_em.ne."CRU")
      end if

      model_arr_stat&exp(iexp) =  var
      weight := new(dimsizes(ind_mod), float)
      weight = 1
      weight@_FillValue = 1e+20
      model_arr_stat(iexp, 0, :) = (/weighted_mean(model_arr_em(ind_mod, :), \
                                    weight)/)
      model_arr_stat(iexp, 1, :) = \
        (/weighted_percentile(model_arr_em(ind_mod, :), weight, 0.05)/)
      model_arr_stat(iexp, 2, :) = \
        (/weighted_percentile(model_arr_em(ind_mod, :), weight, 0.95)/)

      if (iexp .eq. 1 .and. isatt(info_items[0], "reference_dataset")) then
        ind_ref = new(3, integer)
        ind_ref(0) = ind(datasetnames_em .eq. "GHCN")
        ind_ref(1) = ind(datasetnames_em .eq. "GPCP-SG")
        ind_ref(2) = ind(datasetnames_em .eq. "CRU")
        model_arr_ref = (/model_arr_em(ind_ref, :)/)
      end if

      delete(model_arr)
      delete(A_ref)

      ; --------------------------------------------------------------------
      ; Low pass filter
      ; --------------------------------------------------------------------

      ; create filter
      nwt = 9      ; number of weights
      fca = 0.2
      nsigma = 1.
      wgt = filwgts_lanczos(nwt, 0, fca, 0, nsigma)

      ; apply filter
      model_arr_stat(iexp, 0, :) = \
        wgt_runave_n_Wrap(model_arr_stat(iexp, 0, :), wgt, 1, 0)
      model_arr_stat(iexp, 1, :) = \
        wgt_runave_n_Wrap(model_arr_stat(iexp, 1, :), wgt, 1, 0)
      model_arr_stat(iexp, 2, :) = \
        wgt_runave_n_Wrap(model_arr_stat(iexp, 2, :), wgt, 1, 0)
      if (iexp .eq. 0 .and. isatt(info_items[0], "reference_dataset")) then
        model_arr_ref(0, :) = wgt_runave_n_Wrap(model_arr_ref(0, :), wgt, 1, 0)
        model_arr_ref(1, :) = wgt_runave_n_Wrap(model_arr_ref(1, :), wgt, 1, 0)
        model_arr_ref(2, :) = wgt_runave_n_Wrap(model_arr_ref(2, :), wgt, 1, 0)

        model_arr_stat(2, :, :) = (/model_arr_ref(:, :)/)
      end if

      ; --------------------------------------------------------------------
      ; calculate statistics for boxplot of trends
      ; --------------------------------------------------------------------

      model_arr_box&exp(iexp) =  var

      weight := new(dimsizes(ind_mod), float)
      weight = 1
      weight@_FillValue = 1e+20

      model_arr_box(iexp, 2) = \
        weighted_mean(model_arr_trend_em(ind_mod), weight)
      model_arr_box(iexp, 0) = \
        weighted_percentile(model_arr_trend_em(ind_mod), weight, 0.05)
      model_arr_box(iexp, 1) = \
        weighted_percentile(model_arr_trend_em(ind_mod), weight, 0.25)
      model_arr_box(iexp, 3) = \
        weighted_percentile(model_arr_trend_em(ind_mod), weight, 0.75)
      model_arr_box(iexp, 4) = \
        weighted_percentile(model_arr_trend_em(ind_mod), weight, 0.95)

    end do

    model_arr_stat(2, :, :) = (/model_arr_ref(:, :)/)

    ; *************************************
    ; output to NetCDF
    ; *************************************

    out_path1 = out_path + "precip_" + anom + "_" + ivar + ".nc"
    model_arr_stat@ncdf = out_path
    model_arr_stat@diag_script = DIAG_SCRIPT
    model_arr_stat@var = var0
    ncdf_outfile = ncdf_write(model_arr_stat, out_path1)

    ; *************************************
    ; Plotting
    ; *************************************

    cmap = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                              + "ipcc-ar6_damip_line.rgb")

    res@xyDashPatterns    := (/0, 0/)
    res@xyLineThicknesses := (/3, 3/)
    res@xyLineColors := cmap(:, :)

    ymax = max((/max(model_arr_stat), max(model_arr_ref)/))
    ymin = min((/min(model_arr_stat), min(model_arr_ref)/))
    res@trYMaxF := max((/0.2, ymax + 0.05 * (ymax - ymin)/))
    res@trYMinF := min((/-0.2, ymin - 0.05 * (ymax - ymin)/))

    if (isatt(diag_script_info, "header")) then
      res@tiMainString = diag_script_info@header(ivar)
    else
      res@tiMainString = ""
    end if

    if (isatt(diag_script_info, "plot_units")) then
      units = diag_script_info@plot_units
    else
      units = info_items[0]@units
    end if

    varname = "Precipitation"

    if (anom .eq. "anom") then
      res@tiYAxisString = varname + " Anomaly" + " (" + units + ")"
    else
      res@tiYAxisString = varname + " (" + units + ")"
    end if

    plots(ivar) = gsn_csm_xy(wks, time, model_arr_stat(0:1, 0, :), res)

    res@pmLegendDisplayMode = "Never"

    ; ---------------------------------------------------------------------
    ; Add multi model statistics
    ; ---------------------------------------------------------------------
    res_stat = True
    res_stat@gsnXYFillOpacities = 0.2
    cmap = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                              + "ipcc-ar6_damip_shading.rgb")
    copy_VarMeta(res, res_stat)
    delete(res_stat@xyLineColors)
    res_stat@gsnXYFillColors = cmap(0, :)
    res_stat@xyLineColor := cmap(0, :)
    res_stat@xyLineThicknesses := (/1, 1/)
    shading_plot = gsn_csm_xy(wks, time, model_arr_stat(0, 1:2, :), \
                              res_stat)
    overlay(plots(ivar), shading_plot)

    res_stat@gsnXYFillColors = cmap(1, :)
    res_stat@xyLineColor = cmap(1, :)
    shading_plot = gsn_csm_xy(wks, time, model_arr_stat(1, 1:2, :), \
                              res_stat)
    overlay(plots(ivar), shading_plot)

    ; ---------------------------------------------------------------------
    ; add reference datasets
    ; ---------------------------------------------------------------------
    res@xyDashPatterns    := (/0, 1, 2, 3/)
    res@xyLineThicknesses := (/3, 3, 3, 3/)
    res@xyLineColors      := "black"  ; change line color

    ref_p = gsn_csm_xy(wks, time, model_arr_ref, res)  ; create plot
    overlay(plots(ivar), ref_p)

    ; ---------------------------------------------------------------------
    ; Create legend
    ; ---------------------------------------------------------------------
    if (ivar .eq. 0) then
      res_text                    = True          ; text mods desired
      res_text@txFontHeightF      = 0.018         ; change text size
      res_text@txJust             = "CenterLeft"  ; text justification

      res_lines                   = True          ; polyline mods desired
      res_lines@gsLineDashPattern = 0.            ; solid line
      res_lines@gsLineThicknessF  = 4.            ; line thicker
      cmap = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_damip_line.rgb")
      res_lines@gsLineColor = cmap(0, :)
      xx = (/1975, 1981/)
      yy = (/-.150, -.15/)
      dum1 = gsn_add_polyline(wks, plots(ivar), xx, yy, res_lines)
      dum2 = gsn_add_text(wks, plots(ivar), "Anthropogenic + Natural",\
                          1983, -.150, res_text)

      yy = (/-.18, -.18/)
      res_lines@gsLineColor = cmap(1, :)
      dum3 = gsn_add_polyline(wks, plots(ivar), xx, yy, res_lines)
      dum4 = gsn_add_text(wks, plots(ivar), "Natural", 1983, -.18, res_text)

      xx = (/1952, 1960/)
      yy = (/-.12, -.12/)
      res_lines@gsLineColor := "black"  ; change to black
      res_lines@gsLineDashPattern := 0
      dum5 = gsn_add_polyline(wks, plots(ivar), xx, yy, res_lines)
      dum6 = gsn_add_text(wks, plots(ivar), datasetnames(ind_ref(0)), 1962, \
                          -.12, res_text)
      yy = (/-.15, -.15/)
      res_lines@gsLineDashPattern := 1
      dum7 = gsn_add_polyline(wks, plots(ivar), xx, yy, res_lines)
      dum8 = gsn_add_text(wks, plots(ivar), datasetnames(ind_ref(1)), 1962, \
                          -.15, res_text)
      yy = (/-.18, -.18/)
      res_lines@gsLineDashPattern := 2
      dum9 = gsn_add_polyline(wks, plots(ivar), xx, yy, res_lines)
      dum10 = gsn_add_text(wks, plots(ivar), datasetnames(ind_ref(2)), \
                           1962, -.18, res_text)
    end if

    ; ---------------------------------------------------------------------
    ; Add inset with trend boxplot
    ; ---------------------------------------------------------------------

    resb@tiYAxisString = "Precipitation Trend (" + units + "/yr)"

    box = boxplot(wks, (/0, 2/), model_arr_box, opti, resb, llres)

    getvalues box
      "trXMinF"   : bxmin
      "trYMinF"   : bymin
      "trXMaxF"   : bxmax
      "trYMaxF"   : bymax
    end getvalues

    xbox1 = (/bxmin, bxmax, bxmax, bxmin, bxmin/)
    ybox1 = (/bymin, bymin, bymax, bymax, bymin/)

    gnres = True
    gnres@gsLineColor = "black"
    gnres@gsLineThicknessF = 2.0
    gnres@gsFillColor     = "Gray70"
    gnres@gsFillOpacityF = 0.1
    gnres@tfPolyDrawOrder = "PreDraw"
    plots@$unique_string("box")$ = gsn_add_polygon(wks, box, xbox1, ybox1, \
                                                   gnres)

    ; add reference dataset
    gnres@gsLineThicknessF = 4.0
    gnres@gsLineDashPattern := 0
    plots@$unique_string("box")$ = \
      gsn_add_polyline(wks, box, (/bxmin, bxmax/), \
                       (/model_arr_trend(ind_ref(0)), \
                         model_arr_trend(ind_ref(0))/), gnres)
    gnres@gsLineDashPattern := 1
    plots@$unique_string("box")$ = \
      gsn_add_polyline(wks, box, (/bxmin, bxmax/), \
                       (/model_arr_trend(ind_ref(1)), \
                         model_arr_trend(ind_ref(1))/), gnres)
    gnres@gsLineDashPattern := 2
    plots@$unique_string("box")$ = \
      gsn_add_polyline(wks, box, (/bxmin, bxmax/), \
                       (/model_arr_trend(ind_ref(2)), \
                         model_arr_trend(ind_ref(2))/), gnres)

    newplot = gsn_attach_plots(plots(ivar), box, res, resb)

  end do

  ; Plotting panels

  pres = True

  pres@gsnPanelMainString = "Precipitation Anomaly"

  pres@gsnPanelRowSpec = True
  pres@gsnPanelCenter  = False

  pres@gsnPanelFigureStrings = (/"a)", "c)", "d)", "e)", "f)"/)
  pres@gsnPanelFigureStringsPerimOn = False
  pres@gsnPanelFigureStringsBackgroundFillColor = "transparent"
  pres@gsnPanelFigureStringsJust = "TopLeft"
  pres@gsnPanelFigureStringsFontHeightF = 0.016

  gsn_panel(wks, plots, (/1, 2, 2/), pres)

  log_info(" wrote " + outfile + "." + file_type)

  ; ***************************************
  ; add meta data to plot (for reporting)
  ; ***************************************

  out_path1 = out_path + "precip_" + anom + "_0" + ".nc"

  domain = "reg"
  statistics = (/"anomaly", "trend"/)
  caption = "Global and zonal average annual mean precipitation (mm day-1)"
  contrib_authors = (/"bock_lisa"/)

  ; Call provenance logger
  log_provenance(out_path1, \
                 outfile + "." + file_type, \
                 caption, \
                 statistics, \
                 domain, \
                 "times", \
                 contrib_authors, \
                 (/"eyring21ipcc"/), \
                 metadata_att_as_array(input_file_info, "filename"))

  leave_msg(DIAG_SCRIPT, "")

end
