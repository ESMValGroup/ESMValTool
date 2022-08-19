; #############################################################################
; WRAPPER SCRIPT FOR COLLECTING AND PLOTTING TIMESERIES
; Author: Lisa Bock (DLR, Germany)
; IPCC AR6
; #############################################################################
; Description
;    Collects timeseries and their statistics previously calculated by
;    teas_anom.ncl and passes them to a new tsline plot
;
; Required diag_script_info attributes
;    None
;
; Optional diag_script_info attributes
;    blending: if true, then var="gmst" otherwise "gsat"
;    ref_start: start year of reference period for anomalies
;    ref_end: end year of reference period for anomalies
;    region: name of domain
;    plot_units: variable unit for plotting
;    y-min: set min of y-axis
;    y-max: set max of y-axis
;    order: order in which experiments should be plotted
;    stat_shading: if true: shading of statistic range
;    ref_shading: if true: shading of reference period
;
; Caveats
;
; Modification history
;    20190204-A_bock_lisa: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  diag_script_base = basename(DIAG_SCRIPT)

  var_names = metadata_att_as_array(variable_info, "short_name")

  var0 = variable_info[0]@short_name
  info_var = select_metadata_by_name(input_file_info, var0)
  units0 = info_var[0]@units
  datasets = metadata_att_as_array(info_var, "dataset")
  dim_MOD = ListCount(info_var)

  if (diag_script_info@blending) then
    var = "gmst"
    units0 = "degC"
  else
    var = "gsat"
  end if

  ; Create output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)

  ; Check file type
  file_type = config_user_info@output_file_type
  if(ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "region", "Global")
  set_default_att(diag_script_info, "stat_shading", False)
  set_default_att(diag_script_info, "ref_shading", False)

  ; Determine time range
  start_year = diag_script_info@start_year
  end_year = diag_script_info@end_year
  all_years = ispan(start_year, end_year, 1)

  ; Create time coordinate
  ntime = end_year - start_year + 1
  time = new(ntime, integer)
  time = ispan(start_year, end_year, 1)

  ref_start = diag_script_info@ref_start
  ref_end = diag_script_info@ref_end

  ; -------------------------------------------
  ; Read additional reference datasets
  ; -------------------------------------------

  if (dim_MOD .gt. 0) then

    data_ref = new((/dim_MOD, ntime/), double)
    data_ref!0 = "dataset"
    data_ref!1 = "time"
    data_ref&dataset = datasets
    data_ref&time = time

    ; Loop over datasets
    do imod = 0, dim_MOD - 1

      log_info("Process dataset: " + datasets(imod))

      ; Read data
      A0 = read_data(info_var[imod])

      ndim = dimsizes(A0)

      ; calculate anomaly for each grid point
      start_year_0 = info_var[imod]@start_year
      end_year_0 = info_var[imod]@end_year
      if (start_year_0 .gt. diag_script_info@ref_start) then
        ref_start = start_year_0
      else
        ref_start = diag_script_info@ref_start
      end if
      if (ref_start .ne. 1961 .or. ref_end .ne. 1990) then
        tmp = time_operations(A0, ref_start, ref_end, "average", \
                              "monthlyclim", True)
        do i = 0, ndim(0) - 1
          A0(i, :, :) = A0(i, :, :) - tmp(mod(i, 12), :, :)
        end do
      end if

      ; calculate field mean
      A0_mean := area_operations(A0, -90., 90., \
                                 0., 360., "average", True)

      ; yearly average
      diag = time_operations(A0_mean, -1, -1, "average", \
                             "yearly", True)

      ; Match time coordinate
      date = diag&year
      idx1 = ind(date(0).eq.data_ref&time)
      idx2 = ind(date(dimsizes(date) - 1).eq.data_ref&time)
      delete(date)

      data_ref(imod, idx1:idx2)  = (/diag/)

      delete(A0)
      delete(A0_mean)
      delete(diag)

    end do

  end if

  ; -----------------------------------------------
  ; Read pre-calculated multi model statistics
  ; -----------------------------------------------

  input_files = diag_script_info@input_files + "/" + \
    "tsline_gmst_anom" + "_stat.nc"
  input_files := tostring(input_files)

  ; Filter out non-existing files
  input_files := input_files(ind(isfilepresent(input_files)))

  nfiles = dimsizes(input_files)

  dim_data   = (/nfiles, 3, ntime/)
  data_mean  = new(dim_data, double)
  data_mean!0 = "experiment"
  data_mean!1 = "stat"
  data_mean!2 = "time"
  data_mean&experiment = new(nfiles, string, "exp")
  data_mean&stat = (/"mean", "stat1", "stat2"/)
  data_mean&time = time
  data_mean@var = var

  ; Loop over files in the list, read and append data
  do ii = 0, nfiles - 1

    log_info("Read in " + input_files(ii))

    var_mean = var + "_mean"
    var_stddev = var + "_stddev"

    data_temp = ncdf_read(input_files(ii), var_mean)
    data_stat = ncdf_read(input_files(ii), var_stddev)
    data_stat1 = ncdf_read(input_files(ii), var+"_5quantile")
    data_stat2 = ncdf_read(input_files(ii), var+"_95quantile")

    log_info("Experiment: " + data_temp@experiment)

    if(isatt(diag_script_info, "order")) then
      ii_n = ind(data_temp@experiment .eq. diag_script_info@order)
      if (all(ismissing(ii_n))) then
        error_msg("f", DIAG_SCRIPT, "", "Collected experiment " + \
                  data_temp@experiment + " does not appear in " + \
                  "diag_script_info@order")
      elseif (ii_n .gt. nfiles-1)
        error_msg("f", DIAG_SCRIPT, "", "diag_script_info@order not " \
                  + "well defined (less experiments available)")
      end if
    else
      ii_n = ii
    end if

    ; Match time coordinate
    date = data_temp&time
    idx1 = ind(date(0).eq.data_mean&time)
    idx2 = ind(date(dimsizes(date) - 1).eq.data_mean&time)
    delete(date)

    data_mean(ii_n, 0, idx1:idx2)  = (/data_temp(:)/)
    data_mean(ii_n, 1, idx1:idx2)  = (/data_stat1(:)/)
    data_mean(ii_n, 2, idx1:idx2)  = (/data_stat2(:)/)
    data_mean&experiment(ii_n) = data_temp@experiment
    data_mean@units = units0

    delete(data_temp)
    delete(data_stat)
    delete(data_stat1)
    delete(data_stat2)

  end do

  ; -------------------------------------------
  ; NetCDF Output
  ; -------------------------------------------

  do i = 0, dimsizes(data_mean&experiment) - 1
    if (any(data_mean&experiment(i) .eq. (/"CMIP6_historical", \
                                           "CMIP6_historical-ssp245"/))) then
      data_mean&experiment(i) = "CMIP6"
    elseif (any(data_mean&experiment(i) .eq. (/"CMIP5_historical", \
                                               "CMIP5_historical-rcp45"/)))
      data_mean&experiment(i) = "CMIP5"
    end if
  end do

  out_path = config_user_info@work_dir
  system("mkdir -p " + out_path)
  out_path1 = out_path + "tsline_collect_" + var0 + ".nc"
  data_mean@ncdf = out_path
  data_mean@diag_script = DIAG_SCRIPT
  data_mean@var = var
  ncdf_outfile = ncdf_write(data_mean, out_path1)

  if (isatt(diag_script_info, "ref")) then
    out_path1 = out_path + "tsline_collect_" + var0 + "_ref.nc"
    data_ref@ncdf = out_path
    data_ref@diag_script = DIAG_SCRIPT
    data_ref@var = var
    ncdf_outfile_ref = ncdf_write(data_ref, out_path1)
  end if

  ; -------------------------------------------
  ; Plotting
  ; -------------------------------------------

  ; Define workstation
  outfile = config_user_info@plot_dir + data_mean@var + "_" + \
    str_sub_str(diag_script_info@region, " ", "_") + \
    "_multimodel_anom_" + start_year + "-" + end_year
  wks = gsn_open_wks(file_type, outfile)

  ; Select colors and other plotting attributes
  ; (see ./diag_scripts/lib/ncl/style.ncl)
  colors = project_style(info_var, diag_script_info, "colors")
  dashes = project_style(info_var, diag_script_info, "dashes")
  thicks = project_style(info_var, diag_script_info, "thicks")
  annots = project_style(info_var, diag_script_info, "annots")

  ; ************************************************
  ; plotting parameters
  ; ************************************************

  plot = new(1, graphic)

  res           = True
  res@gsnDraw   = False
  res@gsnFrame  = False

  res@vpXF      = 0.05
  res@vpYF      = 0.7
  res@vpHeightF = 0.4
  res@vpWidthF  = 0.7
  res@pmLegendDisplayMode = "Never"
  res@tmYRLabelFontHeightF = 0.016
  res@tmYLLabelFontHeightF = 0.016
  res@tiXAxisFontHeightF = 0.016
  res@tiYAxisFontHeightF = 0.016
  res@tmXBMajorOutwardLengthF = 0.006
  res@tmYLMajorOutwardLengthF = 0.006
  res@tmXBMinorOutwardLengthF = 0.003
  res@tmYLMinorOutwardLengthF = 0.003
  res@tmXBMajorLengthF = 0.006
  res@tmYLMajorLengthF = 0.006
  res@tmXBMinorLengthF = 0.003
  res@tmYLMinorLengthF = 0.003

  if (isatt(diag_script_info, "ref")) then
    min_tmp = (/min(data_mean), min(data_ref)/)
    max_tmp = (/max(data_mean), max(data_ref)/)
  else
    min_tmp = (/min(data_mean), min(data_stat)/)
    max_tmp = (/max(data_mean), max(data_stat)/)
  end if
  res@trYMinF  = min(min_tmp) - 0.05 * (max(max_tmp) - min(min_tmp))
  res@trYMaxF  = max(max_tmp) + 0.05 * (max(max_tmp) - min(min_tmp))
  res@tiYAxisOn = True
  res@tiXAxisString = ""
  res@gsnStringFontHeightF = 0.016

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
  res@tmXBLabelAngleF = 45
  res@tmXBLabelJust     = "CenterRight"
  if (isatt(diag_script_info, "y_min")) then
    res@trYMinF = diag_script_info@y_min
  end if
  if (isatt(diag_script_info, "y_max")) then
    res@trYMaxF = diag_script_info@y_max
  end if

  res@tmXBMode         = "Manual"
  res@tmXBTickSpacingF = 20

  res@tiMainString = ""

  res@tiYAxisString = "Change in Temperature (~F34~0~F~ C)"

  ; ---------------------------------------------------------------------
  ; add multi model mean of different experiments
  ; ---------------------------------------------------------------------
  ; number of different experiments
  nexp     = dim_data(0)

  cmap = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                            + "ipcc-ar6_cmip_line.rgb")
  res@xyDashPatterns    = (/0, 0, 0, 0, 0, 0/)
  res@xyLineThicknesses = (/5, 5, 5, 5, 5, 5/)
  res@xyLineColors = cmap(::-1, :)

  if (nexp .gt. 6) then
    error_msg("w", scriptname, funcname, "Color palette not defined for " \
              + "more than " + nexp + " experiments")
  end if

  res@pmLegendDisplayMode = "Never"
  res@xyExplicitLegendLabels = data_mean&experiment
  res@lgBoxMinorExtentF      = 0.2       ; Shorten the legend lines

  plot(0) = gsn_csm_xy(wks, data_ref&time, data_mean(:, 0, :), res)

  res@pmLegendDisplayMode = "Never"

  ; ---------------------------------------------------------------------
  ; Add multi model statistics (5% and 95% quantile)
  ; ---------------------------------------------------------------------
  if (isatt(diag_script_info, "stat_shading")) then
    if (diag_script_info@stat_shading .ne. False) then
      res_stat = True
      res_stat@gsnXYFillOpacities = 0.2
      cmap = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                                + "ipcc-ar6_cmip_shading.rgb")
      copy_VarMeta(res, res_stat)
      res_stat@gsnXYFillColors = cmap(2, :)
      delete(res_stat@xyLineColors)
      ; We don't want the line, so make it transparent.
      res_stat@xyLineColor = cmap(2, :)
      delete(res_stat@xyLineThicknesses)
      res_stat@xyLineThicknesses = (/1, 1/)
      shading_plot = gsn_csm_xy(wks, data_ref&time, data_mean(0, 1:2, :), \
                                res_stat)
      overlay(plot(0), shading_plot)
      if(nexp .ge. 2) then
        res_stat@gsnXYFillColors = cmap(1, :)
        res_stat@xyLineColor = cmap(1, :)
        shading_plot = gsn_csm_xy(wks, data_ref&time, data_mean(1, 1:2, :), \
                                  res_stat)
        overlay(plot(0), shading_plot)
      end if
      if(nexp .ge. 3) then
        res_stat@gsnXYFillColors = cmap(0, :)
        res_stat@xyLineColor = cmap(0, :)
        shading_plot = gsn_csm_xy(wks, data_ref&time, data_mean(2, 1:2, :), \
                                  res_stat)
        overlay(plot(0), shading_plot)
      end if
      if(nexp .ge. 4) then
        res_stat@gsnXYFillColors = cmap(3, :)
        res_stat@xyLineColor = cmap(3, :)
        shading_plot = gsn_csm_xy(wks, data_ref&time, data_mean(3, 1:2, :), \
                                  res_stat)
        overlay(plot(0), shading_plot)
      end if
      if(nexp .gt. 4) then
        error_msg("w", scriptname, funcname, "Color palette not defined for " \
                  + nexp + " experiments")
      end if
    end if
  end if

  ; ---------------------------------------------------------------------
  ; add reference datasets
  ; ---------------------------------------------------------------------
  if (isatt(diag_script_info, "ref")) then
    delete(res@xyDashPatterns)
    delete(res@xyLineThicknesses)
    delete(res@xyLineColors)
    res@xyDashPatterns    = dashes
    res@xyLineThicknesses = (/3, 3, 3, 3, 3/)   ; make lines thicker
    res@xyLineColors      = colors  ; change line color

    ref_p = gsn_csm_xy(wks, data_ref&time, data_ref, res)  ; create plot
    overlay(plot(0), ref_p)
  end if

  if (diag_script_info@ref_shading) then

    ; ---------------------------------------------------------------------
    ; grey shading of reference period
    ; ---------------------------------------------------------------------
    gsres                = True
    gsres@gsFillColor    = "Grey70"
    gsres@gsFillOpacityF = 0.1
    xbox = (/ref_start, ref_end, ref_end, ref_start, ref_start/)
    ybox = (/res@trYMinF, res@trYMinF, res@trYMaxF, res@trYMaxF, res@trYMinF/)
    newplot00 = gsn_add_polygon(wks, plot(0), xbox, ybox, gsres)

    txtres               = True
    txtres@txFont        = "helvetica-bold"
    txtres@txFontColor   = "Grey40"
    txtres@txFontHeightF = 0.013
    ref_txt = gsn_add_text(wks, plot, "reference period", \
                           0.5*(ref_start + ref_end), \
                           res@trYMaxF - 0.05 * (res@trYMaxF - res@trYMinF), \
                           txtres)

  end if

  ; ---------------------------------------------------------------------
  ; Draw some lines to create a legend
  ; ---------------------------------------------------------------------
  res_lines                   = True      ; polyline mods desired
  res_lines@tfPolyDrawOrder   = "Predraw"
  res_lines@gsLineColor       = "grey"    ; line color
  res_lines@gsLineThicknessF  = 1.        ; line thicker
  res_lines@gsLineDashPattern = 1.        ; dash pattern

  xx = (/res@trXMinF, res@trXMaxF/)
  yy = (/0.0, 0.0/)
  dum0 = gsn_add_polyline(wks, plot(0), xx, yy, res_lines)

  ; ***********************************************
  ; legend resources
  ; ***********************************************

  lgres                    = True
  lgres@lgItemType         = "Lines"   ; show lines only (default)
  lgres@lgLabelFontHeightF = .06       ; set the legend label font thickness
  lgres@vpWidthF           = 0.2      ; width of legend (NDC)
  lgres@vpHeightF          = 0.15      ; height of legend (NDC)
  lgres@lgPerimColor       = "gray"    ; draw the box perimeter in orange
  lgres@lgPerimThicknessF  = 1.0       ; thicken the box perimeter
  lgres@lgPerimFill        = 0
  lgres@lgPerimFillColor   = "white"

  delete(cmap)
  cmap = read_colormap_file("$diag_scripts/shared/plot/rgb/" \
                            + "ipcc-ar6_cmip_line.rgb")
  tmp := tostring(cmap)
  tmp1 = (/"(/" + tmp(2, 0) + "," + tmp(2, 1) + "," + tmp(2, 2) + "/)", \
           "(/" + tmp(1, 0) + "," + tmp(1, 1) + "," + tmp(1, 2) + "/)", \
           "(/" + tmp(0, 0) + "," + tmp(0, 1) + "," + tmp(0, 2) + "/)"/)
  lgres@lgLineColors      = array_append_record(colors, tmp1, 0)
  lgres@lgDashIndexes     = array_append_record(dashes, (/0, 0, 0/), 0)
  lgres@lgLineThicknessF = 3
  labels = array_append_record(annots, data_mean&experiment, 0)
  nitems = dimsizes(labels)
  lgres@lgItemOrder = ispan(nitems - 1, 0, 1)

  ; Create legend
  lbid = gsn_create_legend(wks, nitems, labels, lgres)

  amres = True
  amres@amParallelPosF   = 0.0
  amres@amOrthogonalPosF = - 0.25
  annoid1 = gsn_add_annotation(plot(0), lbid, amres)

  resP = True
  resP@gsnMaximize = True
  resP@gsnPaperOrientation = "portrait"
  resP@gsnPaperMargin =  0.8

  gsn_panel(wks, plot, (/1, 1/), resP)

  log_info(" wrote " + outfile + "." + file_type)

  ; ***************************************
  ; add meta data to plot (for reporting)
  ; ***************************************

  if (diag_script_info@region .eq. "Global") then
    domain = "global"
  else
    domain = "reg"
  end if

  caption = "Multi-model mean time series of the anomaly of near-surface air" \
            + "temperature"

  infiles = metadata_att_as_array(input_file_info, "filename")
  infiles := array_append_record(infiles, input_files, 0)

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 outfile + "." + file_type, \
                 caption, \
                 "anomaly", \
                 domain, \
                 "times", \
                 (/"bock_lisa"/), \
                 (/"eyring21ipcc"/), \
                 infiles)

  leave_msg(DIAG_SCRIPT, "")

end
