; #############################################################################
; PROCEDURE FOR THE ZONAL PLOT OF THE PERFORMANCE METRICS
; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
; ESMVal project
; #############################################################################

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/zonalmean_profile.ncl"

procedure perfmetrics_ptype_script()
begin

  ; Start loop over datasets, make sure reference dataset is calculated first
  modidx = ispan(0, nDatasets - 1, 1)
  modidx := array_append_record(ref_ind, modidx(ind(modidx.ne.ref_ind)), 0)

  annots = project_style(info_items, diag_script_info, "annots")

  ; Loop over datasets
  do ii = 0, dimsizes(modidx) - 1

    imod = modidx(ii)

    log_debug("Processing " + datasetnames(imod))

    ; Check if difference is possible
    if (diag_script_info@plot_diff .and. imod.ne.ref_ind) then
      l_diff = True
    else
      l_diff = False
    end if

    ; Set path for saving processed data
    system("mkdir -p " + config_user_info@work_dir)
    fname_abs = \
      str_join((/"perfmetrics", "zonal", annots(imod), var0, \
                 diag_script_info@time_avg, diag_script_info@region/), "_")
    workpath_abs = config_user_info@work_dir + fname_abs + ".nc"
    plotpath_abs = config_user_info@plot_dir + fname_abs
    fname_diff = \
      str_join((/"perfmetrics", "zonal", \
                 annots(imod) + "-" + annots(ref_ind), var0, \
                 diag_script_info@time_avg, diag_script_info@region/), "_")
    workpath_diff = config_user_info@work_dir + fname_diff + ".nc"
    plotpath_diff = config_user_info@plot_dir + fname_diff

    ; Determine start/end year
    start_year = info_items[imod]@start_year
    end_year = info_items[imod]@end_year

    ; Read data
    var = read_data(info_items[imod])

    ; Calculate zonal mean
    var_zon = dim_avg_Wrap(var)
    delete(var)

    ; Calculate time average and standard deviation w.r.t.
    ; interannual variability
    var_avg = time_operations(var_zon, start_year, end_year, "average", \
                              diag_script_info@time_avg, True)
    var_std = interannual_variability(var_zon, start_year, end_year, \
                                      diag_script_info@time_avg, "None")
    delete(var_zon)

    ; Store in global array
    var_all = new(array_append_record(dimsizes(var_avg), 2, 0), float)
    var_all!2 = "statistic"
    var_all&statistic = (/"mean", "stddev"/)
    var_all(:, :, 0) = var_avg
    var_all(:, :, 1) = var_std
    var_all@var = var0
    var_all@diag_script = (/DIAG_SCRIPT/)
    copy_VarAtts(diag_script_info, var_all)
    var_all@ref_model = variable_info[0]@reference_dataset
    delete(var_avg)
    delete(var_std)

    ; Store reference dataset
    if (imod.eq.ref_ind) then
      var_ref = var_all
    end if

    ; Calculate difference
    if (l_diff) then
      var_diff = var_all(:, :, 0)
      var_diff = var_all(:, :, 0) - var_ref(:, :, 0)
    end if

    ; Write output
    if (config_user_info@write_netcdf) then
      var_all@ncdf = workpath_abs
      ncdf_abs = ncdf_write(var_all, workpath_abs)
      if (l_diff) then
        var_diff@ncdf = workpath_diff
        ncdf_diff = ncdf_write(var_diff, workpath_diff)
      end if
    end if

    ; Plotting
    if (.not.config_user_info@write_plots) then
      continue
    end if

    ; Convert units for plotting (if required)
    if (isatt(diag_script_info, "plot_units")) then
      var_all = convert_units(var_all, diag_script_info@plot_units)
      if (l_diff)
        var_diff = convert_units(var_diff, diag_script_info@plot_units)
      end if
    end if

    ; Absolute plot
    plot_var = var_all(:, :, 0)
    plot_var@res_tiMainString = annots(imod)
    plot_var@res_cnLevelSelectionMode = "ExplicitLevels"
    plot_var@res_cnLevels = diag_script_info@abs_levs
    plot_var@res_trYMinF = diag_script_info@zonal_ymin
    plot_var@res_gsnLeftString = \
      plot_var@long_name + " [" + format_units(plot_var@units) + "]"
    wks = gsn_open_wks(file_type, plotpath_abs)
    gsn_define_colormap(wks, diag_script_info@zonal_cmap)
    plot = zonalmean_profile(wks, plot_var, var0)
    draw(plot)
    frame(wks)
    delete(plot_var@res_cnLevels)
    delete(wks)
    delete(plot)

    ; Call provenance logger
    log_provenance(ncdf_abs, \
                   plotpath_abs + "." + file_type, \
                   "Zonal mean of variable " + var0, \
                   (/"mean", "stddev"/), \
                   diag_script_info@region, \
                   "zonal", \
                   (/"frank_franziska", "righi_mattia", "eyring_veronika"/), \
                   (/"righi15gmd", "gleckler08jgr"/), \
                   info_items[imod]@filename)

    ; Difference plot to the reference
    if (.not.l_diff) then
      continue
    end if
    plot_var = var_diff
    plot_var@res_tiMainString = annots(imod) + " - " + annots(ref_ind)
    if (isatt(diag_script_info, "t_test")) then
      plot_var@res_gsnLeftString = plot_var@res_gsnLeftString + " - " + \
        sprinti("%2i", toint(100 * diag_script_info@conf_level)) + "% c.l."
    end if
    plot_var@res_cnLevelSelectionMode = "ExplicitLevels"
    plot_var@res_cnLevels = diag_script_info@diff_levs
    plot_var@res_cnMissingValFillColor = "gray70"
    plot_var@res_gsnSpreadColorEnd = -2
    wks = gsn_open_wks(file_type, plotpath_diff)
    gsn_merge_colormaps(wks, "temp_19lev", "gray70")
    plot = zonalmean_profile(wks, plot_var, var0)
    delete(plot_var@res_cnLevels)

    ; Apply t-test (assuming different population variances)
    caption = "Difference in zonal mean of variable " + var0
    if (diag_script_info@t_test) then

      x1 = var_all(:, :, 0)
      x2 = var_ref(:, :, 0)
      s1 = var_all(:, :, 1) ^ 2
      s2 = var_ref(:, :, 1) ^ 2
      n1 = info_items[imod]@end_year - info_items[imod]@start_year + 1
      n2 = info_items[ref_ind]@end_year - info_items[ref_ind]@start_year + 1
      prob = ttest(x1, s1, n1, x2, s2, n2, True, False)
      plot_var = 1. - prob
      delete(prob)

      ; Mask non-significant values
      plot_var@res_cnMissingValFillColor = -1
      plot_var@res_gsnDraw = False
      plot_var@res_gsnFrame = False
      plot_var@res_cnLevelSelectionMode = "ExplicitLevels"
      plot_var@res_cnLevels := diag_script_info@conf_level
      if (diag_script_info@stippling) then
        plot_var@res_cnFillColors = (/"transparent", "black"/)
        plot_var@res_cnFillPattern = 17
        caption = caption + " (non-significant points are stippled)"
      else
        plot_var@res_cnFillColors = (/"transparent", "gray70"/)
        caption = caption + " (non-significant points are masked out in gray)"
      end if
      plot_var@res_cnInfoLabelOn = False
      plot_var@res_cnLinesOn = False
      plot_var@res_cnLineLabelsOn = False
      plot_var@res_lbLabelBarOn = False
      plot_var@res_gsnLeftString = ""
      delete(plot_var@long_name)
      plot2 = zonalmean_profile(wks, plot_var, var0)
      overlay(plot, plot2)
      draw(plot)
      frame(wks)
      delete(wks)
      delete(plot)
      delete(plot2)
      delete(plot_var)

    end if

    ; Call provenance logger
    log_provenance(ncdf_diff, \
                   plotpath_diff + "." + file_type, \
                   caption, \
                   (/"mean", "stddev", "diff"/), \
                   diag_script_info@region, \
                   "zonal", \
                   (/"frank_franziska", "righi_mattia", "eyring_veronika"/), \
                   (/"righi15gmd", "gleckler08jgr"/), \
                   info_items[imod]@filename)

  end do

end
