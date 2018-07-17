; #############################################################################
;  Global timeseries for historical and RCPs 1850-2?00
;  Author: Ruth Lorenz (ETH, Switzerland)
;  PROJECT-NAME CRESCENDO
; #############################################################################;
; Caveats
;    Needs lots of memory, in particular if lots of models and scenarios
;    are included
;    Variable number of scenarios? So far hardcoded 4 scenarios (rcp2.6,
;    rcp4.5, rcp6.0, rcp8.5)
;    every model (ensemble) in RCP needs corresponding historical run as
;    reference, if no corresponding model found in scenario log_info
;
; Modification history
;    20182006-A_lore_ru: ported to version2, put plotting back into here
;    20170523-A_lore_ru: separated plotting from calcuation
;    20160625-A_lore_ru: adapted to ESMValTool
;    20130506-A_sedl_ja: written for IPCC AR5.
;
; #############################################################################

load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/plot/aux_plotting.ncl"
begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(diag_script, "")

  ; 'diag_script' and 'variables' are fetched from the
  ; above 'ncl.interface' file
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; 'datasets', 'variables' and 'field_types' are fetched from the
  ; above 'ncl.interface' file
  dim_MOD = dimsizes(dataset_info@dataset)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; ###########################################
  ; # References Information                  #
  ; ###########################################
;  write_references(diag_script, \  ; script name
;                   "A_lore_ru", \  ; authors
;                   "", \  ; contributors
;                   "D_collins13ipcc",    \  ; diag_references
;                   "", \  ; obs_references
;                   "P_crescendo")     ; proj_references

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", "eyears", \
                                           "begin_ref_year", "end_ref_year"/))
end

begin
  ; Get environment variables
  if (isatt(diag_script_info, "plot_type")) then
    file_type = diag_script_info@plot_type
  else
    file_type = "pdf"
  end if

  ; Basename of diag_script
  diag_script_base = basename(diag_script)  ; Needed?

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  ; Output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
end

begin
  ; #############################################
  ; # Get parameters from ./variable_defs/*.ncl #
  ; # passed via the 'info' attribute           #
  ; #############################################
  if(isvar("MyParam")) then
    delete(MyParam)
  end if
  if(isatt(variable_info, "long_name")) then
    MyParam = variable_info@long_name
    log_debug(" MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if
  ; ###########################################
  ; # Get data and average annual and globally#
  ; ###########################################
  ; get data from first dataset for historical period
  imod = 0  ; NCL array indices start from zero
  log_debug("processing " + dataset_info@dataset(imod) + "_" \
            + dataset_info@exp(imod) + "_" + dataset_info@ensemble(imod))

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field_type0)

  ; Check dimensions
  dims = getvardims(A0)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
  end if
  idx = ind(dims .eq. "lat")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lat dimension")
  end if
  idx = ind(dims .eq. "lon")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lon dimension")
  end if

  ; Average over time (function in ./diag_scripts/shared/statistics.ncl)
  data1 = time_operations(A0, -1, -1, "average", "yearly", True)

  dim_file = dimsizes(data1)  ; file size of historical file
  ; Average over globe/region
  ; (function in ./diag_scripts/shared/latlon.ncl)
  if (isatt(diag_script_info, "ts_minlat")) then
    data1_glob = area_operations(data1, diag_script_info@ts_minlat, \
                                 diag_script_info@ts_maxlat, \
                                 diag_script_info@ts_minlon, \
                                 diag_script_info@ts_maxlon, "average", \
                                 True)
  else
    data1_glob = coswgt_areaave(data1)
  end if

  ; How many historical model runs?
  ind_hist = ind(dataset_info@exp .eq. "historical")
  dim_MOD_hist = dimsizes(dataset_info@dataset(ind_hist))
  ref_data2 = new((/dim_MOD_hist/), float, data1_glob@_FillValue)
  data2_rel = new((/dim_MOD_hist, dim_file(0)/), float, \
                  data1_glob@_FillValue)

  ; How many scenarios?
  dim_scen = dimsizes(diag_script_info@scenarios)
  ; start and end years of timeperiods
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
  nperiods = dimsizes(syears)
  nyears = toint(eyears(nperiods - 1)) - toint(syears(0)) + 1
  nyearsp1 = toint(eyears(nperiods - 1)) - toint(eyears(0))
  ; How many rcp2.6 model runs?
  ind_26 = ind(dataset_info@exp .eq. diag_script_info@scenarios(0))
  dim_MOD_26 = dimsizes(dataset_info@dataset(ind_26))
  ; indicate array for scenarios
  data26_rel = new((/dim_MOD_26, nyearsp1/), float, \
                   data1_glob@_FillValue)
  i26 = 0

  if (dim_scen .gt. 1) then
    ; How many rcp4.5 model runs?
    ind_45 = ind(dataset_info@exp .eq. diag_script_info@scenarios(1))
    dim_MOD_45 = dimsizes(dataset_info@dataset(ind_45))
    data45_rel = new((/dim_MOD_45, nyearsp1/), float, \
                     data1_glob@_FillValue)
    i45 = 0
  end if
  if (dim_scen .gt. 2) then
    ; How many rcp6.0 model runs?
    ind_60 = ind(dataset_info@exp .eq. diag_script_info@scenarios(2))
    dim_MOD_60 = dimsizes(dataset_info@dataset(ind_60))
    data60_rel = new((/dim_MOD_60, nyearsp1/), float, \
                     data1_glob@_FillValue)
    i60 = 0
  end if
  if (dim_scen .gt. 3) then
    ; How many rcp8.5 model runs?
    ind_85 = ind(dataset_info@exp .eq. diag_script_info@scenarios(3))
    dim_MOD_85 = dimsizes(dataset_info@dataset(ind_85))
    data85_rel = new((/dim_MOD_85, nyearsp1/), float, \
                     data1_glob@_FillValue)
    i85 = 0
  end if
  i = 0
  do imod = 0, dim_MOD_hist - 1
    if (imod .eq. 0) then
      data2 = data1
      delete(A0)
    else
      A0 = read_data(imod, var0, field_type0)
      data2 = time_operations(A0, -1, -1, "average", "yearly", True)
      delete(A0)
    end if
    ; Average over globe (function in
    ; ./diag_scripts/shared/statistics.ncl)
    data2_glob = coswgt_areaave(data2)
    ind_start = ind(data2&year .eq. diag_script_info@begin_ref_year)
    ind_end = ind(data2&year .eq. diag_script_info@end_ref_year)
    ref_data2(imod) = dim_avg(data2_glob(ind_start : ind_end))
    dim_data2g = dimsizes(data2_glob)
    if (dim_data2g .ne. dim_file(0)) then
      tdiff = dim_file(0) - dim_data2g
      data2_rel(imod, 0 : tdiff - 1) = data2_glob@_FillValue
      data2_rel(imod, tdiff :) =  data2_glob - ref_data2(imod)
    else
      data2_rel(imod, :) = data2_glob - ref_data2(imod)
    end if
    ; clean up
    delete(data2)
    delete(data2_glob)

    ; find all other runs from this model and calculate relative ts
    imod26 = ind(dataset_info@dataset .eq. dataset_info@dataset(imod) .and. \
                 dataset_info@exp .eq. diag_script_info@scenarios(0) .and. \
                 dataset_info@ensemble .eq. dataset_info@ensemble(imod))
    if (.not.all(ismissing(imod26))) then
      A1 = read_data(imod26, var0, field_type0)
      data3 = time_operations(A1, -1, -1, "average", "yearly", True)
      data3_glob = coswgt_areaave(data3)
      dim_data3g = dimsizes(data3_glob)
      if (dim_data3g .ne. nyearsp1) then
        if (dim_data3g .gt. nyearsp1) then
          error_msg("w", diag_script, "", "Length of dataset " + \
                    dataset_info@dataset(imod) + " in recipe is " + \
                    "longer than period given in recipe (" + nyearsp1 + \
                    " years)")
        else
          data26_rel(i26, 0 : dim_data3g - 1) = data3_glob - ref_data2(imod)
          data26_rel(i26, dim_data3g :) =  data3_glob@_FillValue
        end if
      else
        data26_rel(i26, :) = data3_glob - ref_data2(imod)
      end if
      i26 = i26 + 1
      ; clean up
      delete(A1)
      delete(data3)
      delete(data3_glob)
    else
      log_info("No corresponding model for historial run " + \
               dataset_info@dataset(imod) + " found in scenario " + \
               diag_script_info@scenarios(0) + \
               ", cannot calculate relative change.")
    end if
    delete(imod26)

    if (dim_scen .gt. 1) then
      imod45 = ind(dataset_info@dataset .eq. dataset_info@dataset(imod) .and. \
                   dataset_info@exp .eq. diag_script_info@scenarios(1) .and. \
                   dataset_info@ensemble .eq. dataset_info@ensemble(imod))
      if (.not.all(ismissing(imod45))) then
        A1 = read_data(imod45, var0, field_type0)
        data3 = time_operations(A1, -1, -1, "average", "yearly", True)
        data3_glob = coswgt_areaave(data3)
        dim_data3g = dimsizes(data3_glob)
        if (dim_data3g .ne. nyearsp1) then
          if (dim_data3g .gt. nyearsp1) then
            error_msg("w", diag_script, "", "Length of dataset " + \
                      dataset_info@dataset(imod) + " in recipe is " + \
                      "longer than period given in recipe (" + nyearsp1 + \
                      " years)")
          else
            data45_rel(i45, 0:dim_data3g - 1) = data3_glob - ref_data2(imod)
            data45_rel(i45, dim_data3g:) =  data3_glob@_FillValue
          end if
        else
          data45_rel(i45, :) = data3_glob - ref_data2(imod)
        end if
        i45 = i45 + 1
        ; clean up
        delete(A1)
        delete(data3)
        delete(data3_glob)
      else
        log_info("No corresponding model ensemble for historial run " + \
                 dataset_info@dataset(imod) + " found in scenario " + \
                 diag_script_info@scenarios(1) + \
                 ", cannot calculate relative change.")
      end if
      delete(imod45)
    end if

    if (dim_scen .gt. 2) then
      imod60 = ind(dataset_info@dataset .eq. dataset_info@dataset(imod) .and. \
                   dataset_info@exp .eq. diag_script_info@scenarios(2) .and. \
                   dataset_info@ensemble .eq. dataset_info@ensemble(imod))
      if (.not.all(ismissing(imod60))) then
        A1 = read_data(imod60, var0, field_type0)
        data3 = time_operations(A1, -1, -1, "average", "yearly", True)
        data3_glob = coswgt_areaave(data3)
        dim_data3g = dimsizes(data3_glob)
        if (dim_data3g .ne. nyearsp1) then
          if (dim_data3g .gt. nyearsp1) then
            error_msg("w", diag_script, "", "Length of dataset " + \
                      dataset_info@dataset(imod) + " in recipe is " + \
                      "longer than period given in recipe (" + nyearsp1 + \
                      " years)")
          else
            data60_rel(i60, 0 : dim_data3g - 1) = data3_glob - ref_data2(imod)
            data60_rel(i60, dim_data3g :) =  data3_glob@_FillValue
          end if
        else
          data60_rel(i60, :) = data3_glob - ref_data2(imod)
        end if
        i60 = i60 + 1
        ; clean up
        delete(A1)
        delete(data3)
        delete(data3_glob)
      else
        log_info("No corresponding model ensemble for historial run " + \
                 dataset_info@dataset(imod) + " found in scenario " + \
                 diag_script_info@scenarios(2) + \
                 ", cannot calculate relative change.")
      end if
      delete(imod60)
    end if

    if (dim_scen .gt. 3) then
      imod85 = ind(dataset_info@dataset .eq. dataset_info@dataset(imod) .and. \
                   dataset_info@exp .eq. diag_script_info@scenarios(3) .and. \
                   dataset_info@ensemble .eq. dataset_info@ensemble(imod))
      if (.not.all(ismissing(imod85))) then
        A1 = read_data(imod85, var0, field_type0)
        data3 = time_operations(A1, -1, -1, "average", "yearly", True)
        data3_glob = coswgt_areaave(data3)
        dim_data3g = dimsizes(data3_glob)
        if (dim_data3g .ne. nyearsp1) then
          if (dim_data3g .gt. nyearsp1) then
            error_msg("w", diag_script, "", "Length of dataset " + \
                      dataset_info@dataset(imod) + " in recipe is " + \
                      "longer than period given in recipe (" + nyearsp1 + \
                      " years)")
          else
            data85_rel(i85, 0:dim_data3g - 1) = data3_glob - ref_data2(imod)
            data85_rel(i85, dim_data3g:) =  data3_glob@_FillValue
          end if
        else
          data85_rel(i85, :) = data3_glob - ref_data2(imod)
        end if
        i85 = i85 + 1
        ; clean up
        delete(A1)
        delete(data3)
        delete(data3_glob)
      else
        log_info("No corresponding model ensemble for historial run " + \
                 dataset_info@dataset(imod) + " found in scenario " + \
                 diag_script_info@scenarios(3) + \
                 ", cannot calculate relative change, set to missing.")
      end if
      delete(imod85)
    end if
  end do  ; imod

  ; average and standard deviation over all models
  datahist_avg = dim_avg_n_Wrap(data2_rel, 0)
  datahist_std = dim_stddev_n_Wrap(data2_rel, 0)

  ; put all data from 1850-2100 into mean_val and std_val
  mean_val = new((/dim_scen + 1, nyears/), typeof(data2_rel), \
                 data2_rel@_FillValue)
  std_val = new((/dim_scen + 1, nyears/), typeof(data2_rel), \
                data2_rel@_FillValue)
  p0 = toint(eyears(0)) - toint(syears(0))
  p1 = toint(eyears(0)) - toint(syears(0)) + 1
  mean_val(0, 0 : p0) = datahist_avg
  std_val(0, 0 : p0) = datahist_std

  data26_avg = dim_avg_n_Wrap(data26_rel, 0)
  data26_std = dim_stddev_n_Wrap(data26_rel, 0)
  mean_val(1, p1 :) = data26_avg
  std_val(1, p1 :) = data26_std

  if (dim_scen .gt. 1) then
    data45_avg = dim_avg_n_Wrap(data45_rel, 0)
    data45_std = dim_stddev_n_Wrap(data45_rel, 0)
    mean_val(2, p1 :) = data45_avg
    std_val(2, p1 :) = data45_std
  end if
  if (dim_scen .gt. 2) then
    data60_avg = dim_avg_n_Wrap(data60_rel, 0)
    data60_std = dim_stddev_n_Wrap(data60_rel, 0)
    mean_val(3, p1 :) = data60_avg
    std_val(3, p1 :) = data60_std
  end if
  if (dim_scen .gt. 3) then
    data85_avg = dim_avg_n_Wrap(data85_rel, 0)
    data85_std = dim_stddev_n_Wrap(data85_rel, 0)
    mean_val(4, p1 :) = data85_avg
    std_val(4, p1 :) = data85_std
  end if

  ; check number of models for timeperiods in scenarios
  ; if there are less than 5 models set to missing
  runs = new((/1 + dim_scen, nperiods/), integer)
  runs(0, 0) = dimsizes(ind(dataset_info@exp .eq. "historical" .and. \
                            dataset_info@end_year .eq. \
                            diag_script_info@eyears(0)))
  do scen = 0, dim_scen - 1
    do p = 0, nperiods - 1
      tmp = dimsizes(ind(dataset_info@exp .eq. \
                         diag_script_info@scenarios(scen) \
                         .and. toint(dataset_info@end_year) .gt. \
                         toint(diag_script_info@eyears(p))))
      runs(scen + 1, p) = dimsizes(ind(dataset_info@exp .eq. \
                                   diag_script_info@scenarios(scen) \
                                   .and. dataset_info@end_year .eq. \
                                   diag_script_info@eyears(p))) + tmp
      delete(tmp)
      if (runs(scen + 1, p) .le. 4) .and. (p .ne. 0) then
        p2 = toint(eyears(p - 1)) - toint(syears(0))
        p3 = toint(eyears(p)) - toint(syears(0))
        mean_val(scen + 1, p2 : p3) = mean_val@_FillValue
        std_val(scen + 1, p2 : p3) = std_val@_FillValue
        log_info("Scenario " + diag_script_info@scenarios(scen) + \
                 " in period ending " + diag_script_info@eyears(p) + \
                 " has less than 5 models, set to missing.")
      end if
    end do
  end do
  dims_data = dimsizes(mean_val)
  to_plot = new((/3 * dims_data(0), dims_data(1)/), float)

  do j = 0, dims_data(0) - 1
    ;  if variable sw or lw TOA flux: change direction
    if (var0 .eq. "rlut" .or. var0 .eq. "rsut") then
      to_plot(j * 3, :) = - mean_val(j, :)
    else
      to_plot(j * 3, :) = mean_val(j, :)
    end if
    to_plot((j * 3) + 1, :) = to_plot(j * 3, :) + 1.64 * std_val(j, :)
    to_plot((j * 3) + 2, :) = to_plot(j * 3, :) - 1.64 * std_val(j, :)
  end do
  copy_VarMeta(mean_val, to_plot)
  to_plot!0 = "metric"
  if (dim_scen .eq. 1) then
    to_plot&metric = (/"mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std"/) 
  else if (dim_scen .eq. 2) then
    to_plot&metric = (/"mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std"/) 
  else if (dim_scen .eq. 3) then
    to_plot&metric = (/"mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std"/)    
  else if (dim_scen .eq. 4) then
    to_plot&metric = (/"mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std", \
                       "mean", "+1.64std", "-1.64std"/)
  end if
  end if
  end if
  end if
             
  to_plot!1 = "year"
  to_plot&year = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 1)

  ; ###########################################
  ; # Other Metadata: diag_script, var        #
  ; ###########################################
  ; Add to to_plot, as attributes without prefix
  if(isatt(to_plot, "diag_script")) then  ; Add to existing entries
    temp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot@diag_script = (/diag_script/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if(isatt(variable_info, "long_name")) then
    to_plot@var_long_name = variable_info@long_name
  end if
  if(isatt(variable_info, "units")) then
    to_plot@var_units     = variable_info@units
  end if

  runs!0 = "scenario"
  runs&scenario = diag_script_info@label
  runs!1 = "period"
  runs&period = diag_script_info@syears
  if(isatt(runs, "diag_script")) then  ; Add to existing entries
    temp = runs@diag_script
    delete(runs@diag_script)
    runs@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    runs@diag_script = (/diag_script/)
  end if
  runs@var = "nr"
  runs@var_long_name = "number of model runs per scenario and period"
  runs@var_units     = 1

  ; ###########################################
  ; # Output to netCDF                        #
  ; ###########################################
  to_plot@ncdf = var0 + "_ts_line_" + syears(0) + "-" + eyears(nperiods - 1) + \
                 ".nc"
  ncdf_file = work_dir + "/" + to_plot@ncdf
  ncdf_file@existing = "overwrite"
  ; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)
  if(isatt(diag_script_info, "model_nr") .and. \
     diag_script_info@model_nr .eq. True) then
    ncdf_file@existing = "append"
    ncdf_outfile = ncdf_write(runs, ncdf_file)
  end if

  ; ###########################################
  ; *********** PLOTTING ************
  dims_data_plot = (dimsizes(to_plot&metric) - 1) / 3

  ; Create plot variables
  outfile = var0 + "_ts_line_" + syears(0) + "-" + eyears(nperiods - 1)
  wks = get_wks(0, plot_dir, outfile)
  if(isatt(diag_script_info, "colormap")) then
    colormap_file = "./diag_scripts/shared/plot/rgb/" + \
                    diag_script_info@colormap
    colormap = RGBtoCmap(colormap_file)
  else
    colormap = "StepSeq25"
  end if
  gsn_define_colormap(wks, colormap)

  colors_main = (/1, 3, 4, 5, 6/)
  colors_spread = (/2, 3, 4, 5, 6/)

  res                   = True

  res@gsnDraw = False
  res@gsnFrame = False
  res@txFont   = 25
  res@vpHeightF = 0.5                    ; change aspect ratio of plot
  res@vpWidthF = 0.8
  res@vpXF     = 0.15                   ; start plot at x ndc coord

  if (isatt(diag_script_info, "title")) then
    res@tiMainString = diag_script_info@title
  end if
  if(isatt(diag_script_info, "yaxis")) then
    res@tiYAxisString = diag_script_info@yaxis
  end if
  res@tiXAxisString = "Year"

  res@xyDashPatterns     = 0
  res@trXMinF = toint(syears(0))
  res@trXMaxF = toint(eyears(nperiods - 1))
  if(isatt(diag_script_info, "ymin")) then
    res@trYMinF = diag_script_info@ymin
  end if
  if(isatt(diag_script_info, "ymax")) then
    res@trYMaxF = diag_script_info@ymax
  end if

  res@tmXBMode = "Explicit"
  res@tmXBValues = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 50)
  res@tmXBLabels = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 50)
  res@tmXBMinorValues  = ispan(toint(syears(0)), \
                               toint(eyears(nperiods - 1)), 10)

  res@xyLineThicknesses = 3.5
  res@xyLineColor      = colors_main(0)

  ; create plot
  plot  = gsn_csm_xy(wks, to_plot&year, to_plot(0, :), res)

  res@gsnXYFillColors = colors_spread(0)
  res@xyLineColor     = -1
  spread_plot  = gsn_csm_xy(wks, to_plot&year, to_plot(1 : 2, :), res)
  overlay(plot, spread_plot)

  delete(res@gsnXYFillColors)

  do j = 1, dims_data_plot
    if (all(ismissing(to_plot(j * 3, :)))) then
      error_msg("w", diag_script, "", "all missing values for " + \
                "dimension " + j * 3 + " in to_plot, need at least 5" + \
                " models, go to next")
      continue
    end if
    res@xyLineColor      = colors_main(j)
    plot1  = gsn_csm_xy(wks, to_plot&year, to_plot(j * 3, :), res)
    overlay(plot, plot1)

    res@gsnXYFillColors = colors_spread(j)
    res@gsnXYFillOpacities = (/0.2, 0.2, 0.2/)
    res@xyLineColor     = -1
    spread_plot = gsn_csm_xy(wks, to_plot&year,\
                             to_plot((j * 3) + 1 : (j * 3) + 2, :), res)
    overlay(plot, spread_plot)
    delete(res@gsnXYFillColors)
    delete(res@gsnXYFillOpacities)
  end do

  polyres                  = True
  polyres@gsLineThicknessF = 1.5
  ; add polylines
  dum = new(nperiods, graphic)
  do n = 0, nperiods - 1
    dum(n) = gsn_add_polyline(wks, plot, (/toint(eyears(n)), \
                              toint(eyears(n))/), (/-20, 20/), polyres)
  end do
  delete(res@xyLineColor)

  ; Attach a legend
  lgres                    = True
  lgres@lgLineColors       = colors_main    ; (/1, 3, 4, 5, 6/)
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.13           ; width of legend (NDC)
  lgres@vpHeightF          = 0.15           ; height of legend (NDC)
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgres@lgLineThicknessF   = 4
  lgres@lgPerimOn = False
  lgres@lgItemOrder = (/4, 3, 2, 1, 0/)
  if (isatt(diag_script_info, "label")) then
    labels = "  " + diag_script_info@label
    legend = gsn_create_legend(wks, 5, labels, lgres)
    amres                  = True
    amres@amJust           = "TopLeft"    ; Use bottom right corner of box
                                          ; for determining its location.
    amres@amParallelPosF   = -0.45              ; Move legend to right
    amres@amOrthogonalPosF = -0.45              ; Move legend down.

    annoid = gsn_add_annotation(plot, legend, amres)  ; add legend to plot
  end if
  draw(wks)
  frame(wks)

  ;alltags = array_append_record(tags,\
  ;          (/"DM_glob", "PT_times", "ST_mean", "ST_stddev", "ST_anomaly"/), 0)
  ;string1 = "Time series of global annual mean " + variable_info@long_name \
  ;+ " anomalies (relative to " + datasets@start_year(0)+ "-" + \
  ;datasets@end_year(0) + ") from CMIP5 concentration drive " + \
  ;"experiments. Projections are shown for each RCP for the " + \
  ;"multi-dataset mean (solid lines) and the 5 to 95% range " + \
  ;"(+-1.64 standar deviation) across the distribution of " + \
  ;"individual models (shading). "
  ;if (nperiods .gt. 1) then
  ;  string2 = "Discontinuities at 2100 are due to different numbers of" + \
  ;  "models performing the extension runs beyond the 21st " + \
  ;  "century and have no physical meaning. Only one ensemble " + \
  ;  "member is used from each model and numbers in the " + \
  ;  "figure indicate the number of different models " + \
  ;  "contributing to the different time periods. No ranges " + \
  ;  "are given for the RCP6.0 projections beyond 2100 as only" + \
  ;  " 2 models are available."
  ;else
  ;  string2 = ""
  ;end if
  ;caption = string1 + string2
  id = diag_script + "_" + variables(0)

  contrib_authors = (/"A_lore_ru", "A_sedl_ja"/)
  ;ESMValMD(wks@fullname, alltags, \
  ;         caption, id, variables(0), datasets@datasets, \
  ;         climofiles, diag_script, contrib_authors)

  ;delete([/alltags, caption, id, climofiles/])

  log_info(" Wrote " + plot_dir + "/" + outfile + "." + file_type)
  
end
