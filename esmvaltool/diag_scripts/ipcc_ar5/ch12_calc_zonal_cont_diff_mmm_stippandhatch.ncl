; #############################################################################
; Annual mean zonal mean multi model changes
; Author: Ruth Lorenz (ETH Zurich, Switzerland)
; CRESCENDO project
; #############################################################################
;
; Description
;    Multi-model annual or seasonal mean zonal mean changes
;    e.g. 2081-2100 relative to base period (1986-2005 in CMIP5)
;    for multiple RCPs
;
; Required diag_script_info attributes (diagnostics specific)
;    ncdf: filename for netcdf to write
;    scenarios: which scenarios to include
;    periods: list with start years of periods to be included
;    	       e.g. (/"2046", "2081", "2180"/)
;    time_avg: list with seasons or annual to be included e.g (/"annualclim"/)
;
; Optional diag_script_info attributes (diagnostic specific)
;    title: title of whole figure
;    label: labels of individual panels
;    max_hori: maximum number of panels in horizontal
;    max_vert: maximum number of panels in vertical
;    colormap: colortable different from default
;    explicit_cn_levels: range and levels for contour plot
;    base_cn: flag indicating if base period plotted in
;             contour lines (as in ipcc Fig. 12.19 for wind)
;    base_cnLevelSpacing: spacing for contour levels (e.g. 10)
;                         if base_cn is True
;    base_cnMinLevel: min value for base contours
;    base_cnMaxLevel: max level for base contours
;
; Modification history
;    20181024-A_lore_ru: ported to v2, regridding, zonal means in preprocessor
;    20170517-A_lore_ru: changed to regrid to common grid first before
;                        calculating zonal average and added ability to use
;                        irregular grids (e.g. thetao)
;    20170127-A_lore_ru: changed routine for plotting stippling/hatching to
;                        zonalmean_profile to get dots and hatch
;    20161031-A_lore_ru: increased flexibility and add stippling and hatching
;    20160831-A_lore_ru: adapted for ESMValTool
;    20130503-A_sedl_ja: written for IPCC AR5
;
; #############################################################################
; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/ensemble.ncl"

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

  ; 'models', 'variables' and 'field_types' are fetched from the
  ; above 'ncl.interface' file
  dim_MOD = dimsizes(dataset_info@dataset)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; Check required diag_script_info attributes
  req_atts = (/"scenarios", "periods", "time_avg", "label"/)
  exit_if_missing_atts(diag_script_info, req_atts)
end

begin
  ; Get environment variables
  if (isatt(diag_script_info, "plot_type")) then
    file_type = diag_script_info@plot_type
  else
    file_type = "pdf"
  end if

  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  input_dir = diag_script_info@input_files(1)
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
    log_info(" MyParam = " + MyParam)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; ###########################################
  ; # Get data and average time               #
  ; ###########################################
  ; get data from first model
  imod = 0  ; NCL array indicies start from zero
  log_info("processing " + dataset_info@dataset(imod) + "_"\
           + dataset_info@exp(imod) + "_"\
           + dataset_info@ensemble(imod))

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field_type0)

  ; Check dimensions
  dims = getvardims(A0)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", diag_script, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
  end if
  idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no lat dimension")
  end if
  lat_name = dims(idx)

  idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
  if (ismissing(idx)) then
    error_msg("f", diag_script, "", "no level dimension")
  end if
  lev_name = dims(idx)

  if (isatt(diag_script_info, "seasons")) then
    dim_seas = dimsizes(diag_script_info@seasons)
  else
    dim_seas = 1
  end if

  ; check if diag_script_info@label has correct size, otherwise error message
  dim_periods = dimsizes(diag_script_info@periods)
  dim_scenarios = dimsizes(diag_script_info@scenarios)
  dim_label = dim_seas * dim_scenarios * dim_periods
  if (dim_label .ne. dimsizes(diag_script_info@label)) then
    error_msg("f", diag_script, "", "wrong number of labels, needs to be " + \
              "seasons * scenarios * periods")
  end if

  annots = project_style(diag_script_info, "annots")
  Fill = default_fillvalue(typeof(A0))

  var_reg = new((/dim_MOD, dim_seas, dimsizes(A0&$lev_name$), \
                 dimsizes(A0&lat)/), typeof(A0), Fill)

  var_reg!0 = "models"
  var_reg!1 = "season"
  var_reg!3 = "lat"
  var_reg&models = annots
  if (isatt(diag_script_info, "seasons")) then
    var_reg&season = (/ispan(0, dim_seas - 1, 1)/)
  else
    var_reg&season = (/0/)
  end if
  var_reg!2 = lev_name
  var_reg&$lev_name$ = A0&$lev_name$
  var_reg&lat = A0&lat
  if (isatt(diag_script_info, "iavmode")) then
    if (diag_script_info@iavmode .eq. "each") then
      iav = new((/dim_MOD, dim_seas, dimsizes(A0&lev), \
                 dimsizes(A0&lat)/), typeof(A0), Fill)
      copy_VarMeta(var_reg, iav)
    end if
  end if
  ; Average over time (function in ./diag_scripts/lib/ncl/statistics.ncl)
  data1 = time_operations(A0, -1, -1, "average", diag_script_info@time_avg, \
                          True)

  do imod = 0, dim_MOD - 1
    if (imod .ne. 0) then
      A0 = read_data(imod, var0, field_type0)
      data1 = time_operations(A0, -1, -1, "average", "annualclim", True)
      dims = getvardims(A0)
      idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
      if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no lat dimension")
      end if
      lat_name = dims(idx)
      idx = ind(dims .eq. "lev" .or. dims .eq. "plev")
      if (ismissing(idx)) then
        error_msg("f", diag_script, "", "no level dimension")
      end if
      lev_name = dims(idx)
    end if

    if (dim_seas .eq. 1) then
      var_reg(imod, 0, :, :) = data1
    else
      var_reg(imod, :, :, :) = data1
    end if
    if (isatt(diag_script_info, "iavmode")) then
      if (diag_script_info@iavmode .eq. "each") then
        ; #########################################################
        ; # Read natural variability for stippling for each model #
        ; #########################################################
        iav_file = addfile(input_dir + "/" + \
                           "IAV_piControl_" + \
                           dataset_info@dataset(imod) + "_" + \
                           dataset_info@ensemble(imod) + "_" + field_type0 \
                           + "_" + var0 + "_" + \
                           diag_script_info@time_avg + ".nc", "r")
        iav_zon = iav_file->iav
        if (dim_seas .eq. 1) then
          iav(imod, 0, :, :) = iav_zon
        else
          iav(imod, :, :, :) = iav_zon
        end if
        delete(iav_native)
      end if
    end if
    delete(A0)
  end do
  if (isatt(diag_script_info, "iavmode")) then
    if (diag_script_info@iavmode .eq. "each") then
      ; calculate multi-model mean of iav
      tmp_pow = dim_avg_n(iav ^ 2, 0)
      natvar_mmm = rm_single_dims(tmp_pow * sqrt(2.))
    else
      log_info("Error: This iavmode is not know")
    end if
  else
    ; read already calculated mmm iav
    natvar_file = addfile(input_dir + "/" + \
                          "IAV_mmm_piControl_" + field_type0 + \
                          "_" + var0 + "_" + diag_script_info@time_avg + \
                          ".nc", "r")
    natvar_mmm = natvar_file->iav
  end if

  ; ##############################################
  ; # Calculate change from hist to periods      #
  ; ##############################################
  ; loop over rcps and periods
  dim_scen = dimsizes(diag_script_info@scenarios)
  dim_per = dimsizes(diag_script_info@periods)

  to_plot = new((/dim_scen * dim_per * dim_seas, \
                 dimsizes(var_reg&$lev_name$), dimsizes(var_reg&lat)/), \
                typeof(var_reg), Fill)
  if (isatt(diag_script_info, "base_cn")) then
    to_plot_base = new((/dim_scen * dim_per * dim_seas, \
                       dimsizes(var_reg&$lev_name$), \
                       dimsizes(var_reg&lat)/), typeof(var_reg), Fill)
  end if
  to_plot_signif = new((/dim_scen * dim_per * dim_seas, \
                        dimsizes(var_reg&$lev_name$),  \
                        dimsizes(var_reg&lat)/), \
                       typeof(var_reg), Fill)
  to_plot_not_signif = new((/dim_scen * dim_per * dim_seas, \
                            dimsizes(var_reg&$lev_name$),  \
                            dimsizes(var_reg&lat)/), \
                           typeof(var_reg), Fill)
  model_number = new((/dim_scen * dim_per * dim_seas/), integer)
  i = 0
  do seas = 0, dim_seas - 1  ; loop over seasons
    do rcp = 0, dim_scen - 1  ; loop over rcps
      do per = 0, dim_per - 1  ; loop over periods
        idx_rcp = ind((dataset_info@exp .eq. \
                       diag_script_info@scenarios(rcp)) .and. \
                      (dataset_info@start_year .eq. \
                       diag_script_info@periods(per)))
        proj = var_reg(idx_rcp, seas, :, :)
        proj_avg = dim_avg_n_Wrap(proj, 0)
        ; find historical runs from same models as in rcp
        idx_hist = find_model_index(dataset_info@dataset(idx_rcp), \
                                    dataset_info@dataset, "historical",\
                                    dataset_info@exp, \
                                    dataset_info@ensemble(idx_rcp), \
                                    dataset_info@ensemble, 0)
        base = var_reg(idx_hist, seas, :, :)
        base_avg = dim_avg_n_Wrap(base, 0)
        var_diff = proj - base
        delete([/base, proj/])
        to_plot(i, :, :) = proj_avg - base_avg
        if (isatt(diag_script_info, "base_cn")) then
          to_plot_base(i, :, :) = base_avg
        end if
        ; determine sigma and sign for significance
        log_debug("determine sigma and sign for significance")
        var_pos_signif = dim_num_n(where(var_diff .gt. 0., 1., \
                                         to_plot@_FillValue) .eq. 1., 0)
        var_neg_signif = dim_num_n(where(var_diff .lt. 0., 1., \
                                         to_plot@_FillValue) .eq. 1., 0)
        var_pos_signif_tmp = \
          var_pos_signif / where((var_pos_signif + var_neg_signif) \
                                 .ne. 0., var_pos_signif + var_neg_signif, \
                                 to_plot@_FillValue)
        var_neg_signif_tmp = \
          var_neg_signif / where((var_pos_signif + var_neg_signif) \
                                 .ne. 0., var_pos_signif + var_neg_signif, \
                                 to_plot@_FillValue)
        signif_and_sign = where(var_pos_signif_tmp .ge. 0.9 .or. \
                                var_neg_signif_tmp .ge. 0.9, 1., 0.)
        signif_and_sign = where(ismissing(signif_and_sign), 0., \
                                signif_and_sign)
        if (dim_seas .ne. 1) then
          sigma_and_sign = where(abs(to_plot(i, :, :)) .gt. \
                                 abs(2. * natvar_mmm(seas, :, :)), 1., 0.)
          to_plot_not_signif(i, :, :) = where(abs(to_plot(i, :, :)) .lt. \
                                              abs(natvar_mmm(seas, :, :)), \
                                              1., 0.)
        else
          sigma_and_sign = where(abs(to_plot(i, :, :)) .gt. \
                                 abs(2. * natvar_mmm(seas, :, :)), 1., 0.)
          to_plot_not_signif(i, :, :) = where(abs(to_plot(i, :, :)) .lt. \
                                              abs(natvar_mmm(seas, :, :)), \
                                              1., 0.)
        end if
        to_plot_signif(i, :, :) = where(signif_and_sign .eq. 1, \
                                        sigma_and_sign, 0.)
        model_number(i) = dimsizes(idx_rcp)
        i = i + 1
        delete([/idx_rcp, idx_hist, var_diff/])
      end do  ; per
    end do  ; rcp
  end do  ; seas
  to_plot!0 = "panel"
  to_plot&panel = diag_script_info@label

  to_plot!1 = lev_name
  to_plot&$lev_name$ = var_reg&$lev_name$
  if (isatt(diag_script_info, "base_cn")) then
    to_plot_base!1 = lev_name
    to_plot_base&$lev_name$ = var_reg&$lev_name$
  end if

  to_plot!2 = "lat"
  to_plot&lat = var_reg&lat
  if (isatt(diag_script_info, "base_cn")) then
    copy_VarMeta(to_plot, to_plot_base)
  end if
  copy_VarMeta(to_plot, to_plot_signif)
  copy_VarMeta(to_plot, to_plot_not_signif)

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
    to_plot@var_long_name = "change in " + variable_info@long_name
  end if
  if(isatt(variable_info, "units")) then
    to_plot@var_units     = variable_info@units
  end if

  if (isatt(diag_script_info, "base_cn")) then
    if(isatt(to_plot_base, "diag_script")) then  ; Add to existing entries
      temp = to_plot_base@diag_script
      delete(to_plot_base@diag_script)
      to_plot_base@diag_script = array_append_record(temp, (/diag_script/), 0)
      delete(temp)
    else
      to_plot_base@diag_script = (/diag_script/)
    end if
    to_plot_base@var = var0 + "_base"
    to_plot_base@var_long_name = variable_info@long_name
    to_plot_base@var_units     = variable_info@units
  end if
  if(isatt(to_plot_signif, "diag_script")) then  ; Add to existing entries
    temp = to_plot_signif@diag_script
    delete(to_plot_signif@diag_script)
    to_plot_signif@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot_signif@diag_script = (/diag_script/)
  end if
  to_plot_signif@var = "signif"
  to_plot_signif@var_long_name = "significance: mean change larger than" +\
                                 " 2*natvar and 90% of models have the same"\
                                 + " sign"
  to_plot_signif@var_units     = 1

  if(isatt(to_plot_not_signif, "diag_script")) then  ; Add to existing entries
    temp = to_plot_not_signif@diag_script
    delete(to_plot_not_signif@diag_script)
    to_plot_not_signif@diag_script = array_append_record(temp, \
                                                         (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot_not_signif@diag_script = (/diag_script/)
  end if
  to_plot_not_signif@var = "not_signif"
  to_plot_not_signif@var_long_name = "not significant: mean change smaller"\
                                     + " than 1*natvar"
  to_plot_not_signif@var_units     = 1

  model_number!0 = "panel"
  model_number&panel = diag_script_info@label
  if(isatt(model_number, "diag_script")) then  ; Add to existing entries
    temp = model_number@diag_script
    delete(model_number@diag_script)
    model_number@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    model_number@diag_script = (/diag_script/)
  end if
  model_number@var = "model_nr"
  model_number@var_long_name = "number of models in mmm"
  model_number@var_units     = 1

  ; ###########################################
  ; # Output to netCDF                        #
  ; ###########################################
  ; Output (controlled by diag_script_info)
  log_info("Write to netcdf")
  to_plot@ncdf = variable_info@diagnostic + ".nc"
  ncdf_file = work_dir + "/" + to_plot@ncdf
  ncdf_file@existing = "overwrite"
  ; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)
  ncdf_file@existing = "append"
  ncdf_outfile = ncdf_write(to_plot_signif, ncdf_file)
  ncdf_outfile = ncdf_write(to_plot_not_signif, ncdf_file)
  ncdf_outfile = ncdf_write(model_number, ncdf_file)
  if ((isatt(diag_script_info, "base_cn")) .and. \
      (diag_script_info@base_cn .eq. True)) then
    ncdf_outfile = ncdf_write(to_plot_base, ncdf_file)
  end if

end
