; #############################################################################
; Multi-model mean change map with significance
; Author: Ruth Lorenz (ETH, Switzerland)
; CRESCENDO project
; #############################################################################
;
; Caveats
;    Only handles one variable at a time. Appends multiple output variables
;    (variable to_plot, signficance to_plot_signif and not significant
;    to_plot_not_signif and number of model in mmm to model_nr)
;    to the same file, if the file already exists overwrites
;
; Modification history
;    20181307-A_lorenz_ruth: Ported to version 2
;    20170120-A_lorenz_ruth: separated plotting from calculations
;    20161027-A_lorenz_ruth: added stippling and hatching for significance
;    20161024-A_lorenz_ruth: adjusted to plot multiple figures
;                            with one diagnostic script
;    20160621-A_lorenz_ruth: adapted to ESMValTool
;    20130501-A_sedlacek_jan: written for IPCC AR5.
;
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/ensemble.ncl"
load "$diag_scripts/shared/scaling.ncl"

load "$diag_scripts/shared/plot/style.ncl"

begin

  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(DIAG_SCRIPT, "")

  ; 'datasets', 'variables' are fetched from the above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)

  ; Save some dataset attributes as arrays for index-based selection below
  dataset_names = metadata_att_as_array(info_items, "dataset")
  dataset_exps = metadata_att_as_array(info_items, "exp")
  dataset_ensembles = metadata_att_as_array(info_items, "ensemble")
  dataset_startyears = metadata_att_as_array(info_items, "start_year")

  ; Check required diag_script_info attributes
  req_atts = (/"scenarios", "periods", "time_avg", "label"/)
  exit_if_missing_atts(diag_script_info, req_atts)

end

begin

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
  if(isatt(variable_info[0], "long_name")) then
    MyParam = variable_info[0]@long_name
    log_info(" MyParam = " + MyParam)
  else
    error_msg("f", DIAG_SCRIPT, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; ###########################################
  ; # Get data and average time               #
  ; ###########################################

  ; Get data from first dataset
  imod = 0  ; NCL array indicies start from zero
  log_debug("processing " + info_items[imod]@dataset + "_" \
            + info_items[imod]@exp + "_" \
            + info_items[imod]@ensemble)

  ; See ./interface_scripts/data_handling.ncl
  A0 = read_data(info_items[imod])

  ; Check dimensions
  dims = getvardims(A0)
  if (dimsizes(dims) .lt. 2) then
    error_msg("f", DIAG_SCRIPT, "", dimsizes(dims) + " dimensions, " + \
              "need 2 or 3")
  end if
  idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
  end if
  lat_name = dims(idx)
  idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
  end if
  lon_name = dims(idx)

  if (isatt(diag_script_info, "seasons")) then
    dim_seas = dimsizes(diag_script_info@seasons)
  else
    dim_seas = 1
  end if

  ; Check if diag_script_info@label has correct size, otherwise error message
  dim_periods = dimsizes(diag_script_info@periods)
  dim_scenarios = dimsizes(diag_script_info@scenarios)
  dim_label = dim_seas * dim_scenarios * dim_periods
  if (dim_label .ne. dimsizes(diag_script_info@label)) then
    error_msg("f", DIAG_SCRIPT, "", "wrong number of labels, needs to be " + \
              "seasons * scenarios * periods")
  end if

  annots = project_style(info_items, diag_script_info, "annots")
  Fill = default_fillvalue(typeof(A0))
  if (isatt(diag_script_info, "iavmode")) then
    if (diag_script_info@iavmode .eq. "each") then
      iav = new((/dim_MOD, dim_seas, dimsizes(A0&lat), \
                 dimsizes(A0&lon)/), typeof(A0), Fill)
    end if
  end if
  data2 = new((/dim_MOD, dim_seas, dimsizes(A0&lat), \
               dimsizes(A0&lon)/), typeof(A0), Fill)

  do imod = 0, dim_MOD - 1
    log_info(info_items[imod]@dataset)
    ; Average over time
    if (imod .eq. 0) then
      if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
          (dim_seas .ne. 4)) then
        tmp_seas = time_operations(A0, -1, -1, "average", \
                                   diag_script_info@time_avg, True)
        data1 = tmp_seas(diag_script_info@seasons, :, :)
        delete(tmp_seas)
      else
        data1 = time_operations(A0, -1, -1, "average", \
                                diag_script_info@time_avg, True)
      end if
    else
      if isvar((/"A1"/)) then
        delete(A1)
      end if
      A1 = read_data(info_items[imod])
      dims = getvardims(A1)
      idx = ind(dims .eq. "lat" .or. dims .eq. "rlat" .or. dims .eq. "j")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
      end if
      lat_name = dims(idx)
      idx = ind(dims .eq. "lon" .or. dims .eq. "rlon" .or. dims .eq. "i")
      if (ismissing(idx)) then
        error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
      end if
      lon_name = dims(idx)
      ; Average over time
      if ((diag_script_info@time_avg .eq. "seasonalclim") .and. \
          (dim_seas .ne. 4)) then
        tmp_seas = time_operations(A1, -1, -1, "average", \
                                   diag_script_info@time_avg, True)
        data1 = tmp_seas(diag_script_info@seasons, :, :)
        delete(tmp_seas)
      else
        data1 = time_operations(A1, -1, -1, "average", \
                                diag_script_info@time_avg, True)
      end if
    end if
    if (dim_seas .eq. 1) then
      data2(imod, 0, :, :) = data1
    else
      data2(imod, :, :, :) = data1
    end if

    if (isatt(diag_script_info, "iavmode")) then
      if (diag_script_info@iavmode .eq. "each") then
        ; #########################################################
        ; # Read natural variability for stippling for each model #
        ; #########################################################
        log_info("Read natural variability for stippling for each model")
        iav_file = addfile(input_dir + "/" + "IAV_piControl_" + \
                           info_items[imod]@dataset + "_" +  \
                           info_items[imod]@ensemble + "_" + \
                           "_" + var0 + "_" + diag_script_info@time_avg + \
                           ".nc", "r")
        if ((dim_seas .eq. 1) .and. \
            (diag_script_info@time_avg .eq. "annualclim")) then
          iav(imod, 0, :, :) = iav_file->iav
        elseif ((dim_seas .eq. 1) .and. \
                (diag_script_info@time_avg .eq. "seasonalclim")) then
          tmp_iav = iav_file->iav
          iav(imod, 0, :, :) = tmp_iav(diag_script_info@seasons, :, :)
        else
          iav(imod, :, :, :) = iav_file->iav
        end if
      end if
    end if
  end do

  if (isatt(diag_script_info, "iavmode")) then
    if (diag_script_info@iavmode .eq. "each") then
      ; Calculate multi-model mean of iav
      log_info("Calculate multi-model mean of natural variability")
      natvar_mmm = rm_single_dims(dim_avg_n(iav, 0) * sqrt(2.))
    elseif (diag_script_info@iavmode .eq. "mmm") then
      log_info("Read precalculated natural variability for multi-model")
      natvar_file = addfile(input_dir + "/" + "IAV_mmm_piControl_" + \
                            var0 + "_" + \
                            diag_script_info@time_avg + ".nc", "r")
      natvar_mmm = natvar_file->iav
      if ((dim_seas .ne. 4) .and. \
          (diag_script_info@time_avg .ne. "annualclim")) then
        tmp = natvar_mmm
        delete(natvar_mmm)
        natvar_mmm = rm_single_dims(tmp(diag_script_info@seasons, :, :))
      elseif (diag_script_info@time_avg .eq. "annualclim") then
        tmp = natvar_mmm
        delete(natvar_mmm)
        natvar_mmm = rm_single_dims(tmp)
      end if
    else
      error_msg("f", DIAG_SCRIPT, "", "Error: This iavmode is not know")
    end if
  else
    ; Read already calculated mmm iav as default
    log_info("Read precalculated natural variability for multi-model")
    natvar_file = addfile(input_dir + "/" + "IAV_mmm_piControl_" + \
                          var0 + "_" + \
                          diag_script_info@time_avg + ".nc", "r")
    natvar_mmm = natvar_file->iav
    if ((dim_seas .ne. 4) .and. \
        (diag_script_info@time_avg .ne. "annualclim")) then
      tmp = natvar_mmm
      delete(natvar_mmm)
      natvar_mmm = rm_single_dims(tmp(diag_script_info@seasons, :, :))
    elseif (diag_script_info@time_avg .eq. "annualclim") then
      tmp = natvar_mmm
      delete(natvar_mmm)
      natvar_mmm = rm_single_dims(tmp)
    end if
  end if
  data2!0 = "models"
  data2!1 = "season"
  data2!2 = "lat"
  data2!3 = "lon"
  data2&models = annots
  if (isatt(diag_script_info, "seasons")) then
    data2&season = (/ispan(0, dim_seas - 1, 1)/)
  else
    data2&season = (/0/)
  end if

  ; ##############################################
  ; # Calculate change from hist to periods      #
  ; ##############################################
  if (isatt(diag_script_info, "reference_run")) then
    reference_run_name = diag_script_info@reference_run
  else
    log_info("No reference run name given, set to default 'historical'")
    reference_run_name = "historical"
  end if

  ; Loop over rcps and periods
  dim_scen = dimsizes(diag_script_info@scenarios)
  dim_per = dimsizes(diag_script_info@periods)

  to_plot = new((/dim_scen * dim_per * dim_seas, dimsizes(A0&lat), \
                 dimsizes(A0&lon)/), typeof(data1), Fill)
  if (isatt(diag_script_info, "percent")) then
    percent = diag_script_info@percent
    ref = new((/dim_scen * dim_per * dim_seas, dimsizes(A0&lat), \
               dimsizes(A0&lon)/), typeof(data1), Fill)
  else
    percent = 0
  end if
  to_plot_signif = new((/dim_scen * dim_per * dim_seas, dimsizes(A0&lat), \
                        dimsizes(A0&lon)/), typeof(data1), Fill)
  to_plot_not_signif = new((/dim_scen * dim_per * dim_seas, dimsizes(A0&lat), \
                            dimsizes(A0&lon)/), typeof(data1), Fill)
  model_number = new((/dim_scen * dim_per * dim_seas/), integer)
  i = 0
  do seas = 0, dim_seas - 1  ; loop over seasons
    do rcp = 0, dim_scen - 1  ; loop over rcps
      do per = 0, dim_per - 1  ; loop over periods
        idx_rcp = ind(dataset_exps .eq. diag_script_info@scenarios(rcp) .and. \
                      dataset_startyears .eq. diag_script_info@periods(per))
        proj = data2(idx_rcp, seas, :, :)

        ; Cut low values
        proj = where(abs(proj) .gt. 1.e-14, proj, proj@_FillValue)
        proj_avg = dim_avg_n_Wrap(proj, 0)

        ; Find historical runs from same datasets as in rcp
        do jj = 0, dimsizes(idx_rcp) - 1
          tmp_idx = ind(dataset_names .eq. dataset_names(idx_rcp(jj)) .and. \
                        dataset_exps .eq. reference_run_name .and. \
                        dataset_ensembles .eq. dataset_ensembles(idx_rcp(jj)))
          if (isdefined("idx_hist")) then
            idx_hist := array_append_record(idx_hist, tmp_idx, 0)
          else
            idx_hist = tmp_idx
          end if
          delete(tmp_idx)
        end do

        base = data2(idx_hist, seas, :, :)

        ; Cut low values
        base = where(abs(base) .gt. 1.e-14, base, base@_FillValue)
        base_avg = dim_avg_n_Wrap(base, 0)

        var_diff = proj - base
        delete([/base, proj/])
        to_plot(i, :, :) = proj_avg - base_avg
        if (isatt(diag_script_info, "percent")) then
          ref(i, :, :) = base_avg
        end if

        ; Determine sigma and sign for significance
        var_pos_signif = dim_num_n(where(var_diff .gt. 0., 1., \
                                         to_plot@_FillValue) .eq. 1., 0)
        var_neg_signif = dim_num_n(where(var_diff .lt. 0., 1., \
                                         to_plot@_FillValue) .eq. 1., 0)
        var_pos_signif_tmp = \
          var_pos_signif / where((var_pos_signif + var_neg_signif) \
                                 .ne. 0., var_pos_signif + \
                                 var_neg_signif, to_plot@_FillValue)
        var_neg_signif_tmp = \
          var_neg_signif / where((var_pos_signif + var_neg_signif) \
                                 .ne. 0., var_pos_signif + \
                                 var_neg_signif, to_plot@_FillValue)
        signif_and_sign = where(var_pos_signif_tmp .ge. 0.9 .or. \
                                var_neg_signif_tmp .ge. 0.9, 1., 0.)
        signif_and_sign = where(ismissing(signif_and_sign), 0.,\
                                signif_and_sign)
        if (dim_seas .ne. 1) then
          sigma_and_sign = where(abs(to_plot(i, :, :)) .gt. \
                                 abs(2. * natvar_mmm(seas, :, :)), 1., 0.)
          to_plot_not_signif(i, :, :) = where(abs(to_plot(i, :, :)) .lt. \
                                              abs(natvar_mmm(seas, :, :)),\
                                              1., 0.)
        else
          sigma_and_sign = where(abs(to_plot(i, :, :)) .gt. \
                                 abs(2. * natvar_mmm), 1., 0.)
          to_plot_not_signif(i, :, :) = where(abs(to_plot(i, :, :)) .lt. \
                                              abs(natvar_mmm),\
                                              1., 0.)
        end if
        to_plot_signif(i, :, :) = where(signif_and_sign .eq. 1, \
                                        sigma_and_sign, 0.)
        model_number(i) = dimsizes(idx_rcp)
        i = i + 1
        delete([/idx_hist, idx_rcp, var_diff/])
      end do  ; per
    end do  ; rcp
  end do  ; seas

  to_plot!0 = "panel"
  to_plot&panel = diag_script_info@label
  to_plot!1 = "lat"
  to_plot&lat = A0&lat
  to_plot!2 = "lon"
  to_plot&lon = A0&lon
  copy_VarMeta(to_plot, to_plot_signif)
  copy_VarMeta(to_plot, to_plot_not_signif)
  if (isatt(diag_script_info, "percent")) then
    percent = diag_script_info@percent
    if (percent .eq. 1) then
      to_plot = 100 * to_plot / where(ref .ne. 0., ref, ref@_FillValue)
      to_plot@units     = "%"
      delete(ref)
    end if
  end if

  ; ###########################################
  ; # Other Metadata: diag_script, var        #
  ; ###########################################

  ; Add to to_plot, as attributes without prefix
  if (isatt(to_plot, "diag_script")) then  ; add to existing entries
    temp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot@diag_script = (/DIAG_SCRIPT/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if (isatt(variable_info[0], "long_name")) then
    to_plot@var_long_name = variable_info[0]@long_name
  end if

  ; Check units and adjust for plotting if necessary
  if (percent .eq. 0) then
    to_plot@units = variable_info[0]@units
  end if
  ; difference in K is the same as difference in degC, only change label
  if (to_plot@units.eq."K") then
    to_plot@units = "degC"
  end if
  if (isatt(diag_script_info, "plot_units")) then
    to_plot = convert_units(to_plot, diag_script_info@plot_units)
  end if

  if(isatt(to_plot_signif, "diag_script")) then  ; add to existing entries
    temp = to_plot_signif@diag_script
    delete(to_plot_signif@diag_script)
    to_plot_signif@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot_signif@diag_script = (/DIAG_SCRIPT/)
  end if

  to_plot_signif@var = "signif"
  to_plot_signif@var_long_name = "significance: mean change larger than " \
                                 + "2*natvar and 90% of models have the " \
                                 + "same sign"
  to_plot_signif@units     = 1

  if(isatt(to_plot_not_signif, "diag_script")) then
    temp = to_plot_not_signif@diag_script
    delete(to_plot_not_signif@diag_script)
    to_plot_not_signif@diag_script = array_append_record(temp, \
                                                         (/DIAG_SCRIPT/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot_not_signif@diag_script = (/DIAG_SCRIPT/)
  end if
  to_plot_not_signif@var = "not_signif"
  to_plot_not_signif@var_long_name = "not significant: mean change " \
                                     + "smaller than 1*natvar"
  to_plot_not_signif@units     = 1

  model_number!0 = "panel"
  model_number&panel = diag_script_info@label
  if(isatt(model_number, "diag_script")) then  ; add to existing entries
    temp = model_number@diag_script
    delete(model_number@diag_script)
    model_number@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
    delete(temp)
  else  ; Add as new attribute
    model_number@diag_script = (/DIAG_SCRIPT/)
  end if
  model_number@var = "model_nr"
  model_number@var_long_name = "number of models in mmm"
  model_number@units     = 1

  ; ###########################################
  ; # Optional output to netCDF               #
  ; ###########################################
  ; Optional output (controlled by diag_script_info)
  to_plot@ncdf = variable_info[0]@diagnostic + ".nc"
  ncdf_file = work_dir + "/" + to_plot@ncdf
  ncdf_file@existing = "overwrite"
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)
  ncdf_file@existing = "append"
  ncdf_outfile = ncdf_write(to_plot_signif, ncdf_file)
  ncdf_outfile = ncdf_write(to_plot_not_signif, ncdf_file)
  ncdf_outfile = ncdf_write(model_number, ncdf_file)

  leave_msg(DIAG_SCRIPT, "")
end
