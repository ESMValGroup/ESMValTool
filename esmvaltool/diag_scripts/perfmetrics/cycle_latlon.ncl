; #############################################################################
; PROCEDURE FOR THE CYCLE-LATLON PLOT OF THE PERFORMANCE METRICS
; Authors: Mattia Righi (DLR, Germany) and Franziska Winterstein (DLR, Germany)
; ESMVal project
; #############################################################################

load "$diag_scripts/shared/plot/style.ncl"

procedure perfmetrics_ptype_script()
begin

  ; Set dataset indexes, with reference and alternative datasets first
  modidx = ispan(0, nDatasets - 1, 1)
  if (l_altern) then
    modidx := array_append_record( \
      (/ref_ind, alt_ind/), \
      modidx(ind(modidx.ne.ref_ind .and. modidx.ne.alt_ind)), 0)
  else
    modidx := array_append_record(ref_ind, modidx(ind(modidx.ne.ref_ind)), 0)
  end if

  annots = project_style(info_items, diag_script_info, "annots")

  ; Loop over datasets
  do ii_modidx = 0, dimsizes(modidx) - 1

    imod = modidx(ii_modidx)

    log_debug("Processing " + datasetnames(imod))

    ; Set path for saving processed data
    system("mkdir -p " + config_user_info@work_dir)
    fname = \
      str_join((/"perfmetrics", "cycle_latlon", annots(imod), var0, \
                 diag_script_info@time_avg, diag_script_info@region/), "_")
    workpath = config_user_info@work_dir + fname + ".nc"

    ; Determine start/end year
    start_year = info_items[imod]@start_year
    end_year = info_items[imod]@end_year

    ; Read data
    var = read_data(info_items[imod])

    if (var0 .eq. "tos" .and. var@units .eq. "K") then
      var = convert_units(var, "degC")
    end if

    ; Convert data to floats if necessary
    if (typeof(var) .eq. "double") then
      var := dble2flt(var)
    end if

    ; Extract region
    var_reg = area_operations(var, region(0), region(1), \
                              region(2), region(3), "extract", True)
    delete(var)

    ; Calculate time average
    var_avg = time_operations(var_reg, start_year, end_year, "average", \
                              diag_script_info@time_avg, True)
    delete(var_reg)

    ; Standard deviation calculation for this ptype is not yet implemented in
    ; diag_scripts/shared/statistics.ncl, but is anyway not required here

    ; Store in global array
    var_all = new(array_append_record(dimsizes(var_avg), 2, 0), float)
    var_all!3 = "statistic"
    var_all&statistic = (/"mean", "stddev"/)
    var_all(:, :, :, 0) = var_avg
    var_all@var = var0
    var_all@diag_script = (/DIAG_SCRIPT/)
    copy_VarAtts(diag_script_info, var_all)
    var_all@ref_model = variable_info[0]@reference_dataset
    delete(var_avg)

    grading_models = remove_index(datasetnames, ref_inds)

    ; Store all models in one array
    if(ii_modidx .eq. 0) then
      var_all_mod = new(array_append_record(dimsizes(var_all), \
                        dimsizes(grading_models), 0), float)
      var_all_mod!4 = "models"
      var_all_mod&models = remove_index(datasetnames, ref_inds)
      copy_VarAtts(diag_script_info, var_all_mod)
    end if
    if (l_altern) then
      if (imod .ne. ref_ind .and. imod .ne. alt_ind) then
        igrad = ind(datasetnames(imod) .eq. grading_models)
        var_all_mod(:, :, :, :, igrad) = var_all
      end if
    else
      if (imod.ne.ref_ind) then
        igrad = ind(datasetnames(imod) .eq. grading_models)
        var_all_mod(:, :, :, :, igrad) = var_all
      end if
    end if

    ; Store reference dataset
    if (imod.eq.ref_ind) then
      var_ref = var_all
    end if

    ; Store alternative dataset
    if (l_altern) then
      if (imod.eq.alt_ind) then
        var_alt = var_all
      end if
    end if

    ; Write output
    var_all@ncdf = workpath
    ncdf = ncdf_write(var_all, workpath)

    ; Call provenance logger
    log_provenance(ncdf, \
                   "n/a", \
                   "n/a", \
                   (/"mean", "stddev"/), \
                   diag_script_info@region, \
                   "other", \
                   (/"winterstein_franziska", "righi_mattia", \
                     "eyring_veronika"/), \
                   (/"righi15gmd", "gleckler08jgr"/), \
                   info_items[imod]@filename)

    ; Calculate grading
    if (.not.diag_script_info@calc_grading) then
      continue
    end if

    ; Loop over requested metrics
    do met = 0, nmetrics - 1

      if (diag_script_info@metric(met).eq."taylor") then

        locidx = ind(datasetnames(imod).eq.taylor&models)
        if (ismissing(locidx)) then
          continue
        end if
        taylor(0, locidx, 0) = \
          calculate_metric(var_all(:, :, :, 0), var_ref(:, :, :, 0), \
                           "stddev_ratio")
        taylor(0, locidx, 1) = \
          calculate_metric(var_all(:, :, :, 0), var_ref(:, :, :, 0), \
                           "correlation")

      elseif (diag_script_info@metric(met).eq."SMPI") then

        locidx = ind(datasetnames(imod).eq.smpi&models)
        if (ismissing(locidx)) then
          continue
        end if
        smpi(:, locidx) = \
          (/calculate_metric(var_all(:, :, :, 0), var_ref(:, :, :, 0), \
                             "SMPI")/)

      else

        if (.not.any(datasetnames(imod).eq.grading&models)) then
          continue
        end if
        grading(met, 0, igrad, 0) = \
          calculate_metric(var_all(:, :, :, 0), var_ref(:, :, :, 0), \
                           diag_script_info@metric(met))
        if (l_altern) then
          grading(met, 0, igrad, 1) = \
            calculate_metric(var_all(:, :, :, 0), var_alt(:, :, :, 0), \
                             diag_script_info@metric(met))
        end if

      end if

    end do

  end do

  if (diag_script_info@calc_grading .and. \
      any(diag_script_info@metric.eq."RMSD") .and. \
      all(grading_models.ne."MultiModelMean") .and. \
      all(grading_models.ne."MultiModelMedian")) then

    ; *************************************
    ; calculate means over ensemble members
    ; *************************************

    imet = ind(diag_script_info@metric.eq."RMSD")

    datasetnames := grading@models
    projectnames := remove_index(projectnames, ref_inds)

    ndim := dimsizes(datasetnames)

    ; ******************************************************
    ; calculate multi-model mean and median for each project
    ; ******************************************************

    projects = get_unique_values(projectnames)
    n_proj = dimsizes(projects)

    do iproj = 0, n_proj - 1

      log_debug("Processing mean of " + projects(iproj))

      ind_proj := ind(projectnames .eq. projects(iproj))

      tmp := dim_avg_n_Wrap(var_all_mod(:, :, :, :, ind_proj), 4)
      mean_ref = calculate_metric(tmp(:, :, :, 0), var_ref(:, :, :, 0), \
                                  diag_script_info@metric(imet))
      if (l_altern) then
        mean_alt = calculate_metric(tmp(:, :, :, 0), var_alt(:, :, :, 0), \
                                    diag_script_info@metric(imet))
      end if
      tmp = dim_median_n(var_all_mod(:, :, :, :, ind_proj), 4)
      median_ref = calculate_metric(tmp(:, :, :, 0), var_ref(:, :, :, 0), \
                                    diag_script_info@metric(imet))
      if (l_altern) then
        median_alt = calculate_metric(tmp(:, :, :, 0), var_alt(:, :, :, 0), \
                                      diag_script_info@metric(imet))
      end if

      dims = dimsizes(grading)
      tmp := new((/dims(0), dims(1), dims(2) + 2, dims(3)/), float)
      copy_VarAtts(grading, tmp)
      tmp(imet, 0, :dims(2)-1, :) = grading(imet, 0, :, :)
      tmp(imet, 0, dims(2), 0) = (/mean_ref/)
      tmp(imet, 0, dims(2)+1, 0) = (/median_ref/)
      if (l_altern) then
        tmp(imet, 0, dims(2), 1) = (/mean_alt/)
        tmp(imet, 0, dims(2)+1, 1) = (/median_alt/)
      end if

      delete(grading)
      grading = tmp
      grading&models(dims(2)) = projects(iproj)+"_mean"
      grading&models(dims(2)+1) = projects(iproj)+"_median"
      grading@projects = grading@projects + " " + projects(iproj) + \
                                            " " + projects(iproj)
      delete(tmp)
    end do

  end if

end
