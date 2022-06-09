; #############################################################################
; PROCEDURE FOR THE CYCLE-LATLON PLOT OF THE PERFORMANCE METRICS
; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
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

print("ref_ind = " + ref_ind)
print(datasetnames(modidx))

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

    ; Store all models in one array
    if(ii_modidx .eq. 0) then
      var_all_mod = new(array_append_record(dimsizes(var_all), \
                                            dimsizes(grading&models), 0), float)
      var_all_mod!4 = "models"
      var_all_mod&models = grading&models
print(grading&models)
printVarSummary(var_all_mod)
      copy_VarAtts(diag_script_info, var_all_mod)
    end if
    if (l_altern) then
      if (imod .ne. ref_ind .and. imod .ne. alt_ind) then
        var_all_mod(:,:,:,:,imod) = var_all
      end if
    else
print("imod = " + imod)
      if (imod.ne.ref_ind) then
        var_all_mod(:,:,:,:,imod) = var_all
      end if
    end if 
    ;if ((imod.ne.ref_ind) .and. (.not. (l_altern .and. imod.eq.alt_ind))) then
    ;    var_all_mod(:,:,:,:,imod) = var_all
    ;end if

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
                   (/"frank_franziska", "righi_mattia", "eyring_veronika"/), \
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
        grading(met, 0, imod, 0) = \
          calculate_metric(var_all(:, :, :, 0), var_ref(:, :, :, 0), \
                           diag_script_info@metric(met))
print("grading = " + grading(met, 0, imod, 0))
        if (l_altern) then
          grading(met, 0, imod, 1) = \
            calculate_metric(var_all(:, :, :, 0), var_alt(:, :, :, 0), \
                             diag_script_info@metric(met))
        end if

      end if

    end do

  end do

print("grading&models = " + grading&models)
print("grading = " + grading(0, 0, :, 0))


  if (any(diag_script_info@metric.eq."RMSD")) then

    ;; *************************************
    ;; calculate means over ensemble members
    ;; *************************************

    imet = ind(diag_script_info@metric.eq."RMSD")

    datasetnames := grading&models
    projectnames := remove_index(projectnames, ref_inds)

    ndim := dimsizes(datasetnames)
    ;ensemble_assign = new(ndim, integer)

    ;enum = 1
    ;temp = datasetnames(0)
    ;do gg = 0, ndim - 1
    ;  if (temp.eq.datasetnames(gg)) then ; does the model name match what's in temp?
    ;    ensemble_assign(gg) = enum ; if so, assign it the same number
    ;  else
    ;    enum = enum + 1              ; if not, assign it the next number
    ;    ensemble_assign(gg) = enum
    ;    temp = datasetnames(gg)
    ;  end if
    ;end do
    ;ensemble_assign@model = datasetnames
    ;delete(temp)

    ;grading_em = grading(:, :, :max(ensemble_assign)-1, :)
    ;grading_em = grading_em@_FillValue
    ;datasetnames_em = datasetnames(:max(ensemble_assign)-1)
    ;projectnames_em = projectnames(:max(ensemble_assign)-1)

    ;var_all_mod_em = var_all_mod(:, :, :, :, :max(ensemble_assign)-1)

    ;do gg = 1,max(ensemble_assign)         ; calculate ensemble means
    ;  wind := ind(ensemble_assign.eq.gg)
    ;  if (dimsizes(wind).eq.1) then
    ;    grading_em(imet, 0, gg-1, :) = (/ grading(imet, 0, wind, :) /)
    ;    grading_em&models(gg-1) =  datasetnames(wind)
    ;    datasetnames_em(gg-1) = datasetnames(wind)
    ;    projectnames_em(gg-1) = projectnames(wind)
    ;    var_all_mod_em(:, :, :, :, gg-1) = (/ var_all_mod(:, :, :, :, wind) /)
    ;    var_all_mod_em&models = datasetnames(wind)
    ;  else
    ;    grading_em(imet, 0, gg-1, :) = (/ dim_avg_n(grading(imet, 0, wind, :), 0) /)
    ;    grading_em&models(gg-1) =  datasetnames(wind(0))
    ;    datasetnames_em(gg-1) = datasetnames(wind(0))
    ;    projectnames_em(gg-1) = projectnames(wind(0))
    ;    var_all_mod_em(:, :, :, :, gg-1) = (/ dim_avg_n(var_all_mod(:, :, :, :, wind), 4) /)
    ;    var_all_mod_em&models = datasetnames(wind(0) )
    ;  end if
    ;end do
 
    ;delete(grading)
    ;grading = grading_em
    ;grading@projects = str_join(projectnames_em, " ")
    ;delete(grading_em)

    ;delete(var_all_mod)
    ;var_all_mod = var_all_mod_em
    ;delete(var_all_mod_em)


    ; ******************************************************
    ; calculate multi-model mean and median for each project
    ; ******************************************************

    projects = get_unique_values(projectnames)
    n_proj = dimsizes(projects)

    do iproj = 0, n_proj - 1
      ind_proj := ind(projectnames .eq. projects(iproj))

print("project = " + projects(iproj))

      tmp := dim_avg_n_Wrap(var_all_mod(:, :, :, :, ind_proj), 4)
print("mean_ref")
      mean_ref = calculate_metric(tmp(:, :, :, 0), var_ref(:, :, :, 0), \
                             diag_script_info@metric(imet))
      if (l_altern) then
        mean_alt = calculate_metric(tmp(:, :, :, 0), var_alt(:, :, :, 0), \
                               diag_script_info@metric(imet))
      end if
      tmp = dim_median_n(var_all_mod(:, :, :, :, ind_proj), 4)
print("median_ref")
      median_ref = calculate_metric(tmp(:, :, :, 0), var_ref(:, :, :, 0), \
                             diag_script_info@metric(imet))
      if (l_altern) then
        median_alt = calculate_metric(tmp(:, :, :, 0), var_alt(:, :, :, 0), \
                               diag_script_info@metric(imet))
      end if

      dims = dimsizes(grading)
      tmp := new((/dims(0), dims(1), dims(2) +2, dims(3)/), float)
      copy_VarAtts(grading, tmp)
      tmp(imet, 0, :dims(2)-1, :) = grading(imet, 0, :, :)
      ;tmp(imet, 0, :dims(2)-1, :) = grading
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
      ;grading&models(dims(2)+1) = projects(iproj)+"_median"
      ;grading@projects = grading@projects + " " + projects(iproj) + \
      ;                   " " + projects(iproj)
      delete(tmp)
      ;delete(ind_proj)
    end do

  end if

end
