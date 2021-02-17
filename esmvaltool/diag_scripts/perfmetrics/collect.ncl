; #############################################################################
; WRAPPER SCRIPT FOR COLLECTING AND PLOTTING PRECALCULATED METRICS
; Author: Franziska Frank (DLR, Germany)
; ESMVal project
; #############################################################################
; Description
;    Collects metrics previously calculated by main.ncl and passes them to the
;    portrait_plot or to the taylor function
;
; Required diag_script_info attributes
;    metric: selected metric (RMSD, BIAS or taylor)
;    label_bounds (for RMSD and BIAS): min and max of the labelbar
;    label_scale (for RMSD and BIAS): bin width of the labelbar
;    colormap (for RMSD and BIAS): color table of the labelbar
;
; Optional diag_script_info attributes
;    label_lo: add lower triange for values outside range
;    label_hi: add upper triange for values outside range
;    cm_interval: min and max color of the color table
;    cm_reverse: reverse color table
;    sort: sort datasets in alphabetic order (excluding MMM)
;    diag_order: sort diagnostics in the specified order
;    title: plot title
;    scale_font: scaling factor applied to the default font size
;    disp_values: switch on/off the grading values on the plot
;    disp_rankings: switch on/off the rankings on the plot
;    rank_order: display rankings in increasing (1) or decreasing (-1) order
;
; Caveats
;
; Modification history
;    20190315-hassler_birgit: added SMPI case
;    20180508-righi_mattia: renamed and generalized
;    20151027-lauer_axel: moved call to 'write_references' to the beginning
;                         of the code
;    20150325-lauer_axel: modified reference tags used for acknowledgements
;                         (projects, observations, etc)
;    20150109-gottschaldt_klaus-dirk: distinguish RMSD & BIAS in output file
;                                     naming.
;    20140620-gottschaldt_klaus-dirk: modified output file naming (variable
;                                     plotpath) and added option for explicit
;                                     plot title.
;    20140204-frank_franziska: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/portrait_plot.ncl"
load "$diag_scripts/shared/plot/taylor_plot.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")
  diag_script_base = basename(DIAG_SCRIPT)

  ; Check file type
  file_type = config_user_info@output_file_type
  if(ismissing(file_type)) then
    file_type = "ps"
  end if

  ; List of grading files
  input_files = diag_script_info@input_files + "/" + \
    diag_script_info@metric + ".nc"
  input_files := tostring(input_files)

  ; Filter out non-existing files
  input_files := input_files(ind(isfilepresent(input_files)))
  data_files = new(dimsizes(input_files), string)

  ; Loop over files in the list, read and append data
  do ii = 0, dimsizes(input_files) - 1

    f = addfile(input_files(ii), "r")
    data_files(ii) = tostring(f->temp_list)

    if (diag_script_info@metric.eq."SMPI") then
      data_temp_smpi = ncdf_read(data_files(ii), "performance_index")

      if (.not.isdefined("data_all_smpi")) then

        dim_temp_smpi = dimsizes(data_temp_smpi)
        n = dimsizes(dim_temp_smpi)
        dim_data_smpi = new(n + 1, integer)
        dim_data_smpi(0) = dimsizes(input_files)
        dim_data_smpi(1:n) = dim_temp_smpi
        vars = new(dimsizes(input_files), string)
        data_all_smpi = new(dim_data_smpi, float)
        delete(dim_data_smpi)
        data_all_smpi(0, :, :) = data_temp_smpi
        data_all_smpi!0 = "vars"
        data_all_smpi!1 = "bootstrap_member"
        data_all_smpi!2 = "models"
        data_all_smpi&models = data_temp_smpi&models
        data_all_smpi&bootstrap_member = data_temp_smpi&bootstrap_member

      end if

      ; Make sure dataset coordinate is consistent
      consistent = False
      if (dimsizes(data_temp_smpi&models).eq. \
          dimsizes(data_all_smpi&models)) then
        if (all(data_temp_smpi&models.eq.data_all_smpi&models)) then
          consistent = True
        end if
      end if

      ; Append data
      if (consistent) then
        data_all_smpi(ii, :, :) = (/data_temp_smpi/)
        copy_VarCoords(data_temp_smpi, data_all_smpi(ii, :, :))
      else
        ; Loop over datasets in new data entry
        do imod_temp = 0, dimsizes(data_temp_smpi&models) - 1
          ; If current dataset is not already part of the model coordinate
          if (.not.any(data_temp_smpi&models(imod_temp) .eq. \
              data_all_smpi&models)) then
            ; Append record for dataset(imod)
            data_new = extend_var_at(data_all_smpi, 2, \
                                     dimsizes(data_all_smpi&models))
            data_new(ii, :, dimsizes(data_all_smpi&models)) = \
              (/data_temp_smpi(:, imod_temp)/)
            data_new&models(dimsizes(data_all_smpi&models)) = \
              (/data_temp_smpi&models(imod_temp)/)
            delete(data_all_smpi)
            data_all_smpi = data_new
            delete(data_new)
          else
            ; Loop over datasets of data
            do imod = 0, dimsizes(data_all_smpi&models) - 1
              ; if new dataset is identical to current entry, write data entry
              if (data_all_smpi&models(imod).eq. \
                  data_temp_smpi&models(imod_temp)) then
                data_all_smpi(ii, :, imod) = (/data_temp_smpi(:, imod_temp)/)
                copy_VarCoords(data_temp_smpi(:, imod_temp), \
                               data_all_smpi(ii:ii, :, imod))
              end if
            end do
          end if
        end do
      end if

      vars(ii) = data_temp_smpi@invar
      delete(data_temp_smpi)

    else  ; metric not SMPI

      data_temp = ncdf_read(data_files(ii), "grade")

      ; Change to 3 dimensional
      if (dimsizes(dimsizes(data_temp)).eq.2) then
        dim_temp = array_append_record(dimsizes(data_temp), 2, 0)
        temp = new(dim_temp, typeof(data_temp))
        temp(:, :, 0) = data_temp
        data_temp := temp
        delete(temp)
      end if

      ; Create array for collecting all datasets
      if (.not.isdefined("data_all")) then
        dim_temp = dimsizes(data_temp)
        dim_data = (/dimsizes(input_files), dim_temp(1), 2/)
        data_all = new(dim_data, float)
        data_all(0, :, :) = data_temp
        data_all!1 = "models"
        data_all&models = data_temp&models
      end if

      ; Make sure model coordinate is consistent
      consistent = False
      if (dimsizes(data_temp&models).eq.dimsizes(data_all&models)) then
        if (all(data_temp&models.eq.data_all&models)) then
          consistent = True
        end if
      end if

      ; Append data
      if (consistent) then
        data_all(ii, :, :) = (/data_temp/)
        copy_VarCoords(data_temp, data_all(ii:ii, :, :))
      else
        do imod2 = 0, dimsizes(data_temp&models) - 1
          if (.not.any(data_temp&models(imod2).eq.data_all&models)) then
            ; Append record for dataset(imod)
            data_new = extend_var_at(data_all, 1, dimsizes(data_all&models))
            data_new(ii, dimsizes(data_all&models), :) = \
              (/data_temp(0, imod2, :)/)
            data_new&models(dimsizes(data_all&models)) = \
              (/data_temp&models(imod2)/)
            data_all := data_new
            delete(data_new)
          else
            ; Loop over datasets of data
            do imod = 0, dimsizes(data_all&models) - 1
              ; if no dataset is similar to curreny entry, write data entry
              if (data_all&models(imod).eq. data_temp&models(imod2)) then
                data_all(ii, imod, :) = (/data_temp(0, imod2, :)/)
                copy_VarCoords(data_temp(0:0, imod2, :), \
                               data_all(ii:ii, imod, :))
              end if
            end do
          end if
        end do
      end if
      delete(data_temp)
    end if  ; metric distinction end 1
  end do

  if (diag_script_info@metric.eq."SMPI") then

    data_all_smpi&vars = vars
    delete(vars)

    ; Check for incomplete array, SMPI only for datasets with all vars
    if any(ismissing(data_all_smpi)) then
      do ivar = 0, dimsizes(data_all_smpi&vars) - 1
        do ii = 0, dimsizes(data_all_smpi&models) - 1
          if all(ismissing(data_all_smpi(ivar, :, ii))) then
            if (isvar("miss_mods")) then
              newmm = array_append_record(miss_mod, \
                                          data_all_smpi&models(ii), 0)
              delete(miss_mod)
              miss_mod = newmm
              delete(newmm)
            else
              miss_mod = data_all_smpi&models(ii)
            end if
          end if
        end do
        if (isvar("miss_mods")) then
          log_info("Variable " + data_all_smpi&vars(ivar) \
                   + " is missing dataset(s) " \
                   + str_join(miss_mods, ", "))
          delete(miss_mods)
          missing_mods = True
        end if
      end do
      if (isvar("missing_mods")) then
        error_msg("f", DIAG_SCRIPT, "", "SMPI only supported for" \
                  + " datasets containing all variables")
      end if
    end if

    i2mean = dim_avg_n(data_all_smpi, 0)  ; eq. 3
    dims = dimsizes(i2mean)
    nboot = dims(0) - 1
    nmod = dims(1)
    delete(dims)

    data_smpi = new((/dimsizes(data_all_smpi&models), 2/), float)
    data_smpi!0 = "models"
    data_smpi&models = data_all_smpi&models
    data_smpi(:, 0) = (/i2mean(0, :)/)  ; clim. mean

    i5 = toint(0.05 * nboot)
    i95 = toint(0.95 * nboot)

    do imod = 0, nmod - 1
      data_sorted = i2mean(1:nboot, imod)
      qsort(data_sorted)
      data_smpi(imod, 1) = 0.5 * (data_sorted(i95) - data_sorted(i5))
      delete(data_sorted)
    end do

    data_smpi!1 = "statistic"
    data_smpi&statistic = (/"SMPI", "95_conf"/)
    data_smpi@diag_script = DIAG_SCRIPT
    data_smpi@var = "SMPI"
    data_smpi@ensemble_name = data_all_smpi@ensemble_name

    ; Sort datasets in alphabetical order (exclude MMM and set to end)
    pid = sort_alphabetically(data_smpi&models, \
                              ind(data_smpi&models.eq."Mean-model".or. \
                                  data_smpi&models.eq."Median-model"), "end")
    tmp = data_smpi(pid, :)
    tmp&models = data_smpi&models(pid)
    delete(data_smpi)
    data_smpi = tmp
    delete(tmp)
    delete(pid)

    ; Define output filenames
    workpath = config_user_info@work_dir + "SMPI.nc"
    plotpath = config_user_info@plot_dir + "SMPI"

    ; Write output
    if (diag_script_info@write_netcdf) then
      system("mkdir -p " + config_user_info@work_dir)
      ncdf_outfile = ncdf_write(data_smpi, workpath)
    end if

    ; Attach plotting options
    copy_VarAtts(diag_script_info, data_smpi)
    if (isatt(diag_script_info, "title")) then
      data_smpi@res_tiMainString = diag_script_info@title
    end if

  else  ; if metric not SMPI

    ; Reduce size if all entries have only one reference
    if (all(ismissing(data_all(:, :, 1)))) then
      data_new = data_all(:, :, 0)
      delete(data_all)
      data_all = data_new
      delete(data_new)
    end if
    data_all@var = "grade"
    if (isatt(data_all, "reference")) then
      delete(data_all@reference)
    end if

    ; Sort datasets in alphabetical order, excluding multi-model mean/median
    ; which are placed at the beginning
    if (isatt(diag_script_info, "sort")) then
      if (diag_script_info@sort) then
        idx = ind(data_all&models.eq."MultiModelMean" .or. \
                  data_all&models.eq."MultiModelMedian")
        pid = sort_alphabetically(data_all&models, idx, "begin")
        if (dimsizes(dimsizes(data_all)).eq.3) then
          data_all := data_all(:, pid, :)
        else
          data_all := data_all(:, pid)
        end if
        delete(pid)
      end if
    end if

    ; Sort diagnostics in the order specified in the settings
    if (isatt(diag_script_info, "diag_order")) then
      l_ok = True
      if (dimsizes(data_all&diagnostics).ne. \
          dimsizes(diag_script_info@diag_order)) then
        error_msg("w", DIAG_SCRIPT, "", "specified order of diagnostics " + \
                  "cannot be applied, number of diagnostics does not match")
        l_ok = False
      end if
      pid = new(dimsizes(diag_script_info@diag_order), integer)
      do ii = 0, dimsizes(diag_script_info@diag_order) - 1
        tmp = ind(data_all&diagnostics.eq.diag_script_info@diag_order(ii))
        if (any(ismissing(tmp)) .or. dimsizes(tmp).gt.1) then
          error_msg("w", DIAG_SCRIPT, "", "specified order of diagnostics " + \
                    "cannot be applied, invalid entry in diag_order")
          break
        end if
        pid(ii) = tmp
        delete(tmp)
      end do
      if (l_ok) then
        if (dimsizes(dimsizes(data_all)).eq.3) then
          data_all := data_all(pid, :, :)
        else
          data_all := data_all(pid, :)
        end if
      end if
    end if

    ; Define output filenames
    if (dimsizes(data_all&diagnostics).gt.1) then
      diags = data_all&diagnostics(0) + "_to_" + \
        data_all&diagnostics(dimsizes(data_all&diagnostics) - 1)
    else
      diags = str_concat(data_all&diagnostics)
    end if
    workpath = config_user_info@work_dir + diags
    plotpath = config_user_info@plot_dir + diags
    if (isatt(data_all, "metric")) then
      workpath = workpath + "_" + data_all@metric + ".nc"
      plotpath = plotpath + "_" + data_all@metric
    end if
    delete(diags)

    ; Write output
    if (diag_script_info@write_netcdf) then
      system("mkdir -p " + config_user_info@work_dir)
      ncdf_outfile = ncdf_write(data_all, workpath)
    end if

    ; Attach plotting options
    copy_VarAtts(diag_script_info, data_all)
    if (isatt(diag_script_info, "title")) then
      data_all@res_tiMainString = diag_script_info@title
    end if
    if (diag_script_info@metric.eq."taylor") then
      data_all@res_varLabels = data_all&models
      data_all@res_caseLabels = data_all&diagnostics
    end if

  end if

  ; Create workspace
  system("mkdir -p " + config_user_info@plot_dir)
  wks = gsn_open_wks(file_type, plotpath)
  if (diag_script_info@metric.eq."taylor") then
    wks@legendfile = plotpath + "_legend"
  end if

  ; Call plotting function
  if (diag_script_info@metric.eq."taylor") then
    plot = taylor_plot(wks, data_all, "grade")
  elseif (diag_script_info@metric.eq."SMPI") then
    plot = circle_plot(wks, data_smpi, "performance_index", \
                       data_smpi@ensemble_name)
  else
    plot = portrait_plot(wks, data_all, "grade")
  end if

  ; Call provenance logger
  if (diag_script_info@metric.eq."RMSD") then
    caption = "RMSD performance metric"
    statistics = (/"rmsd"/)
    authors = (/"frank_franziska", "righi_mattia", "eyring_veronika"/)
    plottype = "portrait"
    references = (/"righi15gmd", "gleckler08jgr"/)
  elseif (diag_script_info@metric.eq."BIAS") then
    caption = "Bias performance metric"
    statistics = (/"diff"/)
    authors = (/"frank_franziska", "righi_mattia", "eyring_veronika"/)
    plottype = "portrait"
    references = (/"righi15gmd", "gleckler08jgr"/)
  elseif (diag_script_info@metric.eq."taylor") then
    caption = "Taylor diagram"
    statistics = (/"rmsd", "corr"/)
    authors = (/"frank_franziska", "righi_mattia", "eyring_veronika"/)
    plottype = "taylor"
    references = (/"righi15gmd", "gleckler08jgr"/)
  elseif (diag_script_info@metric.eq."SMPI") then
    caption = \
      "Performance index I2. Similar to Figure 1 of Reichler and Kim (2008)."
    statistics = "smpi"
    authors = (/"gier_bettina", "hassler_birgit"/)
    plottype = "circle"
    references = (/"rk2008bams"/)
  else
    error_msg("f", DIAG_SCRIPT, "", "cannot add provenance information " + \
              "for metric " + diag_script_info@metric)
  end if
  domains = (/"global"/)

  ; Call provenance logger
  log_provenance(ncdf_outfile, plotpath, caption, statistics, domains, \
                 plottype, authors, references, data_files)

  leave_msg(DIAG_SCRIPT, "")

end
