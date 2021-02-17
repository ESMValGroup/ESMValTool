; #############################################################################
; DIAGNOSTIC PLOT SCRIPT for reproducing IPCC ch. 9 fig. 9.6
; Author: Bettina Gier (DLR, Germany)
; CRESCENDO project
; #############################################################################
;
; Description
;    Calculated centred pattern correlations for annual mean climatologies
;    and plots them. Like IPCC ch. 9 fig 9.6
;
; Required diag_script_info attributes (diagnostics specific)
;
; Optional diag_script_info attributes (diagnostic specific)
;     diag_script_info@diag_order: give order of plotting variables on the
;                                  x-axis
;
; Required variable_info attributes (variable specific)
;     none
;
; Optional variable_info attributes (variable specific)
;     none
;
; Required variable attributes (defined in namelist)
;    reference_dataset: name of reference data set (observations)
;
; Caveats
;     Effect of different regridding methods not yet determined
;
; Modification history
;    20191011-A_bock_lisa: Add customizable order of variables
;    20190205-A_gier_bettina: Adapted to new ncl structures
;    20181101-A_gier_bettina: moved collect from main script
;
; #############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")

  ; Define file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

; -----------------------------------------------------------------------------
; ------------------- Collecting Data -----------------------------------------
; -----------------------------------------------------------------------------

  ; List of correlation files and project names
  file_list = tostring(diag_script_info@input_files) + "/pattern_cor.nc"
  mp_file_list = tostring(diag_script_info@input_files) + "/modprojnames.txt"

  ; Filter non-existing files (mp files only useful if cor file exists)
  file_list := file_list(ind(isfilepresent(file_list)))
  mp_file_list := mp_file_list(ind(isfilepresent(file_list)))

  ; Set up auxiliary variables
  var_collect = new(dimsizes(file_list), string)
  var_diag = new(dimsizes(file_list), string)
  alt_obs = new(dimsizes(file_list), string)

  ; Loop over files in list, read and append data
  do ii = 0, dimsizes(file_list) - 1
    data_temp = ncdf_read(file_list(ii), "cor")
    var_collect(ii) = data_temp@corvar
    var_diag(ii) = data_temp@diagnostics
    alt_obs(ii) = data_temp@alt_obs

    ; Make 2D array to store all data
    if (.not.isdefined("data_all")) then
      data_all = new((/dimsizes(data_temp), dimsizes(file_list)/), float)
      data_all(:, ii) = data_temp
      data_all!0 = "models"
      data_all&models = data_temp&models

      ; Input file list for provenance
      prov_files = str_split(data_temp@input, ",")
    else
      ; If model coordinates are identical
      if (dimsizes(data_temp&models).eq.dimsizes(data_all&models)) \
        .and. all(data_temp&models.eq.data_all&models) then
        data_all(:, ii) = (/data_temp/)
      else
        ; Loop over models in new data entry
        do imod_temp = 0, dimsizes(data_temp&models) - 1
          ; If current model is not already part of the model coordinate
          if (.not.any(data_temp&models(imod_temp) .eq. data_all&models)) then
            ; Append record for model(imod)
            data_new = extend_var_at(data_all, 0, \
                                     dimsizes(data_all&models))
            data_new(dimsizes(data_all&models), ii) = (/data_temp(imod_temp)/)
            data_new&models(dimsizes(data_all&models)) = \
              (/data_temp&models(imod_temp)/)
            delete(data_all)
            data_all = data_new
            delete(data_new)
          else
            ; Loop over models of data
            do imod = 0, dimsizes(data_all&models)-1
              ; if neq data model is similar to current
              ; entry, write data entry
              if (data_all&models(imod).eq. data_temp&models(imod_temp)) then
                data_all(imod, ii) = (/data_temp(imod_temp)/)
              end if
            end do
          end if
        end do
      end if
      ; Append input file list for provenance
      prov_files := array_append_record(prov_files, \
                                        str_split(data_temp@input, ","), 0)
    end if
    delete(data_temp)
  end do
  data_all!1 = "vars"
  data_all&vars = var_collect
  delete(var_collect)

  ; Get project for models
  projects = new(dimsizes(data_all&models), string)

  ; Loop over model-project files to complete project list
  do ii = 0, dimsizes(mp_file_list) - 1
    modproj = asciiread(mp_file_list(ii), -1, "string")
    mods = modproj(:dimsizes(modproj)/2-1)
    projs = modproj(dimsizes(modproj)/2:)

    overlap_index = get1Dindex(data_all&models, mods)
    projects(overlap_index) = projs
    delete([/modproj, mods, projs, overlap_index/])
  end do

  data_all&models@project = projects
  delete(projects)

  ; Sort diagnostics in the order specified in the settings
    if (isatt(diag_script_info, "diag_order")) then
      l_ok = True
      if (dimsizes(data_all&vars).ne. \
          dimsizes(diag_script_info@diag_order)) then
        error_msg("w", DIAG_SCRIPT, "", "specified order of diagnostics " + \
                  "cannot be applied, number of diagnostics does not match")
        l_ok = False
      end if
      pid = new(dimsizes(diag_script_info@diag_order), integer)
      do ii = 0, dimsizes(diag_script_info@diag_order) - 1
        tmp = ind(var_diag.eq.diag_script_info@diag_order(ii))
        if (any(ismissing(tmp)) .or. dimsizes(tmp).gt.1) then
          error_msg("w", DIAG_SCRIPT, "", "specified order of diagnostics " + \
                    "cannot be applied, invalid entry in diag_order")
          break
        end if
        pid(ii) = tmp
        delete(tmp)
      end do
      if (l_ok) then
        data_all := data_all(:, pid)
        alt_obs  := alt_obs(pid)
      end if
    end if

; -------------------------------------------------------------------------
; ----------------- Interim Functions -------------------------------------
; -------------------------------------------------------------------------

  undef("get_unique_entries")
  function get_unique_entries(array)
  ;
  ; Arguments:
  ;      array: 1D array
  ;
  ; Return value: 1D array of unique entries in array
  ;
  ; Modification history:
  ;    20170406-A_gier_bettina: written.
  local dummy_array, unique_new, new_array, nodupes
  begin
  dummy_array = array
  do while (dimsizes(dummy_array).ne.0)
    if (.not.isdefined("unique")) then
      unique = dummy_array(0)
    else
      unique_new = array_append_record(unique, dummy_array(0), 0)
      delete(unique)
      unique = unique_new
      delete(unique_new)
    end if
    nodupes = ind(dummy_array.ne.dummy_array(0))
    ; Missing value index are dim 1 and would give an error
    if (dimsizes(dummy_array).eq. \
        dimsizes(ind(dummy_array.eq.dummy_array(0)))) then
      break
    end if
    new_array = dummy_array(nodupes)
    delete(nodupes)
    delete(dummy_array)
    dummy_array = new_array
    delete(new_array)
  end do
  return(unique)
  end
; -----------------------------------------------------------------------------
; ---------------------------- Plotting ---------------------------------------
; -----------------------------------------------------------------------------

  ; Calculating necessary values
  ; Number of Projects needed to determine span
  ; For now just CMIP projects
  c_projects = str_match_ic(data_all&models@project, "CMIP")
  projects = get_unique_entries(c_projects)
  n_var = dimsizes(data_all&vars)

  nr_projects = dimsizes(projects)
  x_val = ispan(1, n_var*nr_projects, nr_projects)

  ; Mean and Median of Ensemble - without alt obs
  obs_ind = get1Dindex(data_all&models, alt_obs)
  if all(alt_obs.eq."none") then
    mod_ind = ispan(0, dimsizes(data_all&models)-1, 1)
  else
    ex_ind = obs_ind(ind(.not.ismissing(obs_ind)))
    mods_ind = ispan(0, dimsizes(data_all&models)-1, 1)
    ex_ind@_FillValue = default_fillvalue("integer")
    mods_ind@_FillValue = default_fillvalue("integer")
    mod_ind = get1Dindex_Collapse(mods_ind, ex_ind)
    delete(mods_ind)
    delete(ex_ind)
  end if
  ; Split by project
  means = new((/nr_projects, n_var/), float)
  median = new((/nr_projects, n_var/), float)

  do iproj = 0, nr_projects - 1
    mod_proj_ind = ind(data_all&models@project(mod_ind).eq.projects(iproj))
    means(iproj, :) = dim_avg_n(data_all(mod_proj_ind, :), 0)
    median(iproj, :) = dim_median_n(data_all(mod_proj_ind, :), 0)
    delete(mod_proj_ind)
  end do

  ; Create outfile directory
  system("mkdir -p " + config_user_info@plot_dir)

  ; Plotting preparation
  name = ""
  outfile = config_user_info@plot_dir + name + "patterncor." + file_type
  wks = gsn_open_wks(file_type, outfile)
  wks@fullname = outfile

  ; Calc limits
  y_min = min(data_all)
  y_min := decimalPlaces(y_min-0.05, 1, True)
  x_max = max(x_val) + nr_projects

  ; Set half line length
  l_length = 0.3

  ; Project Colors - TODO: let them be specified in cfg
  fcolors = (/"black", "blue", "red"/)

  res = True
  res@gsnDraw          = False
  res@gsnFrame         = False
  res@vpWidthF = 0.8
  n_square = 16.
  if x_max.le. n_square then
    res@vpHeightF = 0.8
  else
    res@vpHeightF = 0.8*(n_square/x_max)
  end if
  font_height = 0.02/0.6 * res@vpHeightF
  res@gsnMaximize      = True
  res@tiYAxisString = "Correlation"
  res@trYMinF = y_min
  res@trYMaxF = 1
  yspan = res@trYMaxF - res@trYMinF
  res@trXMinF = 0
  res@trXMaxF = x_max
  res@tmXBLabels = data_all&vars
  res@tmXBValues = x_val
  res@tmXBMode = "Explicit"
  res@tmXBLabelAngleF = 90.
  res@tmXBLabelJust = "CenterRight"
  res@tmXBLabelFontHeightF = font_height
  ; Set Marker Size to be half of line_length
  marker_size = res@vpWidthF * l_length / (res@trXMaxF - res@trXMinF)

  ; Resources for model lines
  res_lines                   = True                  ; polyline mods desired
  res_lines@gsLineDashPattern = 0.                    ; solid line
  res_lines@gsLineThicknessF  = 2.5                   ; line thickness
  res_lines@tfPolyDrawOrder = "PreDraw"

  ; Resources for mean lines
  res_mlines                   = True                  ; polyline mods desired
  res_mlines@gsLineDashPattern = 0.                    ; solid line
  res_mlines@gsLineThicknessF  = 4.                    ; line thicker
  res_mlines@tfPolyDrawOrder = "PreDraw"

  ; Resources for obs data markers
  res_circ = True
  res_circ@gsMarkerIndex = 16
  res_circ@gsMarkerColor = "green"
  res_circ@gsMarkerSizeF = marker_size
  res_circ@gsMarkerOpacityF = 0.4

  ; Resources for white markers below median
  res_circw = True
  res_circw@gsMarkerIndex = 16
  res_circw@gsMarkerColor = "white"
  res_circw@gsMarkerSizeF = 0.95*marker_size
  res_circw@tfPolyDrawOrder = "PreDraw"

  ; Resources for median markers if required
  res_circm = True
  res_circm@gsMarkerIndex = 4
  res_circm@gsMarkerSizeF = marker_size
  res_circm@gsMarkerThicknessF = 3.
  res_circm@tfPolyDrawOrder = "Draw"

  ; Resources for legend text
  res_text                    = True                  ; text mods desired
  res_text@txFontHeightF      = font_height           ; change text size
  res_text@txJust             = "CenterLeft"          ; text justification

  ; New x_val according to median!
  x_val_proj = new((/nr_projects, n_var/), float)
  ; space between projects in graph
  if nr_projects.eq.1 then
    d_proj = 1.5  ; offset
  else
    d_proj = 1    ; (nr_projects - 1.)/nr_projects
  end if
  do iproj = 0, nr_projects - 1
    do ivar = 0, n_var - 1
      x_val_proj(iproj, ivar) = ivar*nr_projects - 0.5 + d_proj*(iproj+1)
    end do
  end do
  ; Start with blank plot! gs and xy marker sizes are different..
  plot = gsn_csm_blank_plot(wks, res)

  do iproj = 0, nr_projects - 1
    res_circm@gsMarkerColor = fcolors(iproj)
    plot@$unique_string("dum_median")$ = gsn_add_polymarker( \
      wks, plot, x_val_proj(iproj, :), median(iproj, :), res_circm)
  end do

  ; add lines for individual models
  do ivar = 0, dimsizes(data_all(0, :))-1
    do iproj = 0, dimsizes(projects)-1
      ; Skip Project if no data for it
      proj_mods = ind(data_all&models@project(mod_ind).eq.projects(iproj))
      if .not. all(ismissing(data_all(proj_mods, ivar))) then
        proj_center = x_val_proj(iproj, ivar)
        xx = (/proj_center-l_length, proj_center+l_length/)
        ; Plot lines for mean
        xx_mean = (/proj_center-l_length*1.5, proj_center+l_length*1.5/)
        yy_mean = (/means(iproj, ivar), means(iproj, ivar)/)
        res_mlines@gsLineColor = fcolors(iproj)
        res_lines@gsLineColor = fcolors(iproj)
        plot@$unique_string("dum")$ = gsn_add_polyline( \
          wks, plot, xx_mean, yy_mean, res_mlines)
        do imod = 0, dimsizes(data_all(:, 0)) - 1
          ; Only plot if model in right project
          if data_all&models@project(imod).eq.projects(iproj) then
            ; Don't plot obs as lines
            if (.not.ismissing(data_all(imod, ivar))) then
              if (data_all&models(imod).ne.alt_obs(ivar)) then
                yy = (/data_all(imod, ivar), data_all(imod, ivar)/)
                plot@$unique_string("dum")$ = gsn_add_polyline( \
                  wks, plot, xx, yy, res_lines)
                end if
            end if
          end if
        end do
        plot@$unique_string("dum_ci")$ = gsn_add_polymarker( \
          wks, plot, x_val_proj(iproj, ivar), median(iproj, ivar), res_circw)
      end if
      delete(proj_mods)
    end do
    if (alt_obs(ivar).ne."none") then
      ; Plot obs as circles
      plot@$unique_string("dum_circ")$ = gsn_add_polymarker( \
        wks, plot, x_val(ivar), data_all(obs_ind(ivar), ivar), res_circ)
    end if
  end do

  y_min_label = res@trYMinF + 0.1*yspan
  lb_stride = yspan/res@vpHeightF * font_height * 1.2   ; font_height*3
  plabel = projects
  ; Draw Legend
  do iproj = 0, dimsizes(projects)-1
    res_text@txFontColor = fcolors(iproj)
    ; CMIP5 label has to be reduced to CMIP5 sometimes
    if str_match_ind_ic(plabel, "CMIP5").eq.iproj then
      plabel(iproj) = "CMIP5"
    end if
    plot@$unique_string("dum_l")$ = gsn_add_text(wks, plot, plabel(iproj),\
                                                 x_val(0), y_min_label + \
                                                 lb_stride*(iproj+1), res_text)
  end do
  res_text@txFontColor = "green"
  plot@$unique_string("dum_l")$ = gsn_add_text(wks, plot, "OBS", \
                                               x_val(0), y_min_label, res_text)

  draw(plot)
  frame(wks)
  ; Write output
  if diag_script_info@write_netcdf then
    system("mkdir -p " + config_user_info@work_dir)
    workpath = config_user_info@work_dir + "pattern_cor.nc"
    ncdf_outfile = ncdf_write(data_all, workpath)
  else
    ncdf_outfile = workpath
  end if

  ; collect meta-data and call ESMValMD function
  caption = "Centered pattern correlations between models and observations" \
            + " for the annual mean climatologies " \
            + "(similar to IPCC ch. 9 fig. 9.6)."
  statistics = (/"corr", "clim"/)
  domains = (/"global"/)
  plottype = "other"
  authors = (/"gier_bettina", "bock_lisa"/)
  references = (/"flato13ipcc"/)
  log_provenance(ncdf_outfile, outfile, caption, statistics, domains, \
                 plottype, authors, references, prov_files)
  leave_msg(DIAG_SCRIPT, "")
end
