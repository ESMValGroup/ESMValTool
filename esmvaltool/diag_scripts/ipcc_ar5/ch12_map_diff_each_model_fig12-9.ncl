; #############################################################################
; Mean variable change in one RCP for individual models
; Author: Ruth Lorenz (ETH, Switzerland)
; CRESCENDO project
; #############################################################################
; Caveats
;    All models are on different grids, so not possible to fill in all data
;    into one array, therefore for each model read hist and rcp after each
;    other, draw all panels and panel together at the end. Only include models
;    in namelist which exist for both, historical and rcp.
;
; Modification history
;    20181806-A_lore_ru: ported to version2
;    20171002-A_lore_ru: Added tags for searching
;    20160428-A_lore_ru: adapted to ESMValTool
;    20130501-A_sedl_ja: written for IPCC AR5.
;
; #############################################################################

load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/plot/contour_maps.ncl"
load "./diag_scripts/shared/ensemble.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.yml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(diag_script, "")

  ; 'diag_script' and 'variables' are fetched from the
  ; above 'interface.ncl' file
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; 'datasets', 'variables' and 'field_types' are fetched from the
  ; above 'interface.ncl' file
  dim_MOD = dimsizes(dataset_info@dataset)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ; ###########################################
  ; # References Information                  #
  ; ###########################################
  write_references(diag_script, \  ; script name
                   "A_lore_ru", \  ; authors
                   "", \  ; contributors
                   "D_collins13ipcc",    \  ; diag_references
                   "",          \  ; obs_references
                   "P_crescendo")   ; proj_references

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"experiment", "time_avg"/))

end

begin
  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  ; Create output plot directory
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
  ; # Get data and average time               #
  ; ###########################################
  ; Set dataset indexes
  modidx = ispan(0, dim_MOD / 2 - 1, 1)
  ; get data from first dataset
  imod = modidx(0)  ; NCL array indicies start from zero
  log_debug("processing " + dataset_info@dataset(imod) + "_" \
            + dataset_info@exp(imod) + "_" \
            + dataset_info@ensemble(imod))

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

  ; Average over time (function in ./diag_scripts/lib/ncl/statistics.ncl)
  data1 = time_operations(A0, -1, -1, "average", diag_script_info@time_avg, \
                          True)

  climofiles = new(dim_MOD, string)

  i = 0
  do ii = 0, dimsizes(modidx) - 1
    imod = modidx(ii)
    if (imod .eq. 0) then
      data2 = data1
    else
      A0 = read_data(imod, var0, field_type0)
      data2 = time_operations(A0, -1, -1, "average", \
                              diag_script_info@time_avg, True)
    end if
    climofiles = input_file_info@filename(imod)

    idx_rcp = find_model_index(dataset_info@dataset(imod), \
                               dataset_info@dataset, \
                               diag_script_info@experiment, \
                               dataset_info@exp, dataset_info@ensemble(imod), \
                               dataset_info@ensemble, 0)
    A1 = read_data(idx_rcp, var0, field_type0)
    data3 = time_operations(A1, -1, -1, "average", \
                            diag_script_info@time_avg, True)

    climofiles = input_file_info@filename(idx_rcp)

    ; ##############################################
    ; # Calculate change from period 1 to period 2 #
    ; ##############################################
    to_plot = data3 - data2
    copy_VarMeta(data2, to_plot)

    ; ###########################################
    ; # Other Metadata: diag_script, var        #
    ; ###########################################
    ; Add to to_plot, as attributes without prefix
    if (isatt(to_plot, "diag_script")) then  ; Add to existing entries
      temp = to_plot@diag_script
      delete(to_plot@diag_script)
      to_plot@diag_script = array_append_record(temp, (/diag_script/), 0)
      delete(temp)
    else  ; Add as new attribute
      to_plot@diag_script = (/diag_script/)
    end if
    to_plot@var = var0  ; Overwrite existing entry
    if (isatt(variable_info, "long_name")) then
      to_plot@var_long_name = variable_info@long_name
    end if
    if(isatt(variable_info, "units")) then
      to_plot@var_units     = variable_info@units
    end if

    ; ###########################################
    ; # Optional output to netCDF               #
    ; ###########################################
    to_plot@ncdf = work_dir + "/" + dataset_info@dataset(imod) + "_" + \
      dataset_info@ensemble(imod) + "_" + variable_info@diagnostic + ".nc"
    to_plot@ncdf_existing = "overwrite"
    ; Function in ~/interface_scripts/auxiliary.ncl
    ncdf_outfile = ncdf_write(to_plot, to_plot@ncdf)

    ; ######################################################
    ; # Separate output to a graphics file, Panelling      #
    ; ######################################################
    to_plot@res_gsnDraw           = False  ; Do not draw yet
    to_plot@res_gsnFrame          = False  ; Don't advance frame.
    to_plot@res_cnFillOn          = True            ; Color plot desired
    to_plot@res_cnLineLabelsOn    = False           ; Contour lines
    to_plot@res_cnLinesOn         = False
    if (isatt(diag_script_info, "colormap")) then   ; Set colormap
      col = read_colormap_file(diag_script_info@colormap)
    else
      col = read_colormap_file("./diag_scripts/shared/plot/rgb/" +  \
                               "ipcc_temperature_with_grey.rgb")
    end if
    to_plot@res_cnFillColors     = col
    if (isatt(diag_script_info, "span")) then
      to_plot@res_cnSpanFillPalette     = diag_script_info@span
    else
      to_plot@res_cnSpanFillPalette     = False   ; use full colormap
    end if
    to_plot@res_lbLabelBarOn = False
    ; function in aux_plotting.ncl
    if (imod.eq.0) then
      wks = get_wks("dummy_for_wks", diag_script, variable_info@diagnostic)
      nmap = dim_MOD / 2

      if (isatt(diag_script_info, "max_vert")) then
        ; Maximum allowed number of plots per page (vertical)
        nvert = toint(diag_script_info@max_vert)
      else
        nvert = 10
      end if
      if (isatt(diag_script_info, "max_hori")) then
      ; Maximum allowed number of plots per page (horizontal)
        nhori = toint(diag_script_info@max_hori)
      else
        nhori = 5
      end if
      maps = new(nmap, graphic)  ; collect individual maps in a graphic array
    end if
    to_plot@res_tiMainString = dataset_info@dataset(imod)
    to_plot@res_tiMainFontHeightF = 0.035
    if dataset_info@dataset(imod) .eq. "GISS-E2-R" .and. \
      dataset_info@ensemble(imod) .eq. "r1i1p1" then
      to_plot@res_tiMainString = "GISS-E2-R-P1"
    end if
    if dataset_info@dataset(imod) .eq. "GISS-E2-R" .and. \
      dataset_info@ensemble(imod) .eq. "r1i1p2" then
      to_plot@res_tiMainString = "GISS-E2-R-P2"
    end if
    if dataset_info@dataset(imod) .eq. "GISS-E2-R" .and. \
      dataset_info@ensemble(imod) .eq. "r1i1p3" then
      to_plot@res_tiMainString = "GISS-E2-R-P3"
    end if
    if dataset_info@dataset(imod) .eq. "GISS-E2-H" .and. \
      dataset_info@ensemble(imod) .eq. "r1i1p1" then
      to_plot@res_tiMainString = "GISS-E2-H-P1"
    end if
    if dataset_info@dataset(imod) .eq. "GISS-E2-H" .and. \
      dataset_info@ensemble(imod) .eq. "r1i1p2" then
      to_plot@res_tiMainString = "GISS-E2-H-P2"
    end if
    if dataset_info@dataset(imod) .eq. "GISS-E2-H" .and. \
      dataset_info@ensemble(imod) .eq. "r1i1p3" then
      to_plot@res_tiMainString = "GISS-E2-H-P3"
    end if
    to_plot@res_gsnRightString = ""

    if(isatt(diag_script_info, "projection")) then
      to_plot@res_mpProjection      = diag_script_info@projection
    else
      to_plot@res_mpProjection      = "Robinson"
    end if
    ; Set explicit contour levels
    if(isatt(diag_script_info, "diff_levs")) then
      to_plot@res_cnLevelSelectionMode = "ExplicitLevels"
      to_plot@res_cnLevels = diag_script_info@diff_levs
    end if

    to_plot@res_mpPerimOn         = False
    maps(i) = contour_map(wks, to_plot(:, :), var0)
    i = i + 1
    ; clean up
    delete(A0)
    delete(A1)
    delete(data2)
    delete(data3)
    delete(to_plot)
  end do

  pres                  = True   ; needed to override panelling defaults
  pres@gsnPanelLabelBar = True  ; no general label bar desired here

  a4_height = 29.7  ; in centimeters
  a4_width = 23.0   ; reference is correct
  cm_per_inch = 2.54

  pres@gsnPaperWidth = a4_width / cm_per_inch
  pres@gsnPaperHeight = a4_height / cm_per_inch
  pres@gsnPaperOrientation = "portrait"
  pres@gsnPanelTop = 0.96

  pres@lbTitleString = "(~F35~J~F~C)"
  pres@lbTitleFontHeightF = 0.017
  pres@lbTitleDirection = "across"
  pres@lbTitlePosition = "Right"
  pres@lbTitleJust = "CenterLeft"
  pres@lbLabelFontHeightF = 0.014
  pres@lbLabelJust = "CenterCenter"
  pres@lbLabelAutoStride = False

  pres@txFontHeightF = 0.015
  if (isatt(diag_script_info, "title")) then
    pres@txString  = diag_script_info@title
  end if
  pres@pmLabelBarParallelPosF = 0.06

  outfile = panelling(wks, maps, nvert, nhori, pres)

; add tags
;  alltags = array_append_record(tags, (/"DM_global", "PT_geo", "ST_clim"/), 0)
;  caption = "Surface change in " + variable_info@long_name + " in " + \
;    dataset_info@start_year(idx_rcp) + "-" + \
;    dataset_info@end_year(idx_rcp) + " displayed as anomalies with " + \
;    "respect to " + dataset_info@start_year(0) + "-" + \
;    dataset_info@end_year(0) + " for " + dataset_info@exp(idx_rcp) + \
;    " from one ensemble member of each of the " + \
;    "concentration driven models available in the CMIP5 archive. " + \
;    "Similar to Figure 12.9 of Collins et al. (2013)."
;  id = diag_script + "_" + variables(0) + "_" + diag_script_info@time_avg
;  contrib_authors = (/"A_lore_ru", "A_sedl_ja"/)

;  ESMValMD(outfile, alltags, caption, id, variables(0), \
;           dataset_info@dataset, climofiles, diag_script, contrib_authors)

;  delete([/alltags, caption, id, climofiles/])

  log_info(" Wrote " + outfile)
  leave_msg(diag_script, "")
end
