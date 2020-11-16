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
;    20181806-A_lorenz_ruth: ported to version2
;    20171002-A_lorenz_ruth: Added tags for searching
;    20160428-A_lorenz_ruth: adapted to ESMValTool
;    20130501-A_sedlacek_jan: written for IPCC AR5.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/ensemble.ncl"
load "$diag_scripts/shared/scaling.ncl"

load "$diag_scripts/shared/plot/contour_maps.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.yml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(DIAG_SCRIPT, "")

  ; 'datasets', 'variables' are fetched from the above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"experiment", "time_avg"/))
end

begin
  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)
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
    log_debug(" MyParam = " + MyParam)
  else
    error_msg("f", DIAG_SCRIPT, "", "missing variable attribute " + \
              "'long_name'")
  end if

  ; ###########################################
  ; # Get data and average time               #
  ; ###########################################
  ; Set dataset indexes
  modidx = ispan(0, dim_MOD / 2 - 1, 1)
  ; get data from first dataset
  imod = modidx(0)  ; NCL array indicies start from zero
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
  idx = ind(dims .eq. "lat")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lat dimension")
  end if
  idx = ind(dims .eq. "lon")
  if (ismissing(idx)) then
    error_msg("f", DIAG_SCRIPT, "", "no lon dimension")
  end if

  ; Average over time (function in ./diag_scripts/lib/ncl/statistics.ncl)
  data1 = \
    time_operations(A0, -1, -1, "average", diag_script_info@time_avg, True)

  climofiles = new(dim_MOD, string)

  jj = 0
  do ii = 0, dimsizes(modidx) - 1
    imod = modidx(ii)
    if (imod .eq. 0) then
      data2 = data1
    else
      A0 = read_data(info_items[imod])
      data2 = \
        time_operations(A0, -1, -1, "average", diag_script_info@time_avg, True)
    end if
    climofiles = info_items[imod]@filename

    atts = True
    atts@dataset = info_items[imod]@dataset
    atts@exp = diag_script_info@experiment
    atts@ensemble = info_items[imod]@ensemble
    item_rcp = select_metadata_by_atts(info_items, atts)
    A1 = read_data(item_rcp[0])
    delete(atts)
    delete(item_rcp)
    data3 = \
      time_operations(A1, -1, -1, "average", diag_script_info@time_avg, True)

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
      to_plot@diag_script = array_append_record(temp, (/DIAG_SCRIPT/), 0)
      delete(temp)
    else  ; Add as new attribute
      to_plot@diag_script = (/DIAG_SCRIPT/)
    end if
    to_plot@var = var0  ; Overwrite existing entry
    if (isatt(variable_info[0], "long_name")) then
      to_plot@var_long_name = variable_info[0]@long_name
    end if
    if(isatt(variable_info[0], "units")) then
      to_plot@units = variable_info[0]@units
      ; difference in K is the same as difference in degC, only change label
      if (to_plot@units.eq."K") then
        to_plot@units = "degC"
      end if
      if (isatt(diag_script_info, "plot_units")) then
        to_plot = convert_units(to_plot, diag_script_info@plot_units)
      end if
      plot_units = to_plot@units
    end if

    ; ###########################################
    ; # Output to netCDF               #
    ; ###########################################
    to_plot@ncdf = work_dir + "/" + info_items[imod]@dataset + "_" + \
      info_items[imod]@ensemble + "_" + variable_info[0]@diagnostic + ".nc"
    to_plot@ncdf_existing = "overwrite"
    ; Function in ~/interface_scripts/auxiliary.ncl
    ncdf_outfile = ncdf_write(to_plot, to_plot@ncdf)

    if (.not.isatt(diag_script_info, "styleset")) then
      diag_script_info@styleset = "DEFAULT"
    end if

    ; Select plotting attributes from the project style files
    ; See ./diag_scripts/shared/plot/style.ncl
    annots = project_style(input_file_info, diag_script_info, "annots")

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
      wks = get_wks("dummy_for_wks", DIAG_SCRIPT, variable_info[0]@diagnostic)
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
    annots_split = str_split(annots(imod), "_")
    to_plot@res_tiMainString = annots_split(0) + " " + annots_split(1)
    to_plot@res_gsnRightString = ""
    to_plot@res_tiMainFontHeightF = 0.035

    if(isatt(diag_script_info, "projection")) then
      to_plot@res_mpProjection = diag_script_info@projection
    else
      to_plot@res_mpProjection = "Robinson"
    end if
    ; Set explicit contour levels
    if(isatt(diag_script_info, "diff_levs")) then
      to_plot@res_cnLevelSelectionMode = "ExplicitLevels"
      to_plot@res_cnLevels = diag_script_info@diff_levs
    end if

    to_plot@res_mpPerimOn = False
    maps(jj) = contour_map(wks, to_plot(:, :), var0)
    jj = jj + 1

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

  if (isatt(diag_script_info, "units")) then
    pres@lbTitleString = diag_script_info@units
  else
    unit_string = format_units(plot_units)
    pres@lbTitleString = "(" + unit_string + ")"
  end if
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

  log_info(" Wrote " + outfile)

  ; collect meta-data
  nc_file = ncdf_outfile
  plot_file = outfile
  caption = "Change in " + var0 + " in " + info_items[dim_MOD-1]@start_year + \
    "-" + info_items[dim_MOD - 1]@end_year + \
    " displayed as anomalies with respect to " + info_items[0]@start_year + \
    "-" + info_items[0]@end_year + " for " + \
    diag_script_info@experiment + "."
  statistics = ("mean")
  domains = ("global")
  plot_types = ("geo")
  authors = (/"lorenz_ruth"/)
  references = (/"collins13ipcc"/)
  infiles = metadata_att_as_array(info_items, "filename")
  log_provenance(nc_file, plot_file, caption, statistics, domains, \
                 plot_types, authors, references, infiles)

  leave_msg(DIAG_SCRIPT, "")

end
