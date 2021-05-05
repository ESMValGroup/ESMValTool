; #############################################################################
; Plot global (or regional) timeseries for historical and RCPs 1850-2?00
; Author: Ruth Lorenz (ETH, Switzerland)
; PROJECT-NAME CRESCENDO
; #############################################################################
;
; Description
;    Plot global or regional timeseries in historical and RCPs
;    1850-2?00 for all available models precalculated by
;    ts_line_collins_ipcc_ch12.ncl
;
; Caveats
;    Variable number of scenarios? So far hardcoded 4 scenarios (rcp2.6,
;    rcp4.5, rcp6.0, rcp8.5)
;
; Modification history
;    20171002-A_lorenz_ruth: Added tags for searching
;    20170523-A_lorenz_ruth: separated plotting from calculations
;    20160625-A_lorenz_ruth: adapted to ESMValTool in ts_line_collins_ipcc_ch12
;    20130506-A_sedlacek_jan: written for IPCC AR5.
;
; #############################################################################

; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
begin
  ; ##############################################################
  ; # Fetch general parameters, set in namelist_collins13ipcc.xml#
  ; # passed via environment variables by python code            #
  ; ##############################################################
  enter_msg(DIAG_SCRIPT, "")

  ; 'models', 'variables' are fetched from the above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)
  dataset_names = metadata_att_as_array(info_items, "dataset")
  MOD_list= get_unique_values(dataset_names)
  n_MOD = dimsizes(MOD_list)
  ; Check required diag_script_info attributes
  req_atts = (/"syears", "eyears"/)
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

  netcdf_dir = diag_script_info@input_files(1)

  ; Output netcdf directory
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
  if(isatt(variable_info[0], "long_name")) then
    MyParam = variable_info[0]@long_name
    log_debug(" MyParam = " + MyParam)
  else
    error_msg("f", DIAG_SCRIPT, "", "missing variable attribute " + \
              "'long_name'")
  end if
  if(isatt(diag_script_info, "runave_window")) then
    runave_window = diag_script_info@runave_window
  else
    runave_window= 1
  end if
  ; ###########################################
  ; # Get data                                #
  ; ###########################################
  datapath = netcdf_dir + "/" + variable_info[0]@diagnostic + ".nc"
  to_plot = ncdf_read(datapath, var0)
  ; start and end years of timeperiods
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
  nperiods = dimsizes(syears)

  dims_data_array = dimsizes(to_plot)
  dims_data = dims_data_array(0)
  dims_time = dims_data_array(1)
  ; ###########################################
  ; *********** PLOTTING ************
  ; Create plot variables
  outfile = var0 + "_ts_line_" + syears(0) + "_" + eyears(nperiods - 1)
  wks_type= file_type
  wks_type@wkPaperSize = "A4"
  wks = gsn_open_wks(wks_type, plot_dir + outfile)
  if(isatt(diag_script_info, "colormap")) then
    colormap = RGBtoCmap(diag_script_info@colormap)
  else
    colormap = "StepSeq25"
  end if
  gsn_define_colormap(wks, colormap)

  colors_main = (/ 2, 3, 4, 5, 6, 7/)

  res                   = True
  res@gsnMaximize= True
  res@gsnPaperOrientation= "portrait"
  res@gsnDraw = False
  res@gsnFrame = False
  res@txFont   = 25
  res@vpHeightF = 0.5                    ; change aspect ratio of plot
  res@vpWidthF = 0.80
  res@vpXF     = 0.16                   ; start plot at x ndc coord

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
  res@tmXBValues = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 20)
  res@tmXBLabels = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 20)
  res@tmXBMinorValues  = ispan(toint(syears(0)), \
                               toint(eyears(nperiods - 1)), 10)

  res@xyLineColor      = colors_main(0)
  res@xyLineThicknesses = 2.7
  res@xyLineOpacityF = 0.8
  ; create labels
  labels= new(n_MOD+1,string)
  ; create plot
  plot  = gsn_csm_xy(wks, to_plot&year, to_plot(0, :), res)
  print(to_plot(0, :))
  labels(0) = MOD_list(0)
  delete(res@tiMainString)

  ;display samplings on plot
  if(isatt(diag_script_info, "sampled_model")) then
    labels(n_MOD) = diag_script_info@sampled_model+" random"
    print(labels(n_MOD))
    print(dims_data-1)
    print(n_MOD)
    do j=n_MOD,dims_data-1
      res@xyLineThicknesses = 1.8
      res@xyLineColor      = colors_main(5)
      res@xyLineOpacityF = 0.2
      plot1  = gsn_csm_xy(wks, to_plot&year, to_plot(j, :), res)
      overlay(plot, plot1)
    end do
    res@xyLineOpacityF = 1.0
  end if

  do j = 0, n_MOD-1
    labels(j) = MOD_list(j)
    res@xyLineThicknesses = 2.7
    res@xyLineOpacityF = 0.8
    res@xyLineColor      = colors_main(j)
    plot1  = gsn_csm_xy(wks, to_plot&year, to_plot(j, :), res)
    overlay(plot, plot1)
  end do


  ;display trend if asked
  half_runave_window = round(todouble(runave_window)/2.,3)-1
  if (isatt(diag_script_info, "trend") .and. diag_script_info@trend) then
    polyres= True
    polyres@gsLineThicknessF=res@xyLineThicknesses
    dum = new(n_MOD, graphic)
    do n = 0, n_MOD-1
        rc = regline(to_plot&year,to_plot(n,:))
        ic = rc@yintercept
        y_start = rc*toint(syears(0))+ic
        ;print(y_start)
        y_end= rc*toint(eyears(1))+ic
        ;print(y_end)
        polyres@gsLineColor      = colors_main(n)
        dum(n) = gsn_add_polyline(wks, plot, (/toint(syears(0)), \
        toint(eyears(1))/), (/y_start, y_end/), polyres)
    end do
  end if

  ;delete(res@xyLineColor)
  ; Attach a legend
  lgres                    = True
  lgres@lgLineColors       = colors_main
  lgres@lgItemType         = "Lines"        ; show lines only (default)
  lgres@lgLabelFontHeightF = .08            ; legend label font thickness
  lgres@vpWidthF           = 0.2           ; width of legend (NDC)
  lgres@vpHeightF          = 0.2           ; height of legend (NDC)
  lgres@lgMonoDashIndex    = True
  lgres@lgDashIndex        = 0
  lgres@lgLineThicknessF        = res@xyLineThicknesses
  lgres@lgPerimOn = False
  lgres@lgItemOrder = (/ 5 , 4, 3, 2, 1, 0/)
  legend = gsn_create_legend(wks, n_MOD+1, labels, lgres)
  amres                  = True
  amres@amJust           = "TopRight"
  amres@amParallelPosF= 0.46
  amres@amOrthogonalPosF= -0.52
  annoid = gsn_add_annotation(plot, legend, amres)  ; add legend to plot
  draw(wks)
  frame(wks)

  log_info(" Wrote " + plot_dir + "/" + outfile + "." + file_type)

  ; collect meta-data
  nc_file = datapath
  plot_file = outfile

  caption = "Time series of annual 10 ensembles across-realization standard deviation for TAS global and CMIP6 models."
  statistics = ("stddev")
  domains = ("global")
  plot_types = ("times")
  authors = (/"debeire_kevin"/)
  references = (/"tebaldi et al. 2020"/)
  infiles = metadata_att_as_array(info_items, "filename")
  log_provenance(nc_file, plot_file, caption, statistics, domains, \
                 plot_types, authors, references, infiles)

  leave_msg(DIAG_SCRIPT, "")
end
