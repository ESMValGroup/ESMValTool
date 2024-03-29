; #############################################################################
; Plot timeseries of across-realization Std. dev. for large ensemble experiment
; Author: Kevin Debeire (DLR, Germany)
; #############################################################################
;
;  Description
;     Plot time series of inter-member standard deviation (ex: SSP3-7.0)
;
;  Required diag_script_info attributes (diagnostics specific)
;    ancestors: variable and diagnostics that calculated timeseries to plot
;    scenarios: list of scenarios included in the figure
;    syears: list of start years of historical and SSPs scenarios
;    eyears: list of end years of historical and SSPs scenarios
;    begin_ref_year: start year of reference period
;    end_ref_year: end year of reference period
;    label: list of labels of the scenarios
;    sampled_model: name of dataset on which to sample
;    n_samples: number of samples of size 10 to draw among all the ensembles
;               of sampled_model only
;
;  Optional diag_script_info attributes (diagnostic specific)
;     trend: whether the trend is calculated and displayed
;     runave_window: only used if trend is true, size window used for the
;                    centered running average
;     title: specify plot title
;     yaxis: specify y-axis title
;     ymin: minimim value on y-axis, default calculated from data
;     ymax: maximum value on y-axis
;     colormap: alternative colormap, path to rgb file or ncl name
;
;  Modification history
;     20220803-A_debe_kevin: written for ESMValTool
;
; #############################################################################
; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
begin
  enter_msg(DIAG_SCRIPT, "")

  ; 'models', 'variables' are fetched from the above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)
  dataset_names = metadata_att_as_array(info_items, "dataset")
  MOD_list = get_unique_values(dataset_names)
  n_MOD = dimsizes(MOD_list)
  ; Check required diag_script_info attributes
  req_atts = (/"syears", "eyears"/)
  exit_if_missing_atts(diag_script_info, req_atts)

end


begin
  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  netcdf_dir = diag_script_info@input_files(1)

  ; Output netcdf directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
end

begin
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
    runave_window = 1
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
  wks_type = 0
  wks_type@wkPaperSize = "A4"
  wks = get_wks(wks_type, DIAG_SCRIPT, outfile)
  if(isatt(diag_script_info, "colormap")) then
    colormap = RGBtoCmap(diag_script_info@colormap)
  else
    colormap = "StepSeq25"
  end if
  gsn_define_colormap(wks, colormap)

  colors_main = (/ 2, 3, 4, 5, 6, 7/)

  res                   = True
  res@gsnMaximize = True
  res@gsnPaperOrientation = "portrait"
  res@gsnDraw = False
  res@gsnFrame = False
  res@txFont   = 25
  res@vpHeightF = 0.5                    ; change aspect ratio of plot
  res@vpWidthF = 0.80
  res@vpXF     = 0.16                   ; start plot at x ndc coord

  caption = "Time series of annual ensembles across-realization" + \
            " standard deviation"
  if (isatt(diag_script_info, "title")) then
    res@tiMainString = diag_script_info@title
    caption = diag_script_info@title
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
  labels = new(n_MOD+1, string)
  ; create plot
  plot  = gsn_csm_xy(wks, to_plot&year, to_plot(0, :), res)
  labels(0) = MOD_list(0)
  delete(res@tiMainString)

  ; display samplings on plot
  if(isatt(diag_script_info, "sampled_model")) then
    labels(n_MOD) = diag_script_info@sampled_model+" random"
    do j = n_MOD, dims_data-1
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

  ; display trend if asked
  half_runave_window = round(todouble(runave_window)/2., 3)-1
  if (isatt(diag_script_info, "trend") .and. diag_script_info@trend) then
    polyres = True
    polyres@gsLineThicknessF = res@xyLineThicknesses
    dum = new(n_MOD, graphic)
    do n = 0, n_MOD-1
        rc = regline(to_plot&year, to_plot(n, :))
        ic = rc@yintercept
        y_start = rc*toint(syears(0)) + ic
        y_end = rc * toint(eyears(1)) + ic
        polyres@gsLineColor      = colors_main(n)
        dum(n) = gsn_add_polyline(wks, plot, (/toint(syears(0)), \
                                  toint(eyears(1))/), (/y_start, y_end/), \
                                  polyres)
    end do
  end if

  ; delete(res@xyLineColor)
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
  lgres@lgItemOrder = (/5, 4, 3, 2, 1, 0/)
  legend = gsn_create_legend(wks, n_MOD+1, labels, lgres)
  amres = True
  amres@amJust = "TopRight"
  amres@amParallelPosF = 0.46
  amres@amOrthogonalPosF = -0.52
  annoid = gsn_add_annotation(plot, legend, amres)  ; add legend to plot
  draw(wks)
  frame(wks)

  log_info(" Wrote " + wks@fullname)

  ; collect meta-data
  nc_file = datapath
  statistics = ("stddev")
  domains = ("global")
  plot_types = ("times")
  authors = (/"debeire_kevin"/)
  references = (/"tebaldi21esd"/)
  infiles = metadata_att_as_array(info_items, "filename")
  log_provenance(nc_file, wks@fullname, caption, statistics, domains, \
                 plot_types, authors, references, infiles)

  leave_msg(DIAG_SCRIPT, "")
end
