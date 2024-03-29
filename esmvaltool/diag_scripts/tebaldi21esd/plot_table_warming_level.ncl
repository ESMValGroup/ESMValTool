; #############################################################################
; Plot the table of warming level crossings, as in Table 1 of Tebaldi et al.
; with mean and 5-95 percent confidence intervals
; Author: Debeire Kevin (DLR, Germany)
; #############################################################################
;
;  Description
;     Plot a table of warming level crossing years for specified scenarios
;     (columns) and warming levels (rows):
;
;  Required diag_script_info attributes (diagnostics specific)
;     scenarios: list of scenarios included in the figure
;     warming_levels: list of warming levels
;     syears: list of start years of historical and SSPs scenarios
;     eyears: list of end years of historical and SSPs scenarios
;     begin_ref_year: start year of reference period
;     end_ref_year: end year of reference period
;     label: list of labels of the scenarios
;     offset: offset between reference baseline and 1850-1900
;
;  Modification history
;     20220803-A_debe_kevin: written for ESMValTool
;
; #############################################################################
; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"

load "$diag_scripts/shared/statistics.ncl"

load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
  enter_msg(DIAG_SCRIPT, "")

  ; 'datasets', 'variables' are fetched from the above 'interface.ncl' file
  var0 = variable_info[0]@short_name
  info_items = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info_items)
  dim_VAR = ListCount(variable_info)
  dataset_names = metadata_att_as_array(info_items, "dataset")
  dataset_exps = metadata_att_as_array(info_items, "exp")
  dataset_ensembles = metadata_att_as_array(info_items, "ensemble")
  n_MOD = count_unique_values(dataset_names)
  MOD_list = get_unique_values(dataset_names)
  labels = diag_script_info@label
  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", \
                       "eyears", "warming_levels"/))
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
end

begin
  ; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)
end
begin
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
  ncdf_dir = diag_script_info@input_files(1)
  load_file = ncdf_dir + "/" + variable_info[0]@diagnostic + ".nc"
  to_plot = ncdf_read(load_file, var0)
  mean_ = tostring(to_plot(0, :, :))
  five_perc_ = tostring(to_plot(1, :, :))
  ninetyfive_perc_ = tostring(to_plot(2, :, :))
  num_crossing_models = tostring(to_plot(3, :, :))
  num_models = tostring(to_plot(4, :, :))
  warming_levels = diag_script_info@warming_levels
  n_WAR = dimsizes(warming_levels)
  n_SCEN = dimsizes(diag_script_info@scenarios)
  ; ###########################################
  ; *********** PLOTTING ************
  ; Create plot variables
  outfile = "warming_level_crossings"
  wks = get_wks(0, DIAG_SCRIPT, outfile)
  res        = True
  res2 = True
  res3 = True
  res4 = True
  res@gsLineThicknessF = 2.
  res2@gsLineThicknessF = 2.
  res3@gsLineThicknessF = 2.
  ncr  = (/n_WAR, n_SCEN/)
  ncr2  = (/n_WAR+1, 1/)
  ncr3 = (/1, n_SCEN/)
  mean_@_FillValue = "NA"
  five_perc_@_FillValue = "NA"
  ninetyfive_perc_@_FillValue = "NA"
  text = new((/n_WAR+1, n_SCEN/), string)
  text2 = new((/n_WAR+1, 1/), string)
  text3 = labels
  do i = 0, dimsizes(warming_levels)-1
    text2(i+1, 0) = sprintf("%2.1f", warming_levels(i)) + \
                            "~S~o~N~C warming level"
    do j = 0, n_SCEN - 1
        text(i, j) = mean_(i, j) + "~C~["+ninetyfive_perc_(i, j) + \
                     ", " + five_perc_(i, j) + "]" + \
                     "~C~(" + num_crossing_models(i, j) + "/" + \
                     num_models(i, j) + ")"
    end do
  end do
  res2@txJust = "CenterCenter"
  res@txJust = "CenterCenter"
  res3@txJust = "CenterCenter"

  res@txFontHeightF  = 0.013
  res2@txFontHeightF = 0.02
  res3@txFontHeightF = 0.02

  res2@gsFillColor = "gray"
  res3@gsFillColor = "gray"

  x = (/0.30, 0.95/)
  y = (/0.25, 0.75-0.5/6/)
  x2 = (/0.05, 0.30/)
  y2 = (/0.25, 0.75/)
  x3 = (/0.30, 0.95/)
  y3 = (/0.75-0.5/6, 0.75/)
  ; Title
  x4 = (/0.05, 0.95/)
  y4 = (/0.77, 0.79/)
  ncr4 = (/1, 1/)
  res4@gsLineColor   = "transparent"
  res4@txFontHeightF = 0.02
  if (isatt(diag_script_info, "title")) then
    text4 = diag_script_info@title
  else
    text4 = "Table of warming level crossings"
  end if

  gsn_table(wks, ncr, x, y, text, res)
  gsn_table(wks, ncr2, x2, y2, text2, res2)
  gsn_table(wks, ncr3, x3, y3, text3, res3)
  gsn_table(wks, ncr4, x4, y4, text4, res4)
  draw(wks)
  frame(wks)
  log_info(" Wrote " + wks@fullname)

  ; Collect meta-data
  netcdf_dir = diag_script_info@input_files(1)
  datapath = netcdf_dir + "/" + variable_info[0]@diagnostic + ".nc"
  nc_file = datapath
  caption = text4
  statistics = ("mean")
  domains = ("global")
  plot_types = ("portrait")
  authors = (/"debeire_kevin"/)
  references = (/"tebaldi21esd"/)
  infiles = metadata_att_as_array(info_items, "filename")
  log_provenance(nc_file, wks@fullname, caption, statistics, domains, \
                 plot_types, authors, references, infiles)
  leave_msg(DIAG_SCRIPT, "")

end
