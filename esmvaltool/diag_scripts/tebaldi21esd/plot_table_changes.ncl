; #############################################################################
;  Plots the table of changes between historical and scenario for different periods
;  Author: Debeire Kevin (DLR, Germany)
; #############################################################################

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
  MOD_list= get_unique_values(dataset_names)
  print(MOD_list)
  labels = diag_script_info@label
  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", "eyears"/))
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
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
end
begin
;Plot the table with means and 5-95 confidence intervals of changes
  ; Output netcdf directory
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
  ncdf_dir = diag_script_info@input_files(1)
  load_file = ncdf_dir + "/" + variable_info[0]@diagnostic + ".nc"
  to_plot = ncdf_read(load_file, var0)
  mean_ = sprintf("%3.2f",decimalPlaces(to_plot(0,:,:),2,True))
  five_perc_ = sprintf("%3.2f",decimalPlaces(to_plot(1,:,:),2,True))
  ninetyfive_perc_ = sprintf("%3.2f",decimalPlaces(to_plot(2,:,:),2,True))
  num_models = tostring(toint(to_plot(3,:,:)))
  period_list = to_plot&periods
  n_PER = dimsizes(period_list)
  n_SCEN= dimsizes(diag_script_info@scenarios)
  ; ###########################################
  ; *********** PLOTTING ************
  ; Create plot variables
  outfile = "table_changes_"+var0
  wks = get_wks(0, DIAG_SCRIPT, outfile)
  res        = True
  res2 = True
  res3 = True
  res4 = True
  res@gsLineThicknessF = 2.
  res2@gsLineThicknessF = 2.
  res3@gsLineThicknessF = 2.
  ncr  = (/n_PER,n_SCEN/)
  ncr2  = (/n_PER+1,1/)
  ncr3= (/1,n_SCEN/)
  mean_@_FillValue = "NA"
  five_perc_@_FillValue = "NA"
  ninetyfive_perc_@_FillValue = "NA"
  text= new((/n_SCEN+1,n_SCEN/),string)
  text2= new((/n_SCEN+1,1/),string)
  text3 = labels
  do i = 0 , n_PER-1
    text2(i+1,0) = period_list(i)
    do j = 0, n_SCEN -1
        text(i,j) = mean_(i,j)+"~C~["+five_perc_(i,j)+","+ninetyfive_perc_(i,j)+"]"+"~C~("+num_models(i,j)+")"

    end do
  end do
  res2@txJust = "CenterCenter"
  res@txJust = "CenterCenter"
  res3@txJust = "CenterCenter"

  res@txFontHeightF  = 0.013
  res2@txFontHeightF= 0.02
  res3@txFontHeightF= 0.02

  res2@gsFillColor = "gray"
  res3@gsFillColor = "gray"

  x = (/0.30,0.95/)
  y = (/0.45,0.65/)
  x2 = (/0.05,0.30/)
  y2 = (/0.45,0.75/)
  x3 = (/0.30,0.95/)
  y3 = (/0.65,0.75/)
  ;Title
  x4 = (/0.05,0.95/)
  y4 = (/0.77,0.79/)
  ncr4 = (/1,1/)
  res4@gsLineColor   = "transparent"
  res4@txFontHeightF = 0.02
  if (isatt(diag_script_info, "title")) then
    text4 = diag_script_info@title
  else
    text4 = "Table of change for "+var0
  end if

  gsn_table(wks,ncr,x,y,text,res)
  gsn_table(wks,ncr2,x2,y2,text2,res2)
  gsn_table(wks,ncr3,x3,y3,text3,res3)
  gsn_table(wks,ncr4,x4,y4,text4,res4)
  draw(wks)
  frame(wks)
  log_info(" Wrote " + plot_dir + "/" + outfile + "." + file_type)
end
