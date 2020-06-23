; #############################################################################
; carbon_constraint.ncl
; Author: Sabrina Wenzel (DLR, Germany)
; #############################################################################
;
; Description:
;    Calculates the long term sensitivity GAMMA_LT according to Wenzel et al.
;    2014 with observations from GCP.
;    This plotscript calls a file which contains the results from calculating
;    GAMMA_IAV (carbon_corr_2var.ncl)
;
; Required diag_script_info attributes (diagnostics specific)9
;    * con_latrange   array of latitudes for time series plots if 2D or 3D
;    * con_lonrange   array of pressure levels for time series plots if 2D/3D
;    * gIAV_file      file path for gamma_IAV values calculated with
;                     carbon_corr_2var.ncl
;
; Modification history
;    20200101-zechlau_sabrina: code rewritten for ESMValTool v2.0
;    20170822-bock_lisa: added tags for reporting
;    20151027-lauer_axel: added call to 'write_references'
;    201509-zechlau_sabrina: updated to new ESMValTool structure
;    201303-wenzel_sabrina: written
; #############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/carbon_ec/carbon_aux.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/ensemble.ncl"

load "$diag_scripts/shared/plot/style.ncl"
load "$diag_scripts/shared/plot/scatterplot.ncl"
load "$diag_scripts/shared/plot/xy_line.ncl"
load "$diag_scripts/shared/plot/legends.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"

begin

  enter_msg(DIAG_SCRIPT, "")

  tmp = metadata_att_as_array(input_file_info, "short_name")
  variables = get_unique_values(tmp)
  delete(tmp)

  ; get indices associated with specific variables
  tas_idx = 0
  nbp_idx = 1

  ; Load var_info
  var0 = variable_info[tas_idx]@short_name
  info0 = select_metadata_by_name(input_file_info, var0)
  dim_MOD = ListCount(info0)

  var1 = variable_info[nbp_idx]@short_name
  info = select_metadata_by_name(input_file_info, var1)

  info1 = info[dim_MOD:]    ; 1pctCO2
  info2 = info[:dim_MOD-1]  ; esmFixClim1

  log_info("+++++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + " + " + var1 + ")")
  log_info("+++++++++++++++++++++++++++++++++++++++++++++")

  datasetnames0 = metadata_att_as_array(info0, "dataset")
  datasetnames1 = metadata_att_as_array(info1, "dataset")
  datasetnames2 = metadata_att_as_array(info2, "dataset")

  ; ******************************************************************
  ; Create output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)
  system("mkdir -p " + config_user_info@work_dir)

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Check for required settings
  exit_if_missing_atts(diag_script_info, "styleset")

  ; Set default values for non-required diag_script_info attributes
  set_default_att(diag_script_info, "ts_anom", False)
  set_default_att(diag_script_info, "volcanoes", False)

  if(isatt(diag_script_info, "ec_anom")) then
    anom = diag_script_info@ec_anom
  end if

  if(isatt(diag_script_info, "con_units")) then
    UNITS = diag_script_info@con_units
  end if

  temp_dir = config_user_info@run_dir + diag_script_info@nc_infile
  temp_dir = temp_dir + "gIAV_" \
                      + diag_script_info@gIAV_start + "-" \
                      + diag_script_info@gIAV_end + ".nc"

  ; ==================================================================
  ;  setup data holders for plots
  xMax_year = max(metadata_att_as_array(info0, "end_year"))
  xMin_year = min(metadata_att_as_array(info0, "start_year"))
  nMax = xMax_year - xMin_year + 1

  ; carbon cycle climate feedback
  gamma_LT = new((/dim_MOD/), double)
  gamma_LT_std = new((/dim_MOD/), double)

  ; trend lines
  Yreg = new((/dim_MOD-1/), double)

  ; slope
  arc = new((/dim_MOD/), double)

  ; read uncoupled VAR0
  climo_dir = getenv("ESMValTool_climo_dir")

  ; ==================================================================
  ; loop for models
  do imod = 0, dim_MOD-1

    ; Read data
    A0 = read_data(info0[imod])
    A1 = read_data(info1[imod])
    A2 = read_data(info2[imod])

    idx_mod := imod

    ; average yearly
    A0 := time_operations(A0, -1, -1, "average", "yearly", True)
    A1 := time_operations(A1, -1, -1, "average", "yearly", True)
    A2 := time_operations(A2, -1, -1, "average", "yearly", True)

    ; Convert units if appropriate
    if (isatt(variable_info[tas_idx], "plot_units")) then
      A0 := convert_units(A0, variable_info[tas_idx]@plot_units)
    end if
    if (isatt(variable_info[nbp_idx+1], "plot_units")) then
      A1 := convert_units(A1, variable_info[nbp_idx+1]@plot_units)
    end if
    if (isatt(variable_info[nbp_idx], "plot_units")) then
      A2 := convert_units(A2, variable_info[nbp_idx]@plot_units)
    end if

    ; cummulate nbp over time series
    tmpcA1 = new((/dimsizes(A1)+1/), double)
    tmpcA2 = new((/dimsizes(A2)+1/), double)
    tmpcA1(0) = 0.0
    tmpcA2(0) = 0.0

    do y = 0, dimsizes(A1&year)-1
      tmpcA1(y+1) = tmpcA1(y)+A1(y)
      tmpcA2(y+1) = tmpcA2(y)+A2(y)
    end do
    delete(y)

    ; calculate integral line
    tidx := dimsizes(A0)-1
    Y1data := avg(tmpcA1(tidx-9:tidx))-avg(tmpcA1(:9))
    Y2data := avg(tmpcA2(tidx-9:tidx))-avg(tmpcA2(:9))
    Y3data := avg(A0(tidx-9:tidx))-avg(A0(:9))

    ; calculate gamma_LT
    gamma_LT(imod) = ((Y1data-Y2data)/(Y3data))
    gamma_LT_std(imod) = sqrt((stddev(A1(tidx-9:tidx)) / Y3data) ^ 2 \
                              + (1 / stddev(A0(tidx-9:tidx))) ^ 2 \
                              + (stddev(A2(tidx-9:tidx)) / Y3data) ^ 2)

    delete([/A0, A1, A2, tidx/])
    delete([/Y1data, Y2data, Y3data/])
    delete([/tmpcA1, tmpcA2/])
  end do  ; model loop

  ; read file with GAMMA_IAV values
  IAV = ncdf_read(temp_dir, "gIAV")

  if (dimsizes(IAV(0, :)).eq.dim_MOD+1) then
    Y6data  = IAV(:, 1:dim_MOD)         ; Gamma_IAV from the models
    obs_std = IAV(:, 0)                 ; Gamma_IAV from the OBS
  else
    Y6data = IAV
    obs_std = (/0, 0/)
  end if

  ; exclude models from regression if appropriate
  aX = tofloat(Y6data)
  aY = tofloat(gamma_LT)
  if(isatt(diag_script_info, "reg_models")) then
    ind_MOD = ind(datasetnames0.eq.diag_script_info@reg_models)
    aX(:, ind_MOD) = aX@_FillValue
    aY(ind_MOD) = aY@_FillValue
  end if

  ; Calculate confidence limits on fit
  minx = min(aX(0, :))-1.1*(max(aX(0, :))-min(aX(0, :)))
  maxx = max(aX(0, :))+1.1*(max(aX(0, :))-min(aX(0, :)))
  nfit = fspan(0, 50, 50)

  ; calculate regression line
  rc = regline_stats(aX(0, :), aY(:))
  arc  = rc@r
  prc  = rc@F_pval
  xfit = minx+(maxx-minx)*nfit/dimsizes(nfit)
  yfit = rc*xfit + rc@yintercept

  PDF = Cond_PDF(xfit, yfit, aX(0, :), aY(:), rc@MSE, obs_std)

  if (obs_std(0).ne.0) then
    Py    = PDF@Py
    gamma_mco2 = rc*obs_std(0)+rc@yintercept
    gamma_eco2 = 0.5*cdft_t(0.05/2., 5)*sqrt(rc@rstd ^ 2+obs_std(1) ^ 2)
  end if
  Py_pr = PDF@Py_pr
  y     = PDF@y

  ; this is for plotting the regression line with confidence intervals
  yfit_std = new((/2, dimsizes(yfit)/), typeof(yfit))
  yfit_std(0, :) = yfit+PDF@sigma
  yfit_std(1, :) = yfit-PDF@sigma

  ; ==================================================================
  ; PLOTTING

  ; Get environment variables
  plot_dir = getenv("ESMValTool_plot_dir")
  diag_script_base = DIAG_SCRIPT
  output_dir = config_user_info@plot_dir

  plot_file = "constr_"+var0+"-"+var1+"_"+xMin_year+"-"+xMax_year
  wks = gsn_open_wks(file_type, output_dir + plot_file)

  colors = project_style(info0, diag_script_info, "colors")
  thicks = project_style(info0, diag_script_info, "thicks")
  plot_num = integertochar(ispan(97, 122, 1))

  XStg = "Sensitivity of land+ocean CO~B~2~N~ flux," \
         + "~F8~g~F21~~B~IAV~N~ (GtC/yr/K)"
  YStg = "Carbon Cycle Climate Feedback, ~F8~g~F21~~B~LT~N~ (GtC/K)"

  data_arr = new((/2, dimsizes(gamma_LT)/), float)
  data_arr!0      = "datasets"
  data_arr!1      = "model"
  data_arr&datasets = (/XStg, YStg/)
  data_arr&model = datasetnames0
  data_arr(0, :) = tofloat(Y6data(0, :))
  data_arr(1, :) = tofloat(gamma_LT)
  data_arr@diag_script = DIAG_SCRIPT

  data_arr@res_tiMainString      = ""
  data_arr@res_tiMainFontHeightF = 0.025
  data_arr@res_gsnCenterString   = ""
  data_arr@res_trXMinF           = min(Y6data(0, :))-0.2*abs(min(Y6data(0, :)))
  data_arr@res_trXMaxF           = max(Y6data(0, :))+0.2*abs(max(Y6data(0, :)))
  data_arr@res_trYMinF           = min(gamma_LT)-0.2*abs(min(gamma_LT))
  data_arr@res_trYMaxF           = max(gamma_LT)+0.2*abs(max(gamma_LT))
  plot = scatterplot(wks, data_arr, var0, False, False, input_file_info)

  peres = True
  peres@tfPolyDrawOrder = "PreDraw"
  peres@xyLineColor     = "orange"
  peres@gsFillColor     = "orange"
  xin = array_append_record(xfit, xfit(::-1), 0)
  yin = array_append_record(yfit_std(0, :), yfit_std(1, ::-1), 0)
  add_prediction_error(wks, plot, peres, xin, yin)
  delete([/xin, yin/])

  if (obs_std(0).ne.0) then
      xin = (/obs_std(0)+obs_std(1), obs_std(0)+obs_std(1), \
             obs_std(0)-obs_std(1), obs_std(0)-obs_std(1)/)
      yin = (/500, -500, -500, 500/)
      add_prediction_error(wks, plot, False, xin, yin)
      delete([/xin, yin/])
  end if

  ; add regression line
  lres = True
  lres@gsLineColor       = "red"
  lres@gsLineThicknessF  = 1.5
  fit_line = gsn_add_polyline(wks, plot, xfit, yfit, lres)

  ; add colored markers
  mres = True
  mres@xyMarker      = 0
  do imod = 0, dim_MOD-1
    mres@xyMarkerColor = colors(imod)
    add_markers(wks, plot, mres, Y6data(0, imod), gamma_LT(imod))
  end do

  ; add xy error bars
  eres = True
  eres@gsLineColor = colors
  add_errorbar(wks, plot, eres, Y6data, (/gamma_LT, gamma_LT_std/))
  draw(plot)
  frame(wks)

  ; create line plot with PDFs
  diag_script_info@multi_model_mean = False
  diag_script_info@scatter_log = False
  diag_script_info@xy_line_legend = False
  pres = True
  pres@tiMainString       = plot_num(1)+") "
  pres@tiXAxisString      = XStg
  pres@tiYAxisString      = "Probability Density"
  pres@vpWidthF           = 0.4
  if (obs_std(0).eq.0) then
    pres@xyLineColors   = "black"
    pres@xyDashPatterns = 2
    pres@trYMinF        = 0
    xy_line(wks, Py_pr, y, 0., pres, False)
  else
    pres@xyLineColors   = (/"black", "red"/)
    pres@xyDashPatterns = (/2, 0/)
    pres@trYMinF        = 0
    xy_line(wks, (/Py_pr, Py/), (/y, y/), 0., pres, input_file_info)
  end if

  ; *************************************
  ; output to NetCDF
  out_path = config_user_info@work_dir
  system("mkdir -p " + out_path)
  out_path1 = out_path + "constraint_" + var1 + "_" \
                       + diag_script_info@gIAV_start + "-" \
                       + diag_script_info@gIAV_end + ".nc"
  data_arr@ncdf = out_path
  data_arr@var = var0
  ncdf_outfile  = ncdf_write(data_arr, out_path1)

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 output_dir + "." + file_type, \
                 XStg + " vs " + YStg, \
                 (/"anomaly", "corr", "stddev"/), \
                 (/"trop", "global"/),\
                 (/""/), \
                 (/"wenzel_sabrina"/), \
                 (/"wenzel14jgr"/), \
                 metadata_att_as_array(info0, "filename"))

end
; ###########################################################################
