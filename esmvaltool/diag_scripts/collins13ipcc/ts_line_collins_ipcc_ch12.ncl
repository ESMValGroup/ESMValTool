;;#############################################################################
;; Global timeseries for historical and RCPs 1850-2?00
;; Author: Ruth Lorenz (ETH, Switzerland)
;; PROJECT-NAME CRESCENDO
;;#############################################################################
;;
;; Description
;;    Global or regional timeseries in historical and RCPs
;;    1850-2?00 for all available models
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    scenarios: which scenarios included in figure
;;    syears : start years in time periods
;;             (e.g. start of historical period and rcps)
;;    eyears : end years in time periods
;;             (end year of historical runs and rcps)
;;    begin_ref_year: start year of reference period (e.g. 1986)
;;    end_ref_year: end year of reference period (e.g 2005)
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    ncdf : enable to output to netCDF; either use "default"
;;           or give a full file name
;;    ts_minlat: minimum latitude if not global e.g(/0/)
;;    ts_maxlat: maximum latitude if not global e.g.(/90/)
;;    ts_minlon: minimum longitude if not global e.g.(/0/)
;;    ts_maxlon: maximum longitude if not global e.g.(/360/)
;;    model_nr: boolean, save number of model runs per period and scenario
;;              in netcdf to print in plot?
;;    colormap: choose file to create colormap from
;;    yaxis: y-axis title
;;    ymin: minimim value on y-axis
;;    ymax: maximum value on y-axis
;;    title: title in plot
;;    label: labels to use in legend depending on scenarios
;;
;; Required variable_info attributes (variable specific)
;;
;; Optional variable_info attributes (variable specific)
;;    long_name: long name of variable
;;    units: variable unit
;;
;; Caveats
;;    Needs lots of memory, in particular if lots of models and scenarios
;;    are included
;;    Variable number of scenarios? So far hardcoded 4 scenarios (rcp2.6,
;;    rcp4.5, rcp6.0, rcp8.5)
;;    every model (ensemble) in RCP needs corresponding historical run as
;;    reference
;;
;; Modification history
;;    20170523-A_lore_ru: separated plotting from calcuation
;;    20160625-A_lore_ru: adapted to ESMValTool
;;    20130506-A_sedl_ja: written for IPCC AR5.
;;
;;#############################################################################

load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/statistics.ncl"

begin
  ;;##############################################################
  ;;# Fetch general parameters, set in namelist_collins13ipcc.xml#
  ;;# passed via environment variables by python code            #
  ;;##############################################################
  enter_msg(diag_script, "")

  ;; 'diag_script' and 'variables' are fetched from the
  ;; above 'ncl.interface' file
  info_output("++++++++++++++++++++++++++++++++++++++++++")
  info_output(diag_script + " (var: " + variables(0) + ")")
  info_output("++++++++++++++++++++++++++++++++++++++++++")

  ;; 'models', 'variables' and 'field_types' are fetched from the
  ;; above 'ncl.interface' file
  dim_MOD = dimsizes(models@name)
  dim_VAR = dimsizes(variables)
  var0 = variables(0)
  field_type0 = field_types(0)

  ;;###########################################
  ;;# References Information                  #
  ;;###########################################
  write_references(diag_script, \  ; script name
                   "A_lore_ru", \  ; authors
                   "A_sedl_ja", \  ; contributors
                   "D_collins13ipcc",    \  ; diag_references
                   "",          \  ; obs_references
                   "P_crescendo")   ; proj_references

  ; Check required diag_script_info attributes
  exit_if_missing_atts(diag_script_info, (/"scenarios", "syears", "eyears", \
                                           "begin_ref_year", "end_ref_year"/))
end

; Load variable_info, delete?
begin
  vardeffile = "interface_data/" + var0 + "_info.tmp"
  loadscript(vardeffile)
end

begin
  ;; Get environment variables
  if (isatt(diag_script_info, "plot_type")) then
    file_type = diag_script_info@plot_type
  else
    file_type = "PDF"
  end if

  ;; Basename of diag_script
  diag_script_base = basename(diag_script) ;; Needed?

  ;; Output netcdf directory
  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir) ;; Needed?
  
  ;; Output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir) ;; Needed?
end

begin
  ;;#############################################
  ;;# Get parameters from ./variable_defs/*.ncl #
  ;;# passed via the 'info' attribute           #
  ;;#############################################
  if(isvar("MyParam")) then
    delete(MyParam)
  end if
  if(isatt(variable_info, "long_name")) then
    MyParam = variable_info@long_name
    info_output(" MyParam = " + MyParam, verbosity, 5)
  else
    error_msg("f", diag_script, "", "missing variable attribute " + \
              "'long_name'")
  end if
  ;;###########################################
  ;;# Get data and average annual and globally#
  ;;###########################################
  ; Set model indexes
  modidx = ispan(0, dim_MOD / 2 - 1, 1)
  ;; get data from first model for historical period
  imod = modidx(0)  ; NCL array indices start from zero
  log_debug("processing " + model_info@model(imod) + "_" \
            + model_info@exp(imod) + "_" + model_info@ensemble(imod))

  ;; See ./interface_scripts/data_handling.ncl
  A0 = read_data(imod, var0, field_type0)
  ;log_debug(" ++++++++++++++ Variable summary of original data (A0): " \
  ;          printVarSummary(A0))

  ;; Check dimensions
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

  ;; Average over time (function in ./diag_scripts/shared/statistics.ncl)
  data1 = time_operations(A0, -1, -1, "average", "yearly", True)

  ;log_debug(" ++++++++++++++ Variable summary of time averaged data " + \
  ;              "(data1): " + printVarSummary(data1))

  dim_file = dimsizes(data1)  ; file size of historical file
  ;; Average over globe/region
  ;; (function in ./diag_scripts/shared/latlon.ncl)
  if (isatt(diag_script_info, "ts_minlat")) then
    data1_glob = area_operations(data1, diag_script_info@ts_minlat, \
                                 diag_script_info@ts_maxlat, \
                                 diag_script_info@ts_minlon, \
                                 diag_script_info@ts_maxlon, "average", \
                                 True)
  else
    data1_glob = coswgt_areaave(data1)
  end if
  ;log_debug(" ++++++++++++++ Variable summary of area averaged data " + \
  ;          "(data1_glob): " + printVarSummary(data1_glob))

  ;; How many historical model runs?
  ind_hist = ind(model_info@exp .eq. "historical")
  dim_MOD_hist = dimsizes(model_info@model(ind_hist))
  ref_data2 = new((/dim_MOD_hist/), float, data1_glob@_FillValue)
  data2_rel = new((/dim_MOD_hist, dim_file(0)/), float, \
                  data1_glob@_FillValue)

  ;; How many scenarios?
  dim_scen = dimsizes(diag_script_info@scenarios)
  ;;start and end years of timeperiods
  syears = diag_script_info@syears
  eyears = diag_script_info@eyears
  nperiods = dimsizes(syears)
  nyears = toint(eyears(nperiods - 1)) - toint(syears(0)) + 1
  nyearsp1 =  toint(eyears(nperiods - 1)) - toint(eyears(0))

  ;; How many rcp2.6 model runs?
  ind_26 = ind(model_info@exp .eq. diag_script_info@scenarios(0))
  dim_MOD_26 = dimsizes(model_info@model(ind_26))
  ;; indicate array for scenarios
  data26_rel = new((/dim_MOD_26, nyearsp1/), float, \
                   data1_glob@_FillValue)
  i26 = 0
                   
  ;; How many rcp4.5 model runs?
  ind_45 = ind(model_info@exp .eq. diag_script_info@scenarios(1))
  dim_MOD_45 = dimsizes(model_info@model(ind_45))
  data45_rel = new((/dim_MOD_45, nyearsp1/), float, \
                   data1_glob@_FillValue)
  i45 = 0
                   
  ;; How many rcp6.0 model runs?
  ind_60 = ind(model_info@exp .eq. diag_script_info@scenarios(2))
  dim_MOD_60 = dimsizes(model_info@model(ind_60))
  data60_rel = new((/dim_MOD_60, nyearsp1/), float, \
                   data1_glob@_FillValue)
  i60 = 0
                   
  ;; How many rcp8.5 model runs?
  ind_85 = ind(model_info@exp .eq. diag_script_info@scenarios(3))
  dim_MOD_85 = dimsizes(model_info@model(ind_85))
  data85_rel = new((/dim_MOD_85, nyearsp1/), float, \
                   data1_glob@_FillValue)
  i85 = 0
  i = 0
  do imod = 0, dim_MOD_hist - 1
    if (imod .eq. 0) then
      data2 = data1
      delete(A0)
    else
      A0 = read_data(imod, var0, field_type0)
      data2 = time_operations(A0, -1, -1, "average", "yearly", True)
      delete(A0)
    end if
    ;; Average over globe (function in
    ;; ./diag_scripts/shared/statistics.ncl)
    data2_glob = coswgt_areaave(data2)
    ind_start = ind(data2&year .eq. diag_script_info@begin_ref_year)
    ind_end = ind(data2&year .eq. diag_script_info@end_ref_year)
    ref_data2(imod) = dim_avg(data2_glob(ind_start : ind_end))
    dim_data2g = dimsizes(data2_glob)
    if (dim_data2g .ne. dim_file(0)) then
      tdiff = dim_file(0) - dim_data2g
      data2_rel(imod, 0 : tdiff - 1) = data2_glob@_FillValue
      data2_rel(imod, tdiff :) =  data2_glob - ref_data2(imod)
    else
      data2_rel(imod, :) = data2_glob - ref_data2(imod)
    end if
    ;; clean up
    delete(data2)
    delete(data2_glob)

    ;; find all other runs from this model and calculate relative ts
    imod26 = ind(model_info@model .eq. model_info@model(imod) .and. \
		 model_info@exp .eq. diag_script_info@scenarios(0))
    if (.not.all(ismissing(imod26))) then
      if (dimsizes(imod26) .gt. 1) then
	do ens = 0, dimsizes(imod26) - 1
	  A1 = read_data(imod26(ens), var0, field_type0)
	end do
      else
        A1 = read_data(imod26, var0, field_type0)
      end if
      data3 = time_operations(A1, -1, -1, "average", "yearly", True)
      data3_glob = coswgt_areaave(data3)
      dim_data3g = dimsizes(data3_glob)
      if (dim_data3g .ne. nyearsp1) then
        data26_rel(i26, 0 : dim_data3g - 1) = data3_glob - ref_data2(imod)
	data26_rel(i26, dim_data3g :) =  data3_glob@_FillValue
      else
	data26_rel(i26, :) = data3_glob - ref_data2(imod)
      end if
      i26 = i26 + 1
      ;; clean up
      delete(A1)
      delete(data3)
      delete(data3_glob)
    end if
    delete(imod26)

    imod45 = ind(model_info@model .eq. model_info@model(imod) .and. \
		 model_info@exp .eq. diag_script_info@scenarios(1))
    if (.not.all(ismissing(imod45))) then
      if (dimsizes(imod45) .gt. 1) then
        do ens = 0, dimsizes(imod45) - 1
          A1 = read_data(imod45(ens), var0, field_type0)
        end do
      else
        A1 = read_data(imod45, var0, field_type0)
      end if
      data3 = time_operations(A1, -1, -1, "average", "yearly", True)
      data3_glob = coswgt_areaave(data3)
      dim_data3g = dimsizes(data3_glob)
      if (dim_data3g .ne. nyearsp1) then
        data45_rel(i45, 0:dim_data3g - 1) = data3_glob - ref_data2(imod)
        data45_rel(i45, dim_data3g:) =  data3_glob@_FillValue
      else
        data45_rel(i45, :) = data3_glob - ref_data2(imod)
      end if
      i45 = i45 + 1
      ;; clean up
      delete(A1)
      delete(data3)
      delete(data3_glob)
    end if
    delete(imod45)

    imod60 = ind(model_info@model .eq. model_info@model(imod) .and. \
		 model_info@exp .eq. diag_script_info@scenarios(2))
    if (.not.all(ismissing(imod60))) then
      if (dimsizes(imod60) .gt. 1) then
        do ens = 0, dimsizes(imod60) - 1
          A1 = read_data(imod60(ens), var0, field_type0)
        end do
      else
	A1 = extract_data(imod60, var0, field_type0)
      end if
      data3 = time_operations(A1, -1, -1, "average", "yearly", True)
      data3_glob = coswgt_areaave(data3)
      dim_data3g = dimsizes(data3_glob)
      if (dim_data3g .ne. nyearsp1) then
        data60_rel(i60, 0 : dim_data3g - 1) = data3_glob - ref_data2(imod)
        data60_rel(i60, dim_data3g :) =  data3_glob@_FillValue
      else
	data60_rel(i60, :) = data3_glob - ref_data2(imod)
      end if
      i60 = i60 + 1
      ;; clean up
      delete(A1)
      delete(data3)
      delete(data3_glob)
    end if
    delete(imod60)

    imod85 = ind(model_info@model .eq. model_info@model(imod) .and. \
		 model_info@exp .eq. diag_script_info@scenarios(3))
    if (.not.all(ismissing(imod85))) then
      if (dimsizes(imod85) .gt. 1) then
	    do ens = 0, dimsizes(imod85) - 1
		A1 = extract_data(imod85(ens), data, -1, 0, 0)
	    end do
      else
	    A1 = extract_data(imod85, data, -1, 0, 0)
      end if
      data3 = time_operations(A1, -1, -1, "average", "yearly", True)
      data3_glob = coswgt_areaave(data3)
      dim_data3g = dimsizes(data3_glob)
      if (dim_data3g .ne. nyearsp1) then
        data85_rel(i85, 0:dim_data3g - 1) = data3_glob - ref_data2(imod)
        data85_rel(i85, dim_data3g:) =  data3_glob@_FillValue
      else
        data85_rel(i85, :) = data3_glob - ref_data2(imod)
      end if
      i85 = i85 + 1
      ;; clean up
      delete(A1)
      delete(data3)
      delete(data3_glob)
    end if
    delete(imod85)
  end do  ; imod

  ;; average and standard deviation over all models
  datahist_avg = dim_avg_n_Wrap(data2_rel, 0)
  datahist_std = dim_stddev_n_Wrap(data2_rel, 0)

  data26_avg = dim_avg_n_Wrap(data26_rel, 0)
  data26_std = dim_stddev_n_Wrap(data26_rel, 0)

  data45_avg = dim_avg_n_Wrap(data45_rel, 0)
  data45_std = dim_stddev_n_Wrap(data45_rel, 0)

  data60_avg = dim_avg_n_Wrap(data60_rel, 0)
  data60_std = dim_stddev_n_Wrap(data60_rel, 0)

  data85_avg = dim_avg_n_Wrap(data85_rel, 0)
  data85_std = dim_stddev_n_Wrap(data85_rel, 0)

  ;; put all data from 1850-2100 into mean_val and std_val
  mean_val = new((/dim_scen + 1, nyears/), typeof(data2_rel), \
                 data2_rel@_FillValue)
  std_val = new((/dim_scen + 1, nyears/), typeof(data2_rel), \
                data2_rel@_FillValue)
  p0 = toint(eyears(0)) - toint(syears(0))
  p1 = toint(eyears(0)) - toint(syears(0)) + 1
  mean_val(0, 0 : p0) = datahist_avg
  std_val(0, 0 : p0) = datahist_std

  mean_val(1, p1 :) = data26_avg
  std_val(1, p1 :) = data26_std
  mean_val(2, p1 :) = data45_avg
  std_val(2, p1 :) = data45_std
  mean_val(3, p1 :) = data60_avg
  std_val(3, p1 :) = data60_std
  mean_val(4, p1 :) = data85_avg
  std_val(4, p1 :) = data85_std

  ;; check number of models for timeperiods in scenarios
  ;; if there are less than 5 models set to missing
  runs = new((/1 + dim_scen, nperiods/), integer)
  runs(0, 0) = dimsizes(ind(model_info@exp .eq. "historical" .and. \
                            model_info@end_year .eq. diag_script_info@eyears(0)))
  do scen = 0, dim_scen - 1
    do p = 0, nperiods - 1
      tmp = dimsizes(ind(model_info@exp .eq. \
                         diag_script_info@scenarios(scen) \
                         .and. toint(model_info@end_year) .gt. \
                         toint(diag_script_info@eyears(p))))
      runs(scen + 1, p) = dimsizes(ind(model_info@exp .eq. \
                                   diag_script_info@scenarios(scen) \
                                   .and. model_info@end_year .eq. \
                                   diag_script_info@eyears(p))) + tmp
      delete(tmp)
      if (runs(scen + 1, p) .le. 4) .and. (p .ne. 0) then
        p2 = toint(eyears(p - 1)) - toint(syears(0))
        p3 = toint(eyears(p)) - toint(syears(0))
        mean_val(scen + 1, p2 : p3) = mean_val@_FillValue
        std_val(scen + 1, p2 : p3) = std_val@_FillValue
      end if
    end do
  end do
  dims_data = dimsizes(mean_val)
  to_plot = new((/3 * dims_data(0), dims_data(1)/), float)

  do j = 0, dims_data(0) - 1
    ;;  if variable sw or lw TOA flux: change direction
    if (var0 .eq. "rlut" .or. var0 .eq. "rsut") then
      to_plot(j * 3, :) = - mean_val(j, :)
    else
      to_plot(j * 3, :) = mean_val(j, :)
    end if
    to_plot((j * 3) + 1, :) = to_plot(j * 3, :) + 1.64 * std_val(j, :)
    to_plot((j * 3) + 2, :) = to_plot(j * 3, :) - 1.64 * std_val(j, :)
  end do
  copy_VarMeta(mean_val, to_plot)
  to_plot!0 = "metric"
  to_plot&metric = (/"mean", "+1.64std", "-1.64std", \
                     "mean", "+1.64std", "-1.64std", \
                     "mean", "+1.64std", "-1.64std", \
                     "mean", "+1.64std", "-1.64std", \
                     "mean", "+1.64std", "-1.64std"/)
  to_plot!1 = "year"
  to_plot&year = ispan(toint(syears(0)), toint(eyears(nperiods - 1)), 1)

  ;;###########################################
  ;;# Other Metadata: diag_script, var        #
  ;;###########################################
  ;; Add to to_plot, as attributes without prefix
  if(isatt(to_plot, "diag_script")) then  ; Add to existing entries
    temp = to_plot@diag_script
    delete(to_plot@diag_script)
    to_plot@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    to_plot@diag_script = (/diag_script/)
  end if
  to_plot@var = var0  ; Overwrite existing entry
  if(isatt(variable_info, "long_name")) then
    to_plot@var_long_name = variable_info@long_name
  end if
  if(isatt(variable_info, "units")) then
    to_plot@var_units     = variable_info@units
  end if
  ;log_debug(" ++++++++++++++ Variable summary data to plot" + \
  ;          "(to_plot): " + printVarSummary(to_plot))

  runs!0 = "scenario"
  runs&scenario = diag_script_info@label
  runs!1 = "period"
  runs&period = diag_script_info@syears
  if(isatt(runs, "diag_script")) then  ; Add to existing entries
    temp = runs@diag_script
    delete(runs@diag_script)
    runs@diag_script = array_append_record(temp, (/diag_script/), 0)
    delete(temp)
  else  ; Add as new attribute
    runs@diag_script = (/diag_script/)
  end if
  runs@var = "nr"
  runs@var_long_name = "number of model runs per scenario and period"
  runs@var_units     = 1
  ;log_debug(" ++++++++++++++ Variable summary of runs" + \
  ;          "(runs): " + printVarSummary(runs))
  
  ;;###########################################
  ;;# Output to netCDF                        #
  ;;###########################################
  to_plot@ncdf = diag_script_info@ncdf
  ncdf_file = work_dir + "/" + diag_script_info@ncdf
  ncdf_file@existing = "overwrite"
  ;; Function in ~/interface_scripts/auxiliary.ncl
  ncdf_outfile = ncdf_write(to_plot, ncdf_file)
  if(isatt(diag_script_info, "model_nr") .and. \
    diag_script_info@model_nr .eq. True) then
    ncdf_file@existing = "append"
    ncdf_outfile = ncdf_write(runs, ncdf_file)
  end if

end
