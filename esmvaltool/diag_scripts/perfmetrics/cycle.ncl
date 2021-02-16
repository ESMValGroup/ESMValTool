; #############################################################################
; PROCEDURE FOR THE ANNUAL CYCLE PLOT OF THE PERFORMANCE METRICS
; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
; ESMVal project
; #############################################################################

load "$diag_scripts/shared/plot/xy_line.ncl"

procedure perfmetrics_ptype_script()
begin

  ; Define output array
  if (diag_script_info@time_avg.eq."seasonalclim") then
    var_all = new((/nDatasets, 4, 2/), float)
    var_all!1 = "season"
    var_all&season = (/"DJF", "MAM", "JJA", "SON"/)
  else if (diag_script_info@time_avg.eq."monthlyclim") then
    var_all = new((/nDatasets, 12, 2/), float)
    var_all!1 = "month"
    var_all&month = (/"J", "F", "M", "A", "M", "J",\
                      "J", "A", "S", "O", "N", "D"/)
  else
      error_msg("f", DIAG_SCRIPT, "", "time_avg option " + \
                diag_script_info@time_avg + \
                " not compatible with plot_type cycle")
  end if
  end if
  var_all!0 = "model"
  var_all!2 = "statistic"
  var_all&model = datasetnames
  var_all&statistic = (/"mean", "stddev"/)

  ; Attach attributes
  var_all@var = var0
  var_all@diag_script = (/DIAG_SCRIPT/)
  copy_VarAtts(diag_script_info, var_all)
  var_all@ref_model = variable_info[0]@reference_dataset

  ; Search for level
  f = addfile(info_items[0]@filename, "r")
  if (isfilevar(f, "plev")) then
    if (dimsizes(f->plev).eq.1) then
      level = toint(f->plev/100.)
    end if
  end if

  ; Set path for saving processed data
  system("mkdir -p " + config_user_info@work_dir)
  if (isdefined("level")) then
    vv = var0 + level
  else
    vv = var0
  end if
  fname = str_join((/"perfmetrics", "cycle", vv, \
                     diag_script_info@time_avg, diag_script_info@region/), "_")
  workpath = config_user_info@work_dir + fname + ".nc"
  plotpath = config_user_info@plot_dir + fname

  ; Loop over datasets
  do imod = 0, nDatasets - 1

    log_debug("Processing " + datasetnames(imod))

    ; Determine start/end year
    start_year = info_items[imod]@start_year
    end_year = info_items[imod]@end_year

    ; Read data
    var = read_data(info_items[imod])
    dnames = getVarDimNames(var)

    ; Extract region and average over latitude and longitude
    if (any(dnames.eq."lat") .and. any(dnames.eq."lon")) then
      var_reg = area_operations(var, region(0), region(1), \
                                region(2), region(3), "average", True)
    else
      var_reg = var
    end if
    delete(var)

    ; Calculate time average
    var_avg = time_operations(var_reg, start_year, end_year, "average", \
                              diag_script_info@time_avg, True)

    ; Calculate time standard deviation (with lower/upper bounds)
    if (start_year.lt.end_year) then
      var_std = time_operations(var_reg, start_year, end_year, "stddev", \
                                diag_script_info@time_avg, True)
    else
      var_std = 0.
    end if
    delete(var_reg)

    ; Store in global array
    var_all(imod, :, 0) = var_avg
    var_all(imod, :, 1) = var_std
    delete(var_avg)
    delete(var_std)

  end do

  ; Write output
  if (diag_script_info@write_netcdf) then
    var_all@ncdf = workpath
    ncdf_outfile = ncdf_write(var_all, workpath)
  end if

  ; Plotting
  if (.not.diag_script_info@write_plots) then
    return
  end if

  ; Convert units for plotting (if required)
  if (isatt(diag_script_info, "plot_units")) then
    var_all = convert_units(var_all, diag_script_info@plot_units)
  end if

  ; Annotation and file names
  title = var_all@long_name
  caption = var0
  if (isdefined("level")) then
    title = title + " " + level + " hPa"
    caption = caption + level
    delete(level)
  end if
  title = title + " - " + diag_script_info@region

  ; Draw plot
  wks = gsn_open_wks(file_type, plotpath)
  wks@legendfile = plotpath + "_legend"
  var_all@res_tiMainString = title
  plot = cycle_plot(wks, var_all, var0, info_items)
  draw(plot)
  frame(wks)

  ; Call provenance logger
  log_provenance(ncdf_outfile, \
                 plotpath + "." + file_type, \
                 "Cycle plot of variable " + caption, \
                 (/"mean", "stddev"/), \
                 diag_script_info@region, \
                 "seas", \
                 (/"frank_franziska", "righi_mattia", "eyring_veronika"/), \
                 (/"righi15gmd", "gleckler08jgr"/), \
                 metadata_att_as_array(info_items, "filename"))

end
