; #############################################################################
; carbon_cycle/mvi.ncl
; #############################################################################
; DIAGNOSTIC SCRIPT TO CREATE MVI-IAV-TREND PLOT
; Author: Alessandro Anav (University of Exeter, UK)
; EMBRACE project
; #############################################################################
;
; Description:
;    Calculates the model variability index (MVI), interannual variability
;    (IAV) and mean, and draws them in a 3D scatter plot.
;
; Required diag_script_info attributes:
;     region: region to be averaged.
;     reference_dataset: reference for the MVI calculation.
;
; Optional diag_script_info attributes:
;     mean_time_range: time period over which the mean is calculated (default:
;                      whole time span).
;     trend_time_range: time period over which the trend is calculated
;                       (default: whole time span).
;     mvi_time_range: time period over which the MVI is calculated (default:
;                     whole time span).
;
; Caveats:
;     If CRU is the reference dataset it is important to apply a landmask in
;     the preprocessor section of the recipe.
;
; Modification history:
;     20180814_A-schl_ma: continued porting to v2.0.
;     20180619_A-wenz_sa: adopted script to v2.
;     20170316_A-gier_be: Added cdo regridding, landmask for CRU obs.
;     20151130_A-laue_ax: moved call to write_references to beginning of code.
;     20151104_A-righ_ma: graphic part moved to external plot_scripts.
;     20151102_A-righ_ma: replacing area functions with ESMValTool functions.
;     201507??_A-wenz_sa: adapted to ESMValTool structure.
;     201505??_A-anav_al: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/statistics.ncl"

load "$diag_scripts/shared/plot/scatterplot.ncl"
load "$diag_scripts/shared/plot/style.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")

  ; Variable
  VAR0 = variable_info[0]
  var0 = VAR0@short_name
  DIM_VAR = ListCount(variable_info)
  if (DIM_VAR .gt. 1) then
    error_msg("w", DIAG_SCRIPT, "", "this diagnostic supports only one " + \
              "variable, processing " + var0)
  end if

  ; Input data
  INFO0 = select_metadata_by_name(input_file_info, var0)
  DATASETS = metadata_att_as_array(INFO0, "dataset")
  DIM_DAT = dimsizes(DATASETS)

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Write references (TODO)
  ; write_references(DIAG_SCRIPT, \
  ;                  (/"A-anav_al"/), \
  ;                  (/"A-wenz_sa", "A-righ_ma", "A-schl_ma"/), \
  ;                  (/"D_anav13jclim"/), \
  ;                  (/""/), \
  ;                  (/"P_embrace"/))

end

begin

  ; ---------------------------------------------------------------------------
  ; Read recipe and config data
  ; ---------------------------------------------------------------------------

  ; Plot file type
  file_type = config_user_info@output_file_type
  if (ismissing(file_type)) then
    file_type = "ps"
  end if

  ; Output plot directory
  plot_dir = config_user_info@plot_dir
  system("mkdir -p " + plot_dir)

  ; Check required diag_script_info attributes
  req_atts = (/"region"/)
  exit_if_missing_atts(diag_script_info, req_atts)
  delete(req_atts)

  ; Set region
  region = select_region(diag_script_info@region)

  ; Set time range
  start_year = 9999
  end_year = -9999
  do idat = 0, DIM_DAT - 1
    if (INFO0[idat]@start_year .lt. start_year) then
      start_year = INFO0[idat]@start_year
    end if
    if (INFO0[idat]@end_year .gt. end_year) then
      end_year = INFO0[idat]@end_year
    end if
  end do
  all_years = ispan(start_year, end_year, 1)
  DIM_TIME = dimsizes(all_years)

  ; Specific time range for mean
  if (isatt(diag_script_info, "mean_time_range")) then
    custom_mean_time_range = True
    mean_time_range = diag_script_info@mean_time_range
    mean_time_ind = (/ind(all_years .eq. mean_time_range(0)), \
                     ind(all_years .eq. mean_time_range(1))/)
    log_info("Restricting mean calculation to time period " + \
             mean_time_range(0) + "-" + mean_time_range(1))
  else
    custom_mean_time_range = False
  end if

  ; Specific time ranges for trend
  if (isatt(diag_script_info, "trend_time_range")) then
    custom_trend_time_range = True
    trend_time_range = diag_script_info@trend_time_range
    trend_time_ind = (/ind(all_years .eq. trend_time_range(0)), \
                      ind(all_years .eq. trend_time_range(1))/)
    log_info("Restricting trend calculation to time period " + \
             trend_time_range(0) + "-" + trend_time_range(1))
  else
    custom_trend_time_range = False
  end if

  ; Specific time ranges for MVI
  if (isatt(diag_script_info, "mvi_time_range")) then
    custom_mvi_time_range = True
    mvi_time_range = diag_script_info@mvi_time_range
    mvi_time_ind = (/ind(all_years .eq. mvi_time_range(0)), \
                    ind(all_years .eq. mvi_time_range(1))/)
    log_info("Restricting mvi calculation to time period " + \
             mvi_time_range(0) + "-" + mvi_time_range(1))
  else
    custom_mvi_time_range = False
  end if

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Calculate yearly average
  do idat = 0, DIM_DAT - 1
    data = INFO0[idat]
    dataset_name = DATASETS(idat)
    log_debug("Processing " + dataset_name)

    ; Read data
    data_var = read_data(data)
    if (typeof(data_var) .eq. "double") then
      data_var := dble2flt(data_var)
    end if

    ; Set fill value if necessary
    if (.not. isatt(data_var, "_FillValue")) then
      data_var@_FillValue = 1.0e20
    end if

    ; Get dimensions
    if (isdim(data_var, "lat") .and. isdim(data_var, "lon")) then
      DIM_LAT = dimsizes(data_var&lat)
      DIM_LON = dimsizes(data_var&lon)
    else
      ; TODO
      error_msg("f", DIAG_SCRIPT, "", "Irregular grids are currently not " + \
                "supported")
    end if

    ; Cut off Antarctica if necessary (for land-based observations)
    if ((isatt(VAR0, "reference_dataset") .and. \
         VAR0@reference_dataset .eq. "CRU")) then
      if (dimsizes(data_var&lat) .gt. 1) then
        log_info("Masking Antarctica")
        data_var(:, {:-60}, :) = data_var@_FillValue
      end if
    end if

    ; Computing annual mean
    data_annual = time_operations(data_var, data@start_year, data@end_year, \
                                  "average", "yearly", True)

    ; Collect data in global array
    if (.not. isvar("all_data_yearly")) then
      all_data_yearly = \
        new((/DIM_DAT, DIM_TIME, DIM_LAT, DIM_LON/), float)
      all_data_yearly!0 = "dataset"
      all_data_yearly!1 = "year"
      all_data_yearly!2 = "lat"
      all_data_yearly!3 = "lon"
      all_data_yearly&dataset = DATASETS
      all_data_yearly&year = all_years
      all_data_yearly&lat = data_var&lat
      all_data_yearly&lon = data_var&lon
      all_data_yearly@_FillValue = 1.0e20
    end if

    ; Check lat and lon dimensions
    if (DIM_LAT .ne. dimsizes(all_data_yearly&lat) .or. \
        DIM_LON .ne. dimsizes(all_data_yearly&lon)) then
      error_msg("f", DIAG_SCRIPT, "", "Not all datasets are on the same " + \
                "grid (" + dataset_name + "), select 'regrid' in " + \
                "preprocessor options")
    end if

    ; Save data at correct place
    idx_1 = data@start_year - start_year
    idx_2 = dimsizes(all_years) - (end_year - data@end_year) - 1
    all_data_yearly(idat, idx_1:idx_2, :, :) = data_annual
    copy_VarAtts(data_annual, all_data_yearly)
    copy_VarCoords(data_annual, all_data_yearly(idat, :, :, :))
    delete([/data_var, data_annual, idx_1, idx_2/])
  end do

  ; ---------------------------------------------------------------------------
  ; Compute means, trends and model variability indices (MVI) globally
  ; ---------------------------------------------------------------------------

  ; MVI array
  all_data_mvi = new((/DIM_DAT, DIM_LAT, DIM_LON/), float)
  all_data_mvi!0 = "dataset"
  all_data_mvi!1 = "lat"
  all_data_mvi!2 = "lon"
  all_data_mvi&dataset = DATASETS
  all_data_mvi&lat = all_data_yearly&lat
  all_data_mvi&lon = all_data_yearly&lon
  all_data_mvi@_FillValue = 1.0e20

  ; Trend array
  all_data_trend = new((/DIM_DAT, DIM_LAT, DIM_LON/), float)
  all_data_trend!0 = "dataset"
  all_data_trend!1 = "lat"
  all_data_trend!2 = "lon"
  all_data_trend&dataset = DATASETS
  all_data_trend&lat = all_data_yearly&lat
  all_data_trend&lon = all_data_yearly&lon
  all_data_trend@_FillValue = 1.0e20

  ; Determine reference dataset
  ref_idx = ind(DATASETS .eq. VAR0@reference_dataset)
  ref_idx := ref_idx(0)
  if (ismissing(ref_idx)) then
      error_msg("f", DIAG_SCRIPT, "", "no adequate reference dataset provided")
  end if
  log_info("Reference dataset: " + DATASETS(ref_idx))

  ; Restrict time range of reference dataset if desired
  if (custom_mvi_time_range) then
    data_ref = all_data_yearly(ref_idx, mvi_time_ind(0):mvi_time_ind(1), \
                               :, :)
  else
    data_ref = all_data_yearly(ref_idx, :, :, :)
  end if
  ref = dim_stddev_n(data_ref, 0)
  ref = where(ref .ne. 0, ref, ref@_FillValue)

  ; Loop over datasets
  do idat = 0, DIM_DAT - 1

    ; Datasets for MVI (restrict time range if desired)
    if (custom_mvi_time_range) then
      data_temp = all_data_yearly(idat, mvi_time_ind(0):mvi_time_ind(1), \
                                  :, :)
    else
      data_temp = all_data_yearly(idat, :, :, :)
    end if
    dataset = dim_stddev_n(data_temp, 0)
    dataset = where(dataset .ne. 0, dataset, dataset@_FillValue)

    ; Compute MVI
    all_data_mvi(idat, :, :) = (dataset / ref - ref / dataset) ^ 2
    delete([/dataset, data_temp/])

    ; Compute trend (restrict time range if desired)
    if (custom_trend_time_range) then
      dataset = all_data_yearly(idat, \
                                trend_time_ind(0):trend_time_ind(1), :, :)
      rc = regCoef_n( \
        1.0 * ispan(trend_time_range(0), trend_time_range(1), 1), \
        dataset, 0, 0)
    else
      dataset = all_data_yearly(idat, :, :, :)
      rc = regCoef_n(1.0 * all_years, dataset, 0, 0)
    end if
    rc!0 = "lat"
    rc!1 = "lon"
    rc&lat = dataset&lat
    rc&lon = dataset&lon
    all_data_trend(idat, :, :) = rc
    delete([/dataset, rc/])
  end do
  copy_VarAtts(all_data_yearly, all_data_trend)

  ; Special case for lai:
  ; Since datasets have different land covers, all sparse vegetated points must
  ; be masked
  if (isStrSubset(var0, "lai")) then
    all_data_yearly = where(all_data_yearly .ge. 0.025, \
                            all_data_yearly, all_data_yearly@_FillValue)
  end if

  ; Compute total temporal mean
  if (custom_mean_time_range) then
    all_data_mean = dim_avg_n_Wrap(all_data_yearly(\
      :, mean_time_ind(0):mean_time_ind(1), :, :), 1)
  else
    all_data_mean = dim_avg_n_Wrap(all_data_yearly, 1)
  end if
  all_data_mean!0 = "dataset"
  all_data_mean!1 = "lat"
  all_data_mean!2 = "lon"
  all_data_mean&dataset = DATASETS
  all_data_mean&lat = all_data_yearly&lat
  all_data_mean&lon = all_data_yearly&lon
  all_data_mean@units = all_data_yearly@units
  delete(all_data_yearly)

  ; ---------------------------------------------------------------------------
  ; Compute means, trends and MVI regionally averaged
  ; ---------------------------------------------------------------------------

  ; Setup arrays
  regional_mean = new((/DIM_DAT/), typeof(all_data_mean))
  regional_mean!0 = "dataset"
  regional_mean&dataset = DATASETS
  regional_trend = new((/DIM_DAT/), typeof(all_data_trend))
  regional_trend!0 = "dataset"
  regional_trend&dataset = DATASETS
  regional_mvi = new((/DIM_DAT/), typeof(all_data_mvi))
  regional_mvi!0 = "dataset"
  regional_mvi&dataset = DATASETS

  ; Set type of area operation (flux vs. non-flux variables)
  ; Flux variables need to be integrated and area weighted ("sum"), non-flux
  ; variables need to be averaged, area weighted and normalized ("average")
  if (all_data_mean@units .eq. "kg m-2 s-1" .and. \
      var0 .ne. "pr") then
    sp_opt = "sum"
  else
    sp_opt = "average"
  end if

  ; Regional averages/sums for every dataset
  do idat = 0, DIM_DAT - 1
    regional_mean(idat) = area_operations(all_data_mean(idat, :, :), \
                                          region(0), region(1), region(2), \
                                          region(3), sp_opt, True)
    regional_trend(idat) = area_operations(all_data_trend(idat, :, :), \
                                           region(0), region(1), region(2), \
                                           region(3), sp_opt, True)
    regional_mvi(idat) = area_operations(all_data_mvi(idat, :, :), \
                                         region(0), region(1), region(2), \
                                         region(3), "average", True)
  end do

  ; Diagnostic- and variable-specific units conversions
  regional_mvi@units = "1"
  trend_tmp = regional_trend
  if (isatt(VAR0, "plot_units")) then
    new_units = VAR0@plot_units
    regional_mean = convert_units(regional_mean, new_units)
    if (new_units .eq. "degC") then
      trend_tmp@units = "degC"
    else
      trend_tmp = convert_units(regional_trend, new_units)
    end if
  end if

  ; Units conversion for trend
  if (INFO0[0]@mip .eq. "Lmon") then
    regional_trend = trend_tmp
    temp_unit = "yr-1"
  else
    regional_trend = trend_tmp * 10.0
    temp_unit = "decade-1"
  end if
  if (trend_tmp@units .ne. "1") then
    regional_trend@units = trend_tmp@units + " " + temp_unit
  else
    regional_trend@units = temp_unit
  end if
  delete([/trend_tmp, temp_unit/])

  ; ---------------------------------------------------------------------------
  ; Write NETCDF
  ; ---------------------------------------------------------------------------

  if (config_user_info@write_netcdf) then
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)

    ; Mean
    new_path = work_dir + "mean_" + var0 + ".nc"
    regional_mean@var = var0
    regional_mean@diag_script = DIAG_SCRIPT
    regional_mean@ncdf = new_path
    ncdf_outfile_mean = ncdf_write(regional_mean, new_path)

    ; Trend
    new_path = work_dir + "trend_" + var0 + ".nc"
    regional_trend@var = var0
    regional_trend@diag_script = DIAG_SCRIPT
    regional_trend@ncdf = new_path
    ncdf_outfile_trend = ncdf_write(regional_trend, new_path)

    ; MVI
    new_path = work_dir + "MVI_" + var0 + ".nc"
    regional_mvi@var = var0
    regional_mvi@diag_script = DIAG_SCRIPT
    regional_mvi@ncdf = new_path
    ncdf_outfile_mvi = ncdf_write(regional_mvi, new_path)
  end if

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  if (config_user_info@write_plots) then
    ; Set plot output file
    outfile = plot_dir + var0 + "_" + region@name
    wks = gsn_open_wks(file_type, outfile)
    gsn_define_colormap(wks, "BlAqGrYeOrReVi200")
    gsn_reverse_colormap(wks)

    ; Axis labels
    if (custom_mean_time_range) then
      mean_label = mean_time_range(0) + "-" + mean_time_range(1) + " Mean"
    else
      mean_label = start_year + "-" + end_year + " Mean"
    end if
    if (custom_trend_time_range) then
      trend_label = trend_time_range(0) + "-" + trend_time_range(1) + \
        " Linear trend"
    else
      trend_label = start_year + "-" + end_year + " Linear trend"
    end if
    if (custom_mvi_time_range) then
      mvi_label = mvi_time_range(0) + "-" + mvi_time_range(1) + "~C~MVI"
    else
      mvi_label = start_year + "-" + end_year + "~C~MVI"
    end if

    ; Pack data
    data_arr = new((/3, dimsizes(regional_mean)/), float)
    data_arr(0, :) = (/regional_mean/)
    data_arr(1, :) = (/regional_trend/)
    data_arr(2, :) = (/regional_mvi/)
    data_arr!0 = "statistic"
    data_arr!1 = "datasets"
    data_arr&statistic = (/mean_label, trend_label, mvi_label/)
    data_arr&datasets = DATASETS
    data_arr@units = (/regional_mean@units, regional_trend@units, "1"/)

    ; Set levels for the color coding
    tmp = regional_mvi
    tmp(ref_idx) = tmp@_FillValue
    nlevs = 11
    tmp := nice_mnmxintvl(min(tmp) - 0.01, max(tmp) + 0.01, nlevs, True)
    levs = decimalPlaces(fspan(tmp(0), tmp(1), 11), 2, True)
    delete(tmp)

    ; Draw a 3D scatterplot
    plot_path = outfile + "." + file_type
    data_arr@res_tiMainString = var0 + " - " + diag_script_info@region
    plot = scatterplot3D(wks, data_arr, var0, levs)
    log_info("Wrote " + plot_path)
    draw(plot)
    frame(wks)

    ; Provenance tracking
    caption =  ("Scatterplot for multiyear average " + var0 + " in x " + \
                "axis, its linear trend in y axis, and MVI. Similar to " + \
                "Anav et al. (2013) figure 1 (bottom).")
    domains = diag_script_info@region
    statistics = (/"mean", "trend", "mvi"/)
    authors = (/"anav_alessandro", "wenzel_sabrina", "righi_mattia", \
                "schlund_manuel"/)
    plot_type = "scatter"
    references = (/"anav13jclim"/)
    ancestors = metadata_att_as_array(INFO0, "filename")
    log_provenance(ncdf_outfile_mvi, plot_path, caption, statistics, domains, \
                   plot_type, authors, references, ancestors)
    delete(data_arr)
  end if

  leave_msg(DIAG_SCRIPT, "")

end
