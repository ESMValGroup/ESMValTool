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
; Required recipe attributes:
;     region: region to be averaged.
;     ref_dataset: reference for the MVI calculation.
;
; Optional recipe attributes:
;
; Caveats:
;     The derived variables pr-mmday and tas-degC needs to be used in order
;     to get the correct linear trend. For land variables nbp, gpp and lai, it
;     is recommended to use the masked version (NBP, GPP and LAI, respectively)
;     in order to have a consistent land mask for all datasets.
;
; Modification history:
;   * 20180814_A-schl_ma: continued porting to v2.0
;   * 20180619_A-wenz_sa: adopted script to v2
;   * 20151130_A-laue_ax: moved call to write_references to beginning of code
;   * 20151104_A-righ_ma: graphic part moved to external plot_scripts
;   * 20151102_A-righ_ma: replacing area functions with ESMValTool functions
;   * 201507??_A-wenz_sa: adapted to ESMValTool structure
;   * 201505??_A-anav_al: written
;
; #############################################################################

load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/scaling.ncl"
load "./diag_scripts/shared/statistics.ncl"

load "./diag_scripts/shared/plot/scatterplot.ncl"
load "./diag_scripts/shared/plot/style.ncl"


begin

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (var: " + variables(0) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Variable
  var0 = variables(0)
  field_type0 = field_types(0)

  ; Check field type
  plot_fields = (/"T2Ms"/)
  ind_f = any(field_type0 .eq. plot_fields)
  if (all(ismissing(ind_f))) then
    error_msg("f", diag_script, "", "cannot use field " + field_type0)
  end if
  delete([/plot_fields, ind_f/])

  ; Write references (TODO)
  ; write_references(diag_script, \
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
  start_year = min(stringtoint(input_file_info@start_year))
  end_year = max(stringtoint(input_file_info@end_year))
  all_years = ispan(start_year, end_year, 1)

  ; Array sizes
  DIM_VAR = dimsizes(variables)
  DIM_DAT = dimsizes(input_file_info@dataset)
  DIM_TIME = dimsizes(all_years)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Calculate yearly average
  do idat = 0, DIM_DAT - 1
    dataset_name = input_file_info@dataset(idat)
    log_debug("Processing " + dataset_name)

    ; Read data
    data_var = read_data(idat, var0, field_type0)
    if (typeof(data_var) .eq. "double") then
      data_var := dble2flt(data_var)
    end if

    ; Get dimensions
    if (isdim(data_var, "lat") .and. isdim(data_var, "lon")) then
      DIM_LAT = dimsizes(data_var&lat)
      DIM_LON = dimsizes(data_var&lon)
    else
      ; TODO
      error_msg("f", diag_script, "", "Irregular grids are currently not " + \
                "supported")
    end if

    ; Computing annual mean
    data_annual = time_operations(data_var, \
                                  toint(input_file_info@start_year(idat)), \
                                  toint(input_file_info@end_year(idat)), \
                                  "average", "yearly", True)

    ; Collect data in global array
    if (.not. isvar("all_data_yearly")) then
      all_data_yearly = \
        new((/DIM_DAT, DIM_TIME, DIM_LAT, DIM_LON/), float)
      all_data_yearly!0 = "dataset"
      all_data_yearly!1 = "year"
      all_data_yearly!2 = "lat"
      all_data_yearly!3 = "lon"
      all_data_yearly&dataset = input_file_info@dataset
      all_data_yearly&year = all_years
      all_data_yearly&lat = data_var&lat
      all_data_yearly&lon = data_var&lon
    end if

    ; Check lat and lon dimensions
    if (DIM_LAT .ne. dimsizes(all_data_yearly&lat) .or. \
        DIM_LON .ne. dimsizes(all_data_yearly&lon)) then
      error_msg("f", diag_script, "", "Not all datasets are on the same " + \
                "grid (" + dataset_name + "), select 'regrid' in " + \
                "preprocessor options")
    end if

    ; Save data at correct place
    idx_1 = start_year - toint(input_file_info@start_year(idat))
    idx_2 = dimsizes(all_years) - \
      (end_year - toint(input_file_info@end_year(idat))) - 1
    all_data_yearly(idat, idx_1:idx_2, :, :) = data_annual
    copy_VarAtts(data_annual, all_data_yearly)
    copy_VarCoords(data_annual, all_data_yearly(idat, :, :, :))
    delete([/data_var, data_annual/])
  end do

  ; ---------------------------------------------------------------------------
  ; Compute means, trends and model variability indices (MVI) globally
  ; ---------------------------------------------------------------------------

  ; MVI array
  all_data_mvi = new((/DIM_DAT, DIM_LAT, DIM_LON/), float)
  all_data_mvi!0 = "dataset"
  all_data_mvi!1 = "lat"
  all_data_mvi!2 = "lon"
  all_data_mvi&dataset = input_file_info@dataset
  all_data_mvi&lat = all_data_yearly&lat
  all_data_mvi&lon = all_data_yearly&lon
  all_data_mvi@_FillValue = 1.0e20

  ; Trend array
  all_data_trend = new((/DIM_DAT, DIM_LAT, DIM_LON/), float)
  all_data_trend!0 = "dataset"
  all_data_trend!1 = "lat"
  all_data_trend!2 = "lon"
  all_data_trend&dataset = input_file_info@dataset
  all_data_trend&lat = all_data_yearly&lat
  all_data_trend&lon = all_data_yearly&lon
  all_data_trend@_FillValue = 1.0e20

  ; Determine reference dataset
  ref_idx = get_ref_dataset_idx(input_file_info, \
                                variable_info@reference_dataset)
  if (ismissing(ref_idx)) then
      error_msg("f", diag_script, "", "no adequate reference dataset provided")
  end if
  log_info("Reference dataset: " + input_file_info@dataset(ref_idx))
  data_ref = all_data_yearly(ref_idx, :, :, :)
  ref = dim_stddev_n(data_ref, 0)
  ref = where(ref .ne. 0, ref, all_data_yearly@_FillValue)

  ; Loop over datasets
  do idat = 0, DIM_DAT - 1
    data_temp = all_data_yearly(idat, :, :, :)
    dataset = dim_stddev_n(data_temp, 0)

    ; Prevents a division by 0
    dataset = where(dataset .ne. 0, dataset, all_data_yearly@_FillValue)

    ; Commpute MVI
    all_data_mvi(idat, :, :) = (dataset / ref - ref / dataset) ^ 2
    delete([/dataset, data_temp/])

    ; Compute trend
    dataset = all_data_yearly(idat, :, :, :)
    rc = regCoef_n(1.0 * all_years, dataset, 0, 0)
    rc!0 = "lat"
    rc!1 = "lon"
    rc&lat = dataset&lat
    rc&lon = dataset&lon
    all_data_trend(idat, :, :) = rc
    delete([/dataset, rc/])
  end do
  copy_VarAtts(all_data_yearly, all_data_trend)

  ; Cap MVI
  all_data_mvi = where(all_data_mvi .lt. 10, all_data_mvi, 50)

  ; Special case for LAI:
  ; Since datasets have different land covers, all sparse vegetated points must
  ; be masked
  if (var0 .eq. "lai") then
    all_data_yearly = where(all_data_yearly .ge. 0.025, \
                            all_data_yearly, 1.0e20)
  end if

  ; Compute total temporal mean
  all_data_mean = dim_avg_n_Wrap(all_data_yearly, 1)
  all_data_mean!0 = "dataset"
  all_data_mean!1 = "lat"
  all_data_mean!2 = "lon"
  all_data_mean&dataset = input_file_info@dataset
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
  regional_mean&dataset = input_file_info@dataset
  regional_trend = new((/DIM_DAT/), typeof(all_data_trend))
  regional_trend!0 = "dataset"
  regional_trend&dataset = input_file_info@dataset
  regional_mvi = new((/DIM_DAT/), typeof(all_data_mvi))
  regional_mvi!0 = "dataset"
  regional_mvi&dataset = input_file_info@dataset

  ; Set type of area operation (flux vs. non-flux variables)
  if (all_data_mean@units .eq. "kg m-2 s-1" .and. var0 .ne. "pr") then
    sp_opt = "sum"
  else
    sp_opt = "average"
  end if

  ; Regional averages for every dataset
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
  if (isatt(variable_info, "plot_units")) then
    new_units = variable_info@plot_units
    regional_mean = convert_units(regional_mean, new_units)
    if (new_units .ne. "degC") then
      trend_tmp = convert_units(regional_trend, new_units)
    else
      trend_tmp = regional_trend
      trend_tmp@units = "degC"
    end if
    if (input_file_info@mip(0) .eq. "Lmon") then
      regional_trend = trend_tmp
      regional_trend@units = trend_tmp@units + " y-1"
    else
      regional_trend = trend_tmp * 10.0
      regional_trend@units = trend_tmp@units + " decade-1"
    end if
    delete(trend_tmp)
  end if

  ; ---------------------------------------------------------------------------
  ; Write NETCDF
  ; ---------------------------------------------------------------------------

  if (config_user_info@write_netcdf) then
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)

    ; Mean
    new_path = work_dir + "mean_" + var0 + ".nc"
    regional_mean@var = var0
    regional_mean@diag_script = diag_script
    regional_mean@ncdf = new_path
    ncdf_outfile = ncdf_write(regional_mean, new_path)

    ; Trend
    new_path = work_dir + "trend_" + var0 + ".nc"
    regional_trend@var = var0
    regional_trend@diag_script = diag_script
    regional_trend@ncdf = new_path
    ncdf_outfile = ncdf_write(regional_trend, new_path)

    ; MVI
    new_path = work_dir + "MVI_" + var0 + ".nc"
    regional_mvi@var = var0
    regional_mvi@diag_script = diag_script
    regional_mvi@ncdf = new_path
    ncdf_outfile = ncdf_write(regional_mvi, new_path)
  end if

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  ; Set plot output file
  outfile = plot_dir + var0 + "_" + region@name
  wks = gsn_open_wks(file_type, outfile)
  gsn_define_colormap(wks, "BlAqGrYeOrReVi200")
  gsn_reverse_colormap(wks)

  ; Pack data
  data_arr = new((/3, dimsizes(regional_mean)/), float)
  data_arr(0, :) = (/regional_mean/)
  data_arr(1, :) = (/regional_trend/)
  data_arr(2, :) = (/regional_mvi/)
  data_arr!0 = "statistic"
  data_arr!1 = "datasets"
  data_arr&statistic = (/"mean", "linear trend", "MVI"/)
  data_arr&datasets = input_file_info@dataset
  data_arr@units = (/regional_mean@units, regional_trend@units, "1"/)

  ; Set levels for the color coding
  tmp = regional_mvi
  tmp(ref_idx) = tmp@_FillValue
  nlevs = 11
  tmp := nice_mnmxintvl(min(tmp), max(tmp), nlevs, True)
  levs = decimalPlaces(fspan(tmp(0), tmp(1), 11), 2, True)
  delete(tmp)

  ; Draw a 3D scatterplot
  data_arr@res_tiMainString = var0 + " - " + diag_script_info@region
  plot = scatterplot3D(wks, data_arr, var0, levs)
  log_info("Wrote " + outfile + "." + file_type)
  draw(plot)
  frame(wks)

  ; ; Add meta data (TODO)
  ; climofiles = new(dimsizes(input_file_info@dataset) * DIM_VAR, \
  ;                  string)
  ; alltags = array_append_record(tags, region@DM_tag, 0)
  ; alltags := array_append_record(alltags, \
  ;                                (/"PT_scatter", "PT_zonal", "ST_mean", \
  ;                                 "ST_trend", "ST_clim", "ST_var"/), 0)
  ; caption =  "Scatterplot for multiyear average " + var0 + \
  ;            " in x axis, its linear trend in y axis, and MVI." + \
  ;            " Like Anav et al. Fig 1 bottom"
  ; id = diag_script
  ; do n_var = 0, DIM_VAR - 1
  ;    do idat = 0, dimsizes(input_file_info@dataset) - 1
  ;     climofiles(n_var * dimsizes(input_file_info@dataset) + idat) = \
  ;       interface_get_inpaths(idat) + "/" + \
  ;       interface_get_infile(variables(n_var), field_types(n_var), idat)
  ;   end do
  ; end do
  ; contrib_authors = (/"A-anav_al", "A-wenz_sa", "A-righ_ma", "A-schl_ma"/)

  ; ESMValMD(outfile + "." + file_type, alltags, caption, id, variables, \
  ;          input_file_info@dataset, climofiles, diag_script, contrib_authors)

  ; delete([/alltags, caption, id, climofiles/])
  delete(data_arr)

  leave_msg(diag_script, "")

end
