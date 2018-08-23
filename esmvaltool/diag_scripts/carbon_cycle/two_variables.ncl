; #############################################################################
; carbon_cycle/two_variables.ncl
; #############################################################################
; DIAGNOSTIC SCRIPT FOR SCATTER PLOT WITH TWO VARIABLES
; Author: Alessandro Anav (University of Exeter, UK)
; EMBRACE project
; #############################################################################
;
; Description:
;   Draw a scatter plot with two variables.
;
; Required recipe attributes:
;    region: the region to be averaged.
;
; Optional recipe attributes:
;
; Caveats:
;
; Modification history
;   * 20180823_A-schl_ma: continued porting to v2.0
;   * 20180619_A-wenz_sa: adopted script to v2
;   * 20170511_A-gier_be: added metadata to files
;   * 20151130_A-laue_ax: corrected passing of file type from namelist to
;   *                     plotting functions; moved call to write_references
;   *                     to beginning of code
;   * 20151105_A-righ_ma: graphic part moved to external plot_scripts
;   * 20151104_A-righ_ma: removing regridding, not necessary for this
;   *                     diagnostic
;   * 201507??_A-wenz_sa: adapted to ESMValTool structure
;   * 201505??_A-anav_al: written
;
; #############################################################################

load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/scaling.ncl"
load "./diag_scripts/shared/set_operators.ncl"

load "./diag_scripts/shared/plot/scatterplot.ncl"
load "./diag_scripts/shared/plot/style.ncl"

; load "./diag_scripts/shared/watermark.ncl"


begin

  if (dimsizes(variables) .ne. 2) then
    error_msg("f", diag_script, "", "script needs exactly two variables")
  end if

  enter_msg(diag_script, "")
  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(diag_script + " (variables: " + variables(0) + ", " + \
           variables(1) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

  ; Variables
  var0 = variables(0)
  var1 = variables(1)
  field_type0 = field_types(0)
  field_type1 = field_types(1)

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

  ; Global arrays for collecting all data
  all_data_A = new((/DIM_DAT/), float)
  all_data_A!0 = "dataset"
  all_data_A&dataset = input_file_info@dataset
  all_data_B = new((/DIM_DAT/), float)
  all_data_B!0 = "dataset"
  all_data_B&dataset = input_file_info@dataset

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Loop over all datasets
  do idat = 0, DIM_DAT - 1
    dataset_name = input_file_info@dataset(idat)
    log_debug("Processing " + dataset_name)

    ; First variable
    ; FIXME: var_attr_exclude?
    if (var_attr_exclude(0) .ne. dataset_attr_id(idat)) then
      var_A = read_data(idat, var0, field_type0)

      ; Compute annual mean
      tmp_A = time_operations(var_A, \
                              toint(input_file_info@start_year(idat)), \
                              toint(input_file_info@end_year(idat)), \
                              "average", "annualclim", True)

      ; Compute area-weighted sum
      avg_A = area_operations(tmp_A, region(0), region(1), region(2), \
                              region(3), "sum", True)

      ; Assign to global array
      all_data_A(idat) = avg_A
      delete([/var_A, tmp_A, avg_A/])
    end if

    ; Second variable
    ; FIXME: var_attr_exclude?
    if (var_attr_exclude(1) .ne. dataset_attr_id(idat)) then
      var_B = read_data(idat, var1, field_type1)

      ; Compute annual mean
      tmp_B = time_operations(var_B, \
                              toint(input_file_info@start_year(idat)), \
                              toint(input_file_info@end_year(idat)), \
                              "average", "annualclim", True)
      delete(var_B)

      ; Computing area-weighted sum
      avg_B = area_operations(tmp_B, region(0), region(1), region(2), \
                              region(3), "sum", True)

      ; Assign to global array
      all_data_B(idat) = avg_B
      delete([/var_B, tmp_B, avg_B/])
    end if
  end do

  ; Convert units if appropriate
  if (isatt(variable_info, "plot_units")) then
    all_data_A = convert_units(all_data_A, variable_info@plot_units)
    all_data_B = convert_units(all_data_B, variable_info@plot_units)
  end if

  ; ---------------------------------------------------------------------------
  ; Write NETCDF
  ; ---------------------------------------------------------------------------

  if (config_user_info@write_netcdf) then
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)

    ; First variable
    new_path = work_dir + var0 + ".nc"
    all_data_A@var = var0
    all_data_A@diag_script = diag_script
    all_data_A@ncdf = new_path
    ncdf_outfile = ncdf_write(all_data_A, new_path)

    ; Second variable
    new_path = work_dir + var1 + ".nc"
    all_data_B@var = var1
    all_data_B@diag_script = diag_script
    all_data_B@ncdf = new_path
    ncdf_outfile = ncdf_write(all_data_B, new_path)
  end if

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  ; Get reference dataset index
  idx_ref_A = get_ref_dataset_idx(input_file_info, \
                                  variable_info@reference_dataset)
  idx_ref_A = idx_ref_B

  ; Replace FillValue in the dataset array with the value of the reference
  ; dataset. Since it is not possible to reduce the size of
  ; input_file_info@dataset, the reference dataset is simply duplicated.
  all_data_A = \
    where(ismissing(all_data_A), all_data_A(idx_ref_A), all_data_A)
  all_data_B = \
    where(ismissing(all_data_B), all_data_B(idx_ref_B), all_data_B)

  ; Collect data
  data_arr = new((/2, DIM_DAT/), float)
  data_arr!0 = "statistic"
  data_arr!1 = "dataset"
  data_arr&statistic = (/var0, var1/)
  data_arr&dataset = input_file_info@dataset
  data_arr(0, :) = (/all_data_A/)
  data_arr(1, :) = (/all_data_B/)
  data_arr@units = (/all_data_A@units, all_data_B@units/)
  data_arr@legend_outside = diag_script_info@legend_outside

  ; Draw plot
  outfile = plot_dir "_" + var0 + "-" + var1 + "_scatter_" + region@name
  wks = gsn_open_wks(file_type, outfile)
  wks@legendfile = outfile + "_legend"
  data_arr@res_tiMainString = var0 + "-" + var1 + " - " + \
    diag_script_info@region
  plot = scatterplot_markers(wks, data_arr, var0)
  plot@outfile = outfile + "." + file_type
  draw(plot)
  frame(wks)

  ; ; Add meta data (TODO)
  ; alltags = array_append_record(tags, region@DM_tag, 0)
  ; alltags := array_append_record(alltags, (/"PT_scatter", "ST_mean"/), 0)
  ; caption = "Scatterplot for " + var0 + " and " + var1 + \
  ;   " over the period " + start_year + "-" + end_year + \
  ;   ", similar to Anav et al., 2013, Fig. 12."
  ; id = diag_script
  ; climofiles = new(dimsizes(input_file_info@dataset) * dimsizes(variables), \
  ;                  string)
  ; do n_var = 0, dimsizes(variables) - 1
  ;   do idat = 0, dimsizes(input_file_info@dataset) - 1
  ;     climofiles(n_var * dimsizes(input_file_info@dataset) + idat) = \
  ;       interface_get_inpaths(idat) + "/" + \
  ;       interface_get_infile(variables(n_var), field_types(n_var), idat)
  ;   end do
  ; end do
  ; contrib_authors = (/"A-anav_al", "A-wenz_sa", "A-righ_ma", "A-schl_ma"/)
  ; ESMValMD(plot@outfile, alltags, caption, id, variables, \
  ;          input_file_info@dataset, climofiles, diag_script, contrib_authors)
  ; delete([/alltags, caption, id, climofiles/])
  delete(data_arr)

  leave_msg(diag_script, "")

end
