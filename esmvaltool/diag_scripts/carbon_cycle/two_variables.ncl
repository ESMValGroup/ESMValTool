; #############################################################################
; carbon_cycle/two_variables.ncl
; #############################################################################
; DIAGNOSTIC SCRIPT FOR SCATTER PLOT WITH TWO VARIABLES
; Author: Alessandro Anav (University of Exeter, UK)
; EMBRACE project
; #############################################################################
;
; Description:
;     Draw a scatter plot with two variables.
;
; Required diag_script_info attributes:
;     region: the region to be averaged.
;
; Modification history:
;     20180823_A-schl_ma: continued porting to v2.0.
;     20180619_A-wenz_sa: adopted script to v2.
;     20170511_A-gier_be: added metadata to files.
;     20151130_A-laue_ax: corrected passing of file type from namelist to
;                         plotting functions; moved call to write_references
;                         to beginning of code.
;     20151105_A-righ_ma: graphic part moved to external plot_scripts.
;     20151104_A-righ_ma: removing regridding, not necessary for this
;                         diagnostic.
;     201507??_A-wenz_sa: adapted to ESMValTool structure.
;     201505??_A-anav_al: written.
;
; #############################################################################

load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/scaling.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/statistics.ncl"

load "./diag_scripts/shared/plot/scatterplot.ncl"
load "./diag_scripts/shared/plot/style.ncl"


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
  DIM_VAR = dimsizes(variables)

  ; Check field types
  plot_fields = (/"T2Ms"/)
  do ivar = 0, DIM_VAR - 1
    ind_f := any(field_types(ivar) .eq. plot_fields)
    if (all(ismissing(ind_f))) then
      error_msg("f", diag_script, "", "cannot use field " + field_types(ivar))
    end if
  end do
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

  ; Set time range, array sizes and new arrays
  all_years = NewList("lifo")
  all_data = NewList("lifo")
  start_year = new(DIM_VAR, integer)
  end_year = new(DIM_VAR, integer)
  DIM_DAT = new(DIM_VAR, integer)
  DIM_TIME = new(DIM_VAR, integer)
  do ivar = 0, DIM_VAR - 1
    ; Time
    start_year(ivar) = min(stringtoint(input_file_info[ivar]@start_year))
    end_year(ivar) = min(stringtoint(input_file_info[ivar]@end_year))
    ListAppend(all_years, ispan(start_year(ivar), end_year(ivar), 1))

    ; Array sizes
    DIM_DAT(ivar) = dimsizes(input_file_info[ivar]@dataset)
    DIM_TIME(ivar) = dimsizes(all_years[ivar])

    ; New arrays
    ListAppend(all_data, new((/DIM_DAT(ivar)/), float))
    all_data[ivar]!0 = "dataset"
    all_data[ivar]&dataset = input_file_info[ivar]@dataset
  end do
  if (DIM_DAT(0) .ne. DIM_DAT(1)) then
    error_msg("f", diag_script, "", "the two variables do not cover the " + \
              "same number of datasets, this is not supported yet")
  end if

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Loop over all variables and datasets
  do ivar = 0, DIM_VAR - 1
    do idat = 0, DIM_DAT(ivar) - 1
      dataset_name = input_file_info[ivar]@dataset(idat)
      log_debug("Processing " + variables(ivar) + " of " + dataset_name)

      ; Read variable
      var = read_data(idat, variables(ivar), field_types(ivar))

      ; Compute annual mean
      tmp = time_operations(var, \
                            toint(input_file_info[ivar]@start_year(idat)), \
                            toint(input_file_info[ivar]@end_year(idat)), \
                            "average", "annualclim", True)

      ; Compute area-weighted sum
      area_op = area_operations(tmp, region(0), region(1), region(2), \
                                region(3), "sum", True)

      ; Assign to global array
      all_data[ivar](idat) = area_op
      delete([/var, tmp, area_op/])
    end do
  end do

  ; Convert units if appropriate
  do ivar = 0, DIM_VAR - 1
    if (isatt(variable_info[ivar], "plot_units")) then
      all_data[ivar] = convert_units(all_data[ivar], \
                                     variable_info[ivar]@plot_units)
    end if
  end do

  ; ---------------------------------------------------------------------------
  ; Write NETCDF
  ; ---------------------------------------------------------------------------

  if (config_user_info@write_netcdf) then
    work_dir = config_user_info@work_dir
    system("mkdir -p " + work_dir)

    ; Loop over variables
    do ivar = 0, DIM_VAR - 1
      new_path = work_dir + variables(ivar) + ".nc"
      all_data[ivar]@var = variables(ivar)
      all_data[ivar]@diag_script = diag_script
      all_data[ivar]@ncdf = new_path
      ncdf_outfile = ncdf_write(all_data[ivar], new_path)
    end do
  end if

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  if (config_user_info@write_plots) then
    ; Output array
    data_arr = new((/DIM_VAR, DIM_DAT(0)/), float)
    data_arr!0 = "statistic"
    data_arr!1 = "dataset"
    data_arr&statistic = (/var0, var1/)
    data_arr&dataset = input_file_info[0]@dataset
    data_arr@units = (/all_data[0]@units, all_data[1]@units/)
    data_arr@legend_outside = diag_script_info@legend_outside

    ; Collect data
    do ivar = 0, DIM_VAR - 1
      idx_ref = get_ref_dataset_idx(input_file_info[ivar], \
                                    variable_info[ivar]@reference_dataset)

      ; Replace FillValue in the dataset array with the value of the reference
      ; dataset. Since it is not possible to reduce the size of
      ; input_file_info[ivar]@dataset, the reference dataset is simply
      ; duplicated.
      all_data[ivar] = \
        where(ismissing(all_data[ivar]), all_data[ivar](idx_ref), \
              all_data[ivar])
      data_arr(ivar, :) = (/all_data[ivar]/)
    end do

    ; Draw plot
    outfile = plot_dir + "_" + var0 + "-" + var1 + "_scatter_" + region@name
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
    ; climofiles = new(dimsizes(input_file_info@dataset) * DIM_VAR, \
    ;                  string)
    ; do ivar = 0, DIM_VAR - 1
    ;   do idat = 0, dimsizes(input_file_info@dataset) - 1
    ;     climofiles(ivar * dimsizes(input_file_info@dataset) + idat) = \
    ;       interface_get_inpaths(idat) + "/" + \
    ;       interface_get_infile(variables(ivar), field_types(ivar), idat)
    ;   end do
    ; end do
    ; contrib_authors = (/"A-anav_al", "A-wenz_sa", "A-righ_ma", "A-schl_ma"/)
    ; ESMValMD(plot@outfile, alltags, caption, id, variables, \
    ;          input_file_info@dataset, climofiles, diag_script,
    ;          contrib_authors)
    ; delete([/alltags, caption, id, climofiles/])
    delete(data_arr)
  end if

  leave_msg(diag_script, "")

end
