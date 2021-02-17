; #############################################################################
; carbon_cycle/two_variables.ncl
; #############################################################################
;
; Description:
;     Draw a scatter plot with two variables.
;
; Required diag_script_info attributes:
;     region: the region to be averaged.
;
; Modification history:
;     20191209-schlund_manuel: added new provenance tracking.
;     20180823-schlund_manuel: continued porting to v2.0.
;     20180619-wenzel_sabrina: adopted script to v2.
;     20170511-gier_bettina: added metadata to files.
;     20151130-lauer_axel: corrected passing of file type from namelist to
;                          plotting functions; moved call to write_references
;                          to beginning of code.
;     20151105-righi_mattia: graphic part moved to external plot_scripts.
;     20151104-righi_mattia: removing regridding, not necessary for this
;                            diagnostic.
;     201507??-wenzel_sabrina: adapted to ESMValTool structure.
;     201505??-anav_alessandro: written.
;
; #############################################################################

load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/latlon.ncl"
load "$diag_scripts/shared/scaling.ncl"
load "$diag_scripts/shared/set_operators.ncl"
load "$diag_scripts/shared/statistics.ncl"

load "$diag_scripts/shared/plot/scatterplot.ncl"
load "$diag_scripts/shared/plot/style.ncl"


begin

  enter_msg(DIAG_SCRIPT, "")
  DOMAIN = diag_script_info@region

  ; Variable
  DIM_VAR = ListCount(variable_info)
  if (DIM_VAR .ne. 2) then
    error_msg("f", DIAG_SCRIPT, "", "this diagnostic expects exactly two " + \
              "variables")
  end if
  var = new(DIM_VAR, string)
  var(0) = variable_info[0]@short_name
  var(1) = variable_info[1]@short_name

  ; Input data
  info0 = select_metadata_by_name(input_file_info, var(0))
  info1 = select_metadata_by_name(input_file_info, var(1))
  projects0 = metadata_att_as_array(info0, "project")
  projects1 = metadata_att_as_array(info1, "project")
  DATASETS = NewList("lifo")
  ListAppend(DATASETS, metadata_att_as_array(info0, "dataset"))
  ListAppend(DATASETS, metadata_att_as_array(info1, "dataset"))
  if any(DATASETS[0](ind(projects0 .ne. "OBS")) .ne. \
         DATASETS[1](ind(projects1 .ne. "OBS"))) then
    error_msg("f", DIAG_SCRIPT, "", "the two variables do not cover the " + \
              "same datasets, this is not supported yet")
  end if
  if any(projects0 .ne. projects1) then
    error_msg("f", DIAG_SCRIPT, "", "the two variables do not include the " + \
              "same number of observations")
  end if
  DIM_DAT = dimsizes(DATASETS[0])

  log_info("++++++++++++++++++++++++++++++++++++++++++")
  log_info(DIAG_SCRIPT + " (vars: " + var(0) + ", " + var(1) + ")")
  log_info("++++++++++++++++++++++++++++++++++++++++++")

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
  if (.not. isatt(diag_script_info, "legend_outside")) then
    diag_script_info@legend_outside = False
  end if

  ; Set region
  region = select_region(DOMAIN)

  ; ---------------------------------------------------------------------------
  ; Read and preprocess data
  ; ---------------------------------------------------------------------------

  ; Loop over all variables and datasets
  all_data = NewList("lifo")
  do ivar = 0, DIM_VAR - 1
    INFO = select_metadata_by_name(input_file_info, var(ivar))

    ; Array
    ListAppend(all_data, new((/DIM_DAT/), float))
    all_data[ivar]!0 = "dataset"
    all_data[ivar]&dataset = DATASETS[ivar]

    ; Time range
    start_year = 9999
    end_year = -9999
    do idat = 0, DIM_DAT - 1
      if (INFO[idat]@start_year .lt. start_year) then
        start_year = INFO[idat]@start_year
      end if
      if (INFO[idat]@end_year .gt. end_year) then
        end_year = INFO[idat]@end_year
      end if
    end do

    ; Read data
    do idat = 0, DIM_DAT - 1
      data = INFO[idat]
      log_debug("Processing " + var(ivar) + " of " + DATASETS[ivar](idat))

      ; Read variable
      var_data = read_data(data)

      ; Compute annual mean
      tmp = time_operations(var_data, data@start_year, data@end_year, \
                            "average", "annualclim", True)

      ; Compute area-weighted sum
      area_op = area_operations(tmp, region(0), region(1), region(2), \
                                region(3), "sum", True)

      ; Assign to global array
      all_data[ivar](idat) = area_op
      delete([/var_data, tmp, area_op/])
    end do

    ; Convert units if appropriate
    if ((isatt(variable_info[ivar], "plot_units")) .and. \
        (all_data[ivar]@units .ne. variable_info[ivar]@plot_units)) then
      all_data[ivar] = convert_units(all_data[ivar], \
                                     variable_info[ivar]@plot_units)
    end if
  end do

  ; ---------------------------------------------------------------------------
  ; Write NETCDF
  ; ---------------------------------------------------------------------------

  work_dir = config_user_info@work_dir
  system("mkdir -p " + work_dir)

  ; Save individual variables
  do ivar = 0, DIM_VAR - 1
    new_path = work_dir + var(ivar) + "_" + DOMAIN + ".nc"
    all_data[ivar]@var = var(ivar)
    all_data[ivar]@diag_script = DIAG_SCRIPT
    all_data[ivar]@ncdf = new_path
    ncdf_outfile = ncdf_write(all_data[ivar], new_path)
    delete(ncdf_outfile)
  end do

  ; Save both variables
  new_path = work_dir + var(0) + "_" + var(1) + "_" + DOMAIN + ".nc"
  both_vars = new((/DIM_VAR, DIM_DAT/), typeof(all_data[0]))
  do ivar = 0, DIM_VAR - 1
    both_vars(ivar, :) = all_data[ivar](:)
  end do
  both_vars!0 = "variable"
  both_vars!1 = "dataset"
  both_vars&variable = var
  both_vars&dataset = DATASETS[0]
  both_vars&dataset(ind(projects0 .eq. "OBS")) = "OBS"
  both_vars@var = "scatter"
  both_vars@diag_script = DIAG_SCRIPT
  both_vars@ncdf = new_path
  ncdf_outfile_both_vars = ncdf_write(both_vars, new_path)

  ; ---------------------------------------------------------------------------
  ; Plots
  ; ---------------------------------------------------------------------------

  if (diag_script_info@write_plots) then
    data_arr = new((/DIM_VAR, DIM_DAT(0)/), float)

    ; Collect data
    do ivar = 0, DIM_VAR - 1
      ref_idx = ind(DATASETS[ivar] .eq. variable_info[ivar]@reference_dataset)
      ref_idx := ref_idx(0)
      if (ismissing(ref_idx)) then
          error_msg("f", DIAG_SCRIPT, "", "no adequate reference dataset " + \
                    "for " + var(ivar) + " provided")
      end if
      log_info("Reference dataset for " + var(ivar) + ": " + \
               DATASETS[ivar](ref_idx))

      ; Replace FillValue in the dataset array with the value of the reference
      all_data[ivar] = \
        where(ismissing(all_data[ivar]), all_data[ivar](ref_idx), \
              all_data[ivar])
      data_arr(ivar, :) = (/all_data[ivar]/)
    end do

    ; Setup metadata of output array
    data_arr!0 = "statistic"
    data_arr!1 = "dataset"
    data_arr&statistic = (/var(0), var(1)/)
    data_arr&dataset = DATASETS[0]
    data_arr&dataset(ind(projects0 .eq. "OBS")) = "OBS"
    data_arr@units = (/all_data[0]@units, all_data[1]@units/)
    data_arr@legend_outside = diag_script_info@legend_outside
    data_arr@annots = data_arr&dataset

    ; Draw plot
    outfile = plot_dir + var(0) + "-" + var(1) + "_scatter_" + DOMAIN
    wks = gsn_open_wks(file_type, outfile)
    wks@legendfile = outfile + "_legend"
    data_arr@res_tiMainString = var(0) + "-" + var(1) + " - " + region@name
    plot = scatterplot_markers(wks, data_arr, var(0), info0)
    plot_path = outfile + "." + file_type
    plot@outfile = plot_path
    log_info("Wrote " + plot_path)
    draw(plot)
    frame(wks)

    ; Provenance tracking
    caption = ("Scatterplot for " + var(0) + " and " + var(1) + \
               " over the period " + start_year + "-" + end_year + \
               ". Similar to Anav et al. (2013), Figure 12.")
    statistics = (/"mean"/)
    authors = (/"anav_alessandro", "wenzel_sabrina", "righi_mattia", \
                "schlund_manuel"/)
    plot_type = "scatter"
    references = (/"anav13jclim"/)
    ancestors = metadata_att_as_array(info0, "filename")
    ancestors := array_append_record( \
      ancestors, metadata_att_as_array(info1, "filename"), 0)
    log_provenance(ncdf_outfile_both_vars, plot_path, caption, statistics, \
                   DOMAIN, plot_type, authors, references, ancestors)
    delete(data_arr)
  end if

  leave_msg(DIAG_SCRIPT, "")

end
