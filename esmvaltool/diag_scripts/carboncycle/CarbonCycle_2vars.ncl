;;#############################################################################
;; DIAGNOSTIC SCRIPT FOR SCATTER PLOT WITH TWO VARIABLES
;; Author: Alessandro Anav (University of Exeter, UK)
;; EMBRACE project
;;#############################################################################
;;
;; Description
;;    Draw a scatter plot with two variables.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    region: the region to be averaged
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;
;; Required variable_info attributes (variable specific)
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;
;; Modification history
;;    20180619-A_wenz_sa: adopted script to v2
;;    20170511-A_gier_be: added metadata to files
;;    20151130-A_laue_ax: corrected passing of file type from namelist to
;;                        plotting functions; moved call to write_references
;;                        to beginning of code
;;    20151105-A_righ_ma: graphic part moved to external plot_scripts.
;;    20151104_A_righ_ma: removing regridding, not necessary for this
;;                        diagnostic.
;;    201507??_A_wenz_sa: adapted to ESMValTool structure.
;;    201505??_A_anav_al: written.
;;
;;#############################################################################
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/scaling.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/watermark.ncl"

load "./diag_scripts/shared/plot/scatterplot.ncl"
load "./diag_scripts/shared/plot/style.ncl"

begin

    enter_msg(diag_script, "")
    log_info("++++++++++++++++++++++++++++++++++++++++++")
    log_info(diag_script + " (var: " + variables(0) + ", " + \
                variables(1) + ")")
    log_info("++++++++++++++++++++++++++++++++++++++++++")

    dim_MOD = dimsizes(dataset_info@dataset)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)

    write_references(diag_script,                  \  ; script name
                     (/"A_anav_al"/),              \  ; authors
                     (/"A_wenz_sa", "A_righ_ma"/), \  ; contributors
                     (/"D_anav13jclim"/),          \  ; diag_references
                     (/""/),                       \  ; obs_references
                     (/"P_embrace"/))                 ; proj_references

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = config_user_info@output_file_type
    if(ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Check watermark
    watermark = config_user_info@watermark
    if(ismissing(watermark)) then
        watermark = "False"
    end if
    print(watermark)

    ;; Check required diag_script_info attributes
    req_atts = (/"region"/)
    exit_if_missing_atts(diag_script_info, req_atts)

    ;; Check field type
    plot_fields = (/"T2Ms"/)
    if (field_type0.ne.plot_fields) then
        error_msg("f", diag_script, "", "can't use field " + field_type + \
                  " for this diag_script")
    end if
    delete(plot_fields)

    ;; Set region
    region = select_region(diag_script_info@region)

    ;; Output plot directory
    plot_dir = config_user_info@plot_dir
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; Determine time range
    start_year = min(stringtoint(dataset_info@start_year))
    end_year = max(stringtoint(dataset_info@end_year))
    all_years = ispan(start_year, end_year, 1)

    ;; GLobal array for collecting all data
    all_data_A = new((/dim_MOD/), float)
    all_data_A!0 = "dataset"
    all_data_A&dataset = dataset_info@dataset

    all_data_B = new((/dim_MOD/), float)
    all_data_B!0 = "dataset"
    all_data_B&dataset = dataset_info@dataset

    ;; Loop over datasets
    do imod = 0, dim_MOD - 1

        log_info("Processing " + all_data_A&dataset(imod))

        ;; First variable
        if (var_attr_exclude(0).ne.dataset_attr_id(imod)) then

            ;; Read data
            varA = read_data(imod, var0, field_type0)

            ;; Computing annual mean
            tmpA = time_operations(varA, toint(dataset_info@start_year(imod)), \
                                   toint(dataset_info@end_year(imod)), "average", \
                                   "annualclim", True)
            delete(varA)

            ;; Computing area-weighted sum
            avg_A = area_operations(tmpA, region(0), region(1), region(2), \
                                    region(3), "sum", True)
            delete(tmpA)

            ;; Assign to global array
            all_data_A(imod) = avg_A
            delete(avg_A)

        end if

        ;; Second variable
        if (var_attr_exclude(1).ne.dataset_attr_id(imod)) then

            ;; Read data
            varB = read_data(imod, var1, field_type1)

            ;; Computing annual mean
            tmpB = time_operations(varB, toint(dataset_info@start_year(imod)), \
                                   toint(dataset_info@end_year(imod)), "average", \
                                   "annualclim", True)
            delete(varB)

            ;; Computing area-weighted sum
            avg_B = area_operations(tmpB, region(0), region(1), region(2), \
                                    region(3), "sum", True)
            delete(tmpB)

            ;; Assign to global array
            all_data_B(imod) = avg_B
            delete(avg_B)

        end if

    end do  ; end loop over datasets

    ;; Convert units if appropriate
    if (isatt(variable_info, "plot_units")) then
        all_data_A = convert_units(all_data_A, variable_info@plot_units)
        all_data_B = convert_units(all_data_B, variable_info@plot_units)
    end if

    ;; Optional output to NetCDF
    if (config_user_info@write_netcdf.eq."True") then

        
        new_path = config_user_info@plot_dir + \
                   basename(diag_script) + "_" + var0 + ".nc"
        all_data_A@var = var0
        all_data_A@diag_script = diag_script
        all_data_A@ncdf = new_path
        ncdf_outfile = ncdf_write(all_data_A, new_path)

        new_path = config_user_info@plot_dir + \
                   basename(diag_script) + "_" + var1 + ".nc"
        all_data_B@var = var1
        all_data_B@diag_script = diag_script
        all_data_B@ncdf = new_path
        ncdf_outfile = ncdf_write(all_data_B, new_path)

    end if

    ;; Get reference dataset index
    idx_ref_A = get_ref_dataset_idx(dataset_info, variable_info@reference_dataset)
    idx_ref_B = get_ref_dataset_idx(dataset_info, variable_info@altreference_dataset)

    ;; Replace FillValue in the dataset array with the value of the reference
    ;; dataset. Since it is not possible to reduce the size of dataset_info@dataset, the
    ;; reference dataset is simply duplicated.
    all_data_A = \
        where(ismissing(all_data_A), all_data_A(idx_ref_A), all_data_A)
    all_data_B = \
        where(ismissing(all_data_B), all_data_B(idx_ref_B), all_data_B)

    ;; Collect data
    data_arr = new((/2, dim_MOD/), float)
    data_arr!0 = "statistic"
    data_arr!1 = "dataset"
    data_arr&statistic = (/var0, var1/)
    data_arr&dataset = dataset_info@dataset
    data_arr(0, :) = (/all_data_A/)
    data_arr(1, :) = (/all_data_B/)
    data_arr@units = (/all_data_A@units, all_data_B@units/)
    data_arr@legend_outside = diag_script_info@legend_outside

    ;; Draw plot
    outfile = output_dir + diag_script + "_" + var0 + "-" + var1 + \
        "_scatter_" + region@name
    wks = gsn_open_wks(file_type, outfile)
    wks@legendfile = outfile + "_legend"
    data_arr@res_tiMainString = var0 + "-" + var1 + " - " + \
        diag_script_info@region
    plot = scatterplot_markers(wks, data_arr, var0)
    plot@outfile = outfile + "." + file_type
    draw(plot)
    frame(wks)
    
    ;; Attach Metadata to plot    
    ;alltags = array_append_record(tags, region@DM_tag, 0)
    ;alltags := array_append_record(alltags,\
    ;       (/"PT_scatter", "ST_mean"/), 0)
    ;caption = "Scatterplot for " + variables(0) + " and " + \
    ;          variables(1) + " over the period " + start_year + \
    ;          "-" + end_year + ", similar to Anav et al., 2013, Fig. 12."
    ;id = diag_script
    ;climofiles = new(dimsizes(dataset_info@dataset) * dimsizes(variables), string)
    ;do n_var = 0, dimsizes(variables) - 1
    ;    do imod = 0, dimsizes(dataset_info@dataset) - 1
    ;        climofiles(n_var * dimsizes(dataset_info@dataset) + imod) = \
    ;                    interface_get_inpaths(imod) + "/" + \
    ;                    interface_get_infile(variables(n_var), field_types(n_var), imod)
    ;    end do
    ;end do
    ;contrib_authors = (/"A_anav_al", "A_wenz_sa", "A_righ_ma"/)

    ;;NOTE: not implemented in v2 yet
    ;ESMValMD(plot@outfile, alltags, caption, id, variables, dataset_info@dataset, climofiles, diag_script, contrib_authors)
    
    delete([/alltags, caption, id, climofiles/])

    

    leave_msg(diag_script, "")

end
