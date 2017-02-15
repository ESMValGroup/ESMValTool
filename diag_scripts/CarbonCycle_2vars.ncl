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
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"

load "./diag_scripts/lib/ncl/scaling.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"

load "./plot_scripts/ncl/scatterplot.ncl"

begin

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ", " + \
                variables(1) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    var1 = variables(1)
    field_type0 = field_types(0)
    field_type1 = field_types(1)

end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    write_references(diag_script,                  \  ; script name
                     (/"A_anav_al"/),              \  ; authors
                     (/"A_wenz_sa", "A_righ_ma"/), \  ; contributors
                     (/"D_anav13jclim"/),          \  ; diag_references
                     (/""/),                       \  ; obs_references
                     (/"P_embrace"/))                 ; proj_references

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "ps"
    end if

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
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; Determine time range
    start_year = min(stringtoint(models@start_year))
    end_year = max(stringtoint(models@end_year))
    all_years = ispan(start_year, end_year, 1)

    ;; GLobal array for collecting all data
    all_data_A = new((/dim_MOD/), float)
    all_data_A!0 = "model"
    all_data_A&model = models@name

    all_data_B = new((/dim_MOD/), float)
    all_data_B!0 = "model"
    all_data_B&model = models@name

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        info_output("Processing " + all_data_A&model(imod), verbosity, 1)

        ;; First variable
        if (var_attr_exclude(0).ne.model_attr_id(imod)) then

            ;; Read data
            dataA = read_data(imod, var0, field_type0)
            varA = extract_data_extended(imod, var0, dataA, -1, 0, 0)
            delete(dataA)

            ;; Computing annual mean
            tmpA = time_operations(varA, toint(models@start_year(imod)), \
                                   toint(models@end_year(imod)), "average", \
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
        if (var_attr_exclude(1).ne.model_attr_id(imod)) then

            ;; Read data
            dataB = read_data(imod, var1, field_type1)
            varB = extract_data_extended(imod, var1, dataB, -1, 0, 0)
            delete(dataB)

            ;; Computing annual mean
            tmpB = time_operations(varB, toint(models@start_year(imod)), \
                                   toint(models@end_year(imod)), "average", \
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

    end do  ; end loop over models

    ;; Convert units if appropriate
    if (isatt(variable_info, "plot_units")) then
        all_data_A = convert_units(all_data_A, variable_info@plot_units)
        all_data_B = convert_units(all_data_B, variable_info@plot_units)
    end if

    ;; Optional output to NetCDF
    if (getenv("ESMValTool_write_netcdf").eq."True") then

        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + \
            basename(diag_script) + "_" + var0 + ".nc"
        all_data_A@var = var0
        all_data_A@diag_script = diag_script
        all_data_A@ncdf = new_path
        ncdf_outfile = ncdf_write(all_data_A, new_path)

        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + \
            basename(diag_script) + "_" + var1 + ".nc"
        all_data_B@var = var1
        all_data_B@diag_script = diag_script
        all_data_B@ncdf = new_path
        ncdf_outfile = ncdf_write(all_data_B, new_path)

    end if

    ;; Get reference model index
    idx_ref_A = get_ref_model_idx(models, var_attr_ref(0))
    idx_ref_B = get_ref_model_idx(models, var_attr_ref(1))

    ;; Replace FillValue in the model array with the value of the reference
    ;; model. Since it is not possible to reduce the size of models@name, the
    ;; reference model is simply duplicated.
    all_data_A = \
        where(ismissing(all_data_A), all_data_A(idx_ref_A), all_data_A)
    all_data_B = \
        where(ismissing(all_data_B), all_data_B(idx_ref_B), all_data_B)

    ;; Collect data
    data_arr = new((/2, dim_MOD/), float)
    data_arr!0 = "statistic"
    data_arr!1 = "model"
    data_arr&statistic = (/var0, var1/)
    data_arr&model = models@name
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
    draw(plot)
    frame(wks)

    leave_msg(diag_script, "", 4)

end
