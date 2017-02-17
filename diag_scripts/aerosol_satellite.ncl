;;#############################################################################
;; DIAGNOSTIC SCRIPT FOR AEROSOL SATELLITE DATA
;; Author: Mattia Righi (DLR, Germany)
;; ESMVal project
;;#############################################################################
;;
;; Description
;;    Plots optical depth and compares with satelite data.
;;    Considers monthly mean aerosol optical depth at 550 nm.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    ref_model: name of the reference model.
;;    regrid_method: grid selection for regridding.
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    styleset: project-specific plot style.
;;    range_option: time selection option.
;;
;; Required variable_info attributes (variable specific)
;;    contours_abs: array of contour levels for the absolute value map.
;;    contours_diff: array of contour levels for the difference map.
;;    colortable_abs: color table for the absolute value map.
;;    colortable_diff: color table for the difference value map.
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;
;; Modification history
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150325-A_laue_ax: modified reference tags used for acknowledgements
;;                        (projects, observations, etc)
;;    20140212: written by Mattia Righi.
;;
;;#############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/ensemble.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"

load "./plot_scripts/ncl/contour_maps.ncl"

begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, \
                verbosity, 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; Reference information
    write_references(diag_script,\
                     "A_righ_ma", \  ; authors
                     "A_kais_ch", \  ; contributors
                     "D_0000",    \  ; diag_ref
                     (/"E_modis", "E_misr"/), \  ; obs_ref
                     (/"P_esmval", "P_dlrveu"/) \  ; proj_acknow
                     )

    ;; Check required diag_script_info attributes
    exit_if_missing_atts(diag_script_info, (/"ref_model", "regrid_method"/))

    ;; Check valid field types
    field_type0 = field_types(0)
    if (field_type0 .ne. "T2Ms") then
        error_msg("f", diag_script, "", "type " + field_type0 + \
                  " not supported")
    end if

    ;; Read data
    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do
end

;; Load var_info
begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

;; Main
begin

    ;; File type
    valid_types = (/"PS", "ps", "PDF", "pdf"/)  ; Support for multiple pages
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "PS"
    end if
    if (dimsizes(models@name) .gt. 1.and.all(file_type .ne. valid_types)) then
        info_output(file_type + " format does not support multiple pages, " \
                    + "switching to ps", verbosity, 1)
        file_type = "PS"
    end if

    ;; Check for reference model
    refID = ind(models@name.eq.diag_script_info@ref_model)
    if (ismissing(refID)) then
        error_msg("f", diag_script, "", "cannot determine reference model")
    end if

    ;; Determine time range (years) of all models
    if (.not.isatt(diag_script_info, "range_option")) then
        diag_script_info@range_option = 0  ; every model use its own range
    end if

    ;; Find the reference grid for regridding
    dest_grid = find_destination_grid(ispan(0, dimsizes(models@name) - 1, 1), \
                                      var0, field_type0, \
                                      diag_script_info@regrid_method)

    ;; Define output array
    output = new(array_append_record(dimsizes(models@name), \
                                     dimsizes(dest_grid), 0), float)
    output!0 = "model"
    output&model = unique_labels_min((/"name", "case_name", \
                                     "ensemble", "experiment"/))
    output@var = var0
    output@diag_script = (/diag_script/)

    ;; Loop over model, process data
    do mID = 0, dimsizes(models@name) - 1

        ;; Check time period
        tmp = extract_data(mID, data, - 1, 0, 0)
        start_year = get_start_year(mID)
        end_year = get_end_year(mID)

        if (start_year.gt.end_year) then
            info_output("No overlapping time period between the selected " + \
                        "models: specify a different time period or " + \
                        "select a different range_option", \
                        verbosity, 1)
            delete(tmp)
            continue
        end if

        ;; Calculate time average
        tmp_avg = time_operations(tmp, start_year, end_year, "average", \
                                  "annualclim", True)
        delete(tmp)

        ;; Regrid
        tmp_regr = linint2_Wrap(tmp_avg&lon, tmp_avg&lat, tmp_avg, True, \
                                dest_grid&lon, dest_grid&lat, 0)
        delete(tmp_avg)

        ;; Assign to output array
        output(mID, :, :) = tmp_regr
        delete(tmp_regr)

    end do

    ;; Check data availability
    if (all(ismissing(output))) then
        info_output("Observational data not available in the selected " + \
                    "time period, skipping", verbosity, 1)
        exit
    end if

    ;; Optional output to NetCDF
    if (getenv("ESMValTool_write_netcdf").eq."True") then
        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + basename(diag_script) + "_" + \
            diag_script_info@ref_model + "_" + var0 + ".nc"
        output@ncdf = new_path
        ncdf_outfile = ncdf_write(output, new_path)
    end if

    ;; Plot
    annots = project_style(diag_script_info, "annots")
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, basename(diag_script))

    plot_abs = output_dir + diag_script_info@ref_model + "_absolute_" + var0
    wksabs = gsn_open_wks(file_type, plot_abs)
    gsn_define_colormap(wksabs, variable_info@colortable_abs)
    plot_diff = output_dir + diag_script_info@ref_model + "_difference_" + var0
    wksdiff = gsn_open_wks(file_type, plot_diff)
    gsn_define_colormap(wksdiff, variable_info@colortable_diff)

    do mID = 0, dimsizes(models@name) - 1

        var = output(mID, :, :)
        var@res_cnFillMode = "CellFill"
        var@res_cnLinesOn = False
        var@res_cnLevelSelectionMode = "ExplicitLevels"
        var@res_tiMainString = var0 + " (" + annots(mID) + ")"
        var@res_cnLevels = variable_info@contours_abs
        plot = contour_map_ce(wksabs, var, var0)
        delete(var@res_cnLevels)
        draw(plot)
        frame(wksabs)

        if (mID .ne. refID) then
            var = output(mID, :, :) - output(refID, :, :)
            var@res_cnLevels = variable_info@contours_diff
            var@res_tiMainString = var0 + " (" + annots(mID) + " - " + \
                annots(refID) + ")"
            plot = contour_map_ce(wksdiff, var, var0)
            delete(var@res_cnLevels)
            draw(plot)
            frame(wksdiff)
        end if

    end do
    info_output("gv " + plot_abs + "." + file_type, verbosity, 1)
    info_output("gv " + plot_diff + "." + file_type, verbosity, 1)

    info_output(">>>>>>>> Leaving " + diag_script, \
                verbosity, 2)

end
