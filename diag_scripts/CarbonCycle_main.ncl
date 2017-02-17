;;#############################################################################
;; DIAGNOSTIC SCRIPT FOR ERRORBAR AND SEASONAL CYCLE PLOTS
;; Author: Alessandro Anav (University of Exeter, UK)
;; EMBRACE project
;;#############################################################################
;;
;; Description
;;    Calculates temporal and spatial averages and plots the variable as
;;    error-bar and seasonal cycle plots
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    region: the region to be averaged
;;    legend_outside: plot legend in a separate page (not for errorbar plot)
;;    seasonal_cycle_plot: draw seasonal cycle plot [False/True]
;;    errorbar_plot: draw errorbar plot [False/True]
;;    mean_IAV_plot: draw Mean (x-axis), IAV (y-axis) plot [False/True]
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    sort: sort model in alphabetical order.
;;
;; Required variable_info attributes (variable specific)
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;    Regridding is not required, since spatial averages are calculated for
;;    all plots. Performance metrics calculated with perfmetrics_main could
;;    be inconsistent with the plots produced by this script, since a
;;    consistent missing value mask is not applied here. However, for variable
;;    such NBP, GPP, LAI, FGCO2 for which this script is designed, there
;;    should be no inconsistencies, as the missing value mask is similar for
;;    all models and observations and coincides with the land-sea mask.
;;
;; Modification history
;;    20151130-A_laue_ax: corrected passing of file type from namelist to
;;                        plotting functions; moved call to write_references
;;                        to beginning of code
;;    20151105-A_righ_ma: graphic part moved to external plot_scripts.
;;    20151030-A_righ_ma: regridding not working on ocean grids, removed since
;;                        it is not necessary for these diagnostics.
;;    20151029-A_righ_ma: replacing area functions with ESMValTool functions.
;;    20151021-A_righ_ma: pep8 and header standardization.
;;    201510??-A_wenz_sa: adapted to ESMValTool structure.
;;    201505??-A_anav_al: written.
;;
;;#############################################################################

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./interface_scripts/constants.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"

load "./plot_scripts/ncl/xy_line.ncl"
load "./plot_scripts/ncl/scatterplot.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Read model data
    data = True
    do idx_mod = 0, dim_MOD - 1
        data = read_data(idx_mod, var0, field_type0)
    end do

end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    write_references(diag_script,                  \
                     (/"A_anav_al"/),              \  ; authors
                     (/"A_wenz_sa", "A_righ_ma"/), \  ; contributors
                     (/"D_anav13jclim"/),          \  ; diag_references
                     (/""/),                       \  ; obs_references
                     (/"P_embrace"/))                 ; proj_references

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"region", "legend_outside", "seasonal_cycle_plot", \
                 "errorbar_plot", "mean_IAV_plot"/)
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

    ;; Set time range
    start_year = min(stringtoint(models@start_year))
    end_year = max(stringtoint(models@end_year))

    ;; Global array for collecting all data
    all_yearly = new((/dim_MOD, end_year - start_year + 1, 2/), float)
    all_monthly = new((/dim_MOD, 12, 2/), float)

    all_yearly!0          = "model"
    all_yearly!1          = "year"
    all_yearly!2          = "statistic"
    all_yearly&model      = models@name
    all_yearly&year       = ispan(start_year, end_year, 1)
    all_yearly&statistic  = (/"mean", "stddev"/)

    all_monthly!0         = "model"
    all_monthly!1         = "month"
    all_monthly!2         = "statistic"
    all_monthly&model     = models@name
    all_monthly&month     = ispan(1, 12, 1)
    all_monthly&statistic = (/"mean", "stddev"/)

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        info_output("Processing " + models@name(imod), verbosity, 0)

        ;; Extract model corresponding data
        var  = extract_data_extended(imod, var0, data, -1, 0, 0)

        ;; Set type of area operation (flux vs. non-flux variables)
        if (var@units.eq."kg m-2 s-1") then
            sp_opt = "sum"  ; for flux variable
        else
            sp_opt = "average"  ; for non-flux variable
        end if

        ;; Calculate annual mean
        data_yearly = \
            time_operations(var, toint(models@start_year(imod)), \
                            toint(models@end_year(imod)), \
                            "average", "yearly", True)

        ;; Calculate mean annual cycle
        data_monthly = \
            time_operations(var, toint(models@start_year(imod)), \
                            toint(models@end_year(imod)), "average", \
                            "monthlyclim", True)

        ;; Special case for LAI.
        ;; Since models have different land covers, all sparse vegetated points
        ;; must be masked out
        if (var0.eq."LAI" .or. var0.eq."lai") then
            data_yearly = where(data_yearly.ge.0.025, data_yearly, \
                                data_yearly@_FillValue)
            data_monthly = where(data_monthly.ge.0.025, data_monthly, \
                                 data_monthly@_FillValue)
        end if

        ;; Yearly data - Mean
        ymean  = area_operations(data_yearly, region(0), region(1), \
                                 region(2), region(3), sp_opt, True)

        ;; Yearly data - Standard deviation
;; FIX-ME: results to be checked
;;         area = map_area(var&lat, var&lon)
;;         tmp = data_yearly
;;         tmp = data_yearly * conform(data_yearly, area, (/1, 2/))
;;         subr = area_operations(tmp, region(0), region(1), \
;;                                region(2), region(3), "extract", True)
;;         ystd = dim_stddev_n_Wrap(subr, (/1, 2/))
;;         delete(subr)
;;         delete(tmp)

        ;; Assign to global array
        idx1 = ind(min(ymean&year).eq.all_yearly&year)
        idx2 = ind(max(ymean&year).eq.all_yearly&year)
        all_yearly(imod, idx1:idx2, 0) = (/ymean/)
        all_yearly(imod, idx1:idx2, 1) = 0.  ; (/ystd/)
        copy_VarAtts(ymean, all_yearly)
        delete(ymean)
;;        delete(ystd)
        delete(idx1)
        delete(idx2)

        ;; Monthly data - Mean
        mmean = area_operations(data_monthly, region(0), region(1), \
                                region(2), region(3), sp_opt, True)

        ;; Monthly data - Standard deviation
;; FIX-ME: results to be checked
;;         tmp = data_monthly
;;         tmp = data_monthly * conform(data_monthly, area, (/1, 2/))
;;         subr = area_operations(tmp, region(0), region(1), \
;;                                region(2), region(3), "extract", True)
;;         mstd = dim_stddev_n_Wrap(subr, (/1, 2/))
;;         delete(subr)
;;         delete(tmp)

        ;; Assign to global array
        all_monthly(imod, :, 0) = (/mmean/)
        all_monthly(imod, :, 1) = 0.  ; (/mstd/)
        copy_VarAtts(mmean, all_monthly)
        delete(mmean)
;;        delete(mstd)
        delete(var)
        delete(data_yearly)
        delete(data_monthly)
;;        delete(area)

    end do

    ;; Convert units if appropriate [kg m-2 s-1] --> [PgC y-1]
    if (isatt(variable_info, "plot_units")) then
        all_yearly = convert_units(all_yearly, variable_info@plot_units)
        all_monthly = convert_units(all_monthly, variable_info@plot_units)
    end if

    ;; Optional output to NetCDF
    if (getenv("ESMValTool_write_netcdf").eq."True") then

        ;; Yearly data
        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + \
            basename(diag_script) + "_yearly_" + var0 + ".nc"
        all_yearly@var = var0
        all_yearly@diag_script = diag_script
        all_yearly@ncdf = new_path
        ncdf_outfile = ncdf_write(all_yearly, new_path)

        ;; Monthly data
        new_path = get_output_dir(get_work_dir(), basename(diag_script))
        new_path = new_path + \
            basename(diag_script) + "_monthly_" + var0 + ".nc"
        all_monthly@var = var0
        all_monthly@diag_script = diag_script
        all_monthly@ncdf = new_path
        ncdf_outfile = ncdf_write(all_monthly, new_path)

    end if

    ;; Seasonal cycle plot
    if (diag_script_info@seasonal_cycle_plot) then

        ;; Collect data
        data_arr = all_monthly
        data_arr@legend_outside = diag_script_info@legend_outside

        ;; Draw plot
        outfile = output_dir + diag_script + "_" + var0 + "_cycle_" + \
            region@name
        wks = gsn_open_wks(file_type, outfile)
        wks@legendfile = outfile + "_legend"
        data_arr@res_tiYAxisString = var0 + " [" + \
            format_units(data_arr@units) + "]"
        data_arr@res_tiMainString = diag_script_info@region
        plot = cycle_plot(wks, data_arr, var0)
        draw(plot)
        frame(wks)
        delete(plot)
        delete(wks)
        delete(data_arr)

    end if

    ;; Error-bar plot
    if (diag_script_info@errorbar_plot) then

        annual_mean = dim_avg_Wrap(all_yearly(:, :, 0))
        annual_std = dim_stddev_Wrap(all_yearly(:, :, 0))

        ;; Sort model in alphabetical order (observations at the end)
        if (isatt(diag_script_info, "sort")) then
            if (diag_script_info@sort) then
                pid = sort_alphabetically(annual_mean&model, \
                                          ind(models@project.eq."OBS"), \
                                          "end")
                annual_mean_s = annual_mean(pid)
                annual_mean_s&model = annual_mean&model(pid)
                annual_mean = annual_mean_s
                annual_std_s = annual_std(pid)
                annual_std_s&model = annual_std&model(pid)
                annual_std = annual_std_s
                delete([/annual_mean_s, annual_std_s, pid/])
            end if
        end if

        ;; Collect data
        data_arr = new((/2, dim_MOD/), float)
        data_arr!0 = "statistic"
        data_arr!1 = "model"
        data_arr&statistic = (/"mean", "standard deviation"/)
        data_arr&model = annual_mean&model
        data_arr(0, :) = (/annual_mean/)
        data_arr(1, :) = (/annual_std/)
        data_arr@units = annual_mean@units
        delete(annual_mean)
        delete(annual_std)

        ;; Draw plot
        outfile = output_dir + diag_script + "_" + var0 + "_errorbar_" + \
            region@name
        wks = gsn_open_wks(file_type, outfile)
        data_arr@res_tiMainString = var0 + " - " + diag_script_info@region
        plot = errorbar_plot(wks, data_arr, var0)
        draw(plot)
        frame(wks)
        delete(plot)
        delete(wks)
        delete(data_arr)

    end if

    ;; Mean, IAV plot
    if (diag_script_info@mean_IAV_plot) then

        ;; Prepare data
        annual_mean            = new((/dim_MOD/), float)
        annual_std             = new((/dim_MOD/), float)
        annual_mean@_FillValue = all_yearly@_FillValue
        annual_std@_FillValue  = all_yearly@_FillValue
        annual_mean@units      = all_yearly@units
        annual_std@units       = all_yearly@units

        ;; Loop over models
        do imod = 0, dim_MOD - 1
            ;; Compute the spatial mean over time
            annual_mean(imod) = dim_avg(all_yearly(imod, :, 0))
            ;; Compute the std over time for the region
            annual_std(imod) = dim_stddev(all_yearly(imod, :, 0))
        end do

        ;; Collect data
        data_arr = new((/2, dim_MOD/), float)
        data_arr!0 = "statistic"
        data_arr!1 = "model"
        data_arr&statistic = (/"mean", "interannual variability"/)
        data_arr&model = models@name
        data_arr(0, :) = (/annual_mean/)
        data_arr(1, :) = (/annual_std/)
        data_arr@units = (/annual_mean@units, annual_std@units/)
        data_arr@legend_outside = diag_script_info@legend_outside
        delete(annual_mean)
        delete(annual_std)

        ;; Draw plot
        outfile = output_dir + diag_script + "_" + var0 + "_scatter_" + \
            region@name
        wks = gsn_open_wks(file_type, outfile)
        wks@legendfile = outfile + "_legend"
        data_arr@res_tiMainString = var0 + " - " + diag_script_info@region
        plot = scatterplot_markers(wks, data_arr, var0)
        draw(plot)
        frame(wks)

    end if

    leave_msg(diag_script, "", 2)

end
