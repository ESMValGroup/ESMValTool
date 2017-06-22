;;#############################################################################
;; DIAGNOSTIC SCRIPT for reproducing IPCC ch. 9 fig. 9.36
;; Author: Axel Lauer (DLR, Germany)
;; CRESCENDO project
;;#############################################################################
;;
;; Description
;;    Calculates (anomaly) standard deviation averaged over Nino 3 and
;;    Nino 4 regions.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;     none
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;     none
;;
;; Required variable_info attributes (variable specific)
;;     none
;;
;; Optional variable_info attributes (variable specific)
;;     none
;;
;; Required variable attributes (defined in namelist)
;;    ref_model: name of reference data set (observations)
;;
;; Caveats
;;     none
;;
;; Modification history
;;    20170622-A_laue_ax: added add tags for reporting
;;    20170307-A_laue_ax: written.
;;
;;#############################################################################

load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"

load "./plot_scripts/ncl/xy_line.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ref_model = var_attr_ref(0)  ; produces an error if ref_model is not set in the namelist
    info_output("reference model = " + ref_model, verbosity, 1)
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Read model data
    data = True
    do imod = 0, dim_MOD - 1
        data = read_data(imod, var0, field_type0)
    end do

    climofiles = new(dim_MOD, string)

    do ii = 0, dim_MOD - 1
        climofiles(ii) = \
            interface_get_inpaths(ii) + "/" \
            + interface_get_infile(variables(0), field_types(0), ii)
    end do

    alltags = array_append_record(tags, (/"PT_scatter", "ST_anomaly", "ST_stddev", "DM_reg"/), 0)
    allmodelnames = array_append_record(models@name, (/"multi-model-mean"/), 0)
end

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    write_references(diag_script,                  \
                     (/"A_laue_ax"/),              \  ; authors
                     (/""/),                       \  ; contributors
                     (/"D_kim12grl"/),             \  ; diag_references
                     (/""/),                       \  ; obs_references
                     (/"P_crescendo"/))               ; proj_references

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Check field type
    plot_fields = (/"T2Ms"/)
    if (field_type0.ne.plot_fields) then
        error_msg("f", diag_script, "", "can't use field " + field_type + \
                  " for this diag_script")
    end if
    delete(plot_fields)

    ;; define regions
    region = (/"Nino3", "Nino4"/)
    nreg = dimsizes(region)
    latmin = (/-5.0, -5.0/)
    latmax = (/5.0, 5.0/)
    lonmin = (/210.0, 160.0/)
    lonmax = (/270.0, 210.0/)

    ;; Output plot directory
    plot_dir = getenv("ESMValTool_plot_dir")
    output_dir = get_output_dir(plot_dir, diag_script_base)

    ;; Set time range
    start_year = min(stringtoint(models@start_year))
    end_year = max(stringtoint(models@end_year))
    tsteps = (end_year - start_year + 1) * 12

;    ;; arrays for calculating multi model mean
;    mmm = new((/tsteps, nreg/), float)
;    pts = new((/tsteps, nreg/), integer)
;
;    mmm = 0.0
;    pts = 0

    timeseries = new((/tsteps, nreg, dim_MOD/), float)
    timeseries = timeseries@_FillValue

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        info_output("Processing " + models@name(imod), verbosity, 0)

        yr0 = stringtoint(models@start_year(imod))
        yr1 = stringtoint(models@end_year(imod))

        ;; Extract model corresponding data
        var  = extract_data_extended(imod, var0, data, -1, 0, 0)

        units = var@units

        do ireg = 0, nreg - 1

            ; calculate domain average
            regmean = area_operations(var, latmin(ireg), latmax(ireg), \
                                 lonmin(ireg), lonmax(ireg), "average", True)

            ; remove trend

            dimname = getvardims(regmean)
            regmean_dtrend = dtrend_msg(regmean&$dimname$, regmean, False, False)

            regmean_dtrend!0 = regmean!0
            regmean_dtrend&$dimname$ = regmean&$dimname$

            ; calculate climatological seasonal cycle
            seasclim = time_operations(regmean_dtrend, -1, -1, "average", "monthlyclim", True)

            ; subtract climatological seasonal cycle
            dseas = regmean_dtrend
            do i = 0, 11
                dseas(i::12) = dseas(i::12) - seasclim(i)
            end do

            i0 = (yr0 - start_year) * 12
            i1 = i0 + (yr1 - yr0) * 12 + 11

;            ; calculate multi-model mean
;
;            if (imod.ne.idx_ref) then
;                if (models@project(imod).eq."OBS") then
;                    info_output("+++ WARNING: including " + models@name(imod) + \
;                                " in multi-model mean calculation +++", verbosity, 0)
;                end if
;                mmm(i0:i1, ireg) = mmm(i0:i1, ireg) + dseas(:)
;                pts(i0:i1, ireg) = pts(i0:i1, ireg) + 1
;            end if

            ; save results
            timeseries(i0:i1, ireg, imod) = dseas(:)

            delete(regmean)
            delete(regmean_dtrend)
            delete(seasclim)
            delete(dseas)

        end do  ; ireg-loop

        delete(var)

    end do  ; imod-loop

;    ; finish multi-model mean calculation
;
;    pts = where(pts.eq.0, pts@_FillValue, pts)
;    mmm = mmm / pts

    ; ===================================================================================================

    ; calculate standard deviation of deseasonalized, detrended time series

    reg_std = dim_stddev_n(timeseries, 0)

    do ireg = 0, nreg - 1

        ; sort models in alphabetical order (observations at the end)

        pid = sort_alphabetically(models@name, idx_ref, "begin")
        annual_mean = reg_std(ireg, pid)
        annual_mean!0 = "model"

        annual_mean&model = models@name(pid)
        annual_mean&model(0) = "multi-model-mean"
        ; from IPCC ch. 9 fig. 9.36 it seems that "multi-model mean" is *not* the temporal
        ; standard deviation of the multi-model mean time series (mmm(:, ireg))
        ; but rather simply the average of the temporal standard deviations
        ; from the individual models (annual_mean(1:dim_MOD - 1))
        ;annual_mean(0) = stddev(mmm(:, ireg))
        annual_mean(0) = avg(annual_mean(1:dim_MOD - 1))
        annual_std = annual_mean * 0.0  ; create array of same size as annual_mean
        annual_std(0) = stddev(annual_mean(1:dim_MOD - 1))
        ; errorbar_plot draws "mean +- std/2"
        ; --> multiply std by 2 so error bar shows +-std instead of +-std/2
        annual_std(0) = annual_std(0) * 2.0
        delete(pid)

        ;; Collect data
        data_arr = new((/2, dim_MOD/), float)
        data_arr!0 = "statistic"
        data_arr!1 = "model"
        data_arr&statistic = (/"mean", "standard deviation"/)
        data_arr&model = annual_mean&model
        data_arr(0, :) = (/annual_mean/)
        data_arr(1, :) = (/annual_std/)
        data_arr@units = units
        delete(annual_mean)
        delete(annual_std)

        x = 0.0
        y = reg_std(ireg, idx_ref)

        ;; Draw plot
        outfile = output_dir + diag_script_base + "_" + var0 + "_" + region(ireg)
        wks = gsn_open_wks(file_type, outfile)
        data_arr@res_tiMainString = var0 + " - " + region(ireg)
        diag_script_info@region = region(ireg)
        data_arr@yrefline = y
        plot = errorbar_plot(wks, data_arr, var0)

        res = True
        res@gsLineColor = (/1.0, 0.0, 0.0/)
        res@gsLineThicknessF = 4.0
        res@gsLineDashPattern = 11
        res@tfPolyDrawOrder = "Draw"
        regl = gsn_add_polyline(wks, plot, (/-1, dim_MOD/), \
                                (/reg_std(ireg, idx_ref), reg_std(ireg, idx_ref)/), res)

        txres = True
        txres@txFontColor = (/1.0, 0.0, 0.0/)
        txres@txFontHeightF = 0.02
        txres@txJust = "CenterLeft"
        x_ndc = x
        y_ndc = y

        datatondc(plot, x, y, x_ndc, y_ndc)

        x_ndc = 0.925

        gsn_text_ndc(wks, "Obs.", x_ndc, y_ndc, txres)

        draw(plot)
        frame(wks)

        ; add meta data to plot (for reporting)

        caption = variables(0) + " anomaly standard deviation over " + region(ireg) + " region." \
                  + " Reference data set shown as dashed line: " + models@name(idx_ref) + "." \
                  + " Similar to IPCC AR5, fig. 9.36."
        id = diag_script + "_" + variables(0) + "_" + region(ireg)

        ESMValMD(outfile, alltags, caption, id, variables(0), allmodelnames, climofiles)

        delete(plot)
        delete(wks)
        delete(data_arr)

    end do  ; ireg-loop

    leave_msg(diag_script, "", 2)

end
