;;#############################################################################
;; ipcc-fig-9-30.ncl
;; Author(s): Axel Lauer (DLR, Germany)
;;
;; PROJECT-NAME: CRESCENDO
;;#############################################################################
;; Description
;; This script computes and outputs the mean diurnal cycle of precipitation
;; averaged over a given region and season(s).
;; The diagnostic is based on DiurnalCycle_box.ncl but in contrast to
;; DiurnalCycle_box.ncl, the results calculated by this diagnostic are averaged
;; over all grid cells with the same local time.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;   latrange_seasonal: latitude crop area
;;   lonrange_seasonal: lonitude crop area
;;              season: season in letters, e.g., JJAS
;;            styleset: chose styles set for lines
;;
;; Optional diag_script_info attributes (diagnostics specific)
;;    masking         0 = no masking (default)
;;                    1 = land points only
;;                    2 = ocean points only
;;    show_anomalies  False = plot absolute values (default)
;;                    True  = plot anomalies in % of daily mean
;;    curve_interpolation  False = no interpolation of the calculated
;;                                 diurnal cycle
;;                         True  = interpolation of the calculated
;;                                 diurnal cycle (useful for input data with
;;                                 a coarse temporal resolution such as e.g.
;;                                 3-hourly values)
;;
;; Caveats
;;    1) No common masking of missing values, i.e. missing values in the
;;       observations are NOT applied to the model data.
;;    2) Supported maximum temporal resolution of input data is 1 hour,
;;       minutes/seconds in input times will be ignored.
;;    3) Time coordinate of input data is assumed to be UTC.
;;
;; Modification history
;;    20171005-A_laue_ax: written.
;;
;;#############################################################################

;; Temporary file updated by Python on the fly

load "./interface_data/ncl.interface"

;;;load "./interface_scripts/auxiliary.ncl"
;;;load "./interface_scripts/data_handling.ncl"
;;;load "./interface_scripts/messaging.ncl"
;;;load "./diag_scripts/lib/ncl/set_operators.ncl"
;;;load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/meta_data.ncl"

exact_panel_position_file = "./diag_scripts/aux/SAMonsoon/exact_panel_positions_precip_world.ncl"
load "./plot_scripts/ncl/monsoon_panels.ncl"
;*********************************************************************

begin
    ;; ============================================
    ;; ||                                        ||
    ;; ||         READING CONFIGURATION          ||
    ;; ||                                        ||
    ;; ============================================

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    ;; 'diag_script' if fetched from ncl.interface
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot - " + diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    output_file_type = getenv("ESMValTool_output_file_type")
    if (ismissing(output_file_type)) then
        output_file_type = "ps"  ; default
    end if

    diag_script_base = basename(diag_script)

    ;;-----------------------
    ; References Information
    ;;-----------------------
    authors = "A_laue_ax"
    write_references(diag_script,            \  ; script name
                     authors,                \  ; authors
                     "",                     \  ; contributors
                     "",                     \  ; diag_references
                     "",                     \  ; obs_references
                     "P_embrace")               ; proj_acknowledgements

end

; Load variable_info

begin
    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ;; 'show_debuginfo' if fetched from ncl.interface
    if (.not. isdefined("show_debuginfo")) then
        debuginfo = False
    else
        debuginfo = show_debuginfo
    end if

    ;; Check valid plot fields
    valid_plot_field = (/"T2Is"/)
    ind_f = inlist(field_type0, valid_plot_field)
    if (any(ind_f .eq. False)) then
        error_msg("fatal", diag_script, "", "can't use field " + field_type0)
        status_exit(1)
    end if
    delete(valid_plot_field)
    delete(ind_f)

    ;; -------------------------------------------------
    ;; Check that required 'diag_script_info'-attributes
    ;; (from the DiurnalCycle cfg-file) are present
    ;; -------------------------------------------------
    required_attributes = (/"season",\
                            "styleset",\
                            "latrange",\
                            "lonrange"/)
    exit_if_missing_atts(diag_script_info, required_attributes)

    ;; Abbreviate a little...
    di = diag_script_info

    if (isatt(di, "masking")) then
        masking = di@masking
    else
        masking = 0
    end if

    if (isatt(di, "show_anomalies")) then
        show_anomalies = di@show_anomalies
    else
        show_anomalies = False
    end if

    season = di@season

     ;;------------------------------------------------------
    ;; Define transforms apply and and debug data to display
    ;;------------------------------------------------------
    transforms = True
    transforms@apply_these = (/"lonflip", "season", "latrange", "lonrange"/)
    transforms@season = season
    transforms@latrange = di@latrange
    transforms@lonrange = di@lonrange

    ;;------------------------------------------
    ;; Define meta data to debuginfo onto the plot
    ;;------------------------------------------
    min_lat = di@latrange(0)
    max_lat = di@latrange(1)
    debuginfo@latrange = min_lat + "-" + max_lat
    min_lon = di@lonrange(0)
    max_lon = di@lonrange(1)
    debuginfo@lonrange = min_lon + "-" + max_lon
    latlon_header_string = get_latlon_string(min_lat, max_lat, min_lon, max_lon)

    ;; Reference model index
    idx_ref = get_ref_model_idx(models, var_attr_ref(0))

    ;; Read all models and extract number of years
    nml_data_vault = True
    years = (/-999.0/)
    do idx_mod = 0, dim_MOD - 1
        nml_data_vault = read_data(idx_mod, var0, field_type0)
        tmp_model = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
        years_tmp = set_inclusive_OR(years, fetch_years(tmp_model))
        delete(years)
        years = years_tmp
        delete(years_tmp)

        delete(tmp_model)
    end do
    no_years = dimsizes(years)
    di@years = years

    model_years_consistent = True
    first_start_year = models@start_year(0)
    first_end_year = models@end_year(0)

    ;; ============================================
    ;; ||                                        ||
    ;; ||         COMPUTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================
    ;; -----------------------------------------------
    ;; Data structure for storing the processed fields
    ;; -----------------------------------------------
    storage_vault = True

    ;; Loop through models and compute diagnostics
    do idx_mod = 0, dim_MOD - 1

        ;; Extract and modify data
        A0_full = extract_data(idx_mod, nml_data_vault, -1, 0, 0)
        A0 = apply_transforms(transforms, A0_full)
        ; optionally extract second season and append to data array
        ; (needed for IPCC Ch.9 Fig 9.30)
        if (isatt(di, "season2")) then
            transforms@season = di@season2
            A0B = apply_transforms(transforms, A0_full)
            A0_tmp = array_append_record(A0, A0B, 0)
            delete(A0)
            delete(A0B)
            A0 = A0_tmp
            delete(A0_tmp)
            season = di@season + "+" + di@season2
        end if

        delete(A0_full)

        ;; Compute mean for each hour
        timeh = cd_calendar(A0&time, 0)

        nHR = 0
        hharray = new((/24/), integer)
        dims = dimsizes(A0)
        hrtmp = new((/24, dims(1), dims(2)/), float)

        do hh = 0, 23
            idx = ind(tointeger(timeh(:, 3)).eq.hh)
            if (.not.all(ismissing(idx))) then
                hharray(nHR) = hh
                tmp = dim_avg_n_Wrap(A0(idx, :, :), (/0/))
                hrtmp(nHR, :, :) = tmp
                delete(tmp)
                nHR = nHR + 1
            end if
            delete(idx)
        end do

        delete(A0)
        delete(timeh)

        hrtmp!0 = "hour"
        hrtmp&hour = hharray

        ; discard unused hour entries

        hravg = hrtmp(0:nHR-1, :, :)
        delete(hrtmp)

        ; optionally apply masking

        if (masking.ne.0) then
            f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
            a = f->LSMASK
            sftlf = byte2flt(landsea_mask(a, hravg&lat, hravg&lon))
            delete(a)
        end if

        if (masking.eq.1) then      ; 1 = land points only
            sftlf = where(sftlf.lt.1, sftlf@_FillValue, 0.0)
            hravg = hravg + conform(hravg, sftlf, (/1, 2/))
        else if (masking.eq.2) then ; 2 = ocean points only
            sftlf = where(sftlf.ge.1, sftlf@_FillValue, 0.0)
            hravg = hravg + conform(hravg, sftlf, (/1, 2/))
        end if
        end if

        if (masking.ne.0) then
            delete(sftlf)
        end if

;        ; debug output
;        system ("rm debug.nc")
;        debugfile = addfile("debug.nc","c")
;        debugfile->hravg = hravg
;        exit

        ; calculate local time (0-24h) at each grid cell and each step step

        hrlocal = conform(hravg, hravg&hour, (/0/)) + conform(hravg, hravg&lon, (/2/)) * 24.0 / 360.0
        hrlocal = where(hrlocal.gt.24.0, hrlocal - 24.0, hrlocal)
        hrlocal = where(hrlocal.lt.0.0, hrlocal + 24.0, hrlocal)

;        ; debug output
;        system ("rm debug.nc")
;        debugfile = addfile("debug.nc","c")
;        debugfile->hrlocal = hrlocal
;        exit

        ; ---------------------------------------------------------------------------
        ; calculate mean diurnal cycle by averaging over all grid cells in the domain
        ; with a similar local time, i.e. given hour +-0.5*delta_t
        ; ---------------------------------------------------------------------------

        n = nHR  ; number of points in diurnal cycle (use same number of time steps as
                 ; in original input files to avoid generation of artifacts)
        delta_t = 24.0 / n

        hr = fspan(0.0, n * delta_t, n + 1)
        dc = new((/n+1/), float)
        dc!0 = "hour"
        dc&hour = hr
        delete(hr)

        dims = dimsizes(hrlocal)
        hrlocal1D = ndtooned(hrlocal)  ; convert to 1-dim for using function 'ind'
        delete(hrlocal)

        data1D = ndtooned(hravg)
        dataselected1D = data1D

        do i = 0, n - 1
            dataselected1D = dataselected1D@_FillValue  ; reset all values

            tmin = (i - 0.5) * delta_t
            tmax = (i + 0.5) * delta_t
            idx = ind((hrlocal1D.ge.tmin) .and. (hrlocal1D.le.tmax))

            dataselected1D(idx) = data1D(idx)  ; copy values within local time interval
            delete(idx)

            ; if tmin < 0, then also add data in the interval 24-tmin...24

            if (tmin.lt.0.0) then
                tmin = 24.0 + tmin
                tmax = 24.0
                idx = ind((hrlocal1D.ge.tmin) .and. (hrlocal1D.le.tmax))
                dataselected1D(idx) = data1D(idx)  ; copy values within local time interval
                delete(idx)
            end if

            ; if tmax > 24, then also add data in the interval 0...tmax-24

            if (tmax.gt.24.0) then
                tmin = 0
                tmax = tmax - 24.0
                idx = ind((hrlocal1D.ge.tmin) .and. (hrlocal1D.le.tmax))
                dataselected1D(idx) = data1D(idx)  ; copy values within local time interval
                delete(idx)
            end if

            dataselected = onedtond(dataselected1D, dims)  ; convert back to 3-dim
            copy_VarMeta(hravg, dataselected)  ; restore dimensions and attributes
    
            ; average over time
            tmean = dim_avg_n_Wrap(dataselected, (/0/))
            delete(dataselected)

            ; calculate average over all grid cells for the given local time interval
            ; weighing boxes by surface area (domain has already been cropped --> use lon=-180...180,
            ; lat=-90...90 to average all remaining grid cells)

            dc(i) = area_operations(tmean, -90, 90, -180, 180, "average", True)
            delete(tmean)
        end do

        delete(dims)
        delete(hrlocal1D)
        delete(data1D)
        delete(dataselected1D)
        delete(hravg)

        ; ---------------------------------------------------------------------

;        ; calculate first two 24-h harmonics (diurnal and semi-diurnal)
;
;        finfo = fourier_info (dc, 2, delta_t)
;        ; amplitudes
;        A1 = finfo(0, 0)
;        A2 = finfo(0, 1)
;        ; phases (in rad)
;        sig1 = finfo(1, 0) / 24.0 * 2.0 * 3.1415927
;        sig2 = finfo(1, 1) / 24.0 * 2.0 * 3.1415927
;        ; time (in rad)
;        t = ispan(0, 24, 1)    ; the dirunal cycle reconstructed from the harmonics
;                               ; contains 25 1-hr time steps (0h, 1h, 2h, ..., 24h)
;        trad = t / 24.0 * 2.0 * 3.1415927  ; convert to rad
;        ; composit
;        dc_harmonic = avg(dc) + A1 * cos(trad - sig1) + A2 * cos(2.0 * trad - sig2)
;        dc_harmonic!0 = "hour"
;        dc_harmonic&hour = t
;        copy_VarAtts(dc, dc_harmonic)
;        delete(t)
;        delete(trad)
;
;        diurnalcycle = dc_harmonic

        ; ---------------------------------------------------------------------

        ; copy first value (0h) into extra entry (24h)

        t = dc&hour(n)
        dc(n) = dc(0)
        dc&hour(n) = t
        delete(t)
        diurnalcycle = dc

        ; ---------------------------------------------------------------------

        ; optionally, calculate anomalies in % of daily mean

        if (show_anomalies) then
            mean = avg(dc(0:n-1))
            if (mean.ne.0.0) then
                diurnalcycle = (dc - mean) / mean * 100.0
                diurnalcycle@units = "%"
            end if
        end if

        ; optionally use spline interpolation to generate a smoother diurnal cycle
        ; with a time resolution of 1 hour

        if (isatt(di, "curve_interpolation")) then
            if (di@curve_interpolation) then
                ti = ispan(0, 24, 1)  ; contains 25 1-hr time steps (0h, 1h, 2h, ..., 24h)
                dci = ftcurv(diurnalcycle&hour, diurnalcycle, ti)
                dci@units = diurnalcycle@units
                delete(diurnalcycle)
                diurnalcycle = dci
                diurnalcycle!0 = "hour"
                diurnalcycle&hour = ti
                diurnalcycle@units = dci@units
                delete(ti)
                delete(dci)
            end if
        end if

        diurnalcycle&hour@_FillValue = -1  ; needed to prevent "store_variable" from crashing

        ;; Check if start/end years are consistent across all models
        if ((first_start_year .ne. models@start_year(idx_mod))\
             .or. (first_end_year .ne. models@end_year(idx_mod))) then
            model_years_consistent = False
        end if

        ;; Store diurnal cycle
        diurnalcycle@idx = idx_mod
        diurnalcycle@model_name = models@name(idx_mod)
        storage_variable_name = str_join((/"diurnalcycle", sprinti("%i", idx_mod)/), str_vault_sep)
        store_variable(diurnalcycle, storage_vault, storage_variable_name)

        delete(diurnalcycle)
        delete(dc)

    end do  ; idx_mod = 0, dim_MOD - 1

    ;; Write vault to netCDF
    if (write_netcdf) then
        vault_to_netcdf(storage_vault, diag_script_base, var0, field_type0)
    end if

    ;; Set years_string depending on whether all models include the same years
    if (model_years_consistent) then
        years_string = "~F22~yrs:~F21~ " + models@start_year(0) + "-" + models@end_year(0)
    else
        years_string = "~F22~yrs:~F21~ inconsistent"
    end if
    di@left_plot_subheader = latlon_header_string + ", " + years_string
    di@years_string = years_string

    ;; ============================================
    ;; ||                                        ||
    ;; ||          PLOTTING DIAGNOSTICS          ||
    ;; ||                                        ||
    ;; ============================================

    dictkey_units = interface_get_dictkeys(var0, 0) + "_attributes_units"
    lb_units = nml_data_vault@$dictkey_units$

    res = True

    res@tmXBLabelFontHeightF = 0.008
    res@tmYLLabelFontHeightF = 0.008
    res@cnInfoLabelOn = False
    res@tmYRLabelsOn = False
    res@tmYRMode = "Automatic"
    res@tiMainFontHeightF = 0.015
    res@trXMinF = 0. 
    res@trXMaxF = 24.
    orig_res = res

    plot_settings = True
    plot_settings@x_axis_coord = "hour"
    plot_settings@gsnAddCyclic = True
    plot_settings@diag_script = diag_script
    plot_settings@idx_ref = idx_ref
    plot_settings@lb_title = lb_units
    plot_settings@part_of_header = "Precipitation"

    debuginfo@description_ycoord = 0.8
    debuginfo@description = "~F22~Diag def:~F21~ "\
                            + season\
                            + " average across indicated years"

    plot_settings@gsn_right_string = variable_info@long_name
    plot_settings@plot_pres_hgt = True

    if (model_years_consistent) then
        years_string = "~F22~yrs:~F21~ "\
                       + models@start_year(0)\
                       + "-" + models@end_year(0)
    else
        years_string = "~F22~yrs:~F21~ inconsistent"
    end if

    if (isatt(di, "plotfile_suffix")) then
        plot_settings@aux_filename_info = di@box + "-" + season + "-" + di@plotfile_suffix + "-diurnal-cycle"
        plot_settings@title_string = season + " Diurnal cycle over " + di@box + " (" + di@plotfile_suffix + ")"
    else
        plot_settings@aux_filename_info = di@box + "-" + season + "-diurnal-cycle"
        plot_settings@title_string = season + " Diurnal cycle over " + di@box
    end if
    if (show_anomalies) then
        plot_settings@yaxis_string = "Anomaly " + variable_info@long_name + " (%)"
    else
        plot_settings@yaxis_string = variable_info@long_name
    end if
    plot_settings@xaxis_string = "Local hour"
    plot_settings@type = "Not applicable"
    plot_settings@left_plot_subheader = latlon_header_string + ", " + years_string
    if (.not.isatt(di, "minmax")) then
        if (isatt(variable_info, "minmax")) then
            di@minmax = variable_info@minmax
        else
            ; default values to prevent the plotting routine from crashing
            di@minmax = (/0.0, 1.0/)
        end if
    end if

    xy_line_overlap(storage_vault, di, plot_settings, "diurnalcycle", debuginfo)

    ; -------------------------------------------------------------------------
    ; add metadata to plot

    xml_name = getenv("ESMValTool_xml_name")
    suffix = get_file_suffix(xml_name, 0)
    base_name = str_sub_str(xml_name, suffix, "")
    plot_outfile = get_output_dir(plot_dir, diag_script_base) + \
          interface_get_figure_filename(diag_script_base, var0, field_type0, \ 
          plot_settings@aux_filename_info, -1)
    alltags = array_append_record(tags, (/"PT_diurn", "DM_reg"/), 0)
    caption = "Mean Diurnal Cycle for " + variables(0) 
    id = diag_script + "_" + base_name
    climofiles = new(dimsizes(models@name) * dimsizes(variables), string)
    do imod = 0, dimsizes(models@name) - 1
        climofiles(imod) = \
            interface_get_inpaths(imod) + "/" + \
            interface_get_infile(variables(0), field_types(0), imod)
    end do
    ESMValMD(plot_outfile + "." + output_file_type, alltags, caption, id, variables, \
             models@name, climofiles, diag_script, authors)
    delete([/alltags, caption, id, climofiles/])
    ; -------------------------------------------------------------------------

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
