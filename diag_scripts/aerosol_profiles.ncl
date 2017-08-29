;;#############################################################################
;; DIAGNOSTIC SCRIPT FOR AEROSOL SIZE DISTRIBUTIONS AND VERTICAL PROFILES
;; Author: Mattia Righi (DLR, Germany)
;; ESMVal project
;;#############################################################################
;;
;; Description
;;    Plots aerosol vertical profiles (concentration vs. level) and aerosol
;;    size distributions (concentration vs. diameter), and compares with
;;    observational data.
;;    Path to observational data is provided in the cfg namelists.
;;    Observational files must have either a plev (for vertical profiles) or
;;    a diam (for size distributions) dimension. They can also optionally
;;    have a case dimension, for plotting multiple lines of the same
;;    observation on the same plot (e.g., measurements from different flights).
;;    It works with both monthly-mean, daily-mean and instantaneous data,
;;    although the latter is recommended for having better statistics.
;;    A climatology is constructed from the model data, using the same
;;    months (days) of the given campaign.
;;
;; Required diag_script_info attributes (diagnostic specific)
;;    campaign: name of the campaign for observational data
;;    datapath: path to the observational data (observations must be of type
;;              type T0M, must have a diam dimension and must contain the
;;              global attributes period or month, info, latmin, latmax,
;;              lonmin, lonmax (minplev, maxplev for T3M model data and size
;;              distribution diagnostics). They can optionally have a case
;;              dimension: all cases will be plotted on the same figure).
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    styleset: project-specific styleset.
;;    cc_threshold: maximum cloud cover to be considered for the analysis.
;;                  Gridboxes with a larger cloud cover value will be excluded.
;;                  This applies only to T3* fields.
;;
;; Required variable_info attributes (variable specific)
;;
;; Optional variable_info attributes (variable specific)
;;    profile_min: lower limit of variable's axis
;;    profile_max: upper limit of variable's axis
;;
;; Caveats
;;
;; Modification history
;;    20170705-A_laue_ax: added tags for reporting
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150325-A_laue_ax: modified reference tags used for acknowledgements
;;                        (projects, observations, etc)
;;    201409516-A_righ_ma: written.
;;
;;#############################################################################

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/scaling.ncl"
load "./reformat_scripts/constants.ncl"

load "./plot_scripts/ncl/xy_line.ncl"

load "./diag_scripts/lib/ncl/meta_data.ncl"

begin

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 2)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    do i = 0, dim_VAR - 1
        info_output(diag_script + " (var: " + variables(i) + ")", verbosity, 1)
    end do
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    ;; Reference information
    write_references(diag_script, \
                     "A_righ_ma", \  ; authors
                     "A_kais_ch", \  ; contributors
                     (/"D_lauer05acp", "D_aquila11gmd"/), \  ; diag_ref
                     (/"E_access2", "E_cirrus", "E_concert", "E_cr-ave", \
                       "E_dc3", "E_hippo", "E_inca", "E_lace1", "E_lace2", \
                       "E_oct-ave", "E_saltrace1", "E_saltrace2", "E_samum1", \
                       "E_samum2", "E_tc4", "E_ucnpacific"/), \  ; obs_ref
                     (/"P_esmval", "P_dlrveu", "P_dlrveu2"/))  ; proj_acknow

    ;; Check for valid fields
    do i = 0, dim_VAR - 1
        if (all(field_types(i).ne.\
            (/"T3M", "T3D", "T3I", "T2Ms", "T2Ds", "T2Is"/))) then
            error_msg("f", diag_script, "", "type " + field_types(i) + \
                      "not supported")
        end if
    end do

    if (variables(0).eq."conccnmode".and.variables(1).eq."diamcnmode") then
        l_sizedist = True
        varobs = "sizecn"
    else
        l_sizedist = False
        varobs = variables(0)
    end if

    ;; Check if cloud cover filtering is required
    l_cloud = False
    if (any(field_types(0).eq.(/"T3M", "T3D", "T3I"/))) then
        if (variables(dim_VAR - 1).ne."cl") then
            error_msg("f", diag_script, "", "cl variable (cloud cover) " + \
                      "required")
        end if
        l_cloud = True
    end if

    ;; Check if a summary plot (all flights together) is desired
    l_summary = False
    if (isatt(diag_script_info, "summary_plot")) then
        if (diag_script_info@summary_plot) then
            l_summary = True
            if (isatt(diag_script_info, "exclude_flights")) then
                exclude_flights = diag_script_info@exclude_flights
            end if
            nmax = 10000
        end if
    end if

    ;; Not compatible with size distributions
    if (l_sizedist.and.l_summary) then
        l_summary = False
        info_output("Summary plot available only for vertical profile, " + \
                    "ignoring", verbosity, 1)
    end if

    ;; Check required diag_script_info attributes
    exit_if_missing_atts(diag_script_info, (/"campaign", "datapath"/))
    if (l_cloud) then
        exit_if_missing_atts(diag_script_info, "cc_threshold")
    end if

    ;; Read data
    do idx_mod = 0, dim_MOD - 1
        data0 = read_data(idx_mod, variables(0), field_types(0))
        data1 = read_data(idx_mod, variables(1), field_types(1))
        if (dim_VAR.eq.3) then
            data2 = read_data(idx_mod, variables(2), field_types(2))
        end if
    end do

end

;; Load var_info
begin
    vardeffile = "interface_data/" + variables(0) + "_info.tmp"
    loadscript(vardeffile)
end

;; Main
begin

    ;; Get file type
    file_type = getenv("ESMValTool_output_file_type")

    ;; Check observational data
    obsfiles = systemfunc("ls " + diag_script_info@datapath + \
                          "/OBS_" + diag_script_info@campaign + \
                          "_campaign_*_T0M_" + varobs + "_*.nc")
    if (all(ismissing(obsfiles))) then
        info_output(diag_script_info@campaign + \
                    " observational data for " + varobs + \
                    " not available, skipping", verbosity, 1)
        exit
    end if

    ;; Print summary info
    info_output("Campaign: " + diag_script_info@campaign, verbosity, 2)
    do vID = 0, dim_VAR - 1
        info_output("Variable: " + variables(vID), verbosity, 2)
    end do
    info_output("Models:", verbosity, 2)
    do mID = 0, dimsizes(models@name) - 1
        info_output("  " + models@name(mID), verbosity, 2)
    end do

    info_output("Selected time range: " + \
                min(toint(models@start_year)) + "01-" + \
                max(toint(models@end_year)) + "12", \
                verbosity, 2)

    ;; Define statistical quantities
    stats = (/"N", "min", "max", "mean", "stddev", \
              "5", "10", "25", "median", "75", "90", "95"/)

    ;; Define the specific dimension name
    if (l_sizedist) then
        dimname = "diam"  ;; for size distributions
    else
        dimname = "plev"  ;; for vertical profiles
    end if

    ;; Set number of flights for the summary array
    if (l_summary) then
        nsummary = dimsizes(obsfiles)
        if (isdefined("exclude_flights")) then
            nsummary = nsummary - dimsizes(exclude_flights)
        end if
    end if

    ;; Loop over available observation files
    icp = 0
    do cpID = 0, dimsizes(obsfiles) - 1

        f = addfile(obsfiles(cpID), "r")
        cname = str_get_field(systemfunc("basename " + obsfiles(cpID)), 4, "_")

        info_output("Campaign: " + cname, verbosity, 1)

        ;; Define observations array
        if (isfilevar(f, "case")) then
            ncases = dimsizes(f->case)
            lcase = True
        else
            ncases = 1
            lcase = False
        end if

        obser_arr = new((/ncases, dimsizes(stats), dimsizes(f->$dimname$)/), \
                        float)
        obser_arr!0 = "case"
        obser_arr&case = ispan(1, ncases, 1)
        obser_arr!1 = "statistic"
        obser_arr&statistic = stats

        ;; Define observations summary array
        if (l_summary.and..not.isdefined("obser_all")) then
            obser_all = new((/nsummary, dimsizes(stats), \
                                  dimsizes(f->$dimname$)/), float)
        end if

        do st = 0, dimsizes(stats) - 1
            var = varobs + "_" + stats(st)
            if (isfilevar(f, var)) then
                if (lcase) then
                    obser_arr(:, st, :) = f->$var$
                else
                    obser_arr(0, st, :) = f->$var$
                end if
            end if
        end do

        ;; Add to the summary array
        if (l_summary) then
            if (.not.isdefined("exclude_flights")) then
                obser_all(icp, :, :) = obser_arr(0, :, :)
                icp = icp + 1
            else
                if (all(cname.ne.exclude_flights)) then
                    obser_all(icp, :, :) = obser_arr(0, :, :)
                    icp = icp + 1
                end if
            end if
        end if

        ;; Read the period/month attribute
        if (isatt(f, "period")) then
            start_month = toint(str_get_field( \
                                str_get_field(f@period, 1, "-"), 2, "/"))
            end_month = toint(str_get_field( \
                              str_get_field(f@period, 2, "-"), 2, "/"))
            start_day = toint(str_get_field( \
                              str_get_field(f@period, 1, "-"), 3, "/"))
            end_day = toint(str_get_field( \
                            str_get_field(f@period, 2, "-"), 3, "/"))
            ;; For monthly mean data consider full month
            if (isStrSubset(field_types(0), "M")) then
                start_day = 1
                end_day = 31
            end if
        else if (isatt(f, "month")) then
            months = f@month
        else
            error_msg("f", diag_script, "", "observational data must " + \
                      "contain the global attribute period or month")
        end if
        end if

        ;; Read models
        do mID = 0, dimsizes(models@name) - 1

            model0 = extract_data(mID, data0, -1, 0, 0)
            if (l_sizedist) then
                model1 = extract_data(mID, data1, -1, 0, 0)
                if (l_cloud) then
                    cloud = extract_data(mID, data2, -1, 0, 0)
                end if
            else
                if (l_cloud) then
                    cloud = extract_data(mID, data1, -1, 0, 0)
                end if
            end if

            date = cd_calendar(model0&time, 0)

            ;; Find indexes for the selected time period or given months
            if (isatt(f, "period")) then
                do yy = toint(min(date(:, 0))), toint(max(date(:, 0)))
                    idx1 = min(ind(date(:, 0).eq.yy .and. \
                                   date(:, 1).eq.start_month .and. \
                                   date(:, 2).ge.start_day))
                    idx2 = max(ind(date(:, 0).eq.yy .and. \
                                   date(:, 1).eq.end_month .and. \
                                   date(:, 2).le.end_day))
                    if (isdefined("idx")) then
                        tmp = array_append_record(idx, ispan(idx1, idx2, 1), 0)
                        delete(idx)
                        idx = tmp
                        delete(tmp)
                    else
                        idx = ispan(idx1, idx2, 1)
                    end if
                    delete(idx1)
                    delete(idx2)
                end do
            end if
            if (isatt(f, "month")) then
                do mm = 0, dimsizes(months) - 1
                    if (isdefined("idx")) then
                        tmp = array_append_record( \
                            idx, ind(date(:, 1).eq.months(mm)), 0)
                        delete(idx)
                        idx = tmp
                        delete(tmp)
                    else
                        idx = ind(date(:, 1).eq.months(mm))
                    end if
                end do
            end if
            delete(date)

            ;; Extract time steps
            rank = dimsizes(dimsizes(model0))
            if (l_sizedist) then
                if (rank.eq.5) then
                    nn_tmodel = model0(idx, :, :, :, :)
                    dd_tmodel = model1(idx, :, :, :, :)
                else if (rank.eq.4) then
                    nn_tmodel = model0(idx, :, :, :)
                    dd_tmodel = model1(idx, :, :, :)
                else
                    error_msg("f", diag_script, "", "invalid rank1")
                end if
                end if
                delete(rank)
                delete(model1)
            else
                tmodel = model0(idx, :, :, :)
            end if
            if (l_cloud) then
                tcloud = cloud(idx, :, :, :)  ; always rank 4
                delete(cloud)
            end if
            delete(model0)
            delete(idx)

            ;; Process model data for size distributions
            if (l_sizedist) then

                ;; Extract level (if required)
                if (isStrSubset(field_types(0), "T3")) then

                    ;; Set target levels
                    exit_if_missing_atts(f, (/"plevmin", "plevmax"/))
                    if (f@plevmin.eq.f@plevmax) then
                        outlev = f@plevmax
                    else if (abs(f@plevmin - f@plevmax).lt.1.e4) then
                        outlev = (/f@plevmax, f@plevmin/)
                    else
                        ;; FIXME: implement a reasonable level sampling
                        error_msg("f", diag_script, "", \
                                  "the level range is too large")
                    end if
                    end if

                    ;; Check units consistency
                    nn_pmodel = int2p_n_Wrap(nn_tmodel&plev, nn_tmodel, \
                                             outlev, 1, 1)
                    dd_pmodel = int2p_n_Wrap(dd_tmodel&plev, dd_tmodel, \
                                             outlev, 1, 1)
                    if (l_cloud) then
                        pcloud = int2p_n_Wrap(tcloud&plev, tcloud, \
                                              outlev, 1, 1)
                        delete(tcloud)
                    end if
                    delete(outlev)

                ;; No level interpolation, just rename arrays
                else

                    nn_pmodel = nn_tmodel
                    dd_pmodel = dd_tmodel
                    if (l_cloud) then
                        pcloud = tcloud
                        delete(tcloud)
                    end if

                end if
                delete(nn_tmodel)
                delete(dd_tmodel)

                ;; Loop over modes
                rank = dimsizes(dimsizes(nn_pmodel))
                do imod = 0, dimsizes(nn_pmodel&mode) - 1
                    sigma = nn_pmodel&mode@sigma(imod)

                    ;; Extract mode
                    if (rank.eq.4) then
                        nn = nn_pmodel(:, :, :, imod)
                        dd = dd_pmodel(:, :, :, imod)
                    else if (rank.eq.5) then
                        nn = nn_pmodel(:, :, :, :, imod)
                        dd = dd_pmodel(:, :, :, :, imod)
                    else
                        error_msg("f", diag_script, "", "invalid rank2")
                    end if
                    end if

                    ;; Apply cloud filtering
                    if (l_cloud)
                        nn = where(pcloud.ge.diag_script_info@cc_threshold, \
                                   nn@_FillValue, nn)
                        dd = where(pcloud.ge.diag_script_info@cc_threshold, \
                                   dd@_FillValue, dd)
                    end if

                    ;; Extract given region
                    exit_if_missing_atts(f, (/"latmin", "latmax", \
                                              "lonmin", "lonmax"/))
                    nn_extr = area_operations(nn, f@latmin, f@latmax, \
                                              f@lonmin, f@lonmax, \
                                              "extract", False)
                    delete(nn)
                    nn_amodel = ndtooned(nn_extr)
                    delete(nn_extr)
                    dd_extr = area_operations(dd, f@latmin, f@latmax, \
                                              f@lonmin, f@lonmax, \
                                              "extract", False)
                    delete(dd)
                    dd_amodel = ndtooned(dd_extr)
                    delete(dd_extr)

                    ;; Calculate size distribution
                    if (.not.isdefined("amodel")) then
                        amodel = \
                            lognormal_dist(nn_amodel, dd_amodel, sigma, \
                                           obser_arr&diam)
                    else
                        amodel = amodel + \
                            lognormal_dist(nn_amodel, dd_amodel, sigma, \
                                           obser_arr&diam)
                    end if

                    delete(nn_amodel)
                    delete(dd_amodel)
                end do
                delete(nn_pmodel)
                delete(dd_pmodel)
                delete(rank)
                if (l_cloud)
                    delete(pcloud)
                end if

            ;; Process model data for vertical profiles
            else

                ;; Apply cloud filtering
                if (l_cloud) then
                    tmodel = where(tcloud.ge.diag_script_info@cc_threshold, \
                                   tmodel@_FillValue, tmodel)
                    delete(tcloud)
                end if

                ;; Extract selected region
                amodel = area_operations(tmodel, f@latmin, f@latmax, \
                                         f@lonmin, f@lonmax, "extract", False)
                delete(tmodel)

            end if

            ;; Define array of model data
            if (mID.eq.0) then
                dim2 = amodel&$dimname$
                model_arr = new((/dim_MOD, dimsizes(stats), dimsizes(dim2)/), \
                                float)
                model_arr!0 = "models"
                model_arr!1 = "statistic"
                model_arr!2 = dimname
                model_arr&models = project_style(diag_script_info, "annots")
                model_arr&statistic = stats
                model_arr&$dimname$ = dim2
            end if
            copy_VarAtts(obser_arr, model_arr)
            if (l_summary.and..not.isdefined("model_all")) then
                model_all = model_arr
                model1d_all = new((/dim_MOD, dimsizes(dim2), nmax/), float)
            end if

            ;; Calculate relevant statistics
            rank = dimsizes(dimsizes(amodel))
            dpos = ind(dimname.eq.getvardims(amodel))
            do kk = 0, dimsizes(dim2) - 1
                if (rank.eq.4) then
                    if (dpos.eq.0) then
                        model1d = ndtooned(amodel(kk, :, :, :))
                    end if
                    if (dpos.eq.1) then
                        model1d = ndtooned(amodel(:, kk, :, :))
                    end if
                end if
                if (rank.eq.3) then
                    if (dpos.eq.0) then
                        model1d = ndtooned(amodel(kk, :, :))
                    end if
                    if (dpos.eq.1) then
                        model1d = ndtooned(amodel(:, kk, :))
                    end if
                end if
                if (rank.eq.2) then
                    if (dpos.eq.0)
                        model1d = amodel(kk, :)
                    end if
                    if (dpos.eq.1) then
                        model1d = amodel(:, kk)
                    end if
                end if
                if (rank.eq.1) then
                    model1d = amodel
                end if

                if (.not.isdefined("model1d")) then
                    error_msg("f", "aerosol_profiles.ncl", "", "invalid rank3")
                end if

                ;; Skip if all missing
                if (all(ismissing(model1d))) then
                    delete(model1d)
                    continue
                end if

                ;; Add to the summary array
                if (l_summary) then
                    if (all(ismissing(model1d_all(mID, kk, :)))) then
                        i1 = 0
                    else
                        i1 = max(ind(.not.ismissing(model1d_all(mID, kk, :))))
                        i1 = i1 + 1
                    end if
                    i2 = i1 + dimsizes(model1d) - 1
                    if (i2.gt.nmax - 1) then
                        error_msg("f", "aerosol_profiles.ncl", "", \
                                  " array size too small, increase nmax")
                    end if
                    model1d_all(mID, kk, i1:i2) = model1d
                end if

                ;; Calculate statistics
                do st = 0, dimsizes(stats) - 1
                    model_arr(mID, st, kk) = \
                        distrib_stats(model1d, stats(st))
                end do
                delete(model1d)
            end do
            delete(amodel)

        end do  ; models loop

        if (isdefined("months")) then
            delete(months)
        end if

        ;; Set attributes
        obser_arr@var = varobs
        obser_arr@diag_script = (/diag_script/)
        model_arr@var = varobs
        model_arr@diag_script = (/diag_script/)

        ;; Optional output to NetCDF
        if (getenv("ESMValTool_write_netcdf").eq."True") then
            ;; Model
            new_path = get_output_dir(get_work_dir(), basename(diag_script))
            new_path = new_path + basename(diag_script) + "_" + \
                diag_script_info@campaign + "_" + cname + "_" + \
                varobs + "_model.nc"
            model_arr@ncdf = new_path
            ncdf_outfile = ncdf_write(model_arr, new_path)
            ;; Observations
            new_path = get_output_dir(get_work_dir(), basename(diag_script))
            new_path = new_path + basename(diag_script) + "_" + \
                diag_script_info@campaign + "_" + cname + "_" + \
                varobs + "_obser.nc"
            obser_arr@ncdf = new_path
            ncdf_outfile = ncdf_write(obser_arr, new_path)
        end if

        ;; Convert units for plotting (if required)
        if (isatt(variable_info, "plot_units")) then
            idN = ind(stats.ne."N")
            model_arr(:, idN, :) = \
                convert_units(model_arr(:, idN, :), variable_info@plot_units)
            idN = ind(stats.ne."N")
            obser_arr(:, idN, :) = \
                convert_units(obser_arr(:, idN, :), variable_info@plot_units)
        end if

        if (l_sizedist) then
            ;; Convert diameter units [m] --> [um]
            model_arr&diam = model_arr&diam * 1.e6
            model_arr&diam@units = "~F33~m~F21~m"
            obser_arr&diam = obser_arr&diam * 1.e6
            obser_arr&diam@units = "~F33~m~F21~m"
        else
            ;; Convert level units [Pa] --> [hPa]
            model_arr&plev = model_arr&plev / 100.
            model_arr&plev@units = "hPa"
            obser_arr&plev = obser_arr&plev / 100.
            obser_arr&plev@units = "hPa"
        end if

        ;; Plot
        if (all(ismissing(model_arr))) then
            delete(model_arr)
            delete(obser_arr)
            continue
        end if

        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, basename(diag_script))
        if (l_sizedist) then
            plot_file = output_dir + diag_script_info@campaign + "_" + \
                cname + "_sizedist_" + varobs
        else
            plot_file = output_dir + diag_script_info@campaign + "_" + \
                cname + "_vertprof_" + varobs
        end if
        wks = gsn_open_wks(file_type, plot_file)
        model_arr@res_tiMainString = diag_script_info@campaign
        if (f@info.ne."") then
            model_arr@res_tiMainString = model_arr@res_tiMainString + \
                " - " + f@info
        end if
        if (isatt(variables(0), "profile_min")) then
            if (l_sizedist) then
                model_arr@res_trYMinF = variable_info@profile_min
                if (l_summary) then
                    model_all@res_trYMinF = variable_info@profile_min
                end if
            else
                model_arr@res_trXMinF = variable_info@profile_min
                if (l_summary) then
                    model_all@res_trXMinF = variable_info@profile_min
                end if
            end if
        end if
        if (isatt(variables(0), "profile_max")) then
            if (l_sizedist) then
                model_arr@res_trYMaxF = variable_info@profile_max
                if (l_summary) then
                    model_all@res_trYMaxF = variable_info@profile_max
                end if
            else
                model_arr@res_trXMaxF = variable_info@profile_max
                if (l_summary) then
                    model_all@res_trXMaxF = variable_info@profile_max
                end if
            end if
        end if

        if (l_sizedist) then
            plot = aerosol_sizedist(wks, model_arr, obser_arr, varobs)
        else
            plot = aerosol_profile(wks, model_arr, obser_arr, varobs)
        end if
        delete(model_arr)
        delete(obser_arr)

        ; -------------------------------------------------------------------------------
        ; add meta data to plot (for reporting)

        climofiles = new(dim_MOD + 1, string)
        do ii = 0, dim_MOD - 1
            climofiles(ii) = \
                interface_get_inpaths(ii) + "/" \
                + interface_get_infile(variables(0), field_types(0), ii)
        end do
        climofiles(dim_MOD) = obsfiles(cpID)
        
        ttags = array_append_record(tags, (/"ST_clim", "DM_reg"/), 0)

        if (l_sizedist) then
            alltags = array_append_record(ttags, (/"PT_size"/), 0)
            stype = "aerosol size-distribution"
        else
            alltags = array_append_record(ttags, (/"PT_vert"/), 0)
            stype = "aerosol concentration"
        end if

        caption = "Mean (long dashed) and standard deviation (sort-dashed) " \
                  + stype + " (campaign: " + diag_script_info@campaign + ", " + cname + ")"
        id = diag_script + "_" + variables(0) + "_" + diag_script_info@campaign + "_" + cname

        allmodels = array_append_record(models@name, diag_script_info@campaign, 0)

        contrib_authors = (/"A_righ_ma", "A_kais_ch"/)

        ESMValMD(plot_file + "." + file_type, alltags, caption, id, variables(0), allmodels, climofiles, diag_script, contrib_authors)

        delete(climofiles)
        delete(ttags)
        delete(alltags)
        delete(allmodels)
        ; -------------------------------------------------------------------------------

        info_output("wrote " + plot_file + "." + file_type + " &", verbosity, 1)

    end do  ; observations loop

    ;; Create summary plot
    if (l_summary) then

        info_output("Creating summary plot with all campaigns", verbosity, 1)
        info_output("Number of campaigns: " + nsummary, verbosity, 1)
        if (isdefined("exclude_flights")) then
            info_output("Excluded campaigns: ", verbosity, 1)
            do ii = 0, dimsizes(exclude_flights) - 1
                info_output("  - " + exclude_flights(ii), verbosity, 1)
            end do
        end if

        ;; Calculate statistics
        do kk = 0, dimsizes(dim2) - 1
            do mID = 0, dimsizes(models@name) - 1
                do st = 0, dimsizes(stats) - 1
                    maxid = \
                        max(ind(.not.ismissing(model1d_all(mID, kk, :))))
                    model1d = model1d_all(mID, kk, 0:maxid)
                    model_all(mID, st, kk) = \
                        distrib_stats(model1d, stats(st))
                    delete(model1d)
                end do
            end do
        end do

        ;; Convert level units [Pa] --> [hPa]
        model_all&plev = model_all&plev / 100.
        model_all&plev@units = "hPa"
        obser_all&plev = obser_all&plev / 100.
        obser_all&plev@units = "hPa"

        ;; Convert units for plotting (if required)
        if (isatt(variable_info, "plot_units")) then
            idN = ind(stats.ne."N")
            model_all(:, idN, :) = \
                convert_units(model_all(:, idN, :), variable_info@plot_units)
            idN = ind(stats.ne."N")
            obser_all(:, idN, :) = \
                convert_units(obser_all(:, idN, :), variable_info@plot_units)
        end if

        model_all@res_tiMainString = diag_script_info@campaign
        plot_file = output_dir + diag_script_info@campaign + "_ALL_" + \
            "vertprof_" + varobs
        wks = gsn_open_wks(file_type, plot_file)
        plot = aerosol_profile(wks, model_all, obser_all, varobs)
        delete(model_all)
        delete(obser_all)

        ; -------------------------------------------------------------------------------
        ; add meta data to plot (for reporting)

        modelfiles = new(dim_MOD, string)
        do ii = 0, dim_MOD - 1
            modelfiles(ii) = \
                interface_get_inpaths(ii) + "/" \
                + interface_get_infile(variables(0), field_types(0), ii)
        end do

        climofiles = array_append_record(modelfiles, obsfiles, 0)
        delete(modelfiles)
        
        ttags = array_append_record(tags, (/"ST_clim", "DM_reg"/), 0)

        if (l_sizedist) then
            alltags = array_append_record(ttags, (/"PT_size"/), 0)
            stype = "aerosol size-distribution"
        else
            alltags = array_append_record(ttags, (/"PT_vert"/), 0)
            stype = "aerosol concentration"
        end if

        caption = "Summary plot: " \
                  + stype + " (campaign: " + diag_script_info@campaign + ")"
        id = diag_script + "_" + variables(0) + "_" + diag_script_info@campaign + "_" + cname

        allmodels = array_append_record(models@name, diag_script_info@campaign, 0)

        contrib_authors = (/"A_righ_ma", "A_kais_ch"/)

        ESMValMD(plot_file + "." + file_type, alltags, caption, id, variables(0), allmodels, climofiles, diag_script, contrib_authors)

        delete(climofiles)
        delete(ttags)
        delete(alltags)
        delete(allmodels)
        ; -------------------------------------------------------------------------------


    end if


end
