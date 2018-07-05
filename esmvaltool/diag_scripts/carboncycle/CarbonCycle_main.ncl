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
;;    anav_month : Conversion of y-axis to PgC/month instead of /year,
;;                 following Anav2013 fig 7, 9 (cycle plot)
;;    evolution_plot: time evolution of a variable comparing obs to
;;                    multi-model mean; requires ref_model in nml,
;;                    and at least one model, one obs
;;    fignum: Reference to figure in Anav paper for evolution plot. "1 top"
;;            as default 
;;    evolution_plot_ref_model: reference model for evolution_plot
;;    evolution_plot_anomaly: makes evolution_plot an anomaly plot
;;    evolution_plot_ignore: Models to ignore, obs are excluded automatically
;;    evolution_plot_volcanoes: Turns on/off lines of volcano eruptions
;;    evolution_plot_color: Hue of the contours; default: red = 0
;;
;;    FOR EXTERNAL DATA-SET INPUT
;;    diag_script_info@obsname    : Name of ext dataset to read in
;;    diag_script_info@ref_varname: Varname in attributes of ext file
;;    diag_script_info@obsfile    : Filename of external obs to read in
;;
;;    dst_grid: Destination grid for regridding
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
;;    20180619-A_wenz_sa: adopted script to v2
;;    20170330-A_gier_be: added external data-read in (e.g. GCP)
;;    20170316-A_gier_be: added cdo regridding, landmask for CRU obs
;;    20170308-A_gier_be: added "evolution_plot"
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
load "./interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/scaling.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/regridding.ncl"
load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/watermark.ncl"

load "./diag_scripts/shared/plot/scatterplot.ncl"
load "./diag_scripts/shared/plot/xy_line.ncl"
load "./diag_scripts/shared/plot/style."

begin

    enter_msg(diag_script, "")
    log_info("++++++++++++++++++++++++++++++++++++++++++")
    log_info(diag_script + " (var: " + variables(0) + ")")
    log_info("++++++++++++++++++++++++++++++++++++++++++")

    dim_MOD = dimsizes(model_info@model)
    dim_MOD_org = dim_MOD
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    write_references(diag_script,                  \
                     (/"A_anav_al"/),              \  ; authors
                     (/"A_wenz_sa", "A_righ_ma"/), \  ; contributors
                     (/"D_anav13jclim"/),          \  ; diag_references
                     (/""/),                       \  ; obs_references
                     (/"P_embrace"/))                 ; proj_references

    ;; Basename of diag_script
    diag_script_base = basename(diag_script)

    ;; Define file type
    file_type = config_user_info@output_file_type
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

    ;; Set region - supports 2D array for lat and strings
    if isstring(diag_script_info@region) then
        region = select_region(diag_script_info@region)
    else if dimsizes(diag_script_info@region).eq.2 then
        region = (/diag_script_info@region(0), diag_script_info@region(1), 0., 360./)        
        region@name = flt2string(diag_script_info@region(0)) + "_-_" + \
                      flt2string(diag_script_info@region(1))
        region@DM_tag = "DM_reg"
        else
            error_msg("f", diag_script, "", \
                  diag_script_info@region + " is an invalid region type")
        end if
    end if

    ;; Output plot directory, set path for saving processed data
    output_dir = config_user_info@plot_dir
    if (.not.isfilepresent(config_user_info@plot_dir)) then
      system("mkdir -p " + config_user_info@plot_dir)
    end if

    ;; Set time range
    start_year = min(stringtoint(model_info@start_year))
    end_year = max(stringtoint(model_info@end_year))

    ;; Only evolution plot utilizes timeranges greater then 1986-2005
    mean_minyear = 1986
    mean_maxyear = 2005
    minyear_ind = ind(ispan(start_year, end_year, 1).eq.mean_minyear)
    maxyear_ind = ind(ispan(start_year, end_year, 1).eq.mean_maxyear)

    ;; Global array for collecting all data
    all_yearly = new((/dim_MOD, end_year - start_year + 1, 2/), float)
    all_yearly!0          = "model"
    all_yearly!1          = "year"
    all_yearly!2          = "statistic"
    all_yearly&model      = model_info@model
    all_yearly&year       = ispan(start_year, end_year, 1)
    all_yearly&statistic  = (/"mean", "stddev"/)

    all_monthly = new((/dim_MOD, 12, 2/), float)
    all_monthly!0         = "model"
    all_monthly!1         = "month"
    all_monthly!2         = "statistic"
    all_monthly&model     = model_info@model
    all_monthly&month     = ispan(1, 12, 1)
    all_monthly&statistic = (/"mean", "stddev"/)

    ;; Check for irregular grid
    l_irreg = new(dim_MOD, logical)
    l_irreg = False
    do imod = 0, dim_MOD - 1
        cnt = interface_get_dictkeys(var0, imod) + "_dimNames"
        dims = data@$cnt$
        ndims = dimsizes(dims)
        if (ndims.lt.2) then
            continue
        end if
        if (ndims.eq.2 .and. dims(0).eq."time") then
            continue
        end if
        if (dims(ndims - 1).ne."lon" .and. dims(ndims - 2).ne."lat") then
            l_irreg(imod) = True
            continue
        end if
    end do
    delete([/cnt, dims, ndims/])

    ;; Set destination grid
    grid_resolution = diag_script_info@grid
    if (grid_resolution.eq."2deg") then
        dst_grid = "./diag_scripts/aux/CDO/cdo_dst_grid_g020"
    end if

    ;; Loop over models
    do imod = 0, dim_MOD - 1

        log_info("Processing " + model_info@model(imod))

        ;; Extract model corresponding data
        var  = read_data(imod, var0, field_type0)

             ;; Regridding
            if (l_irreg(imod)) then
                regr_file = basename(diag_script) + "_" + model_info@model(imod)+ "_" + \
                              var0 + "_" + start_year + "-" + end_year+ "_grid_" + \
                              grid_resolution
                              regr_path = wrk_dir + "/CDO/" + regr_file + ".nc"            
                if (isfilepresent(regr_path)) then
                    var_regrid = ncdf_read(regr_path, var0)
                else
                    var_regrid = cdo_remapdis(var0, var, imod, data, \
                             dst_grid, regr_file)                
                end if
                var_regrid := var_regrid(time|:, lat|:, lon|:)
                delete(var)
                var = var_regrid
                delete(var_regrid)
            end if

            ;; CRU is only land based - Need landmask and also cut off Antarctica
            if (isatt(diag_script_info, "ref_model")).and. \
                                variable_info@reference_model.eq."CRU" then
                f = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc", "r")
                a = f->LSMASK
                sftlf = byte2flt(landsea_mask(a, var&lat, var&lon))
                sftlf = where(sftlf.gt.1., 1., sftlf)
                ;; Mask out Antarctica (below 60degS)
                sftlf({:-60}, :) = sftlf@_FillValue
                ;; Apply mask
                var = var * conform(var, sftlf, (/1, 2/))
                delete([/sftlf, a, f/])           
            end if

            ;; Set type of area operation (flux vs. non-flux variables)
            if (var@units.eq."kg m-2 s-1") then
                sp_opt = "sum"  ; for flux variable
            else
                sp_opt = "average"  ; for non-flux variable
            end if

            ;; Calculate annual mean
            data_yearly = \
                   time_operations(var, toint(model_info@start_year(imod)), \
                                   toint(model_info@end_year(imod)), \
                                    "average", "yearly", True)

            ;; Calculate mean annual cycle
            data_monthly = \
                    time_operations(var, max((/mean_minyear, \
                                    toint(model_info@start_year(imod))/)), \
                                    min((/mean_maxyear, \
                                    toint(models@end_year(imod))/)), "average", \
                                     "monthlyclim", True)

            ;; Special case for LAI.
            ;; Since model_info have different land covers, all sparse vegetated points
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
;;          delete(ystd)
            delete(idx1)
            delete(idx2)

             ;; Monthly data - Mean
            mmean = area_operations(data_monthly, region(0), region(1), \
                                region(2), region(3), sp_opt, True)

         ;; Monthly data - Standard deviation
;; FIX-ME: results to be checked
;;          tmp = data_monthly
;;          tmp = data_monthly * conform(data_monthly, area, (/1, 2/))
;;          subr = area_operations(tmp, region(0), region(1), \
;;                                region(2), region(3), "extract", True)
;;          mstd = dim_stddev_n_Wrap(subr, (/1, 2/))
;;          delete(subr)
;;          delete(tmp)

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
    

    ;; Read in external obs-data if necessary
    if (isatt(diag_script_info, "obsfile")) then
        log_info("Reading in additional data from: " \
                    + diag_script_info@obsfile)
        data_dum = addfile(diag_script_info@obsfile, "r")
        data_dum2 = data_dum->$diag_script_info@ref_varname$
        
        ;; find out the overlap of the timeframes
        obs_time = cd_calendar(data_dum2&time, 0)
        min_obs_year = toint(min(obs_time(:,0)))
        max_obs_year = toint(max(obs_time(:,0)))
        obs_years = ispan(min_obs_year, max_obs_year, 1)
        if start_year.lt.min_obs_year then
            l_idx1 = ind(all_yearly&year.eq.min_obs_year)
            l_idx2 = 0
        else
            l_idx1 = 0
            l_idx2 = ind(start_year.eq.obs_years)
        end if
        if end_year.lt.max_obs_year then
            u_idx1 = dimsizes(all_yearly&year)-1
            u_idx2 = ind(end_year.eq.obs_years)
        else
            u_idx1 = ind(all_yearly&year.eq.max_obs_year)
            u_idx2 = dimsizes(obs_years)-1
        end if

        ;; Make an array like all_yearly to easily append data
        all_yearly_add = new((/1, end_year - start_year + 1, 2/), float)
        all_yearly_add!0          = "model"
        all_yearly_add!1          = "year"
        all_yearly_add!2          = "statistic"
        all_yearly_add&model      = diag_script_info@obsname
        all_yearly_add&year       = ispan(start_year, end_year, 1)
        all_yearly_add&statistic  = (/"mean", "stddev"/)

        ;; Data saved in one month per year, also put it in right timerange
        obs_ind = ind(.not.ismissing(data_dum2))
        all_yearly_add(0, l_idx1:u_idx1, 0) = (/data_dum2(obs_ind(l_idx2:u_idx2))/)
        all_yearly_add(0, l_idx1:u_idx1, 1) = 0.  ; (/ystd/)

        ;; Add data to existing arrays, add project, name to models attributes
        all_yearly_new = array_append_record(all_yearly, all_yearly_add, 0)
        delete([/data_dum, data_dum2, min_obs_year, max_obs_year, obs_time, \
                all_yearly_add, obs_years, obs_ind/])
        model_info@model_new = array_append_record(model_info@model, diag_script_info@obsname, 0)
        model_info@project_new = array_append_record(model_info@project, "OBS", 0)
        delete([/model_info@model, model_info@project, all_yearly/])
        model_info@project = model_info@project_new
        model_info@model = model_info@model_new
        all_yearly = all_yearly_new
        delete([/model_info@project_new, model_info@model_new, all_yearly_new/])
        dim_MOD = dimsizes(model_info@model)  ; update number of models
    end if

    ;; Convert units if appropriate [kg m-2 s-1] --> [PgC y-1]
    if (isatt(variable_info, "plot_units")) then
        all_yearly = convert_units(all_yearly, variable_info@plot_units)
        all_monthly = convert_units(all_monthly, variable_info@plot_units)
    end if
    ;; Don't know how to do this in var_defs without losing irreg grid info
    if (var0.eq."tos") then
        all_yearly = all_yearly - 273.15
        all_yearly@units = "degC"
        all_monthly = all_monthly - 273.15
        all_monthly@units = "degC"
    end if

    ;; Optional output to NetCDF
    if (config_user_info@write_netcdf.eq."True") then

      ;; Set path for saving processed data
      if (.not.isfilepresent(config_user_info@work_dir)) then
        system("mkdir -p " + config_user_info@work_dir)
      end if

        ;; Yearly data
        new_path = config_user_info@work_dir + \
                   basename(diag_script) + "_yearly_" + var0 + "_" + \
            region@name + ".nc"
        all_yearly@var = var0
        all_yearly@diag_script = diag_script
        all_yearly@ncdf = new_path
        ncdf_outfile = ncdf_write(all_yearly, new_path)

        ;; Monthly data
        new_path = config_user_info@work_dir + \
                   basename(diag_script) + "_monthly_" + var0 + "_" + \
            region@name + ".nc"
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
    
        ;; Add option for PgC/month y-axis that is used in Anav2013 plots 7, 9
        if (isatt(diag_script_info, "anav_month")).and. \
                          diag_script_info@anav_month then
            data_arr = data_arr/12.
            data_arr@units = "PgC/month"
        end if
        
        ;; Draw plot
        outfile = output_dir + diag_script + "_" + var0 + "_cycle_" + \
            region@name
        wks = gsn_open_wks(file_type, outfile)
        wks@legendfile = outfile + "_legend"
        data_arr@res_tiYAxisString = var0 + " [" + \
            format_units(data_arr@units) + "]"
        data_arr@res_tiMainString = region@name
        plot = cycle_plot(wks, data_arr, var0)
        draw(plot)
        frame(wks)

        ;; Attach Metadata to plot

        ;alltags = array_append_record(tags, region@DM_tag, 0)
        ;alltags := array_append_record(alltags,\
        ;       (/"PT_seas", "ST_mean"/), 0)
        ;caption = "Seasonal Cycle plot for " + variables(0) + \
       ;           " over the period " + start_year + \
       ;           "-" + end_year + ", similar to Anav et al., 2013, Fig. 9."
      ;  id = diag_script
     ;   isize = dim_MOD_org * dimsizes(variables)
      ;  if (isatt(diag_script_info, "obsfile")) then
      ;      isize = isize + 1
      ;;  end if
       ; climofiles = new(isize, string)
       ; do n_var = 0, dimsizes(variables) - 1
       ;     do imod = 0, dim_MOD_org - 1
       ;         climofiles(n_var * dim_MOD_org + imod) = \
       ;                     interface_get_inpaths(imod) + "/" + \
       ;                     interface_get_infile(variables(n_var), field_types(n_var), imod)
       ;     end do
       ; end do
       ; if (isatt(diag_script_info, "obsfile")) then
       ;     climofiles(isize - 1) = diag_script_info@obsfile
       ; end if
      ;  contrib_authors = (/"A_anav_al", "A_wenz_sa", "A_righ_ma"/)

        ;;NOTE: not implemented in v2 yet
    ;ESMValMD(outfile+"."+file_type, alltags, caption, id, variables, model_info@model, climofiles, diag_script, contrib_authors)

        delete([/alltags, caption, id, climofiles/])

        delete(plot)
        delete(wks)
        delete(data_arr)
    end if

    ;; Error-bar plot
    if (diag_script_info@errorbar_plot) then

        annual_mean = dim_avg_Wrap(all_yearly(:, minyear_ind:maxyear_ind, 0))
        annual_std = dim_stddev_Wrap(all_yearly(:, minyear_ind:maxyear_ind, 0))

        ;; Sort model in alphabetical order (observations at the beginning)
        if (isatt(diag_script_info, "sort")) then
            if (diag_script_info@sort) then
                pid = sort_alphabetically(annual_mean&model, \
                                          ind(model_info@project.eq."OBS"), \
                                          "begin")
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
        data_arr = new((/2, dimsizes(annual_mean)/), float)
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
        data_arr@res_tiMainString = var0 + " - " + region@name
        plot = errorbar_plot(wks, data_arr, var0)
        draw(plot)
        frame(wks)

        ;; Attach Metadata to plot

        alltags = array_append_record(tags, region@DM_tag, 0)
        alltags := array_append_record(alltags,\
               (/"PT_errorbar", "ST_stddev", "ST_mean"/), 0)
        caption = "Errorbar plot for " + variables(0) + \
                  " over the period " + start_year + \
                  "-" + end_year + ", similar to Anav et al., 2013, Fig. 6."
        id = diag_script

        isize = dim_MOD_org * dimsizes(variables)
        if (isatt(diag_script_info, "obsfile")) then
            isize = isize + 1
        end if
        climofiles = new(isize, string)
        do n_var = 0, dimsizes(variables) - 1
            do imod = 0, dim_MOD_org - 1
                climofiles(n_var * dim_MOD_org + imod) = \
                            interface_get_inpaths(imod) + "/" + \
                            interface_get_infile(variables(n_var), field_types(n_var), imod)
            end do
        end do
        if (isatt(diag_script_info, "obsfile")) then
            climofiles(isize - 1) = diag_script_info@obsfile
        end if

        contrib_authors = (/"A_anav_al", "A_wenz_sa", "A_righ_ma"/)

        ;;NOTE: not implemented in v2 yet
    ;ESMValMD(outfile+"."+file_type, alltags, caption, id, variables, model_info@model, climofiles, diag_script, contrib_authors)

        delete([/alltags, caption, id, climofiles/])

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
            annual_mean(imod) = dim_avg(all_yearly(imod, minyear_ind:maxyear_ind, 0))
            ;; Compute the std over time for the region
            annual_std(imod) = dim_stddev(all_yearly(imod, minyear_ind:maxyear_ind, 0))
        end do

        ;; Collect data
        data_arr = new((/2, dim_MOD/), float)
        data_arr!0 = "statistic"
        data_arr!1 = "model"
        data_arr&statistic = (/"mean", "interannual variability"/)
        data_arr&model = model_info@model
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
        data_arr@res_tiMainString = var0 + " - " + region@name
        plot = scatterplot_markers(wks, data_arr, var0)
        draw(plot)
        frame(wks)
        delete(plot)
        delete(wks)
        delete(data_arr)

        ;; Attach Metadata to plot

        alltags = array_append_record(tags, region@DM_tag, 0)
        alltags := array_append_record(alltags,\
               (/"PT_scatter", "ST_mean"/), 0)
        caption = "Scatterplot for " + variables(0) + \
                  " over the period " + start_year + \
                  "-" + end_year + ", similar to Anav et al., 2013, Fig. 17."
        id = diag_script

        isize = dim_MOD_org * dimsizes(variables)
        if (isatt(diag_script_info, "obsfile")) then
            isize = isize + 1
        end if
        climofiles = new(isize, string)
        do n_var = 0, dimsizes(variables) - 1
            do imod = 0, dim_MOD_org - 1
                climofiles(n_var * dim_MOD_org + imod) = \
                            interface_get_inpaths(imod) + "/" + \
                            interface_get_infile(variables(n_var), field_types(n_var), imod)
            end do
        end do
        if (isatt(diag_script_info, "obsfile")) then
            climofiles(isize - 1) = diag_script_info@obsfile
        end if

        contrib_authors = (/"A_anav_al", "A_wenz_sa", "A_righ_ma"/)

        ;;NOTE: not implemented in v2 yet
    ;ESMValMD(outfile+"."+file_type, alltags, caption, id, variables, model_info@model, climofiles, diag_script, contrib_authors)

        delete([/alltags, caption, id, climofiles/])

    end if
    
    ;; Time evolution plot

    if (isatt(diag_script_info, "evolution_plot")).and. \
                                 diag_script_info@evolution_plot then
        ;; Check if it's supposed to be an anomaly plot
        if (isatt(diag_script_info, "evolution_plot_anomaly")).and. \ 
                         diag_script_info@evolution_plot_anomaly then
            anomaly = True
        else
            anomaly = False
        end if

        outfile = output_dir + diag_script + "_" + var0 + "_evolution_" + \
           region@name

        ;; Find out ref model index
        obs_ind = ind(model_info@model.eq.diag_script_info@evolution_plot_ref_model)

        ;; Set list of observations/model_info to exclude in the statistics
        ex_mod = ind(model_info@project.eq."OBS")
        if (isatt(diag_script_info, "evolution_plot_ignore")) then
            ex_mod1 = array_append_record(ex_mod, ind(model_info@model.eq. \
                                 diag_script_info@evolution_plot_ignore), 0)
            delete(ex_mod)
            ex_mod = ex_mod1
            delete(ex_mod1)
        end if
        ;; To not get a warning about fillvalues from get1Dindex_Exclude
        ind_mods = ispan(0, dimsizes(model_info@model) -1, 1)
        ind_mods@_FillValue = default_fillvalue("integer")
        mod_ind = get1Dindex_Exclude(ind_mods, ex_mod)
        delete([/ex_mod, ind_mods/])
          
        ;;data_arr = all_yearly                   ;; all_yearly:  model, year, std = 0
        data_arr = new((/6, dimsizes(all_yearly&year)/), float)
        model_mean = dim_avg_n(all_yearly(mod_ind,:,0), 0)
        model_std = new((/dimsizes(model_mean)/), float)
        
        data_arr@ref_name = diag_script_info@evolution_plot_ref_model
        data_arr@units = all_yearly@units
        data_arr@dim_Mod = dimsizes(mod_ind)

        ;; Subtract mean 1901-1930 (coded as till 1930) if anomaly plot
        ind_1930 = ind(all_yearly&year.eq.1930)
        
        if anomaly then
            data_arr(0, :) = (/all_yearly(obs_ind, :, 0)/) - \
                              avg(all_yearly(obs_ind, :ind_1930, 0))
            mod_means = dim_avg_n(all_yearly(mod_ind, :ind_1930, 0), 1)
            data_arr(1, :) =(/model_mean/) - avg(mod_means)
            do time_i = 0, dimsizes(all_yearly&year)-1
                data_arr(2, time_i) = dim_min_n(all_yearly(mod_ind,time_i,0) - mod_means, 0)
                data_arr(3, time_i) = dim_max_n(all_yearly(mod_ind,time_i,0) - mod_means, 0)
                model_std(time_i)  = dim_stddev_n(all_yearly(mod_ind,time_i,0) - mod_means, 0)
            end do
        else
            data_arr(1, :) = (/model_mean/)
            data_arr(2, :) = dim_min_n(all_yearly(mod_ind,:,0), 0)
            data_arr(3, :) = dim_max_n(all_yearly(mod_ind,:,0), 0)
            data_arr(0, :) = (/all_yearly(obs_ind, :, 0)/)
            model_std  = dim_stddev_n(all_yearly(mod_ind,:,0), 0)

        end if

        data_arr(4, :) = (/all_yearly&year/)
        data_arr(5, :) = (/model_std/)
        
        ;; Pass on contour color
        if (isatt(diag_script_info, "evolution_plot_color")) then
            data_arr@contour_color = diag_script_info@evolution_plot_color
        end if

        wks = gsn_open_wks(file_type, outfile)
        wks@legendfile = outfile + "_legend"
        plot1 = evolution_base_plot(wks, data_arr, var0, anomaly)

        res = True
        ; After overlaying the plots, add titles and legend manually
        res_text                    = True                  ; text mods desired
        res_text@txFontHeightF      = 0.01                 ; change text size
        res_text@txJust             = "CenterLeft"          ; text justification

        res_lines                   = True                  ; polyline mods desired
        res_lines@gsLineDashPattern = 0.                    ; solid line
        res_lines@gsLineThicknessF  = 5.                    ; line thicker
        res_lines@gsLineColor       = "gray60"               ; line color
        wks_out = get_wks(wks, diag_script, var0)

        ;; Create lines for Volcanic eruptions

        if (isatt(diag_script_info, "evolution_plot_volcanoes")).and. \
            diag_script_info@evolution_plot_volcanoes then
            ymax = max((/max(data_arr(3,:)), max(data_arr(0, :))/))
            ymin = min((/min(data_arr(2,:)), min(data_arr(0, :))/))
            ytop = ymax + 0.1*(ymax-ymin)
            xx = (/1963., 1963./)
            yy = (/(ymin - 0.05 * (ymax - ymin)),(ymax + 0.15*(ymax-ymin))/)
        
            plot1@$unique_string("dum")$ = gsn_add_polyline(wks_out,plot1,xx,yy,res_lines)
            plot1@$unique_string("dum")$ = gsn_add_text(wks_out,plot1,"Agung",\
                                      (xx(1)+0.5),(ytop),res_text)

            xx = (/1982., 1982./)
            plot1@$unique_string("dum")$ = gsn_add_polyline(wks_out,plot1,xx,yy,res_lines)
            plot1@$unique_string("dum")$ = gsn_add_text(wks_out,plot1,"El Chichon",\
                                      (xx(1)),(ytop),res_text)

            xx = (/1991., 1991./)
            plot1@$unique_string("dum")$ = gsn_add_polyline(wks_out,plot1,xx,yy,res_lines)
            plot1@$unique_string("dum")$ = gsn_add_text(wks_out,plot1,"Pinatubo",\
                                      (xx(1)),(ytop),res_text)

            res_lines@gsLineColor = "orange"
            xx = (/1998, 1998./)
            plot1@$unique_string("dum")$ = gsn_add_polyline(wks_out,plot1,xx,yy,res_lines)
            plot1@$unique_string("dum")$ = gsn_add_text(wks_out,plot1,"El Nino",\
                                      (xx(1)+0.5),(ytop),res_text)

        end if
        res_lines@gsLineColor = "black"
        res_lines@gsLineDashPattern = 16
        res_lines@gsLineThicknessF = 3.
        if .not.anomaly then
            plot1@$unique_string("dum")$ = gsn_add_polyline(wks_out, plot1, \ 
                   (/start_year, end_year/), \
                   (/0., 0./), res_lines) 
        end if

        res@gsnDraw          = False
        res@gsnFrame         = False
        res@gsnMaximize      = True
        res@gsnPaperOrientation = "landscape"
        res@xyLineThicknesses = (/2,2/)
        res@xyLineColors      = (/"black","white"/)
        res@xyDashPatterns	 = (/0.,0./)      ; make all lines solid

        plot3 = gsn_csm_xy(wks_out, data_arr(4, :), data_arr(:1, :), res)

        NhlAddOverlay(plot1, plot3, plot1)

        draw(plot1)
        frame(wks)

        ;; Attach Metadata to plot

        alltags = array_append_record(tags, region@DM_tag, 0)
        alltags := array_append_record(alltags,\
               (/"PT_times"/), 0)
        if (anomaly) then
            alltags := array_append_record(alltags,\
                   (/"ST_anomaly"/), 0)
            tagstr = " (anomalies) "
        else
            alltags := array_append_record(alltags,\
                   (/"ST_mean"/), 0)
            tagstr = " "
        end if
        if isatt(diag_script_info, "fignum") then
            fignum = diag_script_info@fignum
        else
            fignum = "1 top"
        end if
        caption = "Time series of " + variables(0) + tagstr + \
                  ", similar to Anav et al., 2013, Fig." + fignum
        id = diag_script

        isize = dim_MOD_org * dimsizes(variables)
        if (isatt(diag_script_info, "obsfile")) then
            isize = isize + 1
        end if
        climofiles = new(isize, string)
        do n_var = 0, dimsizes(variables) - 1
            do imod = 0, dim_MOD_org - 1
                climofiles(n_var * dim_MOD_org + imod) = \
                            interface_get_inpaths(imod) + "/" + \
                            interface_get_infile(variables(n_var), field_types(n_var), imod)
            end do
        end do
        if (isatt(diag_script_info, "obsfile")) then
            climofiles(isize - 1) = diag_script_info@obsfile
        end if

        contrib_authors = (/"A_anav_al", "A_wenz_sa", "A_righ_ma"/)

        ;;NOTE: not implemented in v2 yet
        ;ESMValMD(outfile+"."+file_type, alltags, caption, id, variables, model_info@model, climofiles, diag_script, contrib_authors)

        delete([/alltags, caption, id, climofiles/])

        delete([/data_arr, anomaly/])

    end if

    leave_msg(diag_script, "")

end

