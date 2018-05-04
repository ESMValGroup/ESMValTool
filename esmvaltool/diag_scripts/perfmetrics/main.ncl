;;#############################################################################
;; MAIN SCRIPT FOR PERFORMANCE METRICS
;; Authors: Mattia Righi (DLR, Germany) and Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################
;; Description
;;    Calculates and (optionally) plots annual/seasonal cycles, zonal means,
;;    lat-lon fields and time-lat-lon fields from input T3M or T2Ms data.
;;    The calculated fields can be also plotted as difference w.r.t. a given
;;    reference model. It also calculates grading and taylor metrics.
;;    Input data have to be regridded to a common grid in the preprocessor.
;;
;; Required diag_script_info attributes
;;    plot_type: cycle (time), zonal (plev, lat), latlon (lat, lon) or
;;               cycle_latlon (time, lat, lon)
;;    time_avg: type of time average (see time_operations in
;;              diag_scripts/shared/statistics.ncl)
;;    region: selected region (see select_region in
;;            diag_scripts/shared/latlon.ncl)
;;    styleset (for cycle): as in diag_scripts/shared/plot/style.ncl functions
;;    plot_stddev (for cycle): plot standard deviation
;;    legend_outside (for cycle): save legend in a separate file
;;
;; Optional diag_script_info attributes
;;    t_test (for zonal and latlon): calculate t-test in difference plots
;;                                   (default: False)
;;    conf_level (for zonal and latlon): confidence level for the t-test
;;                                       (default: False)
;;    range_option: time range selection option (default: 0)
;;    projection: map projection for lat-lon plots (default:
;;                CylindricalEquidistant)
;;    draw_plots: draw plots (default: True)
;;    plot_diff: draw difference plots (default: False)
;;    calc_grading: calculate grading (default: False)
;;    stippling: use stippling to mark stat. significant differences (default:
;;               False = mask out non-significant differences in gray)
;;    show_global_avg: diplay global avaerage as right string on lat-lon plots
;;                     (default: False)
;;    metric: grading metric (if calc_grading is True)
;;    normalization: metric normalization (for RMSD and BIAS metrics)
;;    abs_levs: (min, max, spacing) contour levels for absolute plot
;;    diff_levs: (min, max, spacing) contour levels for difference plot
;;    zonal_cmap (for zonal): color table (default: "amwg_blueyellowred")
;;    zonal_ymin (for zonal): minimum pressure on the plots (default: 5. hPa)
;;    latlon_cmap (for latlon): color table (default: "amwg_blueyellowred")
;;    plot_units: plotting units (if different from standard CMOR units)
;;    aux_info: additional information for stations data
;;    location: additional information for stations data (location)
;;    altitude: additional information for stations data (altitude)
;;
;; Caveats
;;
;; Modification history
;;    20180503-A_righ_ma: completely rewritten and modularized
;;    20171215-A_righ_ma: merged with perfmetrics_grading and
;;                        permetrics_taylor.ncl
;;    20171124-A_righ_ma: completely revised to adapt it to the new backend
;;                        (level selection, regridding and masking now done
;;                        by the python preprocessor)
;;    20161220-A_laue_ax: added option to set map projection for lat-lon plots
;;                        (diag_script_info@projection)
;;                        added option to choose how to plot t-test results:
;;                        stippling or masking out in gray (lat-lon plots only)
;;    20161019-A_laue_ax: changed plotting of t-test results:
;;                        now stippling significant grid cells (old version:
;;                        masking out non-significant values in gray)
;;    20160628-A_righ_ma: moving ref_model specification from cfg- files to
;;                        main namelist
;;    20160628-A_senf_da: added regridding for irregular grids (ESMF_regrid).
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code.
;;    20151013-A_righ_ma: fixed t-test mask in lat-lon difference plots.
;;    20150325-A_laue_ax: modified reference tags used for acknowledgements
;;                        (projects, observations, etc.).
;;    20150119-A-gott_kl: removed "grid", "region" from req_atts
;;                        (for T2Ms vmrco).
;;    20150113-A_gott_kl: reconciled generalised regridding with T1* & T0*
;;    20140905-A_righ_ma: consistent regridding and missing values mask.
;;    20140701-A_gott_kl: Adapted for T1M.
;;    20140630-A_gott_kl: Adapted for T0Ms.
;;    20131203-A_fran_fr: written.
;;
;;#############################################################################

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/latlon.ncl"
load "./diag_scripts/shared/statistics.ncl"
load "./diag_scripts/shared/regridding.ncl"
load "./diag_scripts/shared/ensemble.ncl"
load "./diag_scripts/shared/scaling.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/xy_line.ncl"
load "./diag_scripts/shared/plot/zonalmean_profile.ncl"
load "./diag_scripts/shared/plot/contour_maps.ncl"

begin

    enter_msg(diag_script, "")
    log_info("++++++++++++++++++++++++++++++++++++++++++")
    log_info(diag_script + " (var: " + variables(0) + ")")
    log_info("++++++++++++++++++++++++++++++++++++++++++")

    dim_MOD = dimsizes(model_info@model)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Write references
    write_references(diag_script, \
                     "A_fran_fr", \
                     (/"A_righ_ma", "A_eyri_ve", "A_gott_kl", "A_senf_da"/), \
                     (/"D_righi15gmd", "D_gleckler08jgr"/), \
                     (/"E_ncep", "E_erainterim", "E_airs", \
                       "E_ceresebaf", "E_srb"/), \
                     (/"P_embrace", "P_esmval"/))

    ;; Check required diag_script_info attributes
    exit_if_missing_atts(diag_script_info, \
                         (/"plot_type", "time_avg", "region"/))

    ;; Store required attributes
    ptype = diag_script_info@plot_type

    if (all(ptype.ne.(/"cycle", "zonal", "latlon", "cycle_latlon"/))) then
        error_msg("f", diag_script, "", "plot_type " + ptype + " is not a " + \
                  "supported plot_type in this diagnostic")
    end if

    ;; Check for plot-type specific settings
    if (ptype.eq."cycle") then
        exit_if_missing_atts(diag_script_info, \
                             (/"legend_outside", "styleset", "plot_stddev"/))
    end if

    if ((ptype.eq."zonal" .or. ptype.eq."latlon") .and. \
        diag_script_info@region.ne."Global") then
        error_msg("f", diag_script, "", "plot_type " + ptype + \
                  " implemented only for region='Global'")
    end if

    ;; Check valid field
    if ((ptype.eq."zonal" .and. \
         all(field_type0.ne.(/"T3M", "T2Mz"/))) .or. \
        (ptype.eq."cycle" .and. \
         all(field_type0.ne.(/"T3M", "T2Ms", "T1M", "T0Ms"/))) .or. \
        (ptype.eq."latlon" .and. \
         all(field_type0.ne.(/"T3M", "T2Ms"/))) .or. \
        (ptype.eq."cycle_latlon" .and. \
         all(field_type0.ne.(/"T3M", "T2Ms"/)))) then
        error_msg("f", diag_script, "", "input field " + field_type0 + \
                  " is not compatible with plot_type " + ptype)
    end if

    ;; Set default values for non-required diag_script_info attributes
    set_default_att(diag_script_info, "range_option", 0)
    set_default_att(diag_script_info, "projection", "CylindricalEquidistant")
    set_default_att(diag_script_info, "draw_plots", True)
    set_default_att(diag_script_info, "plot_diff", False)
    set_default_att(diag_script_info, "calc_grading", False)
    set_default_att(diag_script_info, "stippling", False)
    set_default_att(diag_script_info, "t_test", False)
    set_default_att(diag_script_info, "show_global_avg", False)
    set_default_att(diag_script_info, "zonal_ymin", 5.)
    set_default_att(diag_script_info, "zonal_cmap", "amwg_blueyellowred")
    set_default_att(diag_script_info, "latlon_cmap", "amwg_blueyellowred")

    ;; Check consistency of diff plots settings
    if (diag_script_info@t_test .and. .not.diag_script_info@plot_diff) then
        error_msg("f", diag_script, "", \
                  "plot_diff must be True to apply t-test")
    end if
    if (diag_script_info@t_test .and. .not.diag_script_info@conf_level) then
        error_msg("f", diag_script, "", \
                  "conf_level must be specified to apply t-test")
    end if

    ;; Check metric
    if (diag_script_info@calc_grading) then
        exit_if_missing_atts(diag_script_info, (/"metric"/))
        if (diag_script_info@metric.ne."taylor") then
            exit_if_missing_atts(diag_script_info, (/"normalization"/))
        end if
    end if

    ;; Unique names for models
    modelnames = project_style(diag_script_info, "annots")

    ;; Check for reference model definition
    if (variable_info@reference_model.eq."None") then
        error_msg("f", diag_script, "", "no reference model is specified")
    end if

    ;; Set index of the reference (and alternative) model
    ref_inds = get_ref_model_idx(model_info, \
                                 (/variable_info@reference_model, \
                                   variable_info@alternative_model/))
    ref_ind = ref_inds(0)
    if (dimsizes(ref_inds).gt.1) then
        alt_ind = ref_inds(1)
        l_altern = True
    end if

    ;; Create output plot directory
    if (diag_script_info@draw_plots) then
    	plot_dir = config_user_info@plot_dir
        system("mkdir -p " + plot_dir)
    end if

    ;; Plot file type
    file_type = config_user_info@output_file_type
    if (ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Load plot-type-specific script
    print("diag_scripts/perfmetrics/" + ptype + ".ncl")
    loadscript("diag_scripts/perfmetrics/" + ptype + ".ncl")

end

begin

    perfmetrics_ptype_script()
    leave_msg(diag_script, "")

end

;;     ;; ### DEFINE GLOBAL ARRAYS FOR COLLECTING MODEL DATA #####################

;;     ;; Read common grid from the first model
;;     ;; (all models are expected to be regridded to a common grid)
;;     data = read_data(0, var0, field_type0)
;;     rank = dimsizes(dimsizes(data))
;;     if (rank.eq.4) then
;;         grid = data(0, :, :, :)  ; plev, lat, lon
;;     end if
;;     if (rank.eq.3) then
;;         grid = data(0, :, :)  ; lat, lon
;;     end if
;;     delete(data)
;;     delete(rank)


;;     ;; Case cycle_latlon: (model, time, lat, lon, stats)
;;     if (ptype.eq."cycle_latlon") then
;;         sgrid = area_operations(grid, region(0), region(1), region(2), \
;;                                 region(3), "extract", False)
;;         if (time_avg.eq."seasonalclim") then
;;             var_all = new((/dim_MOD, 4, dimsizes(sgrid&lat), \
;;                           dimsizes(sgrid&lon), 2/), float)
;;             var_all!1 = "season"
;;             var_all&season = (/"DJF", "MAM", "JJA", "SON"/)
;;         else if (time_avg.eq."monthlyclim") then
;;             var_all = new((/dim_MOD, 12, dimsizes(sgrid&lat), \
;;                           dimsizes(sgrid&lon), 2/), float)
;;             var_all!1 = "month"
;;             var_all&month = (/"J", "F", "M", "A", "M", "J",\
;;                               "J", "A", "S", "O", "N", "D"/)
;;         else
;;             error_msg("f", diag_script, "", "time_avg option " + time_avg + \
;;                       " not compatible with plot_type " + ptype)
;;         end if
;;         end if
;;         var_all!0 = "model"
;;         var_all!2 = "lat"
;;         var_all!3 = "lon"
;;         var_all!4 = "statistic"
;;         var_all&model = modelnames
;;         var_all&lat = sgrid&lat
;;         var_all&lon = sgrid&lon
;;         var_all&statistic = (/"mean", "stddev"/)
;;         if (diag_script_info@plot_diff) then
;;             var_diff = var_all(:, :, :, :, 0)
;;         end if
;;         delete(sgrid)
;;     end if

;;     ;; Case zonal: (model, level, lat, stats)
;;     if (ptype.eq."zonal") then
;;         var_all = new((/dim_MOD, dimsizes(grid&plev), \
;;                       dimsizes(grid&lat), 2/), float)
;;         var_all!0 = "model"
;;         var_all!1 = "plev"
;;         var_all!2 = "lat"
;;         var_all!3 = "statistic"
;;         var_all&model = modelnames
;;         var_all&plev = grid&plev
;;         var_all&lat = grid&lat
;;         var_all&statistic = (/"mean", "stddev"/)
;;         if (diag_script_info@plot_diff) then
;;             var_diff = var_all(:, :, :, 0)
;;         end if
;;     end if

;;     ;; Case latlon: (model, lat, lon, stats)
;;     if (ptype.eq."latlon") then
;;         var_all = new((/dim_MOD, dimsizes(grid&lat), \
;;                       dimsizes(grid&lon), 2/), float)
;;         var_all!0 = "model"
;;         var_all!1 = "lat"
;;         var_all!2 = "lon"
;;         var_all!3 = "statistic"
;;         var_all&model = modelnames
;;         var_all&lat = grid&lat
;;         var_all&lon = grid&lon
;;         var_all&statistic = (/"mean", "stddev"/)
;;         if (diag_script_info@plot_diff) then
;;             var_diff = var_all(:, :, :, 0)
;;         end if
;;     end if

;;     ; If needed, create array for statistically significant grid cells
;;     if (diag_script_info@t_test) then
;;         var_sig = var_diff
;;     end if

;;     ;; Variable rank: cycle:3, zonal/latlon:4, cycle_latlon:5
;;     irank = dimsizes(dimsizes(var_all))

;;     ;; Attach attributes
;;     var_all@var = var0
;;     var_all@diag_script = (/diag_script/)
;;     copy_VarAtts(diag_script_info, var_all)
;;     var_all@ref_model = variable_info@reference_model

;;     ;; ### MAIN LOOP OVER MODELS AND DATA PROCESSING ##########################

;;     ;; Loop over models
;;     do imod = 0, dim_MOD - 1

;;         log_debug("Processing " + modelnames(imod))

;;         ;; Set path for saving processed data ;; FIX-ME add preproc_id
;;         fullpath = config_user_info@work_dir + "/"
;;         system("mkdir -p " + fullpath)
;;         fname = basename(systemfunc("basename " + model_info@filename(imod)))
;;         fname = fname + "_" + basename(diag_script)
;;         fname = fname + "_" + diag_script_info@plot_type
;;         fname = fname + "_" + diag_script_info@time_avg
;;         if (isatt(diag_script_info, "level")) then
;;             fname = fname + "_" + diag_script_info@level  ; FIX-ME
;;         end if
;;         fname = fname + "_" + str_sub_str(diag_script_info@region, " ", "")
;;         procpath = fullpath + fname + ".nc"

;;         ;; Determine start/end year
;;         start_year = get_start_year(imod)
;;         end_year = get_end_year(imod)

;;         ;; Read data
;;         var = read_data(imod, var0, field_type0)
;;         dnames = getVarDimNames(var)


;;         ;; CASE: cycle_latlon #################################################
;;         if (ptype.eq."cycle_latlon") then

;;             ;; Extract region
;;             var_reg = area_operations(var, region(0), region(1), \
;;                                       region(2), region(3), "extract", True)
;;             delete(var)

;;             ;; Calculate time average
;;             var_avg = time_operations(var_reg, start_year, end_year, \
;;                                       "average", time_avg, True)
;;             delete(var_reg)

;;             ;; Standard deviation calculation for this rank is not implemented
;;             ;; yet in statistics.ncl (but is anyway not required here)

;;             ;; Assign to global array
;;             var_all(imod, :, :, :, 0) = var_avg
;;             delete(var_avg)
;;             var = var_all(imod, :, :, :, :)

;;         end if



;;     ;; ### GRADING ############################################################

;;     if (diag_script_info@calc_grading) then  ; FIX-ME loop over metrics?

;;         ;; Annotation & file names
;;         region_name = ""  ; priority 3
;;         location = ""  ; priority 3
;;         if(isatt(diag_script_info, "aux_info")) then
;;             region_name = "_" + diag_script_info@aux_info  ; priority 2
;;             location = " - " + diag_script_info@aux_info  ; priority 2
;;         end if
;;         if(isatt(diag_script_info, "location")) then
;;             location = " - " + diag_script_info@location  ; priority 1
;;         end if
;;         if(isatt(diag_script_info, "region")) then
;;             region_name = "_" + region@name  ; priority 1
;;             location = " - " + diag_script_info@region
;;         end if
;;         altitude = ""  ; priority 2
;;         if(any(field_type0.eq.(/"T0M", "T0Ms"/))) then
;;             if(isatt(diag_script_info, "altitude")) then
;;                 altitude = " - " + diag_script_info@altitude  ; priority 1
;;             end if
;;         else
;;             if(isatt(diag_script_info, "level")) then
;;                 altitude = " - " + diag_script_info@level + " hPa"  ; priority1
;;             end if
;;         end if

;;         ;; Set variable name
;;         var0_new = var0
;;         var0_new = var0_new + region_name
;;         if(isatt(diag_script_info, "level").and. \
;;             any(field_type0.eq.(/"T3M", "T1M"/))) then
;;             var0_new = var0_new + "-" + diag_script_info@level
;;         end if

;;         ;; Define grading filename
;;         ncdf_dir = config_user_info@work_dir + "/" + \
;;             diag_script_info@metric + "_" + var0_new + ".nc"

;;         ;; Create separate arrays for models and observations
;;         if (diag_script_info@metric.eq."taylor") then
;;             id_remove = ref_ind
;;         else
;;             id_remove = ref_inds
;;         end if
;;         super_array = var_all
;;         models_only = remove_index(super_array, id_remove)
;;         nmodels = dimsizes(models_only&model)
;;         if (l_altern .and. diag_script_info@metric.ne."taylor") then
;;             if (irank.eq.3) then
;;                 obs_only = super_array(id_remove, :, :)
;;             end if
;;             if (irank.eq.4) then
;;                 obs_only = super_array(id_remove, :, :, :)
;;             end if
;;             if (irank.eq.5) then
;;                 obs_only = super_array(id_remove, :, :, :, :)
;;             end if
;;         else  ; avoid dimension reduction
;;             if (irank.eq.3) then
;;                 obs_only = super_array(id_remove:id_remove, :, :)
;;             end if
;;             if (irank.eq.4) then
;;                 obs_only = super_array(id_remove:id_remove, :, :, :)
;;             end if
;;             if (irank.eq.5) then
;;                 obs_only = super_array(id_remove:id_remove, :, :, :, :)
;;             end if
;;         end if
;;         nobs = dimsizes(obs_only&model)
;;         delete(super_array)

;;         ;; Define result variable
;;         ;; A dummy coordinate "diagnostics" (size 1) is added to facilitate
;;         ;; appending different variable in the _collect script
;;         if (diag_script_info@metric.eq."taylor") then
;;             grading = new((/dim_VAR, nmodels, 2/), float)
;;             grading!2 = "statistic"
;;             grading&statistic = (/"stddev_ratio", "correlation"/)
;;         else
;;             grading = new((/dim_VAR, nmodels, nobs/), float)
;;             grading!2 = "reference"
;;             grading&reference = obs_only&model
;;         end if
;;         grading!0 = "diagnostics"
;;         grading!1 = "models"
;;         grading&diagnostics = var0_new
;;         grading&models = models_only&model

;;         grading@ncdf_dir = ncdf_dir

;;         ;; Loop over models, calculate metric
;;         do imod = 0, nmodels - 1
;;             do iobs = 0, nobs - 1
;;                 if (irank.eq.3) then
;;                     model = models_only(imod, :, 0)
;;                     obs = obs_only(iobs, :, 0)
;;                 end if
;;                 if (irank.eq.4) then
;;                     model = models_only(imod, :, :, 0)
;;                     obs = obs_only(iobs, :, :, 0)
;;                 end if
;;                 if (irank.eq.5) then
;;                     model = models_only(imod, :, :, :, 0)
;;                     obs = obs_only(iobs, :, :, :, 0)
;;                 end if
;;                 if (diag_script_info@metric.eq."taylor") then
;;                     grading(0, imod, 0) = \
;;                         calculate_metric(model, obs, "stddev_ratio")
;;                     grading(0, imod, 1) = \
;;                         calculate_metric(model, obs, "correlation")
;;                 else
;;                     grading(0, imod, iobs) = \
;;                         calculate_metric(model, obs, diag_script_info@metric)
;;                 end if
;;             end do
;;         end do

;;         ;; Apply normalization
;;         if (diag_script_info@metric.ne."taylor") then
;;             do iobs = 0, nobs - 1
;;                 grading(:, :, iobs) = \
;;                     normalize_metric(grading(:, :, iobs), \
;;                                      diag_script_info@normalization)
;;             end do
;;         end if

;;         ;; Reduce dimensionality if no alternative model
;;         if (.not.l_altern .and. diag_script_info@metric.ne."taylor") then
;;             grading := grading(:, :, 0)
;;             delete(grading@reference)
;;         end if

;;         ;; Attach attributes to the results
;;         grading@title = diag_script_info@metric + "metric"
;;         grading@long_name = \
;;             "Grading table of metric " + diag_script_info@metric
;;         grading@metric = diag_script_info@metric
;;         grading@diag_script = (/diag_script/)
;;         grading@var = "grade"  ; FIX-ME ?
;;         grading@region = location
;;         grading@num_climofiles = dimsizes(model_info@model)  ; for tagging
;;         do imod = 0, dimsizes(model_info@model) - 1
;;             num_climo = "climofile_" + imod
;;             grading@$num_climo$ = input_file_info@filename(imod)
;;         end do

;;         ;; Write NetCDF output
;;         ncdf_outfile = ncdf_write(grading, grading@ncdf_dir)

;;         ;; Write results of temporary grading list
;;         temp_dir = config_user_info@work_dir + "/"
;;         if (diag_script_info@metric.eq."taylor") then
;;             temp_dir = temp_dir + "taylor.nc"
;;         else
;;             temp_dir = temp_dir + "grading.nc"
;;         end if

;;         if (fileexists(temp_dir)) then
;;             temp_file = addfile(temp_dir, "r")
;;             temp_list = temp_file->temp_list
;;             temp_list_string = tostring(temp_list)
;;             temp_list_new_string = array_append_record(temp_list_string, \
;;                                                        ncdf_dir, 0)
;;             temp_list_new = tochar(temp_list_new_string)
;;             system("rm -f " + temp_dir)
;;         else
;;             ncdf_char = tochar(ncdf_dir)
;;             temp_list_new = new((/1, dimsizes(ncdf_char)/), character)
;;             temp_list_new(0, :) = ncdf_char
;;         end if

;;         ;; Create new file and add list
;;         temp = addfile(temp_dir, "c")
;;         temp->temp_list = temp_list_new

;;     end if


