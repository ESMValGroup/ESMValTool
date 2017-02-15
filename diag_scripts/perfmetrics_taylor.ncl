;;#############################################################################
;; DIAGNOSTIC SCRIPT FOR TAYLOR DIAGRAM METRICS
;; Author: Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################
;; Description
;;    Calculates grades according to a given metric, with different option for
;;    normalization. It requires field precalculation by perfmetrics_main.ncl.
;;    To get consistent results, the same cfg_ files should be used for the
;;    two scripts.
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    metric: selected metric
;;    normalization: type of metric normalization
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    region: selected region
;;    MultiModelMedian: calcualte multi-model median
;;    MultiModelMean: calcualte multi-model mean
;;
;; Required variable_info attributes (variable specific)
;;
;; Optional variable_info attributes (variable specific)
;;
;; Caveats
;;    The reference model must be specified in the main namelsit in the 
;;    <variable> tag as <variable ref_model="XX">. An alternative model
;;    can be also given, separated by a comme <variable ref_model="XX,"YY">.
;;
;; Modification history
;;    20160628_A_righ_ma: moving ref_model specification from cfg- files to
;;                        main namelist, some code restructuring necessary
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20150325-A_laue_ax: modified reference tags used for acknowledgements
;;                        (projects, observations, etc)
;;    20150119-A_gott_kl: removed "region" from req_atts
;;    20140701-A_gott_kl: Adapted for T1M
;;    20140630-A_gott_kl: Adapted for T0Ms
;;    20131210-A_fran_fr: written.
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Write references
    write_references(diag_script, \
                     "A_fran_fr", \
                     (/"A_righ_ma", "A_eyri_ve", "A_gott_kl"/), \
                     (/"D_righi15gmd", "D_gleckler08jgr"/), \
                     (/"E_ncep", "E_erainterim", "E_airs", \
                       "E_ceresebaf", "E_srb"/), \
                     (/"P_embrace", "P_esmval"/))
end

begin

    vardeffile = "interface_data/" + var0 + "_info.tmp"
    loadscript(vardeffile)

end

begin

    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Check field type
    if (.not.any(field_type0.eq.(/"T3M", "T2Ms", "T1M", "T0M", "T0Ms"/))) then
        error_msg("f", diag_script, "", "field " + field_type0 + \
                  " not supported")
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"metric", "normalization"/)
    exit_if_missing_atts(diag_script_info, req_atts)

    ;; Check for reference model definition
    if (var_attr_ref(0).eq."None") then
        error_msg("f", diag_script, "", "no reference model is specified")
    end if

    ;; Diag_script basename
    diag_script_base = basename(diag_script)

    ;; Unique names for models
    modelnames = project_style(diag_script_info, "annots")

    ;; Selected region
    if (isatt(diag_script_info, "region")) then
        region = select_region(diag_script_info@region)
    end if

    ;; Annotation & file names
    region_name = ""  ; priority 3
    location = ""  ; priority 3
    if(isatt(diag_script_info, "aux_info")) then
        region_name = "_" + diag_script_info@aux_info  ; priority 2
        location = " - " + diag_script_info@aux_info  ; priority 2
    end if
    if(isatt(diag_script_info, "location")) then
        location = " - " + diag_script_info@location  ; priority 1
    end if
    if(isatt(diag_script_info, "region")) then
        region_name = "_" + region@name  ; priority 1
        location = " - " + diag_script_info@region
    end if
    altitude = ""  ; priority 2
    if(any(field_type0.eq.(/"T0M", "T0Ms"/))) then
        if(isatt(diag_script_info, "altitude")) then
            altitude = " - " + diag_script_info@altitude  ; priority 1
        end if
    else
        if(isatt(diag_script_info, "level")) then
            altitude = " - " + diag_script_info@level + " hPa"  ; priority 1
        end if
    end if

    ;; Set variable name
    var0_new = var0
    var0_new = var0_new + region_name
    if(isatt(diag_script_info, "level").and. \
        any(field_type0.eq.(/"T3M", "T1M"/))) then
        var0_new = var0_new + "-" + diag_script_info@level
    end if

end

begin

    ;; Define taylor file directory
    work_dir = getenv("ESMValTool_wrk_dir")
    xml_name = getenv("ESMValTool_xml_name")
    suffix = get_file_suffix(xml_name, 0)
    name = str_sub_str(xml_name, suffix, "")
    ncdf_dir = get_output_dir(work_dir, diag_script_base) + \
        name + "_taylor_" + var0_new + ".nc"

    ;; Start calculation (if necessary)
    if (isfilepresent_esmval(ncdf_dir) .and. \
        getenv("ESMValTool_force_taylorcalc").eq."False") then

        info_output("No recalculation!", verbosity, 1)
        info_output("Reading in file = " + ncdf_dir, verbosity, 1)

    else

        ;; Read data precalcualted in perfmetrics_main.ncl
        do imod = 0, dim_MOD - 1
            new_path = clim_diag_file_path(var0, field_type0, imod, \
                                           "perfmetrics_main", \
                                           diag_script_info)

            ;; Skip current model if no data available
            if (.not.isfilepresent_esmval(new_path)) then
                info_output("No data available for model:" + \
                            modelnames(imod), verbosity, 1)
                info_output("Skipping", verbosity, 1)
                continue
            end if

            ;; Store into global array
            info_output("Read in: " + new_path, verbosity, 1)
            tmp = ncdf_read(new_path, var0)
            if (.not.isdefined("super_array")) then
                dim_super_array = \
                    array_append_record(dim_MOD, dimsizes(tmp), 0)
                super_array = new(dim_super_array, float)
                super_array!0 = "model"
                super_array&model = models@name
            end if
            rank = dimsizes(dim_super_array)
            if (rank.eq.3) then
                super_array(imod, :, :) = tmp
            end if
            if (rank.eq.4) then
                super_array(imod, :, :, :) = tmp
            end if
            if (rank.eq.5) then
                super_array(imod, :, :, :, :) = tmp
            end if
            delete(tmp)
        end do

        ;; Get the index of the reference/alternative model
        refs_ind = get_ref_model_idx(models, var_attr_ref(0))
        l_altern = dimsizes(refs_ind).gt.1
        if (l_altern) then
            info_output("More than one reference model is selected," + \
                        " only the first one will be used.", verbosity, 2)
        end if

        ;; Calculate multi-model mean/median
        models_only = remove_index(super_array, refs_ind(0))
        if (isatt(diag_script_info, "MultiModelMedian")) then
            if (diag_script_info@MultiModelMean) then
                if (rank.eq.3) then
                    median = models_only(0:0, :, :)
                    median(0, :, :) = dim_median_n(models_only, 0)
                end if
                if (rank.eq.4) then
                    median = models_only(0:0, :, :, :)
                    median(0, :, :, :) = dim_median_n(models_only, 0)
                end if
                if (rank.eq.5) then
                    median = models_only(0:0, :, :, :, :)
                    median(0, :, :, :, :) = dim_median_n(models_only, 0)
                end if
                median&model = "Median-model"
                super_array := array_append_record(median, super_array, 0)
                delete(median)
                refs_ind = refs_ind + 1  ; update ref model index
            end if
        end if
        if (isatt(diag_script_info, "MultiModelMean")) then
            if (diag_script_info@MultiModelMean) then
                if (rank.eq.3) then
                    mean = models_only(0:0, :, :)
                    mean(0, :, :) = dim_avg_n(models_only, 0)
                end if
                if (rank.eq.4) then
                    mean = models_only(0:0, :, :, :)
                    mean(0, :, :, :) = dim_avg_n(models_only, 0)
                end if
                if (rank.eq.5) then
                    mean = models_only(0:0, :, :, :, :)
                    mean(0, :, :, :, :) = dim_avg_n(models_only, 0)
                end if
                mean&model = "Mean-model"
                super_array := array_append_record(mean, super_array, 0)
                delete(mean)
                refs_ind = refs_ind + 1  ; update ref model index
            end if
        end if
        delete(models_only)

        ;; Create separate arrays for models and observations
        models_only = remove_index(super_array, refs_ind(0))
        nmodels = dimsizes(models_only&model)
        if (rank.eq.3) then
            obs_only = super_array(refs_ind(0), :, :)
        end if
        if (rank.eq.4) then
            obs_only = super_array(refs_ind(0), :, :, :)
        end if
        if (rank.eq.5) then
            obs_only = super_array(refs_ind(0), :, :, :, :)
        end if
        delete(super_array)

        ;; Define result variable
        ;; A dummy coordinate "diagnostics" (size 1) is added to facilitate
        ;; appending different variable in the _collect script
        val = new((/dim_VAR, nmodels, 2/), float)
        val!0 = "diagnostics"
        val!1 = "models"
        val!2 = "statistic"
        val&diagnostics = var0_new
        val&models = models_only&model
        val&statistic = (/"stddev_ratio", "correlation"/)
        val@ncdf_dir = ncdf_dir

        ;; Loop over models, calculate Taylor metrics
        if (rank.eq.3) then
            obs = obs_only(:, 0)
        end if
        if (rank.eq.4) then
            obs = obs_only(:, :, 0)
        end if
        if (rank.eq.5) then
            obs = obs_only(:, :, :, 0)
        end if
        do imod = 0, nmodels - 1
            if (rank.eq.3) then
                model = models_only(imod, :, 0)
            end if
            if (rank.eq.4) then
                model = models_only(imod, :, :, 0)
            end if
            if (rank.eq.5) then
                model = models_only(imod, :, :, :, 0)
            end if
            val(0, imod, 0) = \
                calculate_metric(model, obs, "stddev_ratio")
            val(0, imod, 1) = \
                calculate_metric(model, obs, "correlation")
        end do

        ;; Attach attributes to the results
        val@title = "taylor metrics"
        val@long_name = "Taylor diagram"
        val@diag_script = (/diag_script/)
        val@var = "var"

        ;; Write NetCDF output
        ncdf_outfile = ncdf_write(val, val@ncdf_dir)

    end if  ; if data available

    ;; Write results of temporary grading list
    xml_name = getenv("ESMValTool_xml_name")
    suffix = get_file_suffix(xml_name, 0)
    name = str_sub_str(xml_name, suffix, "")
    temp_dir = get_output_dir(getenv("ESMValTool_wrk_dir"), "temp") + \
        name + "_taylor.nc"

    if (isfilepresent_esmval(temp_dir)) then
        ; Read grading list
        temp_file = addfile(temp_dir, "r")
        temp_list = temp_file->temp_list
        temp_list_string = tostring(temp_list)
        ; Add new entry
        temp_list_new_string = array_append_record(temp_list_string, \
                                                   ncdf_dir, 0)
        temp_list_new = tochar(temp_list_new_string)
        ; Remove previous file
        system("rm -f " + temp_dir)
    else
        ncdf_char = tochar(ncdf_dir)
        temp_list_new = new((/1, dimsizes(ncdf_char)/), character)
        temp_list_new(0, :) = ncdf_char
    end if

    ;; Create new file and add list
    temp = addfile(temp_dir, "c")
    temp->temp_list = temp_list_new

    leave_msg(diag_script, "", 4)

end
