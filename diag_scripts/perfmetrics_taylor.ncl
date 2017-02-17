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
;;    region: selected region
;;    time_avg: type of time average
;;    metric: selected metric
;;    normalization: type of metric normalization
;;    ref_model: reference model used for comparison (optionally two models can
;;               be provided in an array, the second is used as alternative)
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
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output(diag_script + " (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = 1
    var0 = variables(0)
    field_type0 = field_types(0)

    ;; Write references
    write_references(diag_script, \
                     "A_fran_fr", \
                     (/"A_righ_ma", "A_eyri_ve", "A_gott_kl"/), \
                     "D_righi15gmd", \
                     (/"E_kalnay96bams", "E_erainterim", "E_airs", "E_ceresebaf", "E_srb"/), \
                     (/"P_embrace", "P_esmval"/))
end

begin

    vardeffile = "interface_data/" + variables(0) + "_info.tmp"
    loadscript(vardeffile)

end

begin

    ;; Get environment variables
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "PS"
    end if

    ;; Check plot type
    if (.not.isatt(diag_script_info, "plot_type")) then
        error_msg("f", diag_script, "", "plot_type not specified")
    else
        ptype = diag_script_info@plot_type
    end if

    ;; Check field type
    if (.not.any(field_type0.eq.(/"T3M", "T2Ms", "T1M", "T0M", "T0Ms"/))) then
        error_msg("f", diag_script, "", "field " + field_type0 + \
                  " not supported")
    end if

    ;; Check required diag_script_info attributes
    req_atts = (/"time_avg", "metric", "normalization", "ref_model"/)  ; region
    if (.not.all(isatt(diag_script_info, req_atts))) then
        error_msg("f", diag_script, "", "missing required " + \
                  "diag_script_info attribute:" + \
                  req_atts(ind(.not.isatt(diag_script_info, req_atts))))
    end if

    ;; Check for reference model
    ref_model = diag_script_info@ref_model(0)
    if (dimsizes(diag_script_info@ref_model).gt.1) then
        info_output("More than one reference model is selected" + \
                    "Only the first one will be used.", verbosity, 2)
        delete(diag_script_info@ref_model)
        diag_script_info@ref_model = ref_model
    end if

    ;; Diag_script basename
    diag_script_base = basename(diag_script)

    ;; Unique names for models
    modelnames = project_style(diag_script_info, "annots")

    ;; Check time average method
    if (.not.(isatt(diag_script_info, "time_avg"))) then
        error_msg("f", diag_script, "", "no method for time average provided")
    end if

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

    ;; Check for multimodel median
    flag_median = False
    if (isatt(diag_script_info, "MultiModelMedian")) then
        if (diag_script_info@MultiModelMean) then
            flag_median = True
            dim_MOD = dim_MOD + 1
            modelnames_new = array_append_record((/"Median-model"/),\
                                                 modelnames, 0)
            modelsatname = array_append_record((/"Median-model"/), \
                                               models@name, 0)
            delete(models@name)
            delete(modelnames)
            models@name = modelsatname
            modelnames = modelnames_new
            delete(modelnames_new)
            delete(modelsatname)
        end if
    end if

    ;; Check for multimodel mean
    flag_mean = False
    if (isatt(diag_script_info, "MultiModelMean")) then
        if (diag_script_info@MultiModelMean) then
            flag_mean = True
            dim_MOD = dim_MOD + 1
            modelnames_new = array_append_record((/"Mean-model"/), \
                                                 modelnames, 0)
            modelsatname = array_append_record((/"Mean-model"/), \
                                               models@name, 0)
            delete(models@name)
            delete(modelnames)
            models@name = modelsatname
            modelnames = modelnames_new
            delete(modelnames_new)
            delete(modelsatname)
        end if
    end if

    ;; Start calculation (if necessary)
    if (isfilepresent_esmval(ncdf_dir) .and. \
        getenv("ESMValTool_force_taylorcalc").eq."False") then

        info_output("No recalculation!", verbosity, 1)
        info_output("Reading in file = " + ncdf_dir, verbosity, 1)

    else

        ;; Read data
        do imod = 0, dim_MOD - 1

            if (flag_mean.and.flag_median) then
                if (any(imod.eq.(/0, 1/))) then
                    new_path = "none"
                else
                    new_path = \
                        clim_diag_file_path(var0, field_type0, imod - 2, \
                                            "perfmetrics_main", \
                                            diag_script_info)
                end if
            else if (flag_mean.or.flag_median) then
                if (imod.eq.0) then
                    new_path = "none"
                else
                    new_path = \
                        clim_diag_file_path(var0, field_type0, imod - 1, \
                                            "perfmetrics_main", \
                                            diag_script_info)
                end if
            else
                new_path = \
                    clim_diag_file_path(var0, field_type0, imod, \
                                        "perfmetrics_main", \
                                        diag_script_info)
            end if
            end if

            if (isfilepresent_esmval(new_path)) then
                info_output("Read in: " + new_path, verbosity, 1)
                data_temp = ncdf_read(new_path, var0)
                if (.not.isdefined("data")) then
                    dim_data = \
                        array_append_record((/dim_MOD/), \
                                            dimsizes(data_temp), 0)
                    data = new(dim_data, float)
                    data!0 = "model"
                    data&model = modelnames
                end if
                rank = dimsizes(dim_data)

                if (rank.eq.3) then
                    data(imod, :, :) = data_temp
                end if

                if (rank.eq.4) then
                    data(imod, :, :, :) = data_temp
                end if

                if (rank.eq.5) then
                    data(imod, :, :, :, :) = data_temp
                end if
                delete(data_temp)

            else

                ;; Skip this model is no data are available
                info_output("No data available for model:" + \
                            modelnames(imod), verbosity, 1)
                info_output("Skipping", verbosity, 1)

            end if
        end do

        ;; Get the index of the reference model
        ref_ind = ind(data&model.eq.ref_model)
        if (all(ismissing(ref_ind))) then
            error_msg("f", diag_script, "", "no data available for the " + \
                      "reference model")
        end if

        ;; Number of used references
        dim_REF = dimsizes(diag_script_info@ref_model)

        ;; Define result variable
        ;; FIXME: do not need for dim_VAR dimension here, as
        ;;        we always use 1 variable
        val = new((/dim_VAR, dim_MOD - dim_REF, 2/), float)
        val!0 = "diagnostics"
        val!1 = "models"
        val!2 = "statistic"
        val&diagnostics = var0_new
        val&models = modelnames(ind((models@name.ne.ref_model)))
        val&statistic = (/"stddev_ratio", "correlation"/)
        val@ncdf_dir = ncdf_dir

        ;; Extract data
        if (ptype.eq."cycle") then
            ref = data(ref_ind, :, 0)
        end if
        if (any(ptype.eq.(/"latlon", "zonal"/))) then
            ref = data(ref_ind, :, :)
        end if
        if (ptype.eq."cycle_latlon") then
            ref = data(ref_ind, :, :, :, 0)
        end if

        ;; Loop over models
        ;; (second model index, for models only, omit references)
        mod_i = 0
        do imod = 0, dim_MOD - 1

            ;; Skip reference
            if (imod.eq.ref_ind) then
                print(" Refmodel = " + data&model(imod))
                continue
            end if

            ;; Extract data

            ;; Multimodel mean
            if (flag_mean.and.imod.eq.0) then
                model_ind_woref = ind(ref_model.ne.data&model)

                if (ptype.eq."cycle") then
                    var = dim_avg_n_Wrap( \
                        data(model_ind_woref, :, 0), 0)
                end if
                if (any(ptype.eq.(/"latlon", "zonal"/))) then
                    var = dim_avg_n_Wrap( \
                        data(model_ind_woref, :, :), 0)
                end if
                if (ptype.eq."cycle_latlon") then
                    var = dim_avg_n_Wrap( \
                        data(model_ind_woref, :, :, :, 0), 0)
                end if

            ;; Multimodel median
            else if ((flag_mean.and.(flag_median.and.imod.eq.1)).or.\
                     (.not.flag_mean.and.(flag_median.and.imod.eq.0))) then

                model_ind_woref = ind(ref_model.ne.data&model)
                if (ptype.eq."cycle") then
                    var = dim_median_n( \
                        data(model_ind_woref, :, 0), 0)
                        copy_VarMeta(data(0, :, 0), var)
                end if
                if (any(ptype.eq.(/"latlon", "zonal"/))) then
                    var = dim_median_n( \
                        data(model_ind_woref, :, :), 0)
                        copy_VarMeta(data(0, :, :), var)
                end if
                if (ptype.eq."cycle_latlon") then
                    var = dim_median_n( \
                        data(model_ind_woref, :, :, :, 0), 0)
                        copy_VarMeta(data(0, :, :, :, 0), var)
                end if

            else  ; Given model
                if (ptype.eq."cycle") then
                    var = data(imod, :, 0)
                end if
                if (any(ptype.eq.(/"latlon", "zonal"/))) then
                    var = data(imod, :, :)
                end if
                if (ptype.eq."cycle_latlon") then
                    var = data(imod, :, :, :, 0)
                end if
            end if
            end if

            ;; Calculate rate of standard deviations
            val(0, mod_i, 0) = calculate_metric(var, ref, \
                                                "stddev_ratio")

            ;; Calculate pattern correlation
            val(0, mod_i, 1) = calculate_metric(var, ref, \
                                                "correlation")

            delete(var)
            mod_i = mod_i + 1

        end do  ; loop over models

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
        temp_list_new_string = array_append_record( \
            temp_list_string, ncdf_dir, 0)
        temp_list_new = tochar(temp_list_new_string)
        ; Remove previous file
        system("rm " + temp_dir)
    else
        ncdf_char = tochar(ncdf_dir)
        temp_list_new = new((/1, dimsizes(ncdf_char)/), character)
        temp_list_new(0, :) = ncdf_char
    end if

    ;; Create new file and add list
    temp = addfile(temp_dir, "c")
    temp->temp_list = temp_list_new

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)

end
