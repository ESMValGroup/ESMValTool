;;#############################################################################
;; WRAPPER SCRIPT FOR COLLECTING AND PLOTTING PRECALCULATED METRICS
;; Author: Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################
;; Description
;;    Collects metrics previously calculated by perfmetrics_grading.ncl and
;;    passes them to the portrait_plot function
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    label_bounds: min and max of the labelbar
;;    label_scale: bin width of the labelbar
;;    colormap: color table of the labelbar
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    label_lo: add lower triange for values outside range
;;    label_hi: add upper triange for values outside range
;;    cm_interval: min and max color of the color table
;;    cm_reverse: reverse color table
;;    sort: sort models in alphabetic order (excluding MMM)
;;    title: plot title
;;    scale_font: scaling factor applied to the default font size
;;    disp_values: switch on/off the grading values on the plot
;;    disp_rankings: switch on/off the rankings on the plot
;;    rank_order: display rankings in increasing (1) or decreasing (-1) order
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
;;    20150109-A_gott_kl: distinguish RMSD & BIAS in output file naming.
;;    20140620-A_gott_kl: modified output file naming (variable outfile_plot)
;;                     and added option for explicit plot title.
;;    20140204-A_fran_fr: written.
;;
;;#############################################################################


load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/style.ncl"

load "./plot_scripts/ncl/portrait_plot.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 4)
    diag_script_base = basename(diag_script)

    ;; Check file type
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "ps"
    end if

    ;; Write references
    write_references(diag_script, \
                     "A_fran_fr", \
                     (/"A_righ_ma", "A_eyri_ve", "A_gott_kl"/), \
                     (/"D_righi15gmd", "D_gleckler08jgr"/), \
                     (/"E_ncep", "E_erainterim", "E_airs", \
                       "E_ceresebaf", "E_srb"/), \
                     (/"P_embrace", "P_esmval"/))

    ;; Read grading list and add files together
    xml_name = getenv("ESMValTool_xml_name")
    suffix = get_file_suffix(xml_name, 0)
    name = str_sub_str(xml_name, suffix, "")
    temp_dir = get_output_dir(getenv("ESMValTool_wrk_dir"), "temp") + \
        name + "_grading.nc"

    if (isfilepresent_esmval(temp_dir)) then

        ;; Read grading list
        temp_file = addfile(temp_dir, "r")
        temp_list = temp_file->temp_list
        file_list = tostring(temp_list)
        file_list@_FillValue = default_fillvalue(typeof(file_list))
        delete(temp_list)

        ;; Remove duplicates
        do ii = 0, dimsizes(file_list) - 1
            curr = file_list(ii)
            idx = ind(curr.eq.file_list)
            if (dimsizes(idx).gt.1) then
                file_list(idx) = file_list@_FillValue
                file_list(ii) = curr
            end if
            delete(curr)
            delete(idx)
        end do
        tmp = file_list(ind(.not.ismissing(file_list)))
        delete(file_list)
        file_list = tmp
        delete(tmp)

        ;; Loop over files in the list, append data
        do ii = 0,  dimsizes(file_list) - 1
            data_temp = ncdf_read(file_list(ii), "grade")

            ;; Change to 3 dimensional
            if (dimsizes(dimsizes(data_temp)).eq.2) then

                dim_temp = array_append_record(dimsizes(data_temp), 2, 0)
                temp = new(dim_temp, typeof(data_temp))
                temp(:, :, 0) = data_temp
                delete(data_temp)
                data_temp = temp
                delete(temp)

            end if

            if (.not.isdefined("data_all")) then

                dim_temp = dimsizes(data_temp)
                dim_data = (/dimsizes(file_list), dim_temp(1), 2/)
                data_all = new(dim_data, float)
                data_all(0, :, :) = data_temp

                data_all!1 = "models"
                data_all&models = data_temp&models

            else

                ;; If model coordinates are identical
                if (dimsizes(data_temp&models).eq.dimsizes(data_all&models)) \
                    .and. all(data_temp&models.eq.data_all&models) then
                    data_all(ii, :, :) = (/data_temp/)
                    copy_VarCoords(data_temp, data_all(ii:ii, :, :))
                else
                    ;; Loop over models in new data entry
                    do imod_temp = 0, dimsizes(data_temp&models) - 1
                        ;; If current model is not already part
                        ;; of the model coordinate
                        if (.not.any(data_temp&models(imod_temp) .eq. \
                            data_all&models)) then
                            ;; Append record for model(imod)
                            data_new = extend_var_at(data_all, 1, \
                                                     dimsizes(data_all&models))
                            data_new(ii, dimsizes(data_all&models), :) = \
                                (/data_temp(0, imod_temp, :)/)
                            data_new&models(dimsizes(data_all&models)) = \
                                (/data_temp&models(imod_temp)/)
                            delete(data_all)
                            data_all = data_new
                            delete(data_new)
                        else
                            ;; Loop over models of data
                            do imod = 0, dimsizes(data_all&models) - 1
                                ;; if neq data model is similar to current
                                ;; entry, write data entry
                                if (data_all&models(imod).eq. \
                                    data_temp&models(imod_temp)) then
                                    data_all(ii, imod, :) = \
                                        (/data_temp(0, imod_temp, :)/)
                                    copy_VarCoords( \
                                        data_temp(0:0, imod_temp, :), \
                                        data_all(ii:ii, imod, :))
                                end if
                            end do
                        end if
                    end do
                end if
            end if
            delete(data_temp)
        end do

        ;; Reduce size if all entries have only one reference
        if (all(ismissing(data_all(:, :, 1)))) then
            data_new = data_all(:, :, 0)
            delete(data_all)
            data_all = data_new
            delete(data_new)
        end if

    else

        error_msg("f", diag_script, "", "cannot find grading list")

    end if

    delete(data_all@var)
    delete(data_all@title)
    delete(data_all@ncdf_dir)
    if(isatt(data_all, "reference")) then
        delete(data_all@reference)
    end if

    ;; Sort models in alphabetical order (exclude MMM)
    if (isatt(diag_script_info, "sort")) then
        if (diag_script_info@sort) then
            pid = sort_alphabetically(data_all&models, \
                                      ind(data_all&models.eq."Mean-model".or. \
                                          data_all&models.eq."Median-model"), \
                                      "begin")
            if (dimsizes(dimsizes(data_all)).eq.3) then
                tmp = data_all(:, pid, :)
            else
                tmp = data_all(:, pid)
            end if
            tmp&models = data_all&models(pid)
            delete(data_all)
            data_all = tmp
            delete(tmp)
            delete(pid)
        end if
    end if

    ;; Attach plotting options
    copy_VarAtts(diag_script_info, data_all)
    if(isatt(diag_script_info, "title")) then
        data_all@res_tiMainString = diag_script_info@title
    end if

    ;; Plot
    if (getenv("ESMValTool_write_plots").eq.True) then
        xml_name = getenv("ESMValTool_xml_name")
        suffix = get_file_suffix(xml_name, 0)
        name = str_sub_str(xml_name, suffix, "")

        ;; Create outfile directory
        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, "perfmetrics_grading")
        if(dimsizes(data_all&diagnostics).gt.2) then
            diags = data_all&diagnostics(0) + "_to_" + \
                data_all&diagnostics(dimsizes(data_all&diagnostics) - 1)
        else
            diags = str_concat(data_all&diagnostics)
        end if

        ;; Define outfile name
        outfile_plot = output_dir + name + "_" + diags
        if (isatt(data_all, "metric")) then
            outfile_plot = outfile_plot + "_" + data_all@metric
        end if
        outfile_plot = outfile_plot + "_grading"

        ;; Create workspace
        wks = gsn_open_wks(file_type, outfile_plot)

        ;; Call plotting function
        plot = portrait_plot(wks, data_all, "grade")

        info_output(" gv " + outfile_plot + "." + file_type, verbosity, 1)

    end if

    ;; Remove temporary grading list
    system("rm -f " + temp_dir)

    leave_msg(diag_script, "", 4)

end
