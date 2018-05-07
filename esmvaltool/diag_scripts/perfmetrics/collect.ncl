;;#############################################################################
;; WRAPPER SCRIPT FOR COLLECTING AND PLOTTING PRECALCULATED METRICS
;; Author: Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################
;; Description
;;    Collects metrics previously calculated by main.ncl and passes them to the
;;    portrait_plot or to the taylor function
;;
;; Required diag_script_info attributes
;;    metric: selected metric (RMSD, BIAS or taylor)
;;    label_bounds (for RMSD and BIAS): min and max of the labelbar
;;    label_scale (for RMSD and BIAS): bin width of the labelbar
;;    colormap (for RMSD and BIAS): color table of the labelbar
;;
;; Optional diag_script_info attributes
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
;; Caveats
;;
;; Modification history
;;    20180508-A_righ_ma: renamed and generalized
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

load "interface_scripts/interface.ncl"

load "./diag_scripts/shared/plot/style.ncl"
load "./diag_scripts/shared/plot/portrait_plot.ncl"
load "./diag_scripts/shared/plot/taylor_plot.ncl"

begin

    enter_msg(diag_script, "")
    diag_script_base = basename(diag_script)

    ;; Check file type
    file_type = config_user_info@output_file_type
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

    ;; List of grading files
    input_files = diag_script_info@input_files + "/" + \
        diag_script_info@metric +  ".nc"
    input_files := tostring(input_files)

    ;; Filter out non-existing files
    input_files := input_files(ind(isfilepresent(input_files)))

    ;; Loop over files in the list, read and append data
    do ii = 0, dimsizes(input_files) - 1

        f = addfile(input_files(ii), "r")
        curr_file = tostring(f->temp_list)
        data_temp = ncdf_read(curr_file, "grade")

        ;; Change to 3 dimensional
        if (dimsizes(dimsizes(data_temp)).eq.2) then

            dim_temp = array_append_record(dimsizes(data_temp), 2, 0)
            temp = new(dim_temp, typeof(data_temp))
            temp(:, :, 0) = data_temp
            data_temp := temp
            delete(temp)

        end if

        ;; Create array for collecting all models
        if (.not.isdefined("data_all")) then
            dim_temp = dimsizes(data_temp)
            dim_data = (/dimsizes(input_files), dim_temp(1), 2/)
            data_all = new(dim_data, float)
            data_all(0, :, :) = data_temp
            data_all!1 = "models"
            data_all&models = data_temp&models
        end if

        ;; Make sure model coordinate is consistent
        consistent = False
        if (dimsizes(data_temp&models).eq.dimsizes(data_all&models)) then
            if (all(data_temp&models.eq.data_all&models)) then
                consistent = True
            end if
        end if

        ;; Append data
        if (consistent) then
           data_all(ii, :, :) = (/data_temp/)
           copy_VarCoords(data_temp, data_all(ii:ii, :, :))
        else
            do imod2 = 0, dimsizes(data_temp&models) - 1
                if (.not.any(data_temp&models(imod2).eq.data_all&models)) then
                   ;; Append record for model(imod)
                    data_new = extend_var_at(data_all, 1, \
                                             dimsizes(data_all&models))
                    data_new(ii, dimsizes(data_all&models), :) = \
                        (/data_temp(0, imod2, :)/)
                    data_new&models(dimsizes(data_all&models)) = \
                        (/data_temp&models(imod2)/)
                    data_all := data_new
                    delete(data_new)
                else
                    ;; Loop over models of data
                    do imod = 0, dimsizes(data_all&models) - 1
                        ;; if neq data model is similar to current
                        ;; entry, write data entry
                        if (data_all&models(imod).eq. \
                            data_temp&models(imod2)) then
                            data_all(ii, imod, :) = \
                                (/data_temp(0, imod2, :)/)
                            copy_VarCoords(data_temp(0:0, imod2, :), \
                                           data_all(ii:ii, imod, :))
                        end if
                    end do
                end if
            end do
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
    delete(data_all@var)
    delete(data_all@title)
    delete(data_all@ncdf_dir)
    if (isatt(data_all, "reference")) then
        delete(data_all@reference)
    end if

    ;; Sort models in alphabetical order, excluding multi-model 
    ;; mean/median which are placed at the beginning
    if (isatt(diag_script_info, "sort")) then
        if (diag_script_info@sort) then
            idx = ind(data_all&models.eq."MultiModelMean" .or. \
                      data_all&models.eq."MultiModelMedian")
            pid = sort_alphabetically(data_all&models, idx, "begin")
            if (dimsizes(dimsizes(data_all)).eq.3) then
                data_all := data_all(:, pid, :)
            else
                data_all := data_all(:, pid)
            end if
            delete(pid)
        end if
    end if

    ;; Attach plotting options
    copy_VarAtts(diag_script_info, data_all)
    if (isatt(diag_script_info, "title")) then
        data_all@res_tiMainString = diag_script_info@title
    end if
    if (diag_script_info@metric.eq."taylor") then
        data_all@res_varLabels = data_all&models
        data_all@res_caseLabels = data_all&diagnostics
    end if

    ;; Create outfile directory
    system("mkdir -p " + config_user_info@plot_dir)
    if (dimsizes(data_all&diagnostics).gt.1) then
        diags = data_all&diagnostics(0) + "_to_" + \
            data_all&diagnostics(dimsizes(data_all&diagnostics) - 1)
    else
        diags = str_concat(data_all&diagnostics)
    end if

    ;; Define outfile name
    outfile_plot = config_user_info@plot_dir + diags
    if (isatt(data_all, "metric")) then
        outfile_plot = outfile_plot + "_" + data_all@metric
    end if

    ;; Create workspace
    wks = gsn_open_wks(file_type, outfile_plot)
    if (diag_script_info@metric.eq."taylor") then
        wks@legendfile = outfile_plot + "_legend"
    end if

    ;; Call plotting function
    if (diag_script_info@metric.eq."taylor") then
        plot  = taylor_plot(wks, data_all, "grade")
    else
        plot = portrait_plot(wks, data_all, "grade")
    end if

    log_info(" gv " + outfile_plot + "." + file_type)

    leave_msg(diag_script, "")

end
