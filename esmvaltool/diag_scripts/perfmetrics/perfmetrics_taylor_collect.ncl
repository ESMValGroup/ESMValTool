;;#############################################################################
;; WRAPPER SCRIPT FOR COLLECTING AND PLOTTING TAYLOR METRICS
;; Author: Franziska Frank (DLR, Germany)
;; ESMVal project
;;#############################################################################
;; Description
;;    Collects metrics previously calculated by perfmetrics_taylor.ncl and
;;    passes them to the taylor_plot function
;;
;; Required diag_script_info attributes (diagnostics specific)
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
;;    20140620-A_gott_kl: modified output file naming.
;;    20140204-A_fran_fr: written.
;;
;;#############################################################################

load "interface_scripts/interface.ncl"

load "./plot_scripts/ncl/taylor_plot.ncl"

begin

    enter_msg(diag_script, "")
    diag_script_base = basename(diag_script)

    ;; Check file type
    file_type = config_user_info@output_file_type
    if (ismissing(file_type)) then
        file_type = "PS"
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
    input_files = diag_script_info@input_files + "/taylor.nc"
    input_files := tostring(input_files)

    ;; Loop over files in the list, read and append data
    do ii = 0, dimsizes(input_files) - 1

        f = addfile(input_files(ii), "r")
        curr_file = tostring(f->temp_list)
        if (.not.fileexists(curr_file)) then
            error_msg("f", diag_script, "", \
                      "cannot find grading list element:" + temp_dir)
        end if
        data_temp = ncdf_read(curr_file, "grade")

        ;; Collect all models data
        if (.not.isdefined("data_all")) then
            data_all = data_temp
        else
            tmp = array_append_record(data_all, data_temp, 0)
            delete(data_all)
            data_all = tmp
            delete(tmp)
        end if
        delete(data_temp)

    end do

    ;; Attach resources
    data_all@res_varLabels    = data_all&models
    data_all@res_caseLabels   = data_all&diagnostics
    copy_VarAtts(diag_script_info, data_all)

    ;; Plot
    if (config_user_info@write_plots .eq. True) then
        ;; Create outfile directory
        plot_dir = config_user_info@plot_dir
        system("mkdir -p " + plot_dir)
        if(dimsizes(data_all&diagnostics).gt.2) then
            diags = data_all&diagnostics(0) + "_to_" + \
                data_all&diagnostics(dimsizes(data_all&diagnostics) - 1)
        else
            diags = str_concat(data_all&diagnostics)
        end if
        outfile_plot = plot_dir + diags + "_taylor_diagram"

        ;; Create workspace
        wks = gsn_open_wks(file_type, outfile_plot)
        wks@legendfile = outfile_plot + "_legend"

        ;; Plot
        plot  = taylor_plot(wks, data_all, "grade")

        log_info(" gv " + outfile_plot + "." + file_type)
    end if

    leave_msg(diag_script, "")

end
