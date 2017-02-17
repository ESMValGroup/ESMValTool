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



load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./plot_scripts/ncl/taylor_plot.ncl"

begin

    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(diag_script, "", 4)
    diag_script_base = basename(diag_script)

    ;; Check file type
    file_type = getenv("ESMValTool_output_file_type")
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

    ;; Read grading list and add files together
    xml_name = getenv("ESMValTool_xml_name")
    suffix = get_file_suffix(xml_name, 0)
    name = str_sub_str(xml_name, suffix, "")
    temp_dir = get_output_dir(getenv("ESMValTool_wrk_dir"), "temp") + \
        name + "_taylor.nc"

    if (isfilepresent_esmval(temp_dir)) then

        ;; Read grading list
        temp_file = addfile(temp_dir, "r")
        temp_list = temp_file->temp_list
        file_list = tostring(temp_list)
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
        do ii = 0, dimsizes(file_list) - 1

            data_temp = ncdf_read(file_list(ii), "var")
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

    else

        error_msg("f", diag_script, "", "cannot find list")

    end if

    ;; Attach resources
    data_all@res_varLabels    = data_all&models
    data_all@res_caseLabels   = data_all&diagnostics
    copy_VarAtts(diag_script_info, data_all)

    ;; Plot
    if (getenv("ESMValTool_write_plots").eq.True) then
        xml_name = getenv("ESMValTool_xml_name")
        suffix = get_file_suffix(xml_name, 0)
        name = str_sub_str(xml_name, suffix, "")

        ;; Create outfile directory
        plot_dir = getenv("ESMValTool_plot_dir")
        output_dir = get_output_dir(plot_dir, "perfmetrics_taylor")
        if(dimsizes(data_all&diagnostics).gt.2) then
            diags = data_all&diagnostics(0) + "_to_" + \
                data_all&diagnostics(dimsizes(data_all&diagnostics) - 1)
        else
            diags = str_concat(data_all&diagnostics)
        end if
        outfile_plot = output_dir + name + "_" + \
            diags + "_taylor_diagram"

        ;; Create workspace
        wks = gsn_open_wks(file_type, outfile_plot)
        wks@legendfile = outfile_plot + "_legend"

        ;; Plot
        plot  = taylor_plot(wks, data_all, "var")

        info_output(" gv " + outfile_plot + "." + file_type, verbosity, 1)
    end if

    ;; Remove temporary Taylor list
    system("rm -f " + temp_dir)

    leave_msg(diag_script, "", 4)

end
