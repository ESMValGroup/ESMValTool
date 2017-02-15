;;#############################################################################
;; AUXILIARY ROUTINES FOR PERFMETRICS DIAGNOSTICS
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function clim_diag_file_path
;;
;;#############################################################################

load "./interface_data/ncl.interface"
load "./interface_scripts/messaging.ncl"

;;#############################################################################
undef("clim_diag_file_path")
function clim_diag_file_path(var[1]: string,
                             field_type[1]: string,
                             model_ind[1]: integer,
                             diag_name[1]: string,
                             opt[1]: logical)
;;
;; Arguments
;;    var: variable which is considered in the diagnostic.
;;    field_typ: field type of current variable.
;;    model_ind: index of the current model.
;;    opt: options, needs the attribute @time_avg
;;
;; Return value
;;    A string with full file path.
;;
;; Description
;;    Provides a file name for the intermediate files of the perfmetrics
;;    calculations.
;;
;; Caveats
;;
;; References
;;
;; Modification history:
;;    * 20140312-A_righ_ma: changed subdir structure.
;;    * 20140212-A_fran_fr: written.
;;
local funcname, scriptname, verbosity, workpath, fullpath, outname
begin

    funcname = "clim_diag_file_path"
    scriptname = "diag_scripts/aux/perfmetrics/aux_perfmetrics.ncl"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg(scriptname, funcname, 6)

    ;; Workdir
    workpath = getenv("ESMValTool_wrk_dir")

    ;; Path
    fullpath = workpath + "/" + basename(diag_name) + "/"
    if (.not.isfilepresent_esmval(fullpath)) then
        system("mkdir -p " + fullpath)
    end if

    ;; File name
    outname = basename(systemfunc("basename " + \
                       interface_get_fullpath(var, field_type, model_ind)))
    outname = outname + "_" + basename(diag_name)

    ;; Add extra info to the file
    if (isatt(opt, "plot_type")) then
        outname = outname + "_" + opt@plot_type
    end if
    if (isatt(opt, "time_avg")) then
        outname = outname + "_" + opt@time_avg
    end if
    if (isatt(opt, "level")) then
        outname = outname + "_" + opt@level
    end if
    if (isatt(opt, "region")) then
        outname = outname + "_" + str_sub_str(opt@region, " ", "")
    end if
    outname = outname + ".nc"

    out = fullpath + outname

    leave_msg(scriptname, funcname, 6)
    return(out)

end
