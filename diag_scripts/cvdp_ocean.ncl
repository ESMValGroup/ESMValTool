;;#############################################################################
;; CONVERT FILE NAMES
;; Author: Dominik Kunert (DLR, Germany)
;; ESMVal project
;;#############################################################################
;; Description
;;    The CVDP expects the endings of the input data file names in 
;;    the format "YYYYMM-YYYYMM.nc". This script converts
;;    the OImon (Ocean Ice monthly) and the Omon (Ocean monthly) components
;;    of the climo files to the required format.
;;
;; Required variable_info attributes (variable specific)
;;    none
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    none
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    none
;;
;; Caveats
;;    none
;; Modification history
;;    20151028-A_laue_ax: added saving list of processed files to log-file
;;    20151027-A_laue_ax: moved call to 'write_references' to the beginning
;;                        of the code
;;    20151014-A-laue_ax: replaced "mv" commands with "ln -s", added call to
;;                        'write_references'
;;    20150924-A_laue_ax: replaced hard wired directories with directories
;;                        specified in namelist (climo_dir)
;;    201504??-A_kune_dm: written.
;;
;;#############################################################################

;; Temporary file updated by Python on the fly
load "./interface_data/ncl.interface"

load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

load "./diag_scripts/lib/ncl/ensemble.ncl"

load "./diag_scripts/lib/ncl/set_operators.ncl"

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
;*********************************************************************

begin
    ;; ============================================
    ;; ||                                        ||
    ;; ||         READING CONFIGURATION          ||
    ;; ||                                        ||
    ;; ============================================

    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    ;; 'diag_script' if fetched from ncl.interface
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("Ocean Section", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

    dim_MOD = dimsizes(models@name)
    dim_VAR = dimsizes(variables)
    field_type0 = field_types(0)
    var = variables

    diag_script_base = basename(diag_script)

    mods = new(dim_MOD, "string")
    paths = new(dim_MOD, "string")

; A-laue_ax+
    ;; References Information
    write_references(diag_script,      \  ; script name
                     "A_kune_dm",      \  ; authors
                     "",               \  ; contributors
                     "D_0000",         \  ; diag_references
                     "",               \  ; obs_references
                     (/""/))              ; proj_references
; A-laue_ax-

;------------------------start_namelist_creator-------------------

; A-laue_ax+
    climo_dir = getenv("ESMValTool_climo_dir") + "/"
; A-laue_ax-

    do idx_mod = 0, dim_MOD - 1
        do idx_var = 0, dim_VAR - 1
            if (models@project(idx_mod) .eq. "CMIP5_ETHZ") then
; A-laue_ax+
                oldfile = climo_dir + models@project(idx_mod) + "/" + \
                    models@project(idx_mod) + "_" + models@mip(idx_mod) + "_" +\
                    models@experiment(idx_mod) + "_" + models@name(idx_mod) + \
                    "_" + models@ensemble(idx_mod) + "_" + field_type0 + "_" + \
                    variables(idx_var) + "_" + models@start_year(idx_mod) + \
                    "-" + models@end_year(idx_mod) + ".nc"
                newfile = climo_dir + \
                    models@project(idx_mod) + "/" + models@project(idx_mod) + \
                    "_" + models@experiment(idx_mod) + "_" + \
                    models@name(idx_mod) + "_" + variables(idx_var) + "_" \
                    + models@start_year(idx_mod) + "01-" + \
                    models@end_year(idx_mod) + "12.nc"

                if (isfilepresent_esmval(newfile)) then
                    system("rm -f " + newfile)
                end if
; A-laue_ax-
                system("ln -s " + oldfile + " " + newfile)
; A-laue_ax+
            else
                error_msg("f", diag_script, "", "project type CMIP5_ETHZ " \
                    + "supported only. Abort.")
; A-laue_ax-
            end if
; A-laue_ax+
            ; write file info to log-file
            inquire_and_save_fileinfo(oldfile, models@name(idx_mod), variables(idx_var))
; A-laue_ax-
        end do
    end do

;--------------------------end_namelist_creator----------------

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end
