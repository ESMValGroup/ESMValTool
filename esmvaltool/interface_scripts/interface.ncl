;;#############################################################################
;; INTERFACE TO HANDLE THE COMMUNICATION BETWEEN THE PYTHON WORKFLOW AND NCL
;;#############################################################################
;; Load the automatically-generated interface files <variable>_info.ncl and
;; settings.ncl for the current diagnostic and variable with the current
;; preprocessor set. These files contains the NCL logical _info varialbe with
;; all the required information to run the diagnostic script
;;
;; No functions/procedures shall be further added to this script.
;;#############################################################################

;; Load config_user_info and diag_script_info from settings.ncl
begin
    print("INFO    Loading settings from " + getenv("settings"))
    loadscript("$settings")
end

;; Load input_file_info, model_info and variable_info from <variable>_info.ncl
begin
    vardeffiles = \
        str_match_ic_regex(diag_script_info@input_files, ".*_info\.ncl")
    if (.not. ismissing(vardeffiles)) then
        vardeffile = diag_script_info@input_files(0)
        print("INFO    Loading input data description from " + vardeffile)
        loadscript(vardeffile)
    end if
end


begin

    ;; Add trailing slash to paths
    config_user_info@plot_dir = config_user_info@plot_dir + "/"
    config_user_info@run_dir = config_user_info@run_dir + "/"
    config_user_info@work_dir = config_user_info@work_dir + "/"

    ;; Copy some info into legacy names
    diag_script = diag_script_info@script
    if isvar("variable_info") then
        variables = (/variable_info@short_name/)
        field_types = (/variable_info@field/)
    end if

end

;; Load other interface scripts
load "./interface_scripts/logging.ncl"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"

