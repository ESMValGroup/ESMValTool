; TODO: add a description of the functions in this file
;------------------------------------------------------
; A temporary file written by the invoking Python script
; Passes on a number of variables from Python to NCL
begin
    print("info: Loading settings from " + getenv("settings"))
    loadscript("$settings")
end

; Copy some variables into legacy names
begin
    diag_script = diag_script_info@script
    if (str_lower(config_user_info@log_level) .eq. "debug") then
        config_user_info@verbosity = 100
    else
        config_user_info@verbosity = 1
    end if
end

; Load variable_info
begin
    vardeffiles = str_match_ic_regex(diag_script_info@input_files, ".*_info\.ncl")
end

begin
	  if (.not. ismissing(vardeffiles)) then
	  	  vardeffile = diag_script_info@input_files(0)
	  	  print("info: Loading input data description from " + vardeffile)
	  	  loadscript(vardeffile)
    end if
end

; Copy some variable_info into legacy names
begin
    if isvar("variable_info") then
  	  	variables = (/variable_info@short_name(0)/)
  	  	derived_var = (/variable_info@short_name(0)/)
  	  	field_types = (/variable_info@field(0)/)
  	  	derived_field_type = (/variable_info@field(0)/)
  	  	var_attr_mip = variable_info@mip
  	  	var_attr_exp = variable_info@exp
  	  	models = variable_info
  	  	models@name = variable_info@model
  	  	dictkeys = True
  	  	dictkeys@dictkeys = variable_info@filename
  	  	fullpaths = variable_info@filename
    end if
end

load "./interface_scripts/messaging.ncl"

undef("interface_get_idx_var")
function interface_get_idx_var(variable)
local idx_var
begin
    if (derived_var .eq. "Undefined") then
        idx_var = ind(variable .eq. variables)
    else
        idx_var = ind(derived_var .eq. variables)
    end if
    return idx_var
end

undef("interface_get_var_mip")
function interface_get_var_mip(idx [1] : integer,
                               variable [1] : string,
                               array [*] : string)
; Description:
;    Fetches the current MIP table, either from the models@-data
;    structure or, if specified, from the specific variable attribute
local array_local, idx_var
begin
    array_local = array
    ;; Default case (no variable attribute is specified)
    if (variable .eq. "default") then
        if (isatt(models, "mip")) then
            array_local = str_sub_str(array_local, "${MIP}", models@mip(idx))
        end if
    else  ; Variable attribute specified and should replace the models@-entry
        idx_var = interface_get_idx_var(variable)

        if (var_attr_mip(idx_var) .ne. "None") then
            array_local = str_sub_str(array_local, "${MIP}", var_attr_mip(idx_var))
        else
            if (isatt(models, "mip")) then
                array_local = str_sub_str(array_local, "${MIP}", models@mip(idx))
            end if
        end if
    end if
    return array_local
end

undef("interface_get_var_exp")
function interface_get_var_exp(idx [1] : integer,
                               variable [1] : string,
                               array [*] : string)
; Description:
;    Fetches the current experiment type, either from the models@-data
;    structure or, if specified, from the specific variable attribute
local array_local, idx_var
begin
    array_local = array
    ;; Default case (no variable attribute is specified)
    if (variable .eq. "default") then
        if (isatt(models, "experiment")) then
            array_local = str_sub_str(array_local, "${EXP}", models@experiment(idx))
        end if
    else  ; Variable attribute specified and should replace the models@-entry
        idx_var = interface_get_idx_var(variable)

        if (var_attr_exp(idx_var) .ne. "None") then
            array_local = str_sub_str(array_local, "${EXP}", var_attr_exp(idx_var))
        else
            if (isatt(models, "experiment")) then
                array_local = str_sub_str(array_local, "${EXP}", models@experiment(idx))
            end if
        end if
    end if
    return array_local
end

undef("interface_replace_place_holders")
function interface_replace_place_holders(variable [1] : string,
                                         field [1] : string,
                                           idx [1] : integer,
                      array_with_place_holders [*] : string)
;                                   return val [1] : string
; Arguments:
;         @brief Reconstructs the current (idx) input path + filename
;         @param variable  -  Current variable
;         @param idx  -  The index to the current model
;         @param field  -  Current field type
;         @param array_with_place_holders  - array with place holder strings
local verbosity, idx_var, array_local_copy
begin
    verbosity  = config_user_info@verbosity
    info_output("<<<<<<<< Entering interface_replace_place_holders", verbosity, 6)

    array_local_copy = str_sub_str(array_with_place_holders(idx), "${VARIABLE}", variable)
    array_local_copy = str_sub_str(array_local_copy, "${FIELD}", field)

    array_local_copy = interface_get_var_mip(idx, variable, array_local_copy)
    array_local_copy = interface_get_var_exp(idx, variable, array_local_copy)

    info_output(">>>>>>>> Leaving interface_replace_place_holders", verbosity, 6)
    return array_local_copy
end

undef("interface_get_figure_filename")
function interface_get_figure_filename(diag_script_base [1] : string,
                                           variable [1] : string,
                                         field_type [1] : string,
                                           aux_info [1] : string,
                                            idx_mod [1] : integer)
;                                        return val [1] : string
;   Arguments:
;         @brief Construct a figure output file name
;         @param diags_script_base  -  The current running diag script without
;                                      it's suffix
;         @param variable  -  Current variable
;         @param field_type  -  Current field type
;         @param aux_info  -  User supplied info to put in figure filename
;         @param idx_mod  -  Current model number, set to -1 if not applicable
local verbosity, aux_sep, figure_name, fig_file_local_copy, use_this_for_var
begin
    verbosity  = config_user_info@verbosity
    info_output("<<<<<<<< Entering interface_get_figure_filename", verbosity, 4)
    sep = "_"  ; Default separator
    aux_sep = "_"  ; Auxiliary info separator
    if (aux_info .eq. "") then
        aux_sep = ""  ; Auxiliary info separator if no aux-info
    else
        aux_info = str_sub_str(aux_info, "_", "-")
    end if

    if (idx_mod .eq. -1) then
        figure_name = diag_script_base + sep + variable \
                                       + sep + field_type \
                                       + aux_sep + aux_info
    else
        if (any(ismissing(ind(variable .eq. variables)))) then
            use_this_for_var = "default"
        else
            use_this_for_var = variable
        end if

        fig_file_local_copy = interface_replace_place_holders(use_this_for_var,\
                                                              field_type,\
                                                              idx_mod,\
                                                              figfiles_suffix)
        figure_name = diag_script_base + sep + variable \
                                       + sep + field_type \
                                       + sep + aux_info + aux_sep \
                                       + fig_file_local_copy
    end if
    info_output(">>>>>>>> Leaving interface_get_figure_filename", verbosity, 4)
    return figure_name
end

undef("interface_get_fullpath")
function interface_get_fullpath(variable [1] : string,
                                   field [1] : string,
                                     idx [1] : integer)
;                             return val [1] : string
; Arguments:
;         @brief Reconstructs the current (idx) input path + filename
;         @param variable  -  Current variable
;         @param idx  -  The index to the current model
;         @param field  -  Current field type
local verbosity, idx_var, fullpaths_local_copy
begin
    verbosity  = config_user_info@verbosity
    info_output("<<<<<<<< Entering interface_get_fullpath", verbosity, 4)
    fullpaths_local_copy = interface_replace_place_holders(variable, field, idx, fullpaths)
    info_output(">>>>>>>> Leaving interface_get_fullpath", verbosity, 4)
    return fullpaths_local_copy
end

undef("interface_get_infile")
function interface_get_infile(variable [1] : string,
                                 field [1] : string,
                                   idx [1] : integer)
;                           return val [1] : string
; Arguments:
;         @brief Reconstructs the current (idx) input filename
;         @param variable  -  Current variable
;         @param idx  -  The index to the current model
;         @param field  -  Current field type
local verbosity, infiles_local_copy
begin
    verbosity  = config_user_info@verbosity
    info_output("<<<<<<<< Entering interface_get_infile", verbosity, 4)
    infiles_local_copy = interface_replace_place_holders(variable, field, idx, infiles)
    info_output(">>>>>>>> Leaving interface_get_infile", verbosity, 4)
    return infiles_local_copy
end

undef("interface_get_inpaths")
function interface_get_inpaths(idx [1] : integer)
;                       return val [1] : string
; Arguments:
;         @brief Returns the current (idx) path to the input filename
;         @param idx  -  The index to the current model
local verbosity
begin
    verbosity  = config_user_info@verbosity
    info_output("<<<<<<<< Entering interface_get_inpaths", verbosity, 4)
    info_output(">>>>>>>> Leaving interface_get_inpaths", verbosity, 4)
    return infile_paths(idx)
end

undef("interface_get_dictkeys")
function interface_get_dictkeys(variable [1] : string,
                                     idx [1] : integer)
;                             return val [1] : string
; Arguments:
;         @brief Reconstructs the current (idx) dictionary keys
;         @param variable  -  Current variable
;         @param idx  -  The index to the current model
;         @param field  -  Current field type
local verbosity, idx_var, dictkeys_local_copy
begin
    verbosity  = config_user_info@verbosity
    info_output("<<<<<<<< Entering interface_get_dictkeys", verbosity, 4)
    dictkeys_local_copy = interface_replace_place_holders(variable, "NO_FIELD", idx, dictkeys@dictkeys)
    info_output(">>>>>>>> Leaving interface_get_dictkeys", verbosity, 4)
    return dictkeys_local_copy
end

undef("interface_get_dictkeys_no_var")
function interface_get_dictkeys_no_var(idx [1] : integer)
;                             return val [1] : string
; Arguments:
;         @brief Reconstructs the current (idx) dictionary keys
;         @param idx  -  The index to the current model
;         @param field  -  Current field type
local verbosity, idx_var, dictkeys_local_copy
begin
    verbosity  = config_user_info@verbosity
    info_output("<<<<<<<< Entering interface_get_dictkeys", verbosity, 4)
    dictkeys_local_copy = interface_replace_place_holders("default", "NO_FIELD", idx, dictkeys@dictkeys)
    info_output(">>>>>>>> Leaving interface_get_dictkeys", verbosity, 4)
    return dictkeys_local_copy
end

undef("noop")
procedure noop()
; no-operation
begin
end

; == Wrapper definitions used extend/redefine existing NCL routines ==

undef("addfile_wrapper")
function addfile_wrapper(filename [1]: string, operation [1]:string)
;                       return [file]: file_handle
; Description:
;         Logs any file access using the addfile function to a log file.
;         To use the wrapper, run the following line to replace all
;          addfile-occurances
;
; find . -type f -name '*ncl' ! -iregex '.*to_be_checked.*' -exec sed -i 's/addfile/addfile_wrapper/g' {} \;
;
local fhandle, logfile, datestamp, use_addfile_wrapper
begin
    verbosity  = config_user_info@verbosity
    info_output("<<<<<<<< Entering addfile_wrapper", verbosity, 4)

    logfile = "addfile_access.log"
    datestamp = systemfunc("date -u")
    operation = str_lower(operation)
    if (operation .eq. "r") then
        logfile_prefix = "'READ:   '"
    else if (operation .eq. "w") then
        logfile_prefix = "'WRITE:  '"
    else if (operation .eq. "c") then
        logfile_prefix = "'CREATE: '"
    end if
    end if
    end if

    system("echo " + datestamp + "  --  " + logfile_prefix + filename + " >> " + logfile)
    fhandle = addfile(filename, operation)

    info_output(">>>>>>>> Leaving addfile_wrapper", verbosity, 4)
    return fhandle
end
