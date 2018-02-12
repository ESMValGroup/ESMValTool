;;#############################################################################
;; INTERFACE TO HANDLE THE COMMUNICATION BETWEEN THE PYTHON WORKFLOW AND NCL
;;#############################################################################
;; Please consider using of extending existing routines before adding new ones.
;; Check the header of each routine for documentation.
;;
;; Contents:
;;    function interface_get_idx_var
;;    function interface_get_var_mip
;;    function interface_get_var_exp
;;    function interface_replace_place_holders
;;    function interface_get_figure_filename
;;    function interface_get_fullpath
;;    function interface_get_infile
;;    function interface_get_inpaths
;;    function interface_get_dictkeys
;;    function interface_get_dictkeys_no_var
;;    function noop
;;    funciton addfile_wrapper
;;
;;#############################################################################

begin
    print("info: Loading settings from " + getenv("settings"))
    loadscript("$settings")
end

; Copy some variables into legacy names
begin
    diag_script = diag_script_info@script
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
  	  	variables = (/variable_info@short_name/)
  	  	derived_var = (/variable_info@short_name/)
  	  	field_types = (/variable_info@field/)
  	  	derived_field_type = (/variable_info@field/)
  	  	var_attr_mip = input_file_info@mip
  	  	var_attr_exp = input_file_info@exp
  	  	models = input_file_info
  	  	models@name = input_file_info@model
  	  	dictkeys = True
  	  	dictkeys@dictkeys = input_file_info@filename
  	  	fullpaths = input_file_info@filename
    end if
end

