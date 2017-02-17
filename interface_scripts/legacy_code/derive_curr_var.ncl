load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"

; Loading temporary interface file. This file is updated on the fly by Python
; and used for Python->NCL communication
load "./interface_data/ncl.interface"

load "./interface_scripts/extract_data.ncl"
load "./interface_scripts/get_varName_from_file.ncl"
load "./interface_scripts/info_output.ncl"
load "./interface_scripts/modify_data.ncl"
load "./interface_scripts/read_data.ncl"
load "./interface_scripts/write_references.ncl"

; Load environment variables
verbosity        = stringtointeger(getenv("ESMValTool_verbosity"))

begin
    ;; Load variable attributes (a single variable)
    ;; 'curr_var' is fetched from ncl.interface
    loadscript(variable_def_dir + "/" + curr_var + ".ncl")
end

;#######################################################################
; Load and compute the first (and often only) variable
;#######################################################################
begin
    ;; If variable_info@derived = True, a 'calculate' function must be defined
    ;; the variable_defs/var.ncl script
    info_output("<<<<<<<< Entering derived_variables.ncl", verbosity, 2)
    data = True
    if (variable_info .eq. True) then

        ;; Loop over models
        do i = 0, dimsizes(models@name) - 1
            calcData = False
            ;; Check if data conversion is required (variable_info@derived)
            ;; 'curr_field_type' is fetched from ncl.interface
            if (variable_info@derived) then
                varfile = interface_get_fullpath(curr_var, curr_field_type, i)
                ;; Check if file is present, 
                ;; 'force_processing' is fetched from ncl.interface
                if (.not. isfilepresent_esmval(varfile) .or. \
                    (force_processing.eq."True")) then
                    calcData = True
                else
                    fv = addfile(varfile, "r")
                    att_units = fv->$curr_var$@units
                    ;; Check if units are present and equal to expected units
                    if (.not.ismissing(att_units)) then
                        if (.not. (att_units .eq. variable_info@units)) then
                            calcData = True
                        end if
                    end if
                    delete(fv)
                    delete(att_units)
                end if
            end if

            ;; Rewrite data if necessary
            if (calcData) then
                info_output("Processing data for variable = " + curr_var +\
                            " and field = " + curr_field_type, verbosity, 1)
                data = calculate(i, curr_var, curr_field_type)
                delete(data)
            end if
        end do
    end if
end

info_output(">>>>>>>> Leaving derived_variables.ncl", verbosity, 2)
