;;#############################################################################
;; FUNCTION FOR LOADING AND COMPUTING THE FIRST (AND OFTEN ONLY) VARIABLE
;;#############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "./interface_data/ncl.interface"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"

;;#############################################################################

begin
    ;; Load variable attributes (a single variable)
    loadscript(variable_def_dir + "/" + derived_var + ".ncl")
end

begin

    funcname = "derive_var"
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    enter_msg("interface_scripts/derive_var.ncl", funcname, 2)

    ;; If variable_info@derived = True, a 'calculate' function must be defined
    ;; the variable_defs/var.ncl script
    data = True
    curr_var_idx = ind(derived_var .eq. variables)

    if (variable_info .eq. True) then

        ;; Loop over models
        do model_index = 0, dimsizes(models@name) - 1
            calcData = False
            if (var_attr_exclude(curr_var_idx).eq. \
                model_attr_id(model_index)) then
                continue
            end if
            if (model_attr_skip(model_index) .eq. "True") then
                continue
            end if

            ;; Check if data conversion is required (variable_info@derived)
            if (variable_info@derived) then
                varfile = interface_get_fullpath(derived_var, \
                                                 derived_field_type, \
                                                 model_index)
                ;; Check if file is present
                if (.not. isfilepresent_esmval(varfile) .or. \
                    (force_processing.eq."True")) then
                    calcData = True
                else
                    fv = addfile(varfile, "r")
                    att_units = fv->$derived_var$@units
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
                info_output("Processing data for variable = " + derived_var + \
                            " and field = " + derived_field_type, verbosity, 1)
                data = calculate(model_index, derived_var, derived_field_type)
                delete(data)
            end if
        end do
    end if

    leave_msg("interface_scripts/derive_var.ncl", funcname, 2)

end
