;
;  This function takes a filename and return a variable name.
;
;  2010 - 12-07 HS
;
undef("get_varName_from_file")
function get_varName_from_file(filename[1]: string)
local variable, delimiter, nfields, shortfile, fil, varNames, i
begin
    variable = new(1, string)
    delimiter = "/"
    nfields = str_fields_count(filename, delimiter)
    shortfile = str_get_field(filename, nfields, delimiter)
    fil = addfile(filename, "r")
    varNames = getfilevarnames(fil)
    delete(fil)
    do i = 0, dimsizes(varNames) - 1
        if isStrSubset(shortfile, varNames(i)) then
            variable = varNames(i)
        end if
    end do
    delete(varNames)

    if ismissing(variable) then
        print("[get_varName_from_file.ncl] Unstandard filename: " \
              + "does not contains varName & check the diagnostic code... " \
              + filename)
        exit
    end if
    return(variable)
end
