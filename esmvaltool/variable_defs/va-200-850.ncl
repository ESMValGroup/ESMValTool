;
; Requires: va:T3M
; Extract wind at pressure levels for Monsoon indicies computation
load "interface_scripts/data_handling.ncl"

variable_info = True
variable_info@derived = True
variable_info@units = "m s-1"

undef("calculate")
function calculate(index [1] : integer, 
                   variable [1] : string, 
                   field_type [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - Index to current infile defined in the 'interface_data/ncl.interface'-file
;;    variable - Current variable as string
;;    field_type  - String with field type classification
;; Return value:
;;    data_new - logical
;; Description:
;;    Computes the Monsoon Hadley Index for the designated models. See
;;    Goswami B.N., Krishnamurthy V., Annamalai H., 1999, A broad-scale circu-
;;    lation index for the interannual variability of the Indian summer monsoon,
;;    Quarterly Journal of the Royal Meteorological Society, vol. 125, issue 554,
;;    pp. 611-633
local result_va, data_va_200, data_va_850
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering va-200-850.ncl", verbosity, 4)
    result_va = read_data(index, "va", "T3M")

    data_va_200 = extract_data(index, result_va, -1, 20000, 20000)
    data_va_850 = extract_data(index, result_va, -1, 85000, 85000)

    dim_va = dimsizes(data_va_200)
    va_200_850 = conform_dims((/dim_va(0), 2, dim_va(2), dim_va(3)/), 0.0, -1)

    va_200_850(:, 0, :, :) = data_va_200(time|:, {plev|20000}, lat|:, lon|:)
    va_200_850(:, 1, :, :) = data_va_200(time|:, {plev|85000}, lat|:, lon|:)

    ;; Rewrite field specification to T2Ms (2D timeseries on the surface)
   ; c=stringtochar(field_type)
   ; c(1)=integertochar(50) ; "2"
   ; c(3)=integertochar(115) ; "s"
   ; field_type=chartostring(c)
   ; delete(c)

    data_new = True
    add_data_var(index, data_new, va_200_850, variable)
    info_output(">>>>>>>> Leaving va-200-850.ncl", verbosity, 4)
    return(data_new)
end


