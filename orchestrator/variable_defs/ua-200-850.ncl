;
; Requires: ua:T3M
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
local result_ua, data_ua_200, data_ua_850, verbosity
begin
    verbosity = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering ua-200-850.ncl", verbosity, 4)
;    result_ua = read_data(index, "ua", "T3M")
    result_ua = read_data(index, "ua", "*3*")

    data_ua_200 = extract_data(index, result_ua, -1, 20000, 20000)
    data_ua_850 = extract_data(index, result_ua, -1, 85000, 85000)

    dim_ua = dimsizes(data_ua_200)
    ua_200_850 = conform_dims((/dim_ua(0), 2, dim_ua(2), dim_ua(3)/), 0.0, -1)

    ua_200_850(:, 0, :, :) = data_ua_200(time|:, {plev|20000}, lat|:, lon|:)
    ua_200_850(:, 1, :, :) = data_ua_200(time|:, {plev|85000}, lat|:, lon|:)

    ;; Rewrite field specification to T2Ms (2D timeseries on the surface)
   ; c=stringtochar(field_type)
   ; c(1)=integertochar(50) ; "2"
   ; c(3)=integertochar(115) ; "s"
   ; field_type=chartostring(c)
   ; delete(c)

    data_new = True
    add_data_var(index, data_new, ua_200_850, variable)
    info_output(">>>>>>>> Leaving ua-200-850.ncl", verbosity, 4)
    return(data_new)
end


