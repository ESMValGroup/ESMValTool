;
;  Requires: ita:TO2Ms
;
variable_info = True
variable_info@derived = False
variable_info@long_name="depth weighted temperate to 730m"
variable_info@units="K"
load "interface_scripts/data_handling.ncl"

undef("calculate")
;function calculate(models,variable,field_number)
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_number [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'temp_ncl_code/filesinfo'-file
;;    variable - logical with releveant variable as string attribute
;;    field_number  - string with field number classification
;; Return value:
;;    data_new - logical
local tmp, dum, dimension 
begin
  data_new = True
  tmp = read_data(index,"ita","TO2Ms")
   dimension = 0  ;;time
   dum = extract_data(index,tmp,dimension,1985.,2012.0)
;  dum = tmp
  dum@long_name=variable_info@long_name
  add_data_var(index, data_new, dum, variable)
  return(data_new)
end
