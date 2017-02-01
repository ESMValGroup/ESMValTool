;
;  Requires: vorpot840:*2*s
;
info = True
info@derived = True
info@long_name="Equivalent_latitude_of_vorpot840"
info@units="degrees"
;load "./variable_defs/variable_defs_aux/function_eq_lat.ncl"
;load "./plot_type/uniform_coordinates.ncl"
load "./variable_defs/variable_defs_aux/CalcEqLat.ncl"

undef("calculate")
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
local tmp, dum, data, i
begin
  data_new = True

  tmp = read_data(index,"vorpot840","*2*s")
  dum = extract_data(index,tmp,-1,0,0)
  data=CalcEqLat(dum,"T2Ds")
  data@long_name=info@long_name
  add_data(index,data_new,data)
  delete(dum)
  delete(data)
  delete(tmp)

  return(data_new)
end


