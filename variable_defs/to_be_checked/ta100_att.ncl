;
;  Requires: ta:*3*
;
info = True
info@derived = True
info@long_name="air_temperature"
info@units="K"

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
local tmp, dum, i
begin
  data_new = True
  tmp = read_data(index,"ta","*3*")
    dum = extract_data(index,tmp,1,100.,100.)   ;extract_data(models,tmp,i,"plev",100.,100.)

    dum@long_name=info@long_name
    add_data(index,data_new,dum)

    delete(dum)
  delete(tmp)

  return(data_new)
end
