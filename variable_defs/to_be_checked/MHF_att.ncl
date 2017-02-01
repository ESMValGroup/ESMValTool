;
; REQUIRES: ua:*3*,ta:*3*
;
info          = True
info@derived  = True
info@long_name= "meridional heat flux"
info@units    = "??"

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
; calculate meridional heat flux and return 
begin
 ua=read_data(index,"ua",field_number)
 ta=read_data(index,"ta",field_number)
; calc zonal mean for ua and ta
; result=ua'*ta
  print("Calculating MHF")
  result=1
  return(result)
end
