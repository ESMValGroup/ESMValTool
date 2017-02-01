;
; Requires: none
;
info = True
info@derived = True
info@long_name = "HCl"
info@units = "ppbv"

info@fig05_yr_min = 1980
info@fig05_yr_max = 1999
info@fig05_lat_sel = (/80.,0.,-80./)
info@fig05_month_lat_sel = (/"APR","APR","NOV"/)
info@fig05_lev_sel = (/50.,50./)
info@fig05_month_lev_sel = (/"APR","NOV"/)
info@fig05_obs_name = (/"HALOE"/)
;info@fig05_refModel=(/"UMETRAC"/)
info@fig05_obs_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"
info@fig05_obs_stdv_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"

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
;; Description:
;;    Converts hcl units "mol/mol" to "ppbv"
local hcl,hcl_flt
begin
  data_new = read_data(index, variable, field_number)
  hcl = extract_data(index,data_new,-1,0,0)
  hcl = (/hcl * 1.e9/)
  hcl@units = info@units
  if (typeof(hcl).eq."double") then
    hcl_flt = dble2flt(hcl)
  else
    hcl_flt = hcl
  end if
  modify_data(index, data_new, hcl_flt)
  delete(hcl)
  delete(hcl_flt)
  return(data_new)
end
