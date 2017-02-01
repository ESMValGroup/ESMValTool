;
; Requires: none
;
info = True
info@derived = True
info@long_name = "Cly"
info@units = "ppbv"

info@fig12_xmin_year =1980
info@fig12_xmax_year =2010
info@fig12_lat_sel = (/-80./)
info@fig12_month_lat_sel = (/"NOV"/)
info@fig12_lev_sel = (/50./)
info@fig12_month_lev_sel = (/"OCT"/)
info@fig12_ymax_lev = 50
info@fig12_ymin_lev = 10
info@fig12_obs_name = (/"OBS"/)
info@fig12_obs_file = "./plot_type/input_data/OBS/CCMVal1_1992-2005_OBS_1_T2Mz_Cly.nc"
info@fig12_obs_stdv_file = "./plot_type/input_data/OBS/CCMVal1_1992-2005_OBS_1_T2Mz_Cly_std.nc"


info@fig05_yr_min = 1980
info@fig05_yr_max = 1999
info@fig05_lat_sel = (/-80./)
info@fig05_month_lat_sel = (/"NOV"/)
info@fig05_lev_sel = (/50.,50./)
info@fig05_month_lev_sel = (/"APR","NOV"/)
info@fig05_altitude = True

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
;;    Converts cly units "mol/mol" to "ppbv"
local cly,cly_flt
begin
  data_new = read_data(index, variable, field_number)
  cly = extract_data(index,data_new,-1,0,0)
  cly = (/cly * 1.e9/)
  cly@units = info@units
  if (typeof(cly).eq."double") then
    cly_flt = dble2flt(cly)
  else
    cly_flt = cly
  end if
  modify_data(index, data_new, cly_flt)
  delete(cly)
  delete(cly_flt)
  return(data_new)
end

