;
; Requires: none
;
info = True
info@derived = True
info@long_name = "CH4"
info@units = "ppmv"
info@vcp_cnlvl=(/0,1,2,3,4,5,6,7,8,9,10,11,12,15,18/)/10. ; CH4
info@vcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

;info@fig05_yr_min = 1980
;info@fig05_yr_max = 1999
info@fig05_xmax_lat=3.
info@fig05_xmin_lat=0.
info@fig05_ymax_lev=3.
info@fig05_ymin_lev=0.
info@fig05_lat_sel = (/80.,0.,-80./)
info@fig05_month_lat_sel = (/"MAR","MAR","OCT"/)
info@fig05_lev_sel = (/50.,50./)
info@fig05_month_lev_sel = (/"MAR","OCT"/)
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
local ch4
; convert mol/mol to ppmv
begin
  result = read_data(index,variable,field_number)
  ch4 = extract_data(index,result,-1,0,0)
  ch4 = (/ch4*1e6/)       ; convert to ppm
  ch4@units = info@units
  if (typeof(ch4).eq."double") then
    ch4_flt = dble2flt(ch4)
  else
    ch4_flt = ch4
  end if
  modify_data(index,result,ch4_flt)
  delete(ch4)
  delete(ch4_flt)

  return(result)
end
