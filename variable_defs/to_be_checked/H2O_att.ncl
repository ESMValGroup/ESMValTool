;
; Requires: none
;
info = True
info@derived = True
info@long_name="Water Vapor"
info@units="ppmv"
info@vcp_cnlvl=(/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14/)*2 ; H2O
info@vcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
;info@vcp_ptop=30
info@vcp_ptop=195.
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"
;info@refModel = "HALOE"
info@seasons=(/"DJF"/)

info@ts_latrange = (/-20,20/)
info@ts_levrange = (/80,80/)
info@ts_lonrange = (/0,360/)
info@ts_anom = True
info@ts_multimean = True
info@ts_smooth = False

info@ml_latrange = (/-20,20/)
info@ml_levrange = (/80,80/)
info@ml_lonrange = (/0,360/)
info@ml_anom = False
;info@ml_obs_name = (/"HALOE"/)
info@ml_obs_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"
info@ml_obs_stdv_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"
info@ml_stdv_plot = "HALOE"
info@ml_multimean = True

info@fig08_lat_min = -10.0
info@fig08_lat_max = 10.0
info@fig08_obs_name = (/"HALOE"/)
info@fig08_obs_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"

;info@fig05_xmax_lat=8.
;info@fig05_xmin_lat=0.
;info@fig05_ymax_lev=10.
;info@fig05_ymin_lev=0.
info@fig05_lat_sel = (/80.,0.,-80./)
info@fig05_month_lat_sel = (/"MAR","MAR","OCT"/)
info@fig05_lev_sel = (/50.,50./)
info@fig05_month_lev_sel = (/"MAR","OCT"/)
info@fig05_obs_name = (/"HALOE"/)
;info@fig05_refModel=(/"UMETRAC"/)
info@fig05_obs_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"
info@fig05_obs_stdv_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"

info@fig07_lat_min = (/0./)
info@fig07_lat_max = (/0./)
info@fig07_lev_sel = (/80./)
info@fig07_obs_name = (/"HALOE"/)
;info@fig05_refModel=(/"UMETRAC"/)
info@fig07_obs_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"
info@fig07_obs_stdv_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"


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
;;    Converts h2o units "mol/mol" to "ppmv"
local h2o,h2o_flt
begin
  data_new = read_data(index, variable, field_number)
  h2o = extract_data(index,data_new,-1,0,0)
  h2o = (/h2o * 1.e6/)
  h2o@units = info@units
  if (typeof(h2o).eq."double") then
    h2o_flt = dble2flt(h2o)
  else
    h2o_flt = h2o
  end if
  modify_data(index, data_new, h2o_flt)
  delete(h2o)
  delete(h2o_flt)
  return(data_new)
end

