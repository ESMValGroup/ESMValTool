;
; Requires: none
;
info = True
info@derived = True
info@long_name="Ozone"
info@units="ppmv"
info@vcp_cnlvl=(/0,1,2,3,4,5,6,7,8,9,10,11,12,13,14/) ; O3
info@vcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"
;info@refModel ="NIWA-Obs"


info@ts_latrange = (/-20,-20/)
info@ts_levrange = (/100,100/)
info@ts_lonrange = (/0,360/)
;info@ts_ymax = 220
;info@ts_ymin = 180
;info@ts_xmax_year = 2010
;info@ts_xmin_year = 1975
info@ts_anom = False
;info@ts_anom = True
;info@ts_smooth = False
;info@ts_multimean = True

info@profile_stnid = "./plot_type/input_data/station_id"
info@profile_allstn = False
info@profile_raobs = "./plot_type/input_data/RAOBS.nc"
info@fig12C_xmin_year =1960
info@fig12C_xmax_year =2100

info@fig12C_lat_sel_min = (/60.,-90.,60,-90.,-60.,30.,-60.,30.,-10.,-10/)
info@fig12C_lat_sel_max = (/90.,-60,90.,-60.,-30.,60.,-30.,60.,10.,10/)
info@fig12C_lev_sel = (/50.,50.,5.,5.,5.,5.,70.,70.,5.,70./)
info@fig12C_season_lev_sel = (/(/"FEB","MAR","APR"/),(/"SEP","OCT","NOV"/),(/"FEB","MAR","APR"/),(/"SEP","OCT","NOV"/),(/"ANN","ANN","ANN"/),(/"ANN","ANN","ANN"/),(/"ANN","ANN","ANN"/),(/"ANN","ANN","ANN"/),(/"ANN","ANN","ANN"/),(/"ANN","ANN","ANN"/),(/"ANN","ANN","ANN"/)/)

info@fig12C_ymax_lev = 50
info@fig12C_ymin_lev = 10
;info@fig12C_obs_name = (/"OBS"/)
;info@fig12C_obs_file = "./plot_type/input_data/OBS/CCMVal1_1992-2005_OBS_1_T2Mz_Cly.nc"
;info@fig12C_obs_stdv_file = "./plot_type/input_data/OBS/CCMVal1_1992-2005_OBS_1_T2Mz_Cly_std.nc"

info@fig05_yr_min = 1980
info@fig05_yr_max = 1999
info@fig05_lat_sel = (/80.,0.,-80./)
info@fig05_month_lat_sel = (/"MAR","MAR","OCT"/)
info@fig05_lev_sel = (/50.,50./)
info@fig05_month_lev_sel = (/"MAR","OCT"/)
info@fig05_obs_name = (/"HALOE"/)
info@fig05_refModel=(/"SOCOL"/)
info@fig05_obs_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"
info@fig05_obs_stdv_file = "./plot_type/input_data/OBS/haloe_lat_climat.nc"

;/////////////ozonesondes Logan////////////////////////////////////////
info@zone="tropo"  ;"tropo"=troposphere  "strato"=stratosphere
info@period="all"   ;("specific":climatologies over all the years specified in the namelist; "all": climatologies over the same years as in the observation)
;/////////////////////////////////////////////////////////////////////////////////

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
;;    Converts o3 units "mol/mol" to "ppmv"
local o3,o3_flt
begin
  data_new = read_data(index, variable, field_number)
  o3 = extract_data(index,data_new,-1,0,0)
  o3 = (/o3 * 1.e6/)
  o3@units = info@units
  if (typeof(o3).eq."double") then
    o3_flt = dble2flt(o3)
  else
    o3_flt = o3
  end if
  modify_data(index, data_new, o3_flt)
  delete(o3)
  delete(o3_flt)
  return(data_new)
end

