;
; Requires: none
;
info = True
info@derived = False
info@long_name="Temperature"
info@units="K"
info@vcp_cnlvl=(/180,185,190,200,210,220,230,240,250,260,270,280,290,295,300/) ;T(ta)
info@vcp_dcnlvl=(/-9,-7,-5,-4,-3,-2,-1,0,1,2,3,4,5,7,9/)
info@vcp_ptop=30
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"
info@profile_stnid = "./plot_type/input_data/station_id"
info@profile_allstn = False
info@profile_raobs = "./plot_type/input_data/RAOBS.nc"
info@refModel = "ERA40"

info@ts_latrange = (/-10,10/)
info@ts_levrange = (/150,70/)
info@ts_lonrange = (/0,360/)
info@ts_ymax = 220
info@ts_ymin = 180
info@ts_xmax_year = 2005
info@ts_xmin_year = 1980
info@ts_anom = False

info@ml_latrange = (/-20,20/)
info@ml_levrange = (/80,80/)

;info@fig04_trend_year_min = 1980  ;if not specified or if out the range all year of the run would be used
;info@fig04_trend_year_max = 1999  ;to calculated the trend anomalies
info@fig04_trend_year_min = 1990  ;if not specified or if out the range all year of the run would be used
info@fig04_trend_year_max = 1994  ;to calculated the trend anomalies
info@fig04_detrend_year_min = 1980  ;if not specified or if out the range all year of the run would be used
info@fig04_detrend_year_max = 1999  ;to calculated the trend anomalies
info@fig04_lat_min = (/60.,-90.,-90./)
info@fig04_lat_max = (/90.,-60.,90./)
info@fig04_lev_sel = 50.
info@fig04_season = (/(/"FEB","MAR","APR"/),(/"SEP","OCT","NOV"/),(/"ANN","ANN","ANN"/)/)
;info@fig04_xmin = 1960           ; max and min XAxis in the plot
;info@fig04_xmax = 2006           ; if not specified 1960 and 2100 would be used!
info@fig04_xmin = 1980           ; max and min XAxis in the plot
info@fig04_xmax = 1999           ; if not specified 1960 and 2100 would be used!

info@fig07_lat_min = (/0./)
info@fig07_lat_max = (/0./)
info@fig07_lev_sel = (/100./)
info@fig07_refModel = (/"ERA40"/)

info@fig01_lat_max = (/90.,90.,-60.,-60./)
info@fig01_lat_min = (/60.,60.,-90.,-90./)
info@fig01_season = (/"DJF","MAM","JJA","SON"/)
info@fig01_refModel = (/"ERA40"/)          ;(/"ERA-40"/)
info@fig01_climObs  =(/"NCEP","UKMO"/)      ;Climatological Observation file
info@fig01_climObs_file=(/"./plot_type/input_data/OBS/CCMVal1_1980_2000_NCEP_Obs_T2Mz_ta.nc",\
"./plot_type/input_data/OBS/CCMVal2_1992-2001_UKMO_Obs_C2Mz_ta.nc"/)
;C2Mz file

info@fig03_lev_sel   =     (/50./)
info@fig03_lat_min   =   (/60.,-90./)
info@fig03_lat_max   =   (/90.,-60./)
info@fig03_month_sel =     (/(/"FEB","MAR"/),(/"AUG","SEP"/)/)
info@fig03_vs_lev_sel =      (/100./)
info@fig03_vs_lat_min =     (/40.,-80./)
info@fig03_vs_lat_max =     (/80.,-40./)
info@fig03_vs_month_sel =  (/(/"JAN","FEB"/),(/"JUL","AUG"/)/)

info@fig03_vs_var ="vt100"
info@fig03_vs_field="T1Ms"
info@fig03_xmin_val=(/0.,0./)
info@fig03_xmax_val=(/30.,20./)
info@fig03_ymin_val=(/200.,180./)
info@fig03_ymax_val=(/225.,205./)
info@fig03_obs  =(/"ERA40"/)
info@fig03_obs_file=(/"./plot_type/input_data/OBS/CCMVal1_1979_2005_ERA40_Obs_T0As_ta.nc"/)
info@fig03_obs_file_vs=(/"./plot_type/input_data/OBS/CCMVal1_1979_2005_ERA40_Obs_T0As_vt100.nc"/)


;####

