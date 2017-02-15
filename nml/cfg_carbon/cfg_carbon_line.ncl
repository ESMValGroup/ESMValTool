diag_script_info = True
diag_script_info@styleset = "CMIP5"

;info for carbon_tsline.ncl
diag_script_info@ts_minlat = -30
diag_script_info@ts_maxlat = 30
diag_script_info@ts_minlon = 0
diag_script_info@ts_maxlon = 360
diag_script_info@multi_model_mean = False
diag_script_info@ts_maxyear = "1989"
diag_script_info@ts_minyear = "1850"
diag_script_info@time_avg   = "yearly"
diag_script_info@ts_anomaly = "anom"

if (any(variables.eq.(/"nbp","cumnbp"/))) then
    diag_script_info@area_opper = "sum"
end if





