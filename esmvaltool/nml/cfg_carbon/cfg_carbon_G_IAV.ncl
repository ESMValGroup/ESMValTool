diag_script_info = True
diag_script_info@styleset = "CMIP5"

;info for carbon_tsline.ncl
diag_script_info@ts_latrange  = (/-90,90/)  ; Latitude
diag_script_info@ts_lonrange  = (/0,360/)   ; Longitude 
diag_script_info@ts_latrange2 = (/-30,30/)  ; Latitude of second VAR
diag_script_info@ts_lonrange2 = (/0,360/)   ; Longitude of second VAR
diag_script_info@ts_anom      = True        ; (optional) calculate anomalies
diag_script_info@ts_detr      = True        ; de-trend timeseries

;info for carbon_corr_2var.ncl
diag_script_info@ec_v2_units  = "K"         ; x-axes units
diag_script_info@ec_2range    = "tas"       ; speifiy x-axes variable
diag_script_info@ec_units     = "GtC"       ; y-axes units
diag_script_info@ec_timerange = (/1960,2005/) ; y-axes units
diag_script_info@ec_anom      = True        ; (optional) calculate anomalies
diag_script_info@ec_detr      = False       ; (optional) de-trend timeseries
diag_script_info@ec_volc      = (/"1963","1964",\  ;(optional) exclude Mt. Agung +1 year  
                                  "1982","1983",\  ;(optional) exclude El Chichon +1 year  
                                  "1992","1993"/)  ;(optional) exclude Mt. Pinatubo +1 year  
diag_script_info@scatter_log  = False       ; set logaritmic axes in scatterplot.ncl
