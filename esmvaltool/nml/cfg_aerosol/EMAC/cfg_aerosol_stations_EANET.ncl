; Settings used by aerosol_stations diag_script
diag_script_info = True

; Selected network (AERONET, CASTNET, EANET, EMEP, IMPROVE)
diag_script_info@network  = "EANET"

; Time average (monthly or yearly)
diag_script_info@time_avg = "monthly"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to station data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier3/"+\
                            diag_script_info@network

; Use log axes in scatter plot
diag_script_info@scatter_log = True

; Legend inside the plot
diag_script_info@legend_outside = False

; Lat-lon limits for aerosol_stations_map
diag_script_info@minlon =   50.
diag_script_info@maxlon =  170.
diag_script_info@minlat =  -10.
diag_script_info@maxlat =   60.
