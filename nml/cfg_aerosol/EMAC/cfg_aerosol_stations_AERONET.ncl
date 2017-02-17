; Settings used by aerosol_stations diag_script
diag_script_info = True

; Selected network (AERONET, CASTNET, EANET, EMEP, IMPROVE)
diag_script_info@network  = "AERONET"

; Time average (monthly or yearly)
diag_script_info@time_avg = "monthly"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to station data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier2/"+\
                            diag_script_info@network

; Legend inside the plot
diag_script_info@legend_outside = False

; Use log axes in scatter plot
diag_script_info@scatter_log = True
