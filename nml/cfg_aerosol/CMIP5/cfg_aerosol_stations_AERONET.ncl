; Settings used by aerosol_stations diag_script
diag_script_info = True

; Selected network (AERONET, CASTNET, EANET, EMEP, IMPROVE)
diag_script_info@network  = "AERONET"

; Time average (monthly or yearly)
diag_script_info@time_avg = "yearly"

; Style
diag_script_info@styleset = "CMIP5"

; Path to station data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier2/"+\
                            diag_script_info@network

; Use log axes in scatter plot
diag_script_info@scatter_log = True

; Legend outside the plot
diag_script_info@legend_outside = True
