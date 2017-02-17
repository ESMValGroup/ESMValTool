; Settings used by aerosol_stations plot_types
diag_script_info = True

; Selected network (AERONET, CASTNET, EANET, EMEP, IMPROVE)
diag_script_info@network  = "IMPROVE"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to station data
diag_script_info@datapath = "/nobackup/rossby17/sm_maeva/Data/EMAC/"+\
                            diag_script_info@network

; Use log axes in scatter plot
diag_script_info@scatter_log = True

; Lat-lon limits for aerosol_stations_map
diag_script_info@minlon = -170.
diag_script_info@maxlon =  -50.
diag_script_info@minlat =   10.
diag_script_info@maxlat =   70.
