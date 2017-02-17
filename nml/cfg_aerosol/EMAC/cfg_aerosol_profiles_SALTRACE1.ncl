; Settings used by aerosol_vertprof diag_script
diag_script_info = True

; Selected campaign
diag_script_info@campaign = "SALTRACE"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to campaign data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier3/"+\
                            diag_script_info@campaign

; Cloud cover threshold (%)
diag_script_info@cc_threshold = 100.

; Summary plot
diag_script_info@summary_plot = True
diag_script_info@exclude_flights = (/"130610a", "130611a", "130611b", "130617b", "130618a", "130618b", "130620a", "130621a", "130622a", "130626a", "130630a", "130630b", "130701a", "130701b", "130705a", "130708a", "130710a", "130711a", "130711b", "130712a", "130713a", "130713b", "130714a", "130714b", "130715a", "130715b"/)
