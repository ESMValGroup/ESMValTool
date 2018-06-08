; Settings used by aerosol_vertprof diag_script
diag_script_info = True

; Selected campaign
diag_script_info@campaign = "SALTRACE-SP2"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to campaign data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier3/"+\
                            "SALTRACE/"

; Cloud cover threshold (%)
diag_script_info@cc_threshold = 100.

; Summary plot
diag_script_info@summary_plot = True
diag_script_info@exclude_flights = (/"130610a", "130611a", "130611b", "130612a", "130612b", "130614a", "130614b", "130617a", "130617b", "130618a", "130618b", "130712a", "130713a", "130713b", "130714a", "130714b", "130715a", "130715b"/)
