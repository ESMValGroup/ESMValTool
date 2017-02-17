; Settings used by aerosol_vertprof diag_script
diag_script_info = True

; Selected campaign
diag_script_info@campaign = "ACCESS-2"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to campaign data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier3/"+\
                            diag_script_info@campaign

; Cloud cover threshold (%)
diag_script_info@cc_threshold = 100.

; Summary plot
diag_script_info@summary_plot = True
diag_script_info@exclude_flights = (/"140501a", "140501b", "140502a", "140502b", "140503a", "140503b", "140503c", "140517a", "140518a", "140518b", "140519a", "140519b", "140520a"/)


