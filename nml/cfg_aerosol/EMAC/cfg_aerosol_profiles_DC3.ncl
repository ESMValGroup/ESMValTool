; Settings used by aerosol_vertprof diag_script
diag_script_info = True

; Selected campaign
diag_script_info@campaign = "DC3"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to campaign data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier3/"+\
                            diag_script_info@campaign

; Cloud cover threshold (%)
diag_script_info@cc_threshold = 100.

; Summary plot
diag_script_info@summary_plot = True
diag_script_info@exclude_flights = (/"120525a", "120525b", "120525c", "120526a", "120526b", "120527a", "120527b", "120617a", "120617b", "120618a", "120618b", "120619a", "120619b"/)
