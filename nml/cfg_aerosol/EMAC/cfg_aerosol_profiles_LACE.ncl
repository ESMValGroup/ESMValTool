; Settings used by aerosol_vertprof diag_script
diag_script_info = True

; Selected campaign
diag_script_info@campaign = "LACE"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to campaign data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier2/"+\
                            diag_script_info@campaign

; Cloud cover threshold (%)
diag_script_info@cc_threshold = 5.
