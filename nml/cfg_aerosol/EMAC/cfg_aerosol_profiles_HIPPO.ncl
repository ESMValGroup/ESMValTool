; Settings used by aerosol_vertprof diag_script
diag_script_info = True

; Selected campaign
diag_script_info@campaign = "HIPPO"

; Style
diag_script_info@styleset = "DEFAULT"

; Path to campaign data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier3/"+\
                            diag_script_info@campaign

; Cloud cover threshold (%)
diag_script_info@cc_threshold = 5.
