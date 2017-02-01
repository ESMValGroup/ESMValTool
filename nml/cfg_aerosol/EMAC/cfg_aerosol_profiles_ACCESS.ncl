; Settings used by aerosol_vertprof diag_script
diag_script_info = True

; Selected campaign
diag_script_info@campaign = "ACCESS"

; Style
diag_script_info@styleset = "DEFAULT"
diag_script_info@style_obs = "marker"
diag_script_info@color_obs = "grey20"
; Legend location. Syntax: t[lr] ("t" for top, "l" for left, "r" for right;
; bottom not yet implemented)
; diag_script_info@legloc = "tr"

; Path to campaign data
diag_script_info@datapath = getenv("ESMValTool_OBSPATH") + "/Tier3/"+\
                            diag_script_info@campaign
diag_script_info@datapath = "$HOME/tmp"

; Cloud cover threshold (%)
diag_script_info@cc_threshold = 100.

; More descriptive names for the models/runs to be used for the "MODEL"
; placeholder in the custom strings below as well as for the legend of the
; time series plot. Use the following syntax:
; diag_script_info@modelname_MODEL
;     MODEL:
;         a "name" argument from a <model> tag in the main namelist
diag_script_info@modelname_34128 = "EMAC (MADE3)"
diag_script_info@modelname_34129 = "EMAC (MADE3)"
diag_script_info@modelname_34130 = "EMAC (MADE3)"
diag_script_info@modelname_34050 = "EMAC (MADE)"
diag_script_info@modelname_34051 = "EMAC (MADE)"
diag_script_info@modelname_34052 = "EMAC (MADE)"

; Customize plot
; Plot resources can be specified using the following syntax:
; diag_script_info@res_mod_RESOURCE = VALUE
;     RESOURCE:
;         NCL resource that is applicable to a gsn_scm_blank_plot ("map")
;     VALUE:
;         value to be specified for RESOURCE; special placeholders that will be
;         replaced include "VARIABLE", "UNITS", and "MODEL"
; diag_script_info@res_mod_tiMainString = "Vertical aerosol distribution"
; diag_script_info@res_mod_gsnLeftString = ""
; diag_script_info@res_mod_gsnCenterString = \
;     "Multi-annual median monthly means vs. flight medians"
; diag_script_info@res_mod_gsnRightString = ""
diag_script_info@res_mod_tiXAxisString = "VARIABLE  [UNITS]"
diag_script_info@res_mod_tmXMajorGrid = False
diag_script_info@res_mod_tmXMinorGrid = False
diag_script_info@res_mod_tmYMajorGrid = False
diag_script_info@res_mod_tmYMinorGrid = False
diag_script_info@res_mod_tmXTOn = False
diag_script_info@res_mod_tmYROn = False
diag_script_info@res_mod_tmXBMode = "Automatic"
; diag_script_info@res_mod_trXMinF = 1.e0
; diag_script_info@res_mod_trXMaxF = 2.e3
; diag_script_info@res_mod_trYMinF = 200 ; hPa
