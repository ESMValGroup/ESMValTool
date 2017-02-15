diag_script_info = True

 diag_script_info@e13fig02_latrange = (/(/-90, 90/),  \   ;panel1
                                        (/-90, -60/)/)    ;panel2
 diag_script_info@styleset = "CMIP5"

 diag_script_info@e13fig02_season = (/ "ANN", "OCT"/)

 diag_script_info@e13fig02_XMin = (/1960., 1960./)
 diag_script_info@e13fig02_XMax = (/2010., 2010./)
 diag_script_info@e13fig02_YMin = (/270., 150./)
 diag_script_info@e13fig02_YMax = (/340., 500./)
 diag_script_info@e13fig02_legend_MMM = "False"
 diag_script_info@e13fig02_legend = "True"

 diag_script_info@list_chem_mod = (/"All"/)
 diag_script_info@list_nochem_mod = (/"None"/)

 ; ----------------------------------------------------------------------
 ; mask missing values from all data sets before calculating averages?
 diag_script_info@apply_masking = True
 ; if masking is selected, all data must be regridded to a common grid
 ; (coarsest, finest, ref, 1x1; default = 1x1)
 ; if "ref" is selected, reference data set must be specified
 diag_script_info@grid = "1x1"
 diag_script_info@refmodel = "ESACCI-OZONE"
 ; ----------------------------------------------------------------------

