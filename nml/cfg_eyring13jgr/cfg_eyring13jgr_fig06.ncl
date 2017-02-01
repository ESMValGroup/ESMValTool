diag_script_info = True

 diag_script_info@e13fig06_latrange = (/(/-90,90/),  \   ;panel1
                                        (/-90,-60/)/)    ;panel2

 diag_script_info@e13fig06_season = (/"ANN", "SON"/)
 diag_script_info@styleset = "CMIP5"
 diag_script_info@e13fig06_baseline = 1980
 diag_script_info@e13fig06_baseline_adj = "True"
 diag_script_info@e13fig06_mod_plot = "MMT" ; "MMT" = plot of the MultiModel mean of each scenario
                                            ;         and selection "list_chem_mod" and "list_nochem_mod"
                                            ; "IMT" = plot of each single model trend
                                            ; "RAW" = plot of each model as raw data
 diag_script_info@e13fig06_mod_plot_CI = "True" ; if MMT the Confidence Interval can be plotted

 diag_script_info@e13fig06_mod_plot_PI = "True" ; if MMT the Prediction Interval can be plotted

 ;diag_script_info@e13fig06_XMin = (/1960., 1960./)
 ;diag_script_info@e13fig06_XMax = (/2010., 2010./)
 ;diag_script_info@e13fig06_YMin = (/260., 150./)
 ;diag_script_info@e13fig06_YMax = (/330., 500./)
 diag_script_info@e13fig06_labels_exp_esp = "True" ; if "True" experiment name specified in IMT and RAW case
 ;diag_script_info@e13fig06_lev =

 diag_script_info@list_chem_mod = (/"CCSM4", "IPSL-CM5A-LR"/) ;"All" if not selection required
 diag_script_info@list_nochem_mod = (/"GFDL-CM3", "MIROC-ESM-CHEM"/); comment if not selection required

