diag_script_info = True
diag_script_info@styleset = "CMIP5"

;info for vpline.ncl:
diag_script_info@e06fig05_latsel = (/ (/80./),  \     ;panel1
                                    (/0./),   \     ;panel2
                                    (/-80./)  /)    ;panel3
  
diag_script_info@e06fig05_monthsel = (/ "MAR", \  ;panel1
                                        "MAR", \  ;panel2
                                        "OCT"  /) ;panel3

diag_script_info@e06fig05_XMax = (/8., 12., 8./)
diag_script_info@e06fig05_XMin = (/0.,  0., 0./)

diag_script_info@e06fig05_levmax = 200.
diag_script_info@e06fig05_levmin = 1.
diag_script_info@e06fig05_multimean = True

; data sets to be skipped when calculating the multi-model mean
; (= observations)

diag_script_info@models_to_skip = (/"ESACCI-OZONE", "BDBP"/)

