; Settings used by most mjo plot_types
diag_script_info = True

; Settings specific for the precip basic plot_type
diag_script_info@latrange_basic = (/-35, 35/)
diag_script_info@lonrange_basic = (/40, 320/)

diag_script_info@cn_levels_pr = ispan(2,13,1) 
diag_script_info@cn_levels_u850 = 3 

; Settings specific for the precip seasonal plot_type
diag_script_info@latrange_seasonal = (/-35, 35/)
diag_script_info@lonrange_seasonal = (/40, 320/)

; Settings shared by the precip global domain/intensity plot_type
diag_script_info@summer_season = "MJJASO"
diag_script_info@winter_season = "NDJFMA"
diag_script_info@season_names  = (/ "May-Oct", "Nov-Apr"/)

diag_script_info@seasons = (/"MJJASO","NDJFMA"/)

diag_script_info@season_mjo   = (/ (/ 5, 6, 7, 8, 9,10/) \       ; May-Oct [summer]
              , (/ 1, 2, 3, 4,11,12/) /)      ; Nov-Apr [winter]
