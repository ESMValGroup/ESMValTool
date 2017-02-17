; Settings used by most mjo plot_types
diag_script_info = True

diag_script_info@nameSeason  = (/ "winter", "summer", "annual" /) 
diag_script_info@nameVariSeason  = (/ "winter", "summer"/) 

diag_script_info@eio_latrange = (/-10, 5/)
diag_script_info@eio_lonrange = (/75, 100/)

diag_script_info@global_latrange = (/-30, 30/)

diag_script_info@prop_lag_latrange =(/-10, 10/)
diag_script_info@prop_lag_lonrange = (/80, 100/)


; Settings specific for the precip seasonal plot_type
diag_script_info@latrange_seasonal = (/-35, 35/)
diag_script_info@lonrange_seasonal = (/40, 320/)

;set pr 20-100day bandpass filtered variance level
diag_script_info@cnLevels_pr_vari      = (/2,4,8,12,16,24,32,64,96/)

;Setting for EOF analysis
diag_script_info@latrange_eof = (/-30, 30/)
diag_script_info@lonrange_eof = (/40, 320/)
diag_script_info@num_eof        =3
diag_script_info@name_eof       = (/ "EOF1", "EOF2", "EOF3" /)
diag_script_info@cnLevels      = (/-0.8,-0.6,-0.4,-0.2,0,0.2,0.4,0.6,0.8/) 

;Setting for mjo life cycle
diag_script_info@latrange_life = (/-20, 20/)
diag_script_info@lonrange_life = (/40, 320/)
diag_script_info@season_life_cycle =(/ "summer", "winter"/)

;Setting for CEOF analysis
diag_script_info@latrange_ceof = (/-15, 15/)
diag_script_info@num_ceof        =2
diag_script_info@name_ceof       = (/ "EOF1", "EOF2" /)
diag_script_info@trYMaxF_ceof    =0.14
diag_script_info@trYMinF_ceof    =-0.12

;Set specific for cross spectran between olr and u850
diag_script_info@latrange_cross = (/-15, 15/)
diag_script_info@segLen         = 256
diag_script_info@segOverLap     = -50
diag_script_info@spectra_type   =(/ "Symmetric", "Asymmetric"/)

;setting for 20-100day bandpass filter
diag_script_info@prop_filter_min = 20
diag_script_info@prop_filter_max = 100
diag_script_info@prop_filter_mxlag = 25
diag_script_info@ihp = 2                     ; bpf=>band pass filter
diag_script_info@nWgt = 201
diag_script_info@sigma = 1.0                 ; Lanczos sigma

diag_script_info@max_days = 100.
diag_script_info@min_days = 20.
diag_script_info@fca = 1./diag_script_info@max_days
diag_script_info@fcb = 1./diag_script_info@min_days

diag_script_info@isv_latrange = (/-10, 10/)
