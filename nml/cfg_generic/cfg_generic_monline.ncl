; contour level settings for various plot_types
diag_script_info = True

diag_script_info@ml_lonrange = (/0, 111/)  ; Longitude range if 2D or 3D
diag_script_info@ml_latrange = (/0, 31/)   ; Latitude range if 2D or 3D

 ; diag_script_info@ml_levrange      = "Not defined"  ; Pressure level range if 3D
 ; diag_script_info@ml_ymax          = "Not defined"  ; (optional) y axis maximum for time series plots
 ; diag_script_info@ml_ymin          = "Not defined"  ; (optional) y axis minimum for time series plots
 ; diag_script_info@ml_xmax_year     = "Not defined"  ; (optional) start year for time series plots
 ; diag_script_info@ml_xmin_year     = "Not defined"  ; (optional) end year for time series plots
 ; diag_script_info@ml_anom          = "Not defined"  ; (optional)  plot anomalies.
 ; diag_script_info@ml_obs_file      = "Not defined"  ; (optional) plot_type specific observational data
 ; diag_script_info@ml_obs_stdv_file = "Not defined"
 ; diag_script_info@ml_obs_name      = "Not defined"  ; (optional) plot_type specific observational data
 ; diag_script_info@ml_stdv_plot     = "Not defined"

; Reference model (e.g., observations, used by some plot_types)
diag_script_info@refModel = (/"TRMM-3B43-v7-0.25deg"/)
