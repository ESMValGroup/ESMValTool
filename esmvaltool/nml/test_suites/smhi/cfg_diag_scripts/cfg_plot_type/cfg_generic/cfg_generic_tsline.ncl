plot_type_info = True

;plot_type_info@ts_levrange = (/200, 150/)  ; Pressure level range if 3D

plot_type_info@ts_lonrange = (/0, 360/)   ; Longitude range if 2D or 3D
plot_type_info@ts_latrange = (/-20, 40/)  ; Latitude  range if 2D or 3D

plot_type_info@ts_anom      = False          ; (optional) Plot anomalies
; plot_type_info@ts_ymax      = "Not defined"  ; (optional) Y axis maximum
; plot_type_info@ts_ymin      = "Not defined"  ; (optional) Y axis minimum
; plot_type_info@ts_xmax_year = "Not defined"  ; (optional) Start year
; plot_type_info@ts_xmin_year = "Not defined"  ; (optional) End year
; plot_type_info@ts_smooth    = "Not defined"  ; (optional) Smooth...
; plot_type_info@ts_multimean = "Not defined"  ;

; Reference model (typically observations, reanalysis)
plot_type_info@refModel = (/"TRMM-3B43-v7-0.25deg"/)
