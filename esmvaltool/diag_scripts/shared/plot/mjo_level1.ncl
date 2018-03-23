;;#############################################################################
;; GENERAL ROUTINES FOR STATISTICS
;;#############################################################################
;; Please consider using or extending existing routines before adding new
;; ones.
;; Check the header of each routine for documentation.      
;;
;; Contents:
;;    function pr_u850_mean_plot
;;    function mjo_xcor_lag_plot
;;    function mjo_pr_ua_vari_plot
;;    function mjo_unvari_eof_plot
;;
;;#############################################################################

;; Temporary file updated by Python on the fly


load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/logging.ncl"

load "./diag_scripts/shared/set_operators.ncl"
load "./diag_scripts/shared/statistics.ncl"

load "./diag_scripts/aux/SAMonsoon/SAMonsoon.ncl"

undef("pr_u850_mean_plot")
function pr_u850_mean_plot(wks[1]:graphic,\
                       mean_pr[*][*] : numeric,\
                       mean_ua[*][*] : numeric,\
                            di[1]    : logical,\
                           opt[1]    : logical)
; Description:
;   Place the 2D-model fields stored in the 'storage_vault' on a set of paneled
;   figures. First entry is always the reference data set. The various logcials
;   carry switches for the plot routine behaviour.
local res1, res2, plot
begin
    ;; common resources for global contour plots
    res1 = True
    res1@gsnDraw = False                ; don't draw
    res1@gsnFrame = False               ; don't advance frame
    res1@gsnStringFontHeightF = 0.0125  ; make larger than default

    res1@lbLabelBarOn = False
    res1@lbLabelFontHeightF = 0.01           ; make labels larger
    res1@pmLabelBarOrthogonalPosF = -0.025   ; move closer to plot
    res1@lbLabelAutoStride = True            ; optimal label stride
    res1@gsnSpreadColors = True              ; use full range of colors
    res1@mpCenterLonF = 180.                 ; center plot at 180

    res1@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res1@cnLevels = ispan(2, 13, 1)

    res1@gsnAddCyclic = False

    res1@gsnStringFontHeightF = 0.02
    res1@tmYLLabelFontHeightF = 0.02
    res1@tmXBLabelFontHeightF = 0.02

    res1@mpMinLonF = di@lonrange_basic(0)     ; select a subregion
    res1@mpMaxLonF = di@lonrange_basic(1)
    res1@mpMinLatF = di@latrange_basic(0)
    res1@mpMaxLatF = di@latrange_basic(1)

    res1@mpLandFillColor = "background"       ; color of land
    res1@cnFillDrawOrder = "Predraw"

    if (isatt(opt, "gsnLeftString")) then
        res1@gsnLeftString = opt@gsnLeftString
    end if

    if (isatt(opt, "gsnCenterString")) then
        res1@gsnCenterString = opt@gsnCenterString
    end if

    if (isatt(opt, "gsnRightString")) then
        res1@gsnRightString = opt@gsnRightString
    end if

    ;;************************************************
    ;; resource list for second (contour only) data array
    ;;************************************************
    res2 = True
    res2@gsnContourZeroLineThicknessF = 2.   ; doubles thickness of zero contour

    ;; Sets negative contours to dash pattern 1
    res2@gsnContourNegLineDashPattern = 2
    res2@cnLineColor = "Black"               ; color of second contours
    res2@cnLineThicknessF = 1.
    res2@cnLevelSpacingF = 3.                ; SST contour spacing
    res2@cnInfoLabelOn = False               ; do not want
    res2@cnLineLabelsOn = True
    res2@cnLabelMasking = True
    res2@tiMainString = ""

    res2@gsnAddCyclic = False
    plot = gsn_csm_contour_map_overlay(wks, mean_pr, mean_ua, res1, res2)

    return(plot)
end

undef("mjo_xcor_lag_plot")
function mjo_xcor_lag_plot(wks[1]      : graphic,\
                           ccr_a[*][*] : numeric,\
                           ccr_b[*][*] : numeric,\
                             opt[1]    : logical)
local res1, res2, CCR1,\
    CCR2, plot
begin
    res1 = True               ; color precip
    res1@gsnDraw = False
    res1@gsnFrame = False
    res1@gsnMaximize = True
    res1@gsnPaperOrientation = "portrait"

    res1@cnFillOn = True                        ; turn on color
    res1@cnLinesOn = False
    res1@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
    res1@cnMinLevelValF = -1.0                  ; set min contour level
    res1@cnMaxLevelValF = 1.0                   ; set max contour level
    res1@cnLevelSpacingF = 0.1                  ; set contour spacing

    res1@cnLabelBarEndLabelsOn = True
    res1@cnLabelBarEndStyle = "ExcludeOuterBoxes"
    res1@cnInfoLabelOn = False

    res1@lbLabelBarOn = False     ; turn off individual cb's
    res1@vpWidthF = 0.6           ; change aspect ratio of plot
    res1@vpHeightF = 0.4

    res1@tiYAxisString = "lag (days)"
    if (isatt(opt, "tiMainString")) then
        res1@tiMainString = opt@tiMainString
    end if

    if (isatt(opt, "gsnLeftString")) then
        res1@gsnLeftString = opt@gsnLeftString
    end if

    if (isatt(opt, "gsnCenterString")) then
        res1@gsnCenterString = opt@gsnCenterString
    end if

    if (isatt(opt, "gsnRightString")) then
        res1@gsnRightString = opt@gsnRightString
    end if

    ;;************************************************
    ;; resource list for second data array
    ;;************************************************
    res2 = True              ; U
    res2@gsnDraw = False
    res2@gsnFrame = False
    res2@cnLevelSelectionMode = "ManualLevels"    ; set manual contour levels
    res2@cnMinLevelValF = -1.0
    res2@cnMaxLevelValF = 1.0
    res2@cnLevelSpacingF = 0.1
    res2@cnLineLabelsOn = True
    res2@gsnContourZeroLineThicknessF = 0.      ; Eliminate 0 line
    res2@gsnContourNegLineDashPattern = 1       ; negative contours dash pattern
    res2@cnInfoLabelOn = False

    CCR1 = ccr_a             ; possible smooth and delete of attribute
    CCR2 = ccr_b

    if (opt .and. isatt(opt, "smth9") .and. abs(opt@smth9).eq.0.25) then
        CCR1 = smth9(CCR1, 0.50, opt@smth9, False)
        CCR2 = smth9(CCR2, 0.50, opt@smth9, False)
    end if
    delete(CCR1@long_name)
    plot = gsn_csm_contour(wks, CCR1, res1)       ; contour the variable

    delete(CCR2@long_name)
    plt2 = gsn_csm_contour(wks, CCR2, res2)       ; contour the variable

    overlay(plot, plt2)

    delete(CCR1)                                ; size may change
    delete(CCR2)

    return(plot)
end

undef("mjo_pr_ua_vari_plot")
function mjo_pr_ua_vari_plot(wks[1]    : graphic,\
                           ccr_a[*][*] : numeric,\
                           ccr_b[*][*] : numeric,\
                             opt[1]    : logical,\
                              di[1]    : logical)
local res1, res2, plot
begin
    ;;************************************************
    ;; resource list for first (color) data array
    ;;************************************************
    res1 = True
    res1@gsnDraw = False   ; don't draw
    res1@gsnFrame = False  ; don't advance frame

    res1@lbLabelBarOn = False

    res1@lbLabelFontHeightF = 0.01             ; make labels larger
    res1@pmLabelBarOrthogonalPosF = -0.025     ; move closer to plot
    res1@lbLabelAutoStride = True              ; optimal label stride
    res1@gsnSpreadColors = True                ; use full range of colors
    res1@mpCenterLonF = 180.                   ; center plot at 180

    res1@gsnStringFontHeightF = 0.02
    res1@tmYLLabelFontHeightF = 0.02
    res1@tmXBLabelFontHeightF = 0.02

    res1@mpMinLonF = di@lonrange_seasonal(0)   ; select a subregion
    res1@mpMaxLonF = di@lonrange_seasonal(1)
    res1@mpMinLatF = di@latrange_seasonal(0)
    res1@mpMaxLatF = di@latrange_seasonal(1)

    res1@cnLevelSelectionMode = "ExplicitLevels"  ; set explicit contour levels
    res1@cnLevels = di@cnLevels_pr_vari

    res1@gsnAddCyclic = False

    res1@mpLandFillColor = "background"       ; color of land
    res1@cnFillDrawOrder = "Predraw"

    if (isatt(opt, "gsnLeftString")) then
        res1@gsnLeftString = opt@gsnLeftString
    end if

    if (isatt(opt, "gsnCenterString")) then
        res1@gsnCenterString = opt@gsnCenterString
    end if

    if (isatt(opt, "gsnRightString")) then
        res1@gsnRightString = opt@gsnRightString
    end if

    ;;************************************************
    ;; resource list for second (contour only) data array
    ;;************************************************
    res2 = True
    res2@gsnContourZeroLineThicknessF = 2.   ; doubles thickness of zero contour

    ;; Sets negative contours to dash pattern 1
    res2@gsnContourNegLineDashPattern = 2
    res2@cnLineColor = "Black"               ; color of second contours
    res2@cnLineThicknessF = 1.

    res2@gsnAddCyclic = False

    res2@cnLevelSpacingF = 3.      ; SST contour spacing
    res2@cnInfoLabelOn = False     ; do not want
    res2@cnLineLabelsOn = True
    res2@cnLabelMasking = True
    res2@tiMainString = ""

    plot = gsn_csm_contour_map_overlay(wks, ccr_a, ccr_b, res1, res2)
    return(plot)
end

;---------------------------------
undef("mjo_unvari_eof_plot")
function mjo_unvari_eof_plot(wks[1]    : graphic, \
                             eof[*][*] : numeric, \
                             opt[1]    : logical,
                              di[1]    : logical)
local res, lot
begin
    ;; EOF patterns
    res = True
    res@gsnDraw = False                 ; don't draw yet
    res@gsnFrame = False                ; don't advance frame yet
    res@gsnSpreadColors = True          ; spread out color table
    res@gsnStringFontHeightF = 0.015    ; make larger than default

    res@mpFillOn = False                ; turn off map fill
    res@mpMinLatF = di@latrange_eof(0)  ; zoom in on map
    res@mpMaxLatF = di@latrange_eof(1)
    res@mpMinLonF = di@lonrange_eof(0)  ; select a subregion
    res@mpMaxLonF = di@lonrange_eof(1)

    res@mpCenterLonF = 180.

    res@gsnAddCyclic = False

    res@cnFillOn = True          ; turn on color fill
    res@cnLinesOn = False        ; True is default
    res@cnLineLabelsOn = False   ; True is default
    res@lbLabelBarOn = False     ; turn off individual lb's

    res@lbLabelFontHeightF = 0.01   ; make labels larger

    res@cnLevelSelectionMode = "ExplicitLevels"   ; set explicit contour levels
    res@cnLevels = di@cnLevels

    res@gsnStringFontHeightF = 0.018
    res@tmYLLabelFontHeightF = 0.018
    res@tmXBLabelFontHeightF = 0.018

    if (isatt(opt, "gsnLeftString")) then
        res@gsnLeftString = opt@gsnLeftString
    end if

    if (isatt(opt, "gsnCenterString")) then
        res@gsnCenterString = opt@gsnCenterString
    end if

    if (isatt(opt, "gsnRightString")) then
        res@gsnRightString = opt@gsnRightString
    end if

    plot = gsn_csm_contour_map_ce(wks, eof(:, :), res)
    return(plot)
end
