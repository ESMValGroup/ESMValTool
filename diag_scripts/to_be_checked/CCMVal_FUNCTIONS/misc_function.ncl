;##################################################################################
;# This NCL script contains the new functions for the CCMVal Diagnostic Tool      #
;##################################################################################
;# Author:                                                                        #
;#                                                                                #
;# Irene Cionni (i.cionni@niwa.co.nz)                                             #
;##################################################################################
;# The NCL script is called by the NCL shell scripts E06FIGXX.ncl                 #
;##################################################################################
;Function to calculate the weight area using the functionarea_cell
undef("month_sel")
function month_sel(month_names:string)
local month_names, month_n
begin
    if (.not.any(ismissing(ind(month_names.eq."ANN")))) then
        month_n = new(12,"integer")
        month_n = (/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)
    else
        month_n = new(dimsizes(month_names),"integer")
        do i=0, dimsizes(month_names)-1
            if (.not.ismissing(ind(month_names(i).eq."JAN"))) then
                month_n(i)=1
            end if
            if (.not.ismissing(ind(month_names(i).eq."FEB"))) then
                month_n(i)=2
            end if
            if (.not.ismissing(ind(month_names(i).eq."MAR"))) then
                month_n(i)=3
            end if
            if (.not.ismissing(ind(month_names(i).eq."APR"))) then
                month_n(i)=4
            end if
            if (.not.ismissing(ind(month_names(i).eq."MAY"))) then
                month_n(i)=5
            end if
            if (.not.ismissing(ind(month_names(i).eq."JUN"))) then
                month_n(i)=6
            end if
            if (.not.ismissing(ind(month_names(i).eq."JUL"))) then
                month_n(i)=7
            end if
            if (.not.ismissing(ind(month_names(i).eq."AUG"))) then
                month_n(i)=8
            end if
            if (.not.ismissing(ind(month_names(i).eq."SEP"))) then
                month_n(i)=9
            end if
            if (.not.ismissing(ind(month_names(i).eq."OCT"))) then
                month_n(i)=10
            end if
            if (.not.ismissing(ind(month_names(i).eq."NOV"))) then
                month_n(i)=11
            end if
            if (.not.ismissing(ind(month_names(i).eq."DEC"))) then
                month_n(i)=12
            end if
        end do
    end if

    return(month_n)

end

;*************************************************
undef("lat_names")
function lat_names(lat_val:numeric)
local lat_name, lat_val
begin
    lat_name = new((/dimsizes(lat_val)/),"string")
    do i=0, dimsizes(lat_val)-1
        if (lat_val(i).gt.0.) then
            lat_name(i) = lat_val(i)+"N"
        end if
        if (lat_val(i).eq.0.) then
            lat_name(i) = "Equator"
        end if
        if (lat_val(i).lt.0.) then
            lat_name(i) = -lat_val(i)+"S"
        end if
    end do

    return (lat_name)

end

;***********************color*************************************
undef("CCMVal2_levelsColor")
function CCMVal2_levelsColor(num_levels:integer)
local colori, colorilev
begin
    colorilev = new((/num_levels/),"string")

    colori=new(20,"string")

    colori(0) = "(/0.105, 0.67, 0.8/)"
    colori(1) = "(/0.737, 0.459, 0.718/)"
    colori(2) = "(/0.922, 0.722, 0.012/)"
    colori(3) = "(/1., 0., 0./)"
    colori(4) = "(/0.505, 0.004, 0.494/)"
    colori(5) = "(/0.94, 0.467, 0./)"
    colori(6) = "(/0., 0.534, 0.216/)"
    colori(7)=  "(/0.6, 0.811, 0.086/)"
    colori(8) = "(/0.05, 0.219, 0.588/)"
    colori(9) = "(/0.105, 0.67, 0.8/)"
    colori(10) = "(/0.737, 0.459, 0.718/)"
    colori(11) = "(/0.922, 0.722, 0.012/)"
    colori(12) = "Black"
    ;;colori(12) = "(/1., 0., 0./)"
    colori(13) = "(/0.505, 0.004, 0.494/)"
    colori(14) = "(/0.94, 0.467, 0./)"
    colori(15) = "(/0., 0.534, 0.216/)"
    colori(16)= "(/0.6, 0.811, 0.086/)"
    colori(17) = "(/0.05, 0.219, 0.588/)"
    colori(18) = "(/0.105, 0.67, 0.8/)"
    colori(19) = "(/0.737, 0.459, 0.718/)"

    colorilev=colori(0:num_levels-1)
    return(colorilev)

end

;***********************color*************************************
undef("CCMVal2_Colors1")
function CCMVal2_Colors1(MODEL:string)
local colori
begin
    colori = new((/dimsizes(MODEL)/),"string")

    do i=0, dimsizes(MODEL)-1
        colori(i) = "Black"
        if((MODEL(i).eq."AMTRAC3").or.(MODEL(i).eq."AMTRAC")) then
            colori(i) = "(/0.105, 0.67, 0.8/)"
            ;; colori(i) = "(/27., 170., 204./)"
        end if
        if(MODEL(i).eq."CAM3.5") then
            colori(i) = "(/0.737, 0.459, 0.718/)"
            ;; colori(i) = "(/128., 207., 226./)"
        end if
        if(MODEL(i).eq."CCSRNIES") then
            colori(i) = "(/0.922, 0.722, 0.012/)"
            ;; colori(i) = "(/235., 184., 3./)"
        end if
        if(MODEL(i).eq."CMAM") then
            colori(i) = "(/1., 0., 0./)"
            ;; colori(i) = "(/255., 0., 0./)"
        end if
        if(MODEL(i).eq."CNRM-ACM") then
            colori(i) = "(/0.505, 0.004, 0.494/)"
            ;; colori(i) = "(/129., 1., 126./)"
        end if
        if(MODEL(i).eq."EMAC_2000_02") then
            colori(i) = "Red"
        end if
        if(MODEL(i).eq."EMAC_2000_04") then
            colori(i) = "Green"
        end if
        if(MODEL(i).eq."EMAC_2000_06") then
            colori(i) = "Blue"
        end if
        if(MODEL(i).eq."EMAC") then
            colori(i) = "(/0.94, 0.467, 0./)"
        end if
        if((MODEL(i).eq."E39CA").or.(MODEL(i).eq."E39C")) then
            colori(i) = "(/0.737, 0.459, 0.718/)"
            ;; colori(i) = "(/128., 207., 226./)"
        end if
        if(MODEL(i).eq."GEOSCCM") then
            colori(i) = "(/0.505, 0.004, 0.494/)"
            ;; colori(i) = "(/129., 1., 126./)"
        end if
        if((MODEL(i).eq."LMDZrepro").or.(MODEL(i).eq."LMDZzINCA")) then
            colori(i)="(/1., 0., 0./)"
            ;; colori(i) =  "(/255., 0., 0./)"
        end if
        if(MODEL(i).eq."MAECHAM4CHEM") then
            colori(i) = "Green"
        end if
        if(MODEL(i).eq."MRI") then
            colori(i) =  "(/0., 0.534, 0.216/)"
            ;; colori(i)="(/0., 136., 55./)"
        end if
        if(MODEL(i).eq."Niwa_SOCOL") then
            colori(i)="(/0.6, 0.811, 0.086/)"
            ;; colori(i) = "(/153., 207., 22./)"
        end if
        if(MODEL(i).eq."SOCOL") then
            ;;a colori(i) = "(/5., 74., 94./)"
            ;;a colori(i) = "(/0.0196, 0.29, 0.369/)"
            ;;b colori(i) =  "(/0.712, 0.239, 0./)"
            ;;b colori(i) = "(/184., 61., 0./)"
            ;;c colori(i)="(/153, 207, 22/)"
            colori(i)="(/0.6, 0.811, 0.086/)"
        end if
        if(MODEL(i).eq."ULAQ") then
            ;; colori(i) = "(/240., 119., 0./)"
            colori(i) = "(/0.94, 0.467, 0./)"
        end if
        if(MODEL(i).eq."UMETRAC") then
            ;; colori(i) = "(/27., 170., 204./)"
            colori(i) = "(/0.105, 0.67, 0.8/)"
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            ;; colori(i) = "(/240., 3., 127./)"
            colori(i) =  "(/0.94, 0.012, 0.498/)"
        end if
        if(MODEL(i).eq."UMUKCA-METO") then
            ;; colori(i) = "(/0., 136., 55./)"
            colori(i) =  "(/0., 0.534, 0.216/)"
        end if
        if(MODEL(i).eq."UMUKCA-UCAM") then
            ;; colori(i) = "DeepPink"
            ;; colori(i) = "(/13., 56., 150./)"
            colori(i) = "(/0.05, 0.219, 0.588/)"
        end if
        if(MODEL(i).eq."WACCM") then
            ;; colori(i) =  "(/13., 56., 150./)"
            colori(i) = "(/0.05, 0.219, 0.588/)"
            ;; colori(i) =  "(/0., 0.580, 0.494/)"
        end if
        if(MODEL(i).eq."ERA40") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."ERA-40") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."HALOE") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."UKMO") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."OBS") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."MULTIMODEL_MEAN") then
            ;; colori(i) = "(/0.875, 0.875, 0.875/)"
            colori(i) = "brown"
        end if
        if(MODEL(i).eq."ECEARTH") then
            colori(i) = "red"
        end if
        if(MODEL(i).eq."GPCP") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."ERAINT") then
            colori(i) = "green"
        end if

    end do

    return(colori)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

undef("CCMVal2_Color1")
function CCMVal2_Color1(MODEL:string)
local colori, MODEL
begin
    colori = new((/dimsizes(MODEL)/),"string")

    do i=0, dimsizes(MODEL)-1
        colori(i) = "Black"
        if((MODEL(i)).eq."AMTRAC3") then
            colori(i) = "(/0., 0., 1./)"
        end if
        if((MODEL(i)).eq."ECHAM1") then
            colori(i) = "(/0., 0., 1./)"
        end if
        if(MODEL(i).eq."CAM3.5") then
            colori(i) = "(/1., 0.498, 0./)"
        end if
        if(MODEL(i).eq."CCSRNIES") then
            colori(i) = "(/0.549, 0., 0.831/)"
        end if

        if(MODEL(i).eq."CMAM") then
            colori(i) = "(/1., 0., 0./)"
        end if
        if(MODEL(i).eq."CNRM-ACM") then
            colori(i) = "(/0., 1., 0./)"
        end if
        if(MODEL(i).eq."EMAC_2000_02") then
            colori(i) = "Red"
        end if
        if(MODEL(i).eq."EMAC_2000_04") then
            colori(i) = "Green"
        end if
        if(MODEL(i).eq."EMAC_2000_06") then
            colori(i) = "Blue"
        end if

        if(MODEL(i).eq."EMAC") then
            colori(i) = "(/1., 0.078, 0.58/)"
        end if
        if(MODEL(i).eq."E39CA") then
            colori(i) = "(/0.6,.831, 0./)"
        end if
        if(MODEL(i).eq."GEOSCCM") then
            colori(i) =  "(/0., 0.75, 1./)"
        end if
        if(MODEL(i).eq."LMDZrepro") then
            colori(i) = "(/0.6,.831, 0./)"
        end if
        if(MODEL(i).eq."MRI") then
            colori(i) = "(/0., 0., 1./)"
        end if
        if(MODEL(i).eq."Niwa_SOCOL") then
            colori(i) = "(/1., 1., 0./)"
        end if
        if(MODEL(i).eq."SOCOL") then
            colori(i) = "(/1., 0.498, 0./)"
        end if
        if(MODEL(i).eq."ULAQ") then
            colori(i) =  "(/0.549, 0., 0.831/)"
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            colori(i) =  "(/1., 0., 0./)"
        end if
        if(MODEL(i).eq."UMUKCA-METO") then
            colori(i) = "(/0., 1., 0./)"
        end if
        if(MODEL(i).eq."UMUKCA-UCAM") then
                                ; colori(i) = "DeepPink"
            colori(i)="(/1., 0.078, 0.58/)"
        end if
        if(MODEL(i).eq."WACCM") then
            colori(i) =  "(/0., 0.75, 1./)"
        end if
        if(MODEL(i).eq."ERA40") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."ERA-40") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."HALOE") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."OBS") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."MULTIMODEL_MEAN") then
            colori(i) = "brown"
        end if
        if(MODEL(i).eq."GEOSCCM_B2a_SRESA2") then
            colori(i) = "(/1., 0.498, 0./)"
        end if
        if(MODEL(i).eq."GEOSCCM_B2b") then
            colori(i) =  "(/0.549, 0., 0.831/)"
        end if
        if(MODEL(i).eq."GEOSCCM_B2f") then
            colori(i) =  "(/1., 0., 0./)"
        end if

    end do

    return(colori)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

undef("CCMVal2_LevelDash")
function CCMVal2_LevelDash(num_levels:integer)
local dashiall
begin
    dashi = new((/num_levels/),"integer")

    dashiall=new(20,"integer")

    dashiall(0) = 0
    dashiall(1) = 0
    dashiall(2) = 0
    dashiall(3) = 0
    dashiall(4) = 0
    dashiall(5) = 0
    dashiall(6) = 0
    dashiall(7) = 0
    dashiall(8) = 0
    dashiall(9) = 2
    dashiall(10) = 2
    dashiall(11) = 2
    dashiall(12) = 1
    dashiall(13) = 2
    dashiall(14) = 2
    dashiall(15) = 2
    dashiall(16) = 1
    dashiall(17) = 1
    dashiall(18) = 1
    dashiall(19) = 1

    dashi =dashiall(0:num_levels-1)
    return(dashi)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

undef("CCMVal2_Dash1")
function CCMVal2_Dash1(MODEL:string)
begin
    dashi = new((/dimsizes(MODEL)/),"integer")

    do i=0, dimsizes(MODEL)-1
        dashi(i) = 0
        if((MODEL(i).eq."AMTRAC3").or.(MODEL(i).eq."AMTRAC")) then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CAM3.5") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CCSRNIES") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CMAM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CNRM-ACM") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."EMAC") then
            dashi(i) = 0
        end if
        if((MODEL(i).eq."E39CA").or.(MODEL(i).eq."E39C")) then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."GEOSCCM") then
            dashi(i) = 0
        end if
        if((MODEL(i).eq."LMDZrepro").or.(MODEL(i).eq."LMDZzINCA")) then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."MAECHAM4CHEM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."MRI") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."Niwa_SOCOL") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."SOCOL") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."ULAQ") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."UMETRAC") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."UMUKCA-METO") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."UMUKCA-UCAM") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."WACCM") then
            dashi(i) = 0
        end if

        if(MODEL(i).eq."ERA40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."ERA-40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."NCEP") then
            dashi(i) = 16
        end if
        if(MODEL(i).eq."UKMO") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."HALOE") then
            dashi(i) = 16
        end if
        if(MODEL(i).eq."OBS") then
            dashi(i) = 7
        end if
    end do

    return(dashi)

end
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

undef("CCMVal_Dash")
function CCMVal_Dash(MODEL:string)
;local dashi
begin
    dashi = new((/dimsizes(MODEL)/),"integer")

    do i=0, dimsizes(MODEL)-1
        dashi(i) = 0
        if(MODEL(i).eq."AMTRAC3") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CAM3.5") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CCSRNIES") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CMAM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CNRM-ACM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."EMAC") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."E39CA") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."GEOSCCM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."LMDZrepro") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."MRI") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."SOCOL") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."ULAQ") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMUKCA-METO") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."UMUKCA-UCAM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."WACCM") then
            dashi(i) = 1
        end if

        if(MODEL(i).eq."ERA40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."ERA-40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."NCEP") then
            dashi(i) = 16
        end if
        if(MODEL(i).eq."UKMO") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."HALOE") then
            dashi(i) = 16
        end if
        if(MODEL(i).eq."OBS") then
            dashi(i) = 7
        end if
    end do

    return(dashi)

end
;****************************************************************

undef("CCMVal_Dash1")
function CCMVal_Dash1(MODEL:string)
;local dashi
begin
    dashi = new((/dimsizes(MODEL)/),"integer")

    do i=0, dimsizes(MODEL)-1

        if(MODEL(i).eq."AMTRAC") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."AMTRAC3") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."CCSRNIES") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CMAM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."E39CA") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."GEOSCCM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."LMDZrepro") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."MAECHAM4CHEM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."MRI") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."SOCOL") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."Niwa_SOCOL") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."ULAQ") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMETRAC") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."WACCM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CAM3.5") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."EMAC") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMUKCA-UCAM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CNRM-ACM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."AMTRAC3") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CAM3.5") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CCSRNIES") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CMAM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CNRM-ACM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."EMAC") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."E39CA") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."GEOSCCM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."LMDZrepro") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."MRI") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."SOCOL") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."ULAQ") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMUKCA-METO") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."UMUKCA-UCAM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."WACCM") then
            dashi(i) = 1
        end if

        if(MODEL(i).eq."NCEP") then
            dashi(i) = 16
        end if

        if(MODEL(i).eq."UKMO") then
            dashi(i) = 2
        end if
        if(MODEL(i).eq."ERA40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."ERA-40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."HALOE") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."OBS") then
            dashi(i) = 0
        end if
    end do

    return(dashi)

end
;*************************************************************

undef("CCMVal1_Color")
function CCMVal1_Color(MODEL:string)
local colori
begin
    colori = new((/dimsizes(MODEL)/),"string")

    do i=0, dimsizes(MODEL)-1
        colori(i) = "Black"
        if((MODEL(i)).eq."AMTRAC") then
            colori(i) = "MediumBlue"
        end if
        if((MODEL(i)).eq."ECHAM1") then
            colori(i) = "MediumBlue"
        end if
        if(MODEL(i).eq."CCSRNIES") then
            colori(i) = "DarkOrchid4"
        end if
        if(MODEL(i).eq."EMAC_2000_02") then
            colori(i) = "Red"
        end if
        if(MODEL(i).eq."EMAC_2000_04") then
            colori(i) = "Green"
        end if
        if(MODEL(i).eq."EMAC_2000_06") then
            colori(i) = "Blue"
        end if
        if(MODEL(i).eq."EMAC") then
            colori(i) = "Purple"
        end if

        if(MODEL(i).eq."CMAM") then
            colori(i) = "Red"
        end if
        if(MODEL(i).eq."E39C") then
            colori(i) = "MediumPurple1"
        end if
        if(MODEL(i).eq."GEOSCCM") then
            colori(i) = "Purple"
        end if
        if(MODEL(i).eq."LMDZrepro") then
            colori(i) = "Darkseagreen2"
        end if
                                ;if(MODEL(i).eq."LMDZrepro") then
                                ;   colori(i) = "LightSeaGreen"
                                ;end if
        if(MODEL(i).eq."LMDzINCA") then
            colori(i) = "LightSeaGreen"
        end if
        if(MODEL(i).eq."MAECHAM4CHEM") then
            colori(i) = "Green"
        end if
        if(MODEL(i).eq."MRI") then
            colori(i) = "Magenta"
        end if
        if(MODEL(i).eq."SOCOL") then
            colori(i) = "Yellow"
        end if
        if(MODEL(i).eq."ULAQ") then
            colori(i) = "GreenYellow"
        end if
        if(MODEL(i).eq."UMETRAC") then
            colori(i) = "darkturquoise"
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            colori(i) = "DarkOrange"
        end if
        if(MODEL(i).eq."WACCM") then
            colori(i) = "RoyalBlue1"
        end if
        if(MODEL(i).eq."ERA40") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."HALOE") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."NCEP") then
            colori(i) = "black"
        end if
        if(MODEL(i).eq."UKMO") then
            colori(i) = "black"
        end if
    end do

    return(colori)

end
;***************************************************************

undef("CCMVal1_Dash")
function CCMVal1_Dash(MODEL:string)
;local dashi
begin
    dashi = new((/dimsizes(MODEL)/),"integer")

    do i=0, dimsizes(MODEL)-1
        if(MODEL(i).eq."AMTRAC") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."CCSRNIES") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CMAM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."E39C") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."GEOSCCM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."LMDZrepro") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."MAECHAM4CHEM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."MRI") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."SOCOL") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."ULAQ") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMETRAC") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."WACCM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."ERA-40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."HALOE") then
            dashi(i) = 16
        end if
        if(MODEL(i).eq."NCEP") then
            dashi(i) = 16
        end if

        if(MODEL(i).eq."UKMO") then
            dashi(i) = 2
        end if

    end do

    return(dashi)

end
;%%%%%%%%%%%%%%%%%%%%%%%% CCMVal2 dashes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;FOR THE MOMENT CCMVal_Dash1=CCMValDash

undef("CCMVal_Dash1")
function CCMVal_Dash1(MODEL:string)
;local dashi
begin
    dashi = new((/dimsizes(MODEL)/),"integer")

    do i=0, dimsizes(MODEL)-1
        dashi(i) = 0
        if(MODEL(i).eq."OBS") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."MEAN") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CAM3.5") then
            dashi(i) = 0
        end if      
        if(MODEL(i).eq."AMTRAC3") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."WACCM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."UMUKCA-METO") then
            dashi(i) = 0
        end if      
        if(MODEL(i).eq."SOCOL") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CCSRNIES") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."EMAC") then
            dashi(i) = 0
        end if      
        if(MODEL(i).eq."CMAM") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."GEOSCCM") then
            dashi(i) = 0
        end if

        if(MODEL(i).eq."E39CA") then
            dashi(i) = 1
        end if      
        if(MODEL(i).eq."UMETRAC") then
            dashi(i) = 1
        end if      
        if(MODEL(i).eq."UMUKCA-UCAM") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."MRI") then
            dashi(i) = 1
        end if      
        if(MODEL(i).eq."Niwa_SOCOL") then
            dashi(i) = 1
        end if      
        if(MODEL(i).eq."ULAQ") then
            dashi(i) = 1
        end if
        if(MODEL(i).eq."LMDZrepro") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."CNRM-ACM") then
            dashi(i) = 1
        end if   
        
        if(MODEL(i).eq."ERA40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."ERA-40") then
            dashi(i) = 0
        end if
        if(MODEL(i).eq."NCEP") then
            dashi(i) = 16
        end if
        if(MODEL(i).eq."UKMO") then
            dashi(i) = 4
        end if
        if(MODEL(i).eq."JRA25") then
            dashi(i) = 3
        end if
        if(MODEL(i).eq."HALOE") then
            dashi(i) = 16
        end if
        if(MODEL(i).eq."ERAI") then
            dashi(i) = 1
        end if
    end do

    return(dashi)

end


;%%%%%%%%%%%%%%%%%%%%%%%% CCMVal2 THICKNESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
undef("CCMVal_Thick")
function CCMVal_Thick(MODEL:string)
;local thicki
begin
    thicki = new((/dimsizes(MODEL)/),"integer")

    do i=0, dimsizes(MODEL)-1
        thicki(i) = 2
        if(MODEL(i).eq."OBS") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."MEAN") then
            thicki(i) = 6
        end if
        if(MODEL(i).eq."CAM3.5") then
            thicki(i) = 2
        end if      
        if(MODEL(i).eq."AMTRAC3") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."WACCM") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."UMUKCA-METO") then
            thicki(i) = 2
        end if      
        if(MODEL(i).eq."SOCOL") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."CCSRNIES") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."EMAC") then
            thicki(i) = 2
        end if      
        if(MODEL(i).eq."CMAM") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."UMSLIMCAT") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."GEOSCCM") then
            thicki(i) = 2
        end if

        if(MODEL(i).eq."E39CA") then
            thicki(i) = 2
        end if      
        if(MODEL(i).eq."UMETRAC") then
            thicki(i) = 2
        end if      
        if(MODEL(i).eq."UMUKCA-UCAM") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."MRI") then
            thicki(i) = 2
        end if      
        if(MODEL(i).eq."Niwa_SOCOL") then
            thicki(i) = 2
        end if      
        if(MODEL(i).eq."ULAQ") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."LMDZrepro") then
            thicki(i) = 2
        end if
        if(MODEL(i).eq."CNRM-ACM") then
            thicki(i) = 2
        end if   
        
        if(MODEL(i).eq."ERA40") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."ERA-40") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."NCEP") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."NCEP2") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."UKMO") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."JRA25") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."RICH-ERA40") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."HALOE") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."NIWA-Obs") then
            thicki(i) = 4
        end if
        if(MODEL(i).eq."ERAI") then
            thicki(i) = 4
        end if

    end do

    return(thicki)

end
;********************end of CCMVal2 additions*************************

undef("CCMVal2_ColorTable")
function CCMVal2_ColorTable(MODEL:string)
begin
    return(CCMVal2_Colors1(MODEL))
end

undef("CCMVal1_ColorTable")
function CCMVal1_ColorTable(MODEL:string)
begin
                                ; verion3: use function ccmval2
    return(CCMVal2_Colors1(MODEL))
end

undef("CCMVal2_DashTable")
function CCMVal2_DashTable(MODEL:string)
begin
    return(CCMVal2_Dash1(MODEL))
end

undef("CCMVal1_DashTable")
function CCMVal1_DashTable(MODEL:string)
begin
                                ; verion3: use function ccmval2
    return(CCMVal2_Dash1(MODEL))
end

;*******ADDLINE PROCEDURE***************************************
undef("add_line")
procedure add_line(wks, plot, x, yval, line_color, line_dash)
local plres, str, y, plot, x, yval, line_color, line_dash
begin
    plres = True
    plres@gsLineColor = line_color
    plres@gsLineDashPattern = line_dash
    plres@gsLineThicknessF=2.5
                                ;y=(/yval, yval/)
    str = unique_string("polyline")
    plot@$str$ = gsn_add_polyline(wks, plot, x, yval, plres)
end

;*******ADDSCATT PROCEDURE***************************************
undef("add_scatt")
procedure add_scatt(wks, plot, x, yval, mark_color, mark_ind)
local plres, str, y, plot, x, yval, line_color, line_dash
begin
    plres = True
    plres@gsMarkerColor = mark_color
    plres@gsMarkerIndex = mark_ind
                                ;y=(/yval, yval/)
    str = unique_string("polymarker")
    plot@$str$ = gsn_add_polymarker(wks, plot, x, yval, plres)
end

;********LEGEND PROCEDURE*************************************
undef("add_legenda")
procedure add_legenda(wks, MODEL:string, colors:string, dashs:numeric, place:string)
local MODEL, dim_mod, colors, dashs, edge_x, edge_y, edge_res, tx1res, dim_lb, xleg, xl, yleg, yl, xtxt, ytxt, place, dim_lines, space_lines, lin_sp
begin
                                ;***************Create Legend***********************************
    dim_mod = dimsizes(MODEL)
    dim_lines = floattoint(dim_mod/3)+1
    if ((mod(dim_mod, 3)).eq.0 ) then
        dim_lines = floattoint(dim_mod/3)+1
    else
        dim_lines = floattoint(dim_mod/3)+2
    end if
    lin_sp = 0.035
    space_lines = (dim_lines+1)*lin_sp

    if (place.eq."bottom")
        top_line = 0.28
    end if
    if (place.eq."middle")
        top_line = 0.42
    end if
    if (place.eq."top")
        top_line = 0.99
    end if

    bot_line = top_line-space_lines
    edge_x = (/ 0.12, 0.92, 0.92, 0.12, 0.12/)
    edge_y = (/ space_lines+bot_line, space_lines+bot_line, bot_line, bot_line, space_lines+bot_line/)
    edge_res = True             ; Indicate you want to set some resources.

                                ;gsn_polyline_ndc(wks, edge_x, edge_y, edge_res)

                                ;****************************************************************

    dim_lb = dimsizes(MODEL)

    gs1res               = True
    tx1res               = True
    tx1res@txFontHeightF = lin_sp/5*2

    xleg = new((/3*dim_lines, 4/),"float")
    xl = new((/3, 4/),"float")

    do i=0, 2
        xl(0,:) = fspan((edge_x(0)+0.02),(edge_x(0)+0.07), 4)
        xl(1,:) = fspan(((edge_x(2)+edge_x(0))/2-0.09),((edge_x(2)+edge_x(0))/2-0.04), 4)
        xl(2,:) = fspan((edge_x(2)-0.20),(edge_x(2)-0.15), 4)
    end do

    nt = 0
    do  j=0, dim_lines-1
        do i=0, 2
            xleg(nt,:) = (/xl(i,:)/)
            nt = nt+1
        end do
    end do

    yleg = new((/3*dim_lines, 4/),"float")
    yl = new(dim_lines,"float")

    do i=0, dim_lines-1
        yl(dim_lines-1-i) = edge_y(3)+(((edge_y(0) - (edge_y(3)))/dim_lines)*(i+1))-(lin_sp)
    end do

    nt = 0
    do j=0, dim_lines-1
        do i=0, 2
            yleg(nt,:) = (/yl(j)/)
            nt = nt+1
        end do
    end do

    xtxt = new((/dim_lines*3/),"float")
    ytxt = new((/dim_lines*3/),"float")

    nt = 0
    do j=0, dim_lines-1
        do i=0, 2
            xtxt(nt) = xl(i, 3)+lin_sp/5.
            nt = nt+1
        end do
    end do
    nt = 0
    do i=0, dim_lines-1
        do j=0, 2
            ytxt(nt) = (/yl(i)/)
            nt = nt+1
        end do
    end do

    do i=0,(dimsizes(MODEL)-1)
        gs1res@gsLineColor = colors(i)
        gs1res@gsLineDashPattern = dashs(i)
        gs1res@gsLineThicknessF=2.5
        tx1res@txFontColor = colors(i)
        tx1res@txJust      = "CenterLeft"             ; change justification
        if (((MODEL(i).eq."HALOE").and.(dashs(i).eq.16)).or.((MODEL(i).eq."NCEP").and.(dashs(i).eq.16)).or.((MODEL(i).eq."UKMO").and.(dashs(i).eq.2)))then
            gs1res@gsMarkerColor = colors(i)
            gs1res@gsMarkerIndex = dashs(i)
            gsn_polymarker_ndc(wks, xleg(i, 2), yleg(i, 2), gs1res)
        else

            gsn_polyline_ndc(wks, xleg(i,:), yleg(i,:), gs1res)
        end if

        gsn_text_ndc(wks, MODEL(i), xtxt(i), ytxt(i), tx1res)
    end do

end

;********LEGEND PROCEDURE*************************************
undef("add_legenda_mark")
procedure add_legenda_mark(wks, MODEL:string, colors:string, indexe:numeric)
local MODEL, colors, dashs, edge_x, edge_y, edge_res, tx1res, dim_lb, xleg, x1, yleg, y1, xtxt, ytxt
begin
                                ;***************Create Legend***********************************
    edge_x = (/ 0.1, 0.9, 0.9, 0.1, 0.1/)
    edge_y = (/ .18, .18, 0.001, 0.001, 0.18/)
    edge_res = False             ; Indicate you want to set some resources.
    gsn_polyline_ndc(wks, edge_x, edge_y, edge_res)

                                ;****************************************************************
    dim_lb = dimsizes(MODEL)

    gs1res               = True
    tx1res               = True
    tx1res@txFontHeightF = 0.013

    xleg = new((/18, 6/),"float")
    xl = new((/3, 6/),"float")
    do i=0, 2
        xl(0,:) = fspan((edge_x(0)+0.02),(edge_x(0)+0.09), 6)
        xl(1,:) = fspan(((edge_x(2)-edge_x(0))/2+0.02),((edge_x(2)-edge_x(0))/2+0.09), 6)
        xl(2,:) = fspan((edge_x(2)-0.20),(edge_x(2)-0.15), 6)
    end do
    a=(/xl(0,:), xl(0,:), xl(0,:), xl(0,:), xl(0,:), xl(0,:), xl(1,:), xl(1,:), xl(1,:), xl(1,:), xl(1,:), xl(1,:), xl(2,:), xl(2,:), xl(2,:), xl(2,:), xl(2,:), xl(2,:)/)
    xleg(:,:) = (/xl(0,:), xl(0,:), xl(0,:), xl(0,:), xl(0,:), xl(0,:), xl(1,:), xl(1,:), xl(1,:), xl(1,:), xl(1,:), xl(1,:), xl(2,:), xl(2,:), xl(2,:), xl(2,:), xl(2,:), xl(2,:)/)

    yleg = new((/18, 6/),"float")
    yl = new(6,"float")
    do i=0, 5
        yl(5-i) = (((edge_y(0) - (edge_y(3)))/6)*i)+0.017
    end do

    yleg(:, 0)=(/yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5)/)
    yleg(:, 1)=(/yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5)/)
    yleg(:, 2)=(/yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5)/)
    yleg(:, 3)=(/yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5)/)
    yleg(:, 4)=(/yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5)/)

    xtxt =(/(xl(0, 5)+0.015),(xl(0, 5)+0.015),(xl(0, 5)+0.015),(xl(0, 5)+0.015),(xl(0, 5)+0.015),(xl(0, 5)+0.015),(xl(1, 5)+0.015),(xl(1, 5)+0.015),(xl(1, 5)+0.015),(xl(1, 5)+0.015),(xl(1, 5)+0.015),(xl(1, 5)+0.015),(xl(2, 5)+0.015),(xl(2, 5)+0.015),(xl(2, 5)+0.015),(xl(2, 5)+0.015),(xl(2, 5)+0.015),(xl(2, 5)+0.015)/)
    ytxt = (/yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5), yl(0), yl(1), yl(2), yl(3), yl(4), yl(5)/) ; strings.

    do i=0,(dimsizes(MODEL)-1)
        gs1res@gsMarkerColor = colors(i)
        gs1res@gsMarkerIndex = indexe(i)
        tx1res@txFontColor   = colors(i)
        tx1res@txJust        = "CenterLeft"             ; change justification
        gsn_polymarker_ndc(wks, xleg(i, 2), yleg(i, 2), gs1res)
        gsn_text_ndc(wks, MODEL(i), xtxt(i), ytxt(i), tx1res)
    end do

end
;*********************************************************************

;********LEGEND PROCEDURE with text*************************************
;**********************************************************************

undef("add_legenda_valore")
procedure add_legenda_valore(wks, MODEL:string, colors:string, dashs:numeric, valore:numeric)
local MODEL, colors, dashs, edge_x, edge_y, edge_res, tx1res, dim_lb, xleg, x1, yleg, y1, xtxt, ytxt, valore
begin
                                ;***************Create Legend***********************************
    edge_x = (/ 0.1, 0.9, 0.9, 0.1, 0.1/)
    edge_y = (/ .18, .18, 0.001, 0.001, 0.18/)
    edge_res = False             ; Indicate you want to set some resources.
    gsn_polyline_ndc(wks, edge_x, edge_y, edge_res)

                                ;****************************************************************
    dim_lb = dimsizes(MODEL)

    gs1res               = True
    tx1res               = True
    tx1res@txFontHeightF = 0.013

    xleg = new((/15, 5/),"float")
    xl = new((/3, 5/),"float")
    do i=0, 2
        xl(0,:) = fspan((edge_x(0)+0.02),(edge_x(0)+0.07), 5)
        xl(1,:) = fspan(((edge_x(2)-edge_x(0))/2+0.02),((edge_x(2)-edge_x(0))/2+0.07), 5)
        xl(2,:) = fspan((edge_x(2)-0.20),(edge_x(2)-0.15), 5)
    end do

    xleg(:,:) = (/xl(0,:), xl(0,:), xl(0,:), xl(0,:), xl(0,:), xl(1,:), xl(1,:), xl(1,:), xl(1,:), xl(1,:), xl(2,:), xl(2,:), xl(2,:), xl(2,:), xl(2,:)/)

    yleg = new((/15, 5/),"float")
    yl = new(5,"float")
    do i=0, 4
        yl(4-i) = (((edge_y(0) - (edge_y(3)))/5)*i)+0.017
    end do

    yleg(:, 0) = (/yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4)/)
    yleg(:, 1) = (/yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4)/)
    yleg(:, 2) = (/yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4)/)
    yleg(:, 3) = (/yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4)/)
    yleg(:, 4) = (/yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4)/)

    xtxt = (/(xl(0, 4)+0.015),(xl(0, 4)+0.015),(xl(0, 4)+0.015),(xl(0, 4)+0.015),(xl(0, 4)+0.015),(xl(1, 4)+0.015),(xl(1, 4)+0.015),(xl(1, 4)+0.015),(xl(1, 4)+0.015),(xl(1, 4)+0.015),(xl(2, 4)+0.015),(xl(2, 4)+0.015),(xl(2, 4)+0.015),(xl(2, 4)+0.015),(xl(2, 4)+0.015)/)
    ytxt = (/yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4), yl(0), yl(1), yl(2), yl(3), yl(4)/) ; strings.

    stringa = new(dimsizes(MODEL),"string")
    do i=0,(dimsizes(MODEL)-1)
        stringa(i) = MODEL(i)+" : "+valore(i)
    end do
                                ;print(stringa)

    do i=0,(dimsizes(MODEL)-1)
        gs1res@gsLineColor = colors(i)
        gs1res@gsLineDashPattern = dashs(i)
        tx1res@txFontColor = colors(i)
        tx1res@txJust      = "CenterLeft"             ; change justification
        gsn_polyline_ndc(wks, xleg(i,:), yleg(i,:), gs1res)
        gsn_text_ndc(wks, stringa(i), xtxt(i), ytxt(i), tx1res)
    end do

end
;##################################################################################
