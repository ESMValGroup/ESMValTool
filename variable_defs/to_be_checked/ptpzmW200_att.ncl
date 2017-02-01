;
; Requires: ta:*2*
;
info = True
info@derived = True   
;OLD Requires: ta:T2*
;info@derived = False 
info@long_name="Width of Tropical Tropopause"
info@units="deg"

info@scp_cnlvl=(/-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80/) 
info@pcp_cnlvl=(/-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80/)
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@sdp_cnlvl = (/-30,-25,-20,-15,-10,-5,-2,0,2,5,10,15,20,25,30/)
info@acp_cnlvl=(/-80,-70,-60,-50,-40,-30,-20,-10,10,20,30,40,50,60,70,80/)
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@seasons = (/"ANN","DJF","MAM","JJA","SON"/)
info@regions = (/"GLBL"/)
;info@refModel = "WACCM"

info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

info@ts_latrange = (/60,90/)
info@ts_lonrange = (/0,360/)
info@ts_anom = True
info@ts_smooth = False

info@ml_latrange = (/-20,-20/)
info@ml_lonrange = (/0,360/)
info@ml_anom = False

;FIND ZONAL MEAN PTP.... (from ptpzm)

undef("calculate")
function calculate(index [1] : integer,
                   variable [1] : string,
                   field_number [1] : string)
;;                 return_val [1] : logical
;; Arguments:
;;    index    - index to current infile defined in the 'temp_ncl_code/filesinfo'-file
;;    variable - logical with releveant variable as string attribute
;;    field_number  - string with field number classification
;; Return value:
;;    data_new - logical
local result, i, pr, T, pin, tin, ptp, lt, latout, ltN, tdatN, \
      lttmpN, knots, xlat, ltS, tdatS, lttmpS, t
begin
  data_new=True
  result = read_data(index,"ta","*2Mz")
;  missing=-999.0
    T = extract_data(index,result,-1,0,0)

;    printVarSummary(T)
;    print(max(T))
;    print(min(T))

;exit

;    dims=dimsizes(T)
;    ptp=new((/dims(0),dims(2),dims(3)/),"float")

;check ordering....needs to be top to bottom (lowest pressure first)

    pr=T&plev
    if (pr(0).eq.max(pr)) then
	pin =  pr(::-1)
	tin =  T({time | :},{lat | :},{plev | ::-1})
    else
	pin = pr
	tin = T({time | :},{lat | :},{plev | :})
    end if	 

    if(max(pin).gt.2000.) then
	pin=pin/100.
    end if

;print(pin)
;printVarSummary(tin)
;print(max(tin))
;print(min(tin))
;exit

     ptp = trop_wmo(pin, tin, 0, False)

; original that did not check for ordering
;    ptp = trop_wmo(pr(::-1), T({time | :},{lat | :},{lon | :},{plev | ::-1}), 0, False) 
 
     delete(pin)
     delete(tin)

;   for outputing altitude (km) instead
;    ztp = -7.500*log(ptp/1013.)

    ptp!0=T!0
    ptp!1=T!2
;    ptp!2=T!3
    ptp&time=T&time
    ptp&lat=T&lat
;    ptp&lon=T&lon
    
;    ptp@_FillValue = T@_FillValue
    ptp@_FillValue = -999.

    printVarSummary(ptp)
    print(max(ptp))
    print(min(ptp))

;Set up latitude indicies

    lt=ptp&lat
    lt&lat=ptp&lat

    lttmpN=lt({20:50.})
    lttmpS=lt({-20:-50})

    latout=new((/dimsizes(ptp&time)/),"float")
    latout!0="time"
    latout&time=ptp&time
    latout@long_name=info@long_name
    latout@units=info@units

    xlat=200.

    knots=4

;now find latitude where ptp = xlat in each hemisphere
; each time, and save
    do t=0,dimsizes(ptp&time)-1
   
      tdatN=ptp(t,{20.:50})

;      print(tdatN)
;      print(lttmpN)
	
;could be better (requires monotonicity)
;      ltN=linint1(tdatN,lttmpN,False,xlat,0)	
       ltN=csa1(tdatN,lttmpN,knots,xlat)

;      print(ltN)

;      jumx=max(tdatN)
;      jj=ind(tdatN .eq. jumx)


      tdatS=ptp(t,{-20.:-50})
;      jumxS=max(tdatS)
;      jjS=ind(tdatS .eq. jumxS)

;      ltS=linint1(tdatS,lttmpS,False,xlat,0)
       ltS=csa1(tdatS,lttmpS,knots,xlat)

;print(ltS)

;      print(jumx)
;      print(jj)
;      print(tdat)
;      print(lttmp(jj))

      if ((ltS.lt.0.and.ltN.gt.0).and.(ltS.gt.-50.and.ltN.lt.50.)) then
         if(typeof(lttmpN).eq."double") then
           latout(t)=doubletofloat(ltN-ltS)
         else
           latout(t)=(/ltN-ltS/)
         end if
      end if	

      delete(ltS)
      delete(ltN)

      delete(tdatN)
;      delete(jumx)
;      delete(jj)
      delete(tdatS)
;      delete(jumxS)
;      delete(jjS)

    end do

;    print(latout)

    add_data(index,data_new,latout)
    delete(T)
    delete(ptp)
    delete(pr)

    delete(lttmpN)
    delete(lttmpS)
    delete(lt)
    delete(latout)


  return(data_new)
end

