;
;  Requires: nufl:*2*s
;
info = True
info@derived = True
info@long_name="Tropical OLR Width"
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
info@refModel = "WACCM"

info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

;info@ts_latrange = (/-25,25/)
;info@ts_lonrange = (/0,360/)
;info@ts_ymax = 220
;info@ts_ymin = 180
;info@ts_xmax_year = 2010
;info@ts_xmin_year = 1950
info@ts_anom = False

;following hu and fu: "most poleward latitude where OLR > 250 W/m2
;could also separate NH and SH.
; Could also interpolate between these latitudes and one more.

tpcolr = 250.  ;threshold value

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
local tmp, dum, data, i, j, ds, rank, cn, cnt, attNames, dimNames, rk_dum
begin
  data_new = True
  tmp = read_data(index,"nufl","*2*s")

  dum = extract_data(index,tmp,-1,0,0)
  rk_dum = dimsizes(dimsizes(dum))

;make zonal mean 
  data = dum
  tt = data
  delete(data)
  data = dim_avg_Wrap(tt)      
  delete(tt)

  printVarSummary(data)

  latout=new((/dimsizes(data&time)/),"float")
  latout!0="time"
  latout&time=data&time
  latout@long_name=info@long_name
  latout@units=info@units

;  printVarSummary(latout)

  lati=data&lat
  lati&lat = data&lat

;  print(lati)
  lattmp=lati({-90.:90})

;now find maximum and save latitude
  do t=0,dimsizes(data&time)-1

    tdat=data(t,{-90.:90})
    jj=ind(tdat .gt. tpcolr)
    
;    print(lattmp(jj))
;    print(tdat)
;    exit

    if(typeof(lati).eq."double") then
      latout(t)=doubletofloat(lattmp(max(jj))-lattmp(min(jj)))
    else
      latout(t)=(/lattmp(max(jj))-lattmp(min(jj))/)
    end if


    delete(jj)
    delete(tdat)
  end do

;  print(latout)
  
  add_data(index,data_new,latout)

  delete(dum)
  delete(data)
  delete(lati)
  delete(lattmp)
  delete(latout)

  delete(tmp) 

  return(data_new)
end
