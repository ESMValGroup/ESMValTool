;
;  Requires: ua:*2*z
;
info = True
info@derived = True
info@long_name="eastward_wind_max_latitude_S"
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
info@ts_anom = True
info@ts_smooth = True

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
local tmp, dum, data, i, rk_dum, lt, lttmp, latout, t, tdat, jumx, j0, \
      tmpi, x1, x2, p1, p2, pt
begin
  data_new = True
  tmp = read_data(index,"ua","*2*z")

    ; extract 200 hPa
    data = extract_data(index,tmp,1,200.,200.)   ;extract_data(models,tmp,i,"plev",200.,200.)
    rk_dum = dimsizes(dimsizes(data))

    ; make zonal mean if required
    if ((rk_dum .eq. 3) .and. (isdim(data,"lon"))) then
      dum = dim_avg_Wrap(data)
      delete(data)
      data = dum
      delete(dum)
    end if

    lt=data&lat
    lt&lat=data&lat

    lttmp=lt({0:-80})
;    printVarSummary(data)

    latout=new((/dimsizes(data&time)/),"float")
    latout!0="time"
    latout&time=data&time
    latout@long_name=info@long_name
    latout@units=info@units

;now find maximum and save latitude
    do t=0,dimsizes(data&time)-1

      tdat=data(t,{0.:-80})
      jumx=max(tdat)

;      jj=ind(tdat .eq. jumx)
;;      if (jumx gt 0)
;      if(typeof(lttmp).eq."double") then
;        latout(t)=doubletofloat(lttmp(jj))
;      else
;        latout(t)=(/lttmp(jj)/)
;      end if

; modified by ZL (2011-04-15)

      j0=ind(tdat .eq. jumx)

      if (j0.eq.0.or.j0.eq.dimsizes(lttmp)-1) then
        tmpi= lttmp(j0)
        continue
      end if

      x1= (lttmp(j0-1)+lttmp(j0))*0.5
      x2= (lttmp(j0)+lttmp(j0+1))*0.5

      p1= (tdat(j0)-tdat(j0-1))/(lttmp(j0)-lttmp(j0-1))
      p2= (tdat(j0+1)-tdat(j0))/(lttmp(j0+1)-lttmp(j0))

      pt= 0.

      if(p2.ne.p1) then
        tmpi= ((pt-p1)*x2+(p2-pt)*x1)/(p2-p1)
      else
        tmpi= (p1+p2)*0.5
      end if

      if(typeof(tmpi).eq."double") then
        latout(t)= doubletofloat(tmpi)
      else
        latout(t)= tmpi
      end if

      delete(tmpi)
      delete(tdat)
      delete(jumx)
      delete(j0)
      delete(x1)
      delete(x2)
      delete(p1)
      delete(p2)
      delete(pt)

    end do

    add_data(index,data_new,latout)

    delete(rk_dum)
    delete(data)
    delete(lttmp)
    delete(lt)
    delete(latout)

  delete(tmp) 

  return(data_new)
end
