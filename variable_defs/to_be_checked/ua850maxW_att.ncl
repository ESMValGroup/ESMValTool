;
;  Requires: ua:*2*z
;
info = True
info@derived = True
info@long_name="eastward_wind_max_width"
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
local tmp, dum, data, i, rk_dum, lt, lttmpN, lttmpS, latout, t, \
      tdatN, jumx, j0, tdatS, jumxS, j1, x1, x2, p1, p2, pt, tmp1, \
      lttmpNN, lttmpSS
begin
  data_new = True
  tmp = read_data(index,"ua","*2*z")

    ; extract 850 hPa
    data = extract_data(index,tmp,1,850.,850.)   ;extract_data(models,tmp,i,"plev",850.,850.)
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

    lttmpN=lt({0:90.})
    lttmpS=lt({0:-80})
;    printVarSummary(data)

    latout=new((/dimsizes(data&time)/),"float")
    latout!0="time"
    latout&time=data&time
    latout@long_name=info@long_name
    latout@units=info@units

; modified by ZL (2011-04-15)
;now find maximum and save latitude
    do t=0,dimsizes(data&time)-1

      do nn=0,0
        tdatN=data(t,{0.:90})
        jumx=max(tdatN)
        j0=ind(tdatN .eq. jumx)

        if (dimsizes(j0).ne.1) then
          tmp1=j0(0)
          delete(j0)
          j0= tmp1
        end if

        if (ismissing(j0)) then
          lttmpNN= -999
          break
        end if

        if (j0.eq.0.or.j0.eq.dimsizes(lttmpN)-1) then
          lttmpNN= lttmpN(j0)
          break
        end if

        x1= (lttmpN(j0-1)+lttmpN(j0))*0.5
        x2= (lttmpN(j0)+lttmpN(j0+1))*0.5

        p1= (tdatN(j0)-tdatN(j0-1))/(lttmpN(j0)-lttmpN(j0-1))
        p2= (tdatN(j0+1)-tdatN(j0))/(lttmpN(j0+1)-lttmpN(j0))

        pt= 0.

        if (p2.ne.p1) then
          lttmpNN= ((pt-p1)*x2+(p2-pt)*x1)/(p2-p1)
        else
          lttmpNN= (p1+p2)*0.5
        end if
      end do

      do ns=0,0
        tdatS=data(t,{0.:-80})
        jumxS=max(tdatS)
        j1=ind(tdatS .eq. jumxS)

        if (dimsizes(j1).ne.1) then
          tmp1=j1(0)
          delete(j1)
          j1= tmp1
        end if

        if (ismissing(j1)) then
          lttmpSS= -999
          break
        end if

        if (j1.eq.0.or.j1.eq.dimsizes(lttmpS)-1) then
          lttmpSS= lttmpS(j1)
          break
        end if

        x1= (lttmpS(j1-1)+lttmpS(j1))*0.5
        x2= (lttmpS(j1)+lttmpS(j1+1))*0.5

        p1= (tdatS(j1)-tdatS(j1-1))/(lttmpS(j1)-lttmpS(j1-1))
        p2= (tdatS(j1+1)-tdatS(j1))/(lttmpS(j1+1)-lttmpS(j1))

        pt= 0.

        if (p2.ne.p1) then
          lttmpSS= ((pt-p1)*x2+(p2-pt)*x1)/(p2-p1)
        else
          lttmpSS= (p1+p2)*0.5
        end if
      end do

      if (typeof(lttmpNN).eq."double") then
        latout(t)=doubletofloat(lttmpNN-lttmpSS)
      else
        latout(t)=(/lttmpNN-lttmpSS/)
      end if

;      tdatN=data(t,{0.:90})
;      jumx=max(tdatN)
;      jj=ind(tdatN .eq. jumx)
;      tdatS=data(t,{0.:-80})
;      jumxS=max(tdatS)
;      jjS=ind(tdatS .eq. jumxS)
;      if(typeof(lttmpN).eq."double") then
;        latout(t)=doubletofloat(lttmpN(jj)-lttmpS(jjS))
;      else
;        latout(t)=(/lttmpN(jj)-lttmpS(jjS)/)
;      end if

      delete(tdatN)
      delete(jumx)
      delete(j0)
      delete(tdatS)
      delete(jumxS)
      delete(j1)
      delete(x1)
      delete(x2)
      delete(p1)
      delete(p2)
      delete(pt)
      delete(tmp1)
      delete(lttmpNN)
      delete(lttmpSS)

    end do

    add_data(index,data_new,latout)

    delete(rk_dum)
    delete(data)
    delete(lttmpN)
    delete(lttmpS)
    delete(lt)
    delete(latout)

  delete(tmp)

  return(data_new)
end
