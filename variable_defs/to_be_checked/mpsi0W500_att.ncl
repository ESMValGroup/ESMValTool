;
; Requires: va:*3*,ps:*2*s
;

info = True
info@derived = True

info@long_name="Zonal Mean Meridional Streamfunction Zero Width"
info@units="deg"
info@vcp_cnlvl=(/-20,-15,-10,-5,0,5,10,15,20,25,30,40,50,60,70/)*2.e9 ; U(ua)
info@vcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*2.e9
info@vcp_ptop=6e-4
info@vcp_LogScale = True

info@scp_cnlvl=(/-20,-15,-10,-5,0,5,10,15,20,25,30,40,50,60,70/)*2.e9 ; U(ua)
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*2.e9

info@pcp_cnlvl=(/-20,-15,-10,-5,0,5,10,15,20,25,30,40,50,60,70/)*2.e9 ; U(ua)
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*2.e9

info@acp_cnlvl=(/-20,-15,-10,-5,0,5,10,15,20,25,30,40,50,60,70/)*2.e9 ; U(ua)
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*2.e9

info@sdp_cnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)*2.e9
info@rgb_file = "red-blue.rgb"

;info@ts_latrange = (/-25,25/)
;info@ts_levrange = (/1000,200/)
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
local tmpV, tmpPS, i, dumV, dumPS, vt2b, pPa, p2Pa, psPa, lat, zmpsi
begin
  data_new = True
  tmpV = read_data(index,"va","*3*")
  tmpPS = read_data(index,"ps","*2*s")

  dumV = extract_data(index,tmpV,-1,0,0)
  dumPS = extract_data(index,tmpPS,-1,0,0)

;  printVarSummary(dumV)
;  printVarSummary(dumPS)

  pPa   = dumV&plev*100.     ; make Pa
  pPa!0 = "plev"
  pPa&plev = dumV&plev
  if (pPa(0).lt.pPa(1)) then
    vt2b = dumV
    p2Pa = pPa
  else
    vt2b = dumV(:,::-1,:,:)    ; reverse the "lev" order to top-to-bottom
    p2Pa = pPa(::-1)
  end if

;  vt2b  = dumV(:,{10:1000},:,:)    ; clip from 10-1000 hpa
;  pPa   = dumV&plev({10:1000})*100.     ; reverse and make Pa

  psPa  = dumPS*100.          ; make Pa
  lat = dumV&lat

;  print(p2Pa)

;rerse vertical levels (so top - bottom)
  zmpsi = zonal_mpsi_Wrap(vt2b(:,{50:1000},:,:),lat,p2Pa({50:1000}),psPa)   ; zmpsi(ntim,klev,nlat)

;  printVarSummary(zmpsi)



;n strip out 500 hPa level, N of 20deg.
  zmN = zmpsi(:,{500.},{20:90})
  lt=zmN&lat
  lt&lat=zmN&lat

  zmS = zmpsi(:,{500.},{-90:-20})
  ltS=zmS&lat
  ltS&lat=zmS&lat

;  printVarSummary(zmN)
;  print(lt)

  latout=new((/dimsizes(zmN&time)/),"float")
  latout!0="time"
  latout&time=zmN&time
  latout@long_name=info@long_name
  latout@units=info@units
  latout@_FillValue=zmpsi@_FillValue

;n find first positive and save latitude
;ne, could linearly intepolate to zero.

;  int=(/0.,.1/)
  do t=0,dimsizes(zmN&time)-1

;Nwant lat > 20 first lat with zmpsi < 0.
    jj=ind(zmN(t,:) .lt. 0)

    if (jj(0) .gt. 0.) then
;lear interpolation to zero point...
;    print(jj)
;    print(zmN(t,jj(0)-1:jj(0)))
;    print(lt(jj(0)-1:jj(0)))

    p1=zmN(t,jj(0)-1)
    p2=zmN(t,jj(0))
    x1=lt(jj(0)-1)
    x2=lt(jj(0))
    pt=0.
    tmpi=((pt-p1)*x2 + (p2-pt)*x1)/(p2-p1)

;    print(tmpi)


    if(typeof(tmpi).eq."double") then
;      nlat=doubletofloat(lt(jj(0)))
      nlat=doubletofloat(tmpi)
     else
      nlat=tmpi
;      nlat=(/lt(jj(0))/)
    end if



    delete(tmpi)

     else
       nlat=zmpsi@_FillValue
     end if


    delete(jj)

;Swant lat < -20 first lat with zmpsi > 0.
    jj=ind(zmS(t,:) .gt. 0)


;    if (jj(0) .gt.0) then 
;lear interpolation to zero point...
;    print(jj)
;    print(zmS(t,max(jj)-1:max(jj)))
;    print(ltS(max(jj)-1:max(jj)))

      p1=zmS(t,max(jj)-1)
      p2=zmS(t,max(jj))
      x1=ltS(max(jj)-1)
      x2=ltS(max(jj))
      pt=0.
      tmpi=((pt-p1)*x2 + (p2-pt)*x1)/(p2-p1)

;    print(tmpi)


;    if(typeof(ltS).eq."double") then
    if(typeof(tmpi).eq."double") then
;      slat=doubletofloat(ltS(max(jj)))
      slat=doubletofloat(tmpi)
    else
;      slat=(/ltS(max(jj))/)
      slat=tmpi
     end if

;     else 
;       slat=zmpsi@_FillValue
;     end if 

    delete(jj)
    delete(tmpi)

    delete(p1)
    delete(p2)
    delete(x1)
    delete(x2)

;    print(slat)
;    print(nlat)
    if (nlat .ne.zmpsi@_FillValue .and. slat .ne.zmpsi@_FillValue) then
      latout(t)=nlat-slat
    else
      latout(t)=zmpsi@_FillValue
    end if

    delete(nlat)
    delete(slat)
  end do

;  print(latout)
  
;
  add_data(index,data_new,latout)

  delete(zmpsi)
  delete(zmN)
  delete(zmS)
  delete(lt)
  delete(ltS)
  delete(lat)
  delete(psPa)
  delete(pPa)
  delete(p2Pa)
  delete(dumV)
  delete(dumPS)
  delete(vt2b)
  delete(latout)
  delete(tmpV)
  delete(tmpPS)

  return(data_new)
end

