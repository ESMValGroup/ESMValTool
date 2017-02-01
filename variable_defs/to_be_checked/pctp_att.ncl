;
; Requires: ta:*3*
;
info = True
info@derived = True
info@long_name="Cold Point Tropopause Air Pressure"
info@units="hPa"
info@scp_cnlvl=(/50,150,200,250,300,350,400,450,500,550,600/) ; ptp
info@scp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@pcp_cnlvl=(/50,150,200,250,300,350,400,450,500,550,600/) ; ptp
info@pcp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@sdp_cnlvl=(/50,150,200,250,300,350,400,450,500,550,600/) ; ptp
info@acp_dcnlvl=(/-18,-15,-12,-9,-6,-3,-1,0,1,3,6,9,12,15,18/)
info@refModel = "CMAM"

info@vcp_ptop=6e-6
info@vcp_LogScale = True
info@rgb_file = "amwg.rgb"

info@ts_latrange = (/-25,25/)
info@ts_lonrange = (/0,360/)
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
local tc, pr, T, missing, pctp, xx, ng, z, itime, ilat, tc_flt
begin
  data_new=True
  result = read_data(index,"ta","*3*")
  missing=-999.0
  T = extract_data(index,result,-1,0,0)
  dims=dimsizes(T)
  pctp=new((/dims(0),dims(2),dims(3)/),"float")

  z = -7500*log(T&plev/1013)
  do itime = 0,dimsizes(T&time)-1
    do ilat  = 0,dimsizes(T&lat)-1
      do ilon  = 0,dimsizes(T&lon)-1

       tc=T(itime,:,ilat,ilon)

        if (typeof(tc).eq."double") then
          tc_flt = doubletofloat(tc)
          delete(tc)
          tc=tc_flt
          delete(tc_flt)
        end if

       pr=tc&plev
       ptp1=pr

;-----------------------------------------------------------------------------
;find cold point trop
;-----------------------------------------------------------------------------

        xx=ind(ptp1.lt.500..and.ptp1.gt.60.)
        ng=dimsizes(xx)
        if (ng .gt. 0) then 
            cpt=min(tc(xx))
            pp=ind(tc(xx) .eq. cpt)
             printVarSummary(pr)
             printVarSummary(z)
             printVarSummary(pr(xx(pp(0))))
             print(pr(xx(pp(0))))
             printVarSummary(pctp)
            if (typeof(ptp1).eq."double") then
              tmp=doubletofloat(ptp1(xx(pp(0))))
            else
              tmp=ptp1(xx(pp(0)))
            end if
            pctp(itime,ilat,ilon)=tmp
             print(pctp(itime,ilat,ilon))
            delete(pp)
        else 
            pctp(itime,ilat,ilon)=missing
        end if

     end do  ;lon
   end do    ;lat
 end do      ;time


  pctp!0=T!0
  pctp!1=T!2
  pctp!2=T!3
  pctp&time=T&time
  pctp&lat=T&lat
  pctp&lon=T&lon
  pctp@long_name= info@long_name
  pctp@units= info@units
  add_data(index,data_new,pctp)

  delete(T)
  delete(pctp)
  delete(ptp1)
  delete(tc)
  delete(z)
  delete(pr)
  delete(xx)
  delete(tmp)

  c=stringtochar(field_number)
  c(1)=integertochar(50) ; "2"
  c(3)=integertochar(115) ; "s"
  field_number=chartostring(c)
  print(field_number)
  delete(c)
  return(data_new)
end

