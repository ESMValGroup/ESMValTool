;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Name Surname (Affiliation)                                  #
;#######################################################################
;# PROJECT:     Project name (CCMVal1, CCMVal2, CMIP5...)              #
;# ENSEMBLE:    Affected ensemble(s) or [all]                          #
;# MODEL:       Model name                                             #
;# VARIABLES:   Name(s) of affected variable(s) or [all]               #
;# FIELDS:      Affected field(s) or [all]
;# DESCRIPTION: Brief description of the error fixed by this routine   #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)

local ret, tmp, ti, tt
begin
    verbosity  = config_user_info@verbosity

    ret = 1

    ;if (isdim(var,"time").and.isatt(var&time,"calendar")) then
    ;  var&time@calendar = "365_day"
;      ret = 0
;    end if
    
    if (name .eq. "co2") then
        ;fils0 = systemfunc("ls /gpfs02/ESMVal/model/ETHZ_CMIP5/historical/Amon/clw/HadGEM2-ES/r1i1p1/clw_Amon_HadGEM2-ES_historical_r1i1p1_185912-188411.nc")
        ;tmp0 = addfile(fils0, "r")
        ;b    = tmp0->b
        ;orog = tmp0->orog

        ;z = new((/dimsizes(b),dimsizes(orog(:,0)),dimsizes(orog(0,:))/),typeof(b))
        ;do t = 0,dimsizes(var&time)-1
        ;    do k = 0, dimsizes(b)-1
        ;        do i = 0,dimsizes(orog(:,0))-1
        ;            do j = 0,dimsizes(orog(0,:))-1
        ;                z(k,i,j) = var&plev(k)+b(k)*orog(i,j)
        ;            end do
        ;        end do
        ;    end do
        ;end do
        ;    z!1 = "lev"
        ;    z&lev = z(0,:,0,0)
        ;var&plev = z(:,0,0)
        ;delete([/b,orog,i,k,j/])
   
        ;fils1 = systemfunc("ls /data/ESMValTool/model/ETHZ_CMIP5/historical/Amon/zg/HadGEM2-ES/r1i1p1/zg_Amon_HadGEM2-ES_historical_r1i1p1_*.nc")
        ;tmp1 = addfiles(fils1, "cat")
        ;zg   = tmp1[:]->zg

        ;zg_xo = linint1_n_Wrap(zg&lat,zg,False,var&lat,0,2)
        ;indmin = ind(min(var&time).eq.zg&time) 
        ;indmax = ind(max(var&time).eq.zg&time)  

        ;weight = new((/dimsizes(var&time),dimsizes(var&plev),dimsizes(var&lat),dimsizes(var&lon)/),typeof(var))
        ;tmp_var= new((/dimsizes(var&time),dimsizes(var&plev),dimsizes(var&lat),dimsizes(var&lon)/),typeof(var))
        ;do t = 0,dimsizes(var&time)-1
        ;    do k = 1, dimsizes(var&plev)-1
        ;        do i = 0,dimsizes(var&lat)-1
        ;            do j = 0,dimsizes(var&lon)-1
        ;                weight = (z(k,i,j)-zg_xo(indmin+t,zgk,i,j)) / (z(k,i,j) - z(k-1,i,j))
        ;                tmp = weight*var(t,k,i,j) + (1-weight)*var(t,k,i,j)
        ;                tmp_var(t,k,i,j) = tofloat(tmp)
        ;                delete(weight)
        ;            end do
        ;        end do
        ;    end do
        ;end do

        ;tmp_var = linint1_n_Wrap(z&lev,z,False,zg&plev,0,1)
        ;printVarSummary(tmp_var)
        ;var&plev = tmp_var;mp_p
        ;delete([/tmp_var,weight,zg,zg_xo,indmin/])

      if (isdim(var,"plev")) then
        var&plev = var&plev(::-1)
        ret = 0
      end if
    end if

     ; The return value is 0 if the fix has been applied, and 1 otherwise
  return(ret)

end

