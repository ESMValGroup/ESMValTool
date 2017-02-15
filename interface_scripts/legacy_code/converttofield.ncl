
;
;  Check the different field numbers to see if anything needs to be done
;  to them.  For 3D fields, just return
;            For surface, make sure is lat,lon then return
;            For zonal, try to calculate the zonal means
;
;  2008/11/26  CAF  Created
;


undef("convertotfield")
function converttofield (var, field, variable, orig_variable)
begin
  dim=dimsizes(var)
  rank=dimsizes(dim)
  dimNames=getvardims(var)
;
; 3D fields
;
  if (field.eq."T3M".or.field.eq."T3D".or.field.eq."T3I") then
    if (rank.eq.4) then
      return(var)
    end if    
  end if
;
; Zonal means
;
  if (field.eq."T2Mz".or.field.eq."T2Dz".or.field.eq."T2Iz") then
    if (rank.eq.3) then
      if ((ind(dimNames.eq."time").ge.0).and.(ind(dimNames.eq."lat").ge.0).and.(ind(dimNames.eq."lev").ge.0)) then
        return(var)
      end if
    end if
    if (rank.eq.4) then
      if ((ind(dimNames.eq."time").ge.0).and.(ind(dimNames.eq."lat").ge.0).and.(ind(dimNames.eq."lon").ge.0).and.(ind(dimNames.eq."lev").ge.0)) then
        result=dim_avg_Wrap(var)
        return(result)
      end if
    end if
  end if
;
; Surface  lat,lon.  Not calulated, already needs to be in lat,lon
; Do something special for derived surface lat,lon.  Like ua @ 200
;
  if (field.eq."T2Ms".or.field.eq."T2Ds".or.field.eq."T2Is") then
    if (rank.eq.3) then
      if ((ind(dimNames.eq."time").ge.0).and.(ind(dimNames.eq."lat").ge.0).and.(ind(dimNames.eq."lon").ge.0)) then
        return(var)
      end if
    end if
    if (orig_variable.ne.variable) then
      if (rank.eq.4) then
        return(var)
      end if
    end if
  end if
  print("converttofield.ncl: ERROR!  Problem with field convert to "+field+".  There is most likely a problem with the coordinates.")  
  exit
  return(result)
end

