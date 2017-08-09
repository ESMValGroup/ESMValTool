;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR-Germany)                                  #
;#######################################################################
;# PROJECT:     CMIP5                                                  #
;# ENSEMBLE:    [all]                                                  #
;# MODEL:       BNU-ESM                                                #
;# VARIABLES:   fgco2, sftlf                                           #
;# FIELDS:      [all]                                                  #
;# DESCRIPTION: inconsistent units (fgco2), wrong latitude (sftlf)     #
;#              inconsistent units (ch4,co2)                           #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
local f, files, tmp
begin

  ret = 1

  if (name.eq."fgco2") then
      var =  var * 12. / 44.
      ret = 0
  end if

  ;; Converting ch4 units from kg/kg to ppmb
  if (name.eq."ch4") then
      var =  var * 29/16 * 1.e9
      ret = 0
  end if

  ;; Converting co2 units from kg/kg to ppmb
  if (name.eq."co2") then
      var =  var * 29/44 * 1.e6
      ret = 0
  end if

  if (name.eq."spco2") then
      idx = ind(models@name.eq.MODEL)
      mip = models@mip(idx)
      if (mip.eq."Omon") then
          var =  var * 1.e6
          ret = 0
      end if
  end if

  if (name.eq."sftlf") then
      files = systemfunc("ls " + INFILE)
      f=addfile(files(0), "r")
      tmp=f->lat
      var&lat = tmp
      delete(tmp)
      delete(f)
      ret = 0
  end if

  return(ret)

end
