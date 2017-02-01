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
