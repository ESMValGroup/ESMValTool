;#######################################################################
;# Model-specific fixes for reformatDEFAULT_main                       #
;# Author: Mattia Righi (DLR, Germany)                                 #
;#######################################################################
;# PROJECT:     CCMVal2                                                #
;# ENSEMBLE:    REF-B1                                                 #
;# MODEL:       CCMVal2_UMUKCA-UCAM_fix.ncl                            #
;# VARIABLES:   va                                                     #
;# FIELDS:      T2Mz                                                   #
;# DESCRIPTION: badly defined latitude coordinate                      #
;#######################################################################

undef("fixing_func")
function fixing_func(var:numeric,name:string)
local f,tmp
begin

  ret = 1

  if (name.eq."va".and.FIELD.eq."T2Mz") then
      f=addfile(INFILE,"r")
      tmp=f->latitude
      var&lat = tmp
      delete(tmp)
      delete(f)
      ret = 0
  end if

  return(ret)

end

