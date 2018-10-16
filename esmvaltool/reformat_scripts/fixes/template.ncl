; TEMPLATE FOR MODEL-SPECIFIC FIXES
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
local ...
begin

  ret = 1

  ; Use if-statements if the fix is to be applied only to certain
  ; variable(s) and/or ensemble(s) and/or field(s)
  if (name.eq...) then
      ...
      ...
      ret = 0
  end if

  if (ENSEMBLE.eq...) then
      ...
      ...
      ret = 0
  end if

  if (FIELD.eq...) then
      ...
      ...
      ret = 0
  end if

  ; The return value is 0 if the fix has been applied, and 1 otherwise
  return(ret)

end
