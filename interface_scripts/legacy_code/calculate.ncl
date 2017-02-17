;
; THIS FILE SHOULD NOT BE CHANGED!!!
;
; Defines the calculate function in order to avoid errors when
; there are no derived variables.  If there is a derived variable
; then a new calculate function needs to be defined in the
; variable_defs/*.ncl file that is assocated with the variable.
;
; 2008-06-27  CAF
;

undef("calculate")
function calculate(dummy1, dummy2, dummy3)
begin
    return(0)
end
