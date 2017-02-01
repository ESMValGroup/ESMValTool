;********************************************************
; S. Murphy
; converts a list of comma deliminated variables passed from
; csh to NCL to an actual array of strings.
; e.g. set vars = ,T,U,V,TS,PRECC, (in csh)
;      setenv nclvars $vars        (to pass to NCL must be and env)
;      vs = getenv("nclvars")      (reads csh variable into NCL)
;      vars = cshstringtolist(vs)  (does conversion)
;      vars = (/"T","U","V","TS","PRECC"/)

; note the string "list" is now an NCL keyword, so "list" has been
; shortened to "lst"

undef("cshstringtolist")
function cshstringtolist(cshlist:string)
local lst, breaks, nwords, str, i
begin

 lst=stringtochar(cshlist)            ;convert the string to characters
 breaks = ind(lst .eq. inttochar(44)) ;locate the deliminators (a space is 32)
 nwords=dimsizes(breaks)              ;count the deliminators

 str=new((/nwords-1/),string)
 do i=0,nwords-2                      ;pull out the separate strings
    str(i) = chartostring (lst(breaks(i)+1:breaks(i+1)-1))
 end do
 return(str)
end


