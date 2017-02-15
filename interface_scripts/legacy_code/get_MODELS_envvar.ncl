;
;  This function retrieves the MODELS enviroment variable
;  and parces out the information into a variable and returns it. 
;
; 2008-06-27  CAF
; 2008-11-06  CAF    added ensemble id
; 

undef("get_MODELS_envvar")
function get_MODELS_envvar()
begin
  colon=inttochar(58)
  lst=cshstringtolist(getenv("MODELS"))
  len=stringtoint(lst(0))
  models=True
  models@name=new(len,string)
  models@case_name=new(len,string)
  models@ensemble=new(len,string)
  models@start_year=new(len,string)
  models@end_year=new(len,string)
  models@dir=new(len,string)
  ii=0
  do i=0,dimsizes(lst)-1,1
    index=ind(stringtochar(lst(i)) .eq. colon(0))
    if (ismissing(index) .eq. False) then
      i=i+1
      if (lst(i-1) .eq. "name:") then
        models@name(ii)=lst(i)
      end if
      if (lst(i-1) .eq. "case:") then
        models@case_name(ii)=lst(i)
      end if
      if (lst(i-1) .eq. "ensemble:") then
        models@ensemble(ii)=lst(i)
      end if
      if (lst(i-1) .eq. "start_year:") then
        models@start_year(ii)=lst(i)
        ii=ii+1
      end if
      if (lst(i-1) .eq. "end_year:") then
        models@end_year(ii)=lst(i)
      end if
      if (lst(i-1) .eq. "dir:") then
        models@dir(ii)=lst(i)
      end if
    end if
  end do
  return(models)
end
