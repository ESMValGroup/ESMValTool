; get_field_number.ncl
; fields supported T3M T3D T2Mz T2Dz T2Ms T2Ds T1Ms T1Ds

undef("get_field_number")
function get_field_number(A:numeric)
local fieldno, ds, rank, time, days
begin
   fieldno = ""

   ds = dimsizes(A)
   rank = dimsizes(ds)

   time = A&time
   days = time(1)-time(0)

   ;T1Ms T1Ds
   if (rank .eq. 2) then
      if (days .lt. 2) then
         fieldno = "T1Ds"
      else
         fieldno = "T1Ms"
      end if
   end if

   ;T2Mz T2Dz T2Ms T2Ds
   if (rank .eq. 3) then
      if (isdim(A,"lon")) then
         ;T2Ms T2Ds - "lon"
         if (days .lt. 2) then
            fieldno = "T2Ds"
         else
            fieldno = "T2Ms"
         end if
      else
         ;T2Mz T2Dz - "plev"
         if (days .lt. 2) then
            fieldno = "T2Dz"
         else
            fieldno = "T2Mz"
         end if
      end if
   end if

   ;T3M T3D
   if (rank .eq. 4) then
      if (days .lt. 2) then
         fieldno = "T3D"
      else
         fieldno = "T3M"
      end if
   end if

   if (fieldno .eq. "") then
     print("fatal:function 'get_field_number' only support the following fields T3M T3D T2Mz T2Dz T2Ms T2Ds T1Ms T1Ds")
     status_exit(1)
   end if

   return(fieldno)

end
