;
;   Read in CCSM output and save into a time series file
;
;   2008-12-01  CAF
;   2008-12-30  CAF  added interpolation to pressure levels for monthly data.
;   2009-01-16  CAF  fixed error with pressure level interpolation
;   2009-01-27  CAF  Changed order of level interp and zonal averaging
;   2009-02-06  CAF  Improved the creating of a new netcdf file
;   2009-03-04  CAF  add hyai, hybi, hybm, hybi, and P0 to files
;                    temp fix for 365 day years
;   2009-03-12 CAF  Fixed problem with hard coding year to get around no leap year model runs
;   2009-03-20 CAF   Modified code to save every 3 years of the instant. field for before 1990 and
;                    after 2005
;   2009-03-24 CAF   Fixed error with instanteous derived fields not being saved
;   2009-05-22 CAF   Fixed error with climo fields being skipped
;   2009-08-07 CAF   Subtract 1 from time for monthly fields.
;   2009-08-11 CAF   Subtract 2 from time for monthly fields.
;   2009-08-12 CAF   Subtract 1 from time for monthly fields.
;   2009-08-13 CAF   add T2Dz to level interp.
;   2011-04-14 ZL    change 'lev'->'plev' and type of start_year.
;   2011-05-16 HS    Fixed error with output of time variable.
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

load "./interface_scripts/array_append.ncl"
load "./interface_scripts/converttofield.ncl"

output_file=getenv("OUTPUT_FILE")
variable=getenv("VARIABLE")
old_variable=getenv("OLD_VARIABLE")
orig_variable=getenv("ORIG_VARIABLE")
field_number=getenv("FIELD_NUMBER")
field_number_orig=getenv("ORIG_FIELD_NUMBER")
model_dir=getenv("MODEL_DIR")
start_year=stringtointeger(getenv("START_DATE"))
end_year  =stringtointeger(getenv("END_DATE"))
project_name=getenv("PROJECT_NAME")

;
; Define pressure levels to interp to
;
plevo=(/1000,850,700,500,400,300,250,200,170,150,130,115,100,90,80,70,50,30,20,15,10,7,5,3,2,1.5,1,0.5,0.3,0.2,0.1/)

;
; Determine files to be read in
;
;all_files=systemfunc("ls "+model_dir+"/*cam2.h?.*.nc")
; all_files=systemfunc("cd "+model_dir+" ; ls *cam2.h?.*.nc")
 all_files=systemfunc("cd "+model_dir+" ; ls *.nc")

start_file=-999
end_file=-999

;
;start file
;
do i=0,dimsizes(all_files)-1,1
  fin=addfile(model_dir+"/"+all_files(i),"r")
  time=fin->time
;  print(time)
  if ((field_number.eq."T3M".or.field_number.eq."T2Mz").or.field_number.eq."T2Ms") then
    time=time-1     
  end if
;caf  time@calendar="standard"
  ds=stringtocharacter(time@units) 
  daysfromyear=(charactertointeger(ds(11))-48)*1000+(charactertointeger(ds(12))-48)*100+(charactertointeger(ds(13))-48)*10+(charactertointeger(ds(14))-48)
  tim=floor(time/365)+daysfromyear   
;caf  tim=cd_calendar(time,3)
   
;  print(tim)
  do j=0,dimsizes(tim)-1,1
;caf    if (round(tim(j)/1000000,0).eq.start_year) then

    if (tim(j).eq.start_year) then
      start_time_index=j
      start_file=i
      break
    end if
  end do
  delete(time)
  delete(tim)
  if (start_file.ne.-999) then 
    break
  end if
end do
;
;end file
;
do i=0,dimsizes(all_files)-1,1
  fin=addfile(model_dir+"/"+all_files(i),"r")
  time=fin->time
  if ((field_number.eq."T3M".or.field_number.eq."T2Mz").or.field_number.eq."T2Ms") then
    time=time-1
  end if
;caf  time@calendar="standard"
;caf  tim=cd_calendar(time,3)
  ds=stringtocharacter(time@units) 
  daysfromyear=(charactertointeger(ds(11))-48)*1000+(charactertointeger(ds(12))-48)*100+(charactertointeger(ds(13))-48)*10+(charactertointeger(ds(14))-48)
  tim=floor(time/365)+daysfromyear
  do j=0,dimsizes(tim)-1,1
;caf    if (round(tim(j)/1000000,0).eq.end_year+1) then
    if (tim(j).eq.end_year+1) then
      if (j.eq.0) then
        end_file=i-1
        end_time_index=j_last
      else
        print("End filename="+all_files(i))
        end_time_index=j-1
        end_file=i
        end if
      break
    end if
    j_last=j
  end do
  delete(time)
  delete(tim)
  if (end_file.ne.-999) then 
    break
  end if
end do
if (end_file.eq.-999) then
  end_time_index=j_last
  end_file=dimsizes(all_files)-1
end if

; change me please, BAD hard coding
end_file=end_file+1
if (end_file.ge.dimsizes(all_files)) then 
  end_file=dimsizes(all_files)-1
end if

;print(all_files)
print("start_file="+all_files(start_file))
print("end_file="+all_files(end_file))

;
;-----------------------------------------------
;

;caf print("Output filenmame="+output_file)
;
; Create output file
;
error=systemfunc("rm -rf "+output_file)
fout = addfile(output_file,"c")
dim_names = "time" ; create "time" as unlimited dimension
dim_sizes = -1
dimUnlim = True
filedimdef (fout, dim_names, dim_sizes, dimUnlim) 
current_index=0

fatt = True

;
; Main loop to read in data and append to an array
;
first_write=True
file_pointer=0
do i=start_file,end_file,1
  fin=addfile(model_dir+"/"+all_files(i),"r")
  print("ccsm.ncl: Reading file "+model_dir+"/"+all_files(i))
  vNames=getfilevarnames(fin)
  if (sum(where(vNames.eq.old_variable,1,0)).eq.0) then
    print(" ")
    print("ccsm.ncl: ERROR! "+old_variable+" is not in file "+model_dir+"/"+all_files(i))
    print(" ")
    system("rm "+output_file)
    exit
  end if
  temp = fin->$old_variable$
  if ((field_number.eq."T3M".or.field_number.eq."T2Mz").or.field_number.eq."T2Ms") then
    temp&time=temp&time-1
  end if
  
  if ((field_number.eq."T3M".or.field_number.eq."T2Mz".or.field_number.eq."T2Dz").and.project_name.eq."CCMVal2") then
    ps=fin->PS
    hyam=fin->hyam
    hybm=fin->hybm
    p0=0.01*fin->P0
    var=vinth2p(temp,hyam,hybm,plevo,ps,1,p0,1,True)
    var!1="plev"
    delete(temp)
    temp=var
    delete(var)
    delete(ps)
  end if
  var=converttofield(temp,field_number,variable,orig_variable)
  delete(temp)
  time = fin->time
  if ((field_number.eq."T3M".or.field_number.eq."T2Mz").or.field_number.eq."T2Ms") then
    time=time-1
  end if
;caf  time@calendar="standard"
  s=dimsizes(var)
  rank=dimsizes(s)
  if (i.eq.start_file.and.project_name.ne."CCMVal2") then
    fout->hyam=(fin->hyam)
    fout->hybm=(fin->hybm)
    fout->hyai=(fin->hyai)
    fout->hybi=(fin->hybi)
    fout->p0=(fin->P0)
  end if

  do j=0,s(0)-1
    current_year=floor(time(j)/365+daysfromyear)
    if current_year.ge.start_year.and.current_year.le.end_year then
      if (current_year.ge.1990.and.current_year.lt.2006).or. \
         field_number_orig.ne."T3I".or. \
         (mod(current_year-1,3).eq.0) then

        if first_write then
          atts=getvaratts(fin)
          do  k=0,dimsizes(atts)-1,1
            fatt@$atts(k)$=fin@$atts(k)$
          end do
          fileattdef( fout, fatt )

          ; define time, add attributes and copy value only
          filevardef(fout, "time",typeof(var&time),"time")
          filevarattdef( fout, "time", var&time )
          fout->time=(/var&time(j)/)

          if (rank.eq.4)
            fout->$variable$=var(j:j,:,:,:)
          end if
          if (rank.eq.3)
            fout->$variable$=var(j:j,:,:)
          end if
          if (rank.eq.2)
            fout->$variable$=var(j:j,:)
          end if
          if (rank.eq.1)
            fout->$variable$=var(j:j)
          end if
          first_write=False
        else
          fout->time(file_pointer)=(/var&time(j)/)
          if (rank.eq.4)
            fout->$variable$(file_pointer,:,:,:)=(/var(j,:,:,:)/)
          end if
          if (rank.eq.3)
            fout->$variable$(file_pointer,:,:)=(/var(j,:,:)/)
          end if
          if (rank.eq.2)
            fout->$variable$(file_pointer,:)=(/var(j,:)/)
          end if
          if (rank.eq.1)
            fout->$variable$(file_pointer)=(/var(j)/)
          end if
        end if
        file_pointer=file_pointer+1
      end if
    end if
  end do
  delete(var)
  delete(time)
end do
;
; Save data to output file
;

return
