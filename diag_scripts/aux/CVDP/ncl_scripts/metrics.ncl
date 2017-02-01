; This script takes all the metrics created by the various scripts and placed
; in metrics_orig.txt, calculates the total scores, reorganizes the data, 
; and writes out a new metrics.txt file. 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  print("Starting: metrics.ncl")
  nsim = numAsciiRow("namelist")                ; retrieve simulation names from namelist
  na = asciiread("namelist",(/nsim/),"string")  ; (It is not done via metrics.txt as there might be a space in names.                                                                             
  blankrow = ind(na.eq."")
  if (.not.any(ismissing(blankrow))) then
;     print("Extra blank row(s) detected in namelist. Check the namelist and make sure there is only one blank line at the end of namelist. Attempting to remove extra blank rows.")
;     print("-----Original namelist---")
;     print(""+na)
;     print("-------------------------")
     goodrows = ind(na.ne."")
     na2 = na(goodrows)
     delete(na)
     na = na2
     delete(na2)
     nsim = dimsizes(na)
;     print("-----Corrected namelist--")
;     print(""+na)
;     print("-------------------------")
  end if      
  
  names = new(nsim,"string")                                                 
  delim = "|"
  do gg = 0,nsim-1
     names(gg) = str_strip(str_get_field(na(gg),1,delim))
  end do
  names_nchar = max(dimsizes(tochar(names)))
  spacer = ""
  do gg = 0,names_nchar
     spacer = spacer+" "
  end do
  delete(names_nchar)

  patcor_rms = new((/nsim,11/),"string")   ; 11 metrics
  pc_score = new(nsim,"string")
  rms_score = new(nsim,"string")

  nrow = numAsciiRow(getenv("OUTDIR")+"metrics_orig.txt")
  if (nrow.ne.16) then
     print("Unexpected number of rows encountered in metric_orig.txt. Expected 16, got "+nrow)
     print("Not calculating total scores for each model, and copying metrics_orig.txt to metrics.txt")
     system("cp "+getenv("OUTDIR")+"metrics_orig.txt "+getenv("OUTDIR")+"metrics.txt")
     exit
  end if
  a = asciiread(getenv("OUTDIR")+"metrics_orig.txt",(/-1/),"string")    ; a(3) = row headers, a(5:15) = various metrics
;  print(a)
  do gg = 0,10   ; 11 metrics
     test = tochar(a(gg+5))
     patcor_rms(:,gg) = str_split(tostring(test(18:))," ")
     delete(test)
  end do

;  print(patcor_rms)
;  exit

  do gg = 0,nsim-1     ; strip out pattern correlations, and calculated score for each model
     pc = new(11,float,9.99)
     rms = pc
     do hh = 0,10   ; 11 metrics
        n1 = str_split(patcor_rms(gg,hh),"/")
;        print(n1)
        pc(hh)  = tofloat(n1(0))   ; strip out pattern correlations. 9.99 = missing.
        rms(hh) = tofloat(n1(1))   ; strip out rms. 9.99 = missing.
        delete(n1)
     end do
     if (any(ismissing(rms))) then
        rms_score(gg) = "----"
     else
        rms_score(gg) = sprintf("%4.2f",avg(rms))
     end if
     delete(rms)

;     total_score(gg) = ""+avg(pc)
;     print("Simple average = "+avg(pc))

     pc_z = pc
     pc_z = pc_z@_FillValue
     if (any(ismissing(pc))) then
;        print("Missing Values detected")
;        print(pc)
        pc_score(gg) = "----"
     else     
        do ii = 0,10   ; use Fisher's z-transformation to translate r->z
           pc_z(ii) = 0.5*(log( (1+pc(ii)) / (1-pc(ii)) ))
        end do
        zavg = avg(pc_z)   ; compute average of z
        delete(pc_z)
        
        pc_score(gg) = sprintf("%4.2f",((2.71828^(2*zavg))-1)/ ((2.71828^(2*zavg))+1))  ; reverse process and convert z-avg -> r.
;        print("average of Z-tranformed correlations = "+pc_score(gg))                   ; r = (e^2Z - 1)/(e^2Z+1)  ; e = 2.71828
        delete(zavg)
     end if
     delete(pc)
  end do
;  print(pc_score)
;  print(rms_score)
  
  header = (/"","Pattern Correlations/RMS Differences   Observations vs. Model(s)",""/)
  write_table(getenv("OUTDIR")+"metrics.txt","w",[/header/],"%s")
  column_header1 = spacer+"   ENSO TAS   ENSO PSL   El Nino    La Nina      AMO        PDO        NAM        SAM     SST sigma  PSL sigma   PR sigma    Mean  "
  column_header2 = spacer+"    (DJF+1)    (DJF+1)     Hov        Hov     (Monthly)  (Monthly)    (DJF)      (DJF)      (Ann)      (Ann)       (Ann)     Score  "
  column_header3 = spacer+"  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------"
  write_table(getenv("OUTDIR")+"metrics.txt","a",[/column_header1/],"%s")
  write_table(getenv("OUTDIR")+"metrics.txt","a",[/column_header2/],"%s")
  write_table(getenv("OUTDIR")+"metrics.txt","a",[/column_header3/],"%s")
  
  patcor_rms = where(patcor_rms.eq."9.99/9.99","----/----",patcor_rms)
  spacer_char = tochar(spacer)
  do gg = 0,nsim-1
     spacer_char1 = spacer_char
     mname_char = tochar(names(gg))
     dimC = dimsizes(mname_char)
     spacer_char1(:dimC-1) = mname_char
     srow = tostring(spacer_char1)
;     print(srow)
     do hh = 0,10
        n1 = str_split(patcor_rms(gg,hh),"/")
;        print("n1 = "+n1)
        if (n1(0).eq."----") then
           srow = srow+"  "+patcor_rms(gg,hh)
        else
           if (tofloat(n1(0)).ge.0) then
              srow = srow+"  "+patcor_rms(gg,hh)
           else
              srow = srow+" "+patcor_rms(gg,hh)
           end if   
        end if
        delete(n1)
     end do
     srow = srow+"  "+pc_score(gg)+"/"+rms_score(gg)     
     write_table(getenv("OUTDIR")+"metrics.txt","a",[/srow/],"%s")
     delete([/spacer_char1,dimC,mname_char,srow/])
  end do
  delete([/patcor_rms,spacer_char,pc_score,rms_score/])  
  print("Finished: metrics.ncl")
end
