; This script takes all the metrics created by the various scripts and placed
; in metrics_orig.txt, calculates the total scores, reorganizes the data, 
; and writes out a new metrics.txt file. 

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

; A_laue_ax+
load "$CVDP_SCRIPTS/functions.ncl"
; A_laue_ax-

begin
  print("Starting: metrics.ncl")
  
  OUTDIR = getenv("OUTDIR")
  
  nsim = numAsciiRow("namelist_byvar/namelist_trefht")    ; retrieve total number of observational + models (all namelist_byvar/namelist have same # of rows)
  na = asciiread("namelist_byvar/namelist_trefht",(/nsim/),"string")  ; (It is not done via metrics.txt as there might be a space in names.                                                                             
  blankrow = ind(na.eq."")
  if (.not.any(ismissing(blankrow))) then
     goodrows = ind(na.ne."")
     na2 = na(goodrows)
     delete(na)
     na = na2
     delete(na2)
     nsim = dimsizes(na)
  end if      
  nsim = nsim - 1    ; first listed dataset is what all others are compared to, thus, output metrics table has nsim-1 column  

  files = (/"sst.indices.1","sst.indices.2","amo","pdo","psl.nam_nao","psl.sam_psa",\
            "sst.mean_stddev","psl.mean_stddev","pr.mean_stddev"/)
  files = "metrics."+files
  files = files+".txt"
  
  do gg = 0,dimsizes(files)-1
; A_laue_ax+
     ; function fileexists does not exist in NCL versions earlier than 6.2.1
     ; in NCL 6.2.1, the behavior of the function isfilepresent changed
     ; ---> function isfilepresent2 (defined in functions.ncl) mimics the old
     ; function isfilepresent for later versions of NCL and can be used universally
     ; with both, earlier and later versions of NCL than 6.2.1
;     if (.not.fileexists(OUTDIR+files(gg))) then
     if (.not.isfilepresent2(OUTDIR+files(gg))) then
; A_laue_ax-
        print("1 or more metrics files missing, exiting metrics.ncl")
        exit
     end if
  end do

  ch = new((/nsim,dimsizes(files)/),"string")   ; hold obs/simulation names
  patcor_rms = new((/nsim,11/),"string")   ; 11 metrics  
  cntr = 0

  do gg = 0,dimsizes(files)-1
     nrow = numAsciiRow(OUTDIR+files(gg))
     a = asciiread(OUTDIR+files(gg),(/-1/),"string")
     
     t0 = tochar(a(3))
     sti0 = str_index_of_substr(a(4), " -",0)   ; read in individual column headers from each metrics file
     do hh = 0,dimsizes(sti0)-1
        if (hh.eq.(dimsizes(sti0)-1)) then
           ch(hh,gg) = str_strip(tostring(t0(sti0(hh):)))
        else 
           ch(hh,gg) = str_strip(tostring(t0(sti0(hh):sti0(hh+1))))
        end if
     end do
     delete([/sti0,t0/])

     test = tochar(a(5:))
     if (dimsizes(dimsizes(test)).eq.2) then
        patcor_rms(:,cntr) = str_split(tostring(test(0,18:))," ")
        patcor_rms(:,cntr+1) = str_split(tostring(test(1,18:))," ")
        cntr = cntr+2
     else
        patcor_rms(:,cntr) = str_split(tostring(test(18:))," ")
        cntr = cntr+1
     end if
     delete(a)
     delete([/test/])
  end do
  delete(cntr)

  do gg = 0,dimsizes(files)-1  ; remove individual metrics files
     system("rm "+OUTDIR+files(gg))
  end do

  names = ch(:,0)
  do gg = 0,nsim-1   ; check to see if data is observations or models by seeing if every name matches
     if (all(ch(gg,0).eq.ch(gg,1:))) then
        names(gg) = ch(gg,0)
     else
        names(gg) = "OBS "+(gg+2)
     end if
  end do

  names_nchar = max(dimsizes(tochar(names)))
  spacer = ""
  do gg = 0,names_nchar
     spacer = spacer+" "
  end do
  delete(names_nchar)     
  pc_score = new(nsim,"string")
  rms_score = new(nsim,"string")
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
           if (pc(ii).eq.1.0) then
              pc_z(ii) = 0.5*(log( (1+1.001) / (1-1.001) ))    ; needed when pattern correlation = 1.0
           else
              pc_z(ii) = 0.5*(log( (1+pc(ii)) / (1-pc(ii)) ))
           end if
        end do
        zavg = avg(pc_z)   ; compute average of z
        delete(pc_z)
        
        pc_score(gg) = sprintf("%4.2f",((2.71828^(2*zavg))-1)/ ((2.71828^(2*zavg))+1))  ; reverse process and convert z-avg -> r.
;        print("average of Z-tranformed correlations = "+pc_score(gg))                   ; r = (e^2Z - 1)/(e^2Z+1)  ; e = 2.71828
        delete(zavg)
     end if
     delete(pc)
  end do
  pc_score = where(pc_score.eq." nan","----",pc_score)  ; needed for when the nan's come out of the z-transform (likey due to numerous pattern correlations = 1)
  
  header = (/"","Pattern Correlations/RMS Differences   Observations vs. Model(s)",""/)
  write_table(OUTDIR+"metrics.txt","w",[/header/],"%s")
  column_header1 = spacer+"   ENSO TAS   ENSO PSL   El Nino    La Nina      AMO        PDO        NAM        SAM     SST sigma  PSL sigma   PR sigma    Mean  "
  column_header2 = spacer+"    (DJF+1)    (DJF+1)     Hov        Hov     (Monthly)  (Monthly)    (DJF)      (DJF)      (Ann)      (Ann)       (Ann)     Score  "
  column_header3 = spacer+"  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------  ---------"
  write_table(OUTDIR+"metrics.txt","a",[/column_header1/],"%s")
  write_table(OUTDIR+"metrics.txt","a",[/column_header2/],"%s")
  write_table(OUTDIR+"metrics.txt","a",[/column_header3/],"%s")
  
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
     write_table(OUTDIR+"metrics.txt","a",[/srow/],"%s")
     delete([/spacer_char1,dimC,mname_char,srow/])
  end do
  delete([/patcor_rms,spacer_char,pc_score,rms_score/])  
  print("Finished: metrics.ncl")
end
