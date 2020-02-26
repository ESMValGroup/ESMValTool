;EYRING06JGR_FIG08
;; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
;; PROJECT-NAME CCMI
;;#############################################################################
;; Description
;;  Time-height sections of water vapor mixing ratio shown as the deviation (in parts per million by volume)
;; from the time mean profile, averaged between min latitude (10°S) and max latitude (10°)
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    @e06fig08_multimean: calculate multi-model mean (i.e. False/True)
;;                         bias calculation (e.g. (/"ERA-Interim"/))
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    @e06fig08_cnMax:   max contour
;;    @e06fig08_cnMmin:   min contour
;;    @e06fig08_cnSpacing: contour interval
;;    @e06fig08_latrange: name of the lat range as in the "Equator"
;; Required variable_info attributes (variable specific)
;;
;; Caveats
;;
;; Modification history
;;    20180320-A_laue_ax: added tags for reporting and call to write_references
;;    20170315-A_cion_ir: written
;;#############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"

begin
   enter_msg(DIAG_SCRIPT, "")
   var0 = variable_info[0]@short_name
   info0 = select_metadata_by_name(input_file_info, var0)
   dim_MOD1 = ListCount(info0)
   names = metadata_att_as_array(info0, "dataset")
   projects = metadata_att_as_array(info0, "project")
   exps = metadata_att_as_array(info0, "exp")
   ensembles = metadata_att_as_array(info0, "ensemble")
   y1 =(/ metadata_att_as_array(info0, "start_year")/)
   y2 =(/ metadata_att_as_array(info0, "end_year") /)
   
   iref=0
   if (isatt(variable_info[0], "reference_dataset")) then
     refname = variable_info[0]@reference_dataset
     do imod = 0, dimsizes(names) - 1
         if (names(imod) .eq. refname) then
            iref = imod
          end if
      end do       
   end if

   if(isatt(variable_info, "units")) then
        UNITS = variable_info@units
   else
        UNITS= "K"
   end if
   
   if(isatt(variable_info, "long_name")) then
        LONG_NAME = variable_info@long_name
   else
        LONG_NAME = var0
   end if

   log_info("++++++++++++++++++++++++++++++++++++++++++")
   log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
   log_info("++++++++++++++++++++++++++++++++++++++++++")
  
   ; Set default values for optional diag_script_info attributes

   set_default_att(diag_script_info, "styleset", "CMIP5")
   set_default_att(diag_script_info, "e06fig08_cnMin", -1.)
   set_default_att(diag_script_info, "e06fig08_cnMax", 1.)
   set_default_att(diag_script_info, "e06fig08_cnSpacing", 0.2)
   set_default_att(diag_script_info, "e06fig08_multimean", False)
   
   set_default_att(diag_script_info, "e06fig08_latrange", "10S-10N")


   flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
   index_mod = ind(flag_mod.gt.0)
   index_obs = ind(flag_mod.eq.0)
   dim_MOD = dimsizes(index_mod)
   dim_OBS = 0

   if (.not.all(ismissing(index_obs))) then
     dim_OBS = dimsizes(index_obs)
     list_obs_mod = names(index_obs)
   end if

         
   cnMax = diag_script_info@e06fig08_cnMax
   cnMin = diag_script_info@e06fig08_cnMin
   cnSpacing =  diag_script_info@e06fig08_cnSpacing
   latrange=diag_script_info@e06fig08_latrange
   multimean = diag_script_info@e06fig08_multimean
   if (multimean)
       dimTot=dim_MOD+dim_OBS+1
   else
       dimTot=dim_MOD+dim_OBS
   end if

  ; make sure path for (mandatory) netcdf output exists

   work_dir = config_user_info@work_dir + "/"
  ; Create work dir
   system("mkdir -p " + work_dir)

end
begin
  aux_title_info = ""
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f08_" + var0)
  output_filename = wks@fullname
  nPanel = dimTot
  plot = new(nPanel, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))
  rad = 4.0 * atan(1.0) / 180.0
  
  A0 = read_data(info0[iref])
  ;printVarSummary(A0)
  lev0 = A0&plev
  dimLev = dimsizes(lev0)
  delete(A0)

 
  tmpData = new((/dim_MOD+dim_OBS,25,dimsizes(lev0),2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "month"
  tmpData!2 = "lev"
  tmpData!3 = "stat"

  Xdata = new((/dimTot, 25,dimsizes(lev0),2/), "double")
  Xdata!0 = "line"
  Xdata!1= "month"
  Xdata!2 = "lev"
  Xdata!3 = "stat"



  plotMod=new(dimTot,typeof(names))
  
  ik=0


;******************************************************************
  do imod = 0, dim_MOD1 - 1
     print(imod)
     ; Extract data
     A0 = read_data(info0[imod])
     if (var0.eq."hus") then
           A0=A0/0.62
           UNITS="ppmv"
           A0=A0*10^6
     end if

     dummy1 = time_operations(A0,y1(imod),y2(imod),"average","monthlyclim",True)
     dummy2_b=tofloat(A0)
     date= yyyymm_time(y1(imod),y2(imod),"integer") 
     yyyy  = date/100 
     mm  = date-yyyy*100                ; mmdd = yyyymmdd%10000
     weights = days_in_month(yyyy, mm)
     dummy2 = dummy2_b(0:11,: )  ; Copy metadata
     do ilev=0,dimsizes(lev0)-1 
        do im = 0, 11
;            print(dummy2_b(im::12,ilev))
            dummy2(im,ilev ) = dim_stddev_Wrap(dummy2_b(im::12,ilev))
             ;dummy2(im,ilev ) = \
             ;  dim_stddev_wgt_Wrap(dummy2_b(im::12,ilev), weights(im::12), 1)
        end do
     end do
     
     dummy3=dummy1
     dummy3=(/dummy1-conform(dummy1,dim_avg_n(dummy1,0),1)/)
     tmpData(imod,0:11,:,0) = (/dummy3/)
     tmpData(imod,0:11,:,1) = (/dummy2/)
     tmpData(imod,12:23,:,0) = (/dummy3/)
     tmpData(imod,12:23,:,1) = (/dummy2/)
     tmpData(imod,24,:,0) = (/dummy3(0,:)/)
     tmpData(imod,24,:,1) = (/dummy2(0,:)/)
     delete(A0)
     delete(weights)
     delete(date)
     delete(yyyy)
     delete(mm)
     delete(dummy2_b)
     delete(dummy1)
     delete(dummy2)
     delete(dummy3)
 end do
 Xdata(:dim_MOD1-1,:,:,:)=(/tmpData/)
 plotMod(:dim_MOD1-1)=names
         

  ;°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  if (multimean) then
             Xdata(dimTot-1,:,:,0)=(/dim_avg_n(Xdata(index_mod,:,:,0),0)/)
             plotMod(dimTot-1)="MMM"
             X_MMM_std=(/dim_stddev_n(Xdata(index_mod,:,:,:),0)/)

  end if
; print(Xdata(line|:,season|:,day|:,stat|0))

  XXData=ispan(1,25,1)
  XXData!0="month"
  XXData&month=XXData

  Ydata= lev0/100.
  Ydata!0="lev"
  Ydata&lev=lev0/100.
  Ydata@units="hPa"   
  plot_d=new(dimTot,"graphic")
;*************************plotting*******************************************************
  do imod=0, dimTot-1
      res                   = True
      res@tiMainString         =plotMod(imod)     ; title
      res@gsnDraw              = False
      res@gsnFrame             = False
      res@lbLabelBarOn         = False
      res@cnLineLabelsOn       = False
      res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
      res@cnLevelSpacingF      = cnSpacing                   ; contour interval
      res@cnMinLevelValF       = cnMin                 ; min level
      res@cnMaxLevelValF       = cnMax                   ; max level
      res@lbLabelStrings       =(/"-1","-0.8","-0.6","-0.4","-0.2","0.","0.2","0.4","0.6","0.8","1."/)
      res@cnFillOn             = True                  ; turn on color fill
      res@cnInfoLabelOn        = False                 ; turn off contour label
      res@gsnSpreadColors      = True                  ; use full range of colors
      res@gsnSpreadColorStart  = 10
      res@lbLabelAutoStride    = True                  ; optimal labels
      res@gsnLeftString        = ""
      res@gsnRightString       = ""
      res@gsnLeftStringFontHeightF = 0.015
      res@tmXBMode             = "Explicit"            ; explicit labels
      ; location of labels
      res@gsnAddCyclic         = True
      res@tmXBValues = (/1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25/)
      res@tmXBLabels = (/"J","F","M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A","S","O","N","D","J"/)                ; labels themselves
      res@trYMaxF              = 100
      res@trYMinF              = 10

      res@trXMinF           = 0 
      res@trXMaxF           = 25
      res@vpHeightF            = 0.4                    ; change aspect ratio of plot
      res@vpWidthF             = 0.75
      res@lbLabelFontHeightF   = 0.011
      res@tiYAxisString        = "Pressure [hPa]"    ; y-axis label
      Xdata!2 = "lev"
      Xdata&lev=Ydata
      plot_d(imod) = gsn_csm_pres_hgt (wks,Xdata(line|imod,lev|:,month|:,stat|0),res)
      delete(res)
   end do
   ;; Call legend-plotting routine
   Pres = True
   Pres@gsnFrame       = False
   Pres@gsnPanelBottom = 0.2
   Pres@gsnPanelFigureStrings=plot_num(0:dimTot-1) 
   Pres@gsnPanelFigureStringsFontHeightF=0.015
   Pres@gsnPanelFigureStringsPerimOn="False"
   Pres@amJust="TopLeft"
   Pres@lbLabelAutoStride    = True 
   Pres@gsnPanelLabelBar = True
;; Call panel routine
   gsn_panel(wks,plot_d,(/1+dimTot/3,3/),Pres) 
 frame(wks)
;**********************common output file**************************
     nc_filename = work_dir + "eyring06f08_" + var0 + ".nc"
     Xdata@var = var0
     Xdata@diag_script="eyring06f08"
     ;printVarSummary(XXdata(line|:,season|ipanel,month|:,stat|0))
     printVarSummary(Xdata)
     Xdata&stat=(/"mean","avg"/)
     Xdata&line=plotMod
     month_names=(/"J","F","M","A","M","J","J","A","S","O","N","D","J","F","M","A","M","J","J","A","S","O","N","D","J"/)
     Xdata&month=month_names
     printVarSummary(Xdata)
     nc_outfile = ncdf_write(Xdata, nc_filename)

end

