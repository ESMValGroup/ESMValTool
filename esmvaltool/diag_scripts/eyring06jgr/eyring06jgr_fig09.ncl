;EYRING06JGR_FIG09
;; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
;; PROJECT-NAME CCMI
;;#############################################################################
;; Description
;; Vertical variation of (a) amplitude and (b) phase lag of annual cycle of water
;; vapor averaged between 10°S and 10°N. 
;;
;; Required diag_script_info attributes (diagnostics specific)
;;    @e06fig09_multimean: calculate multi-model mean (i.e. False/True)
;;                         bias calculation (e.g. (/"ERA-Interim"/))
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    @e06fig09_bMax:   max contour
;;    @e06fig09_bMmin:   min contour
;;    @e06fig09_latrange: name of the lat range as in the "Equator"
;;    @e06fig09a_start_year: year when start the climatology calculation
;;    @e06fig09a_start_year: year when end the climatology calculation
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
   set_default_att(diag_script_info, "e06fig09_aMin", 0.)
   set_default_att(diag_script_info, "e06fig09_aMax", 1.)
   set_default_att(diag_script_info, "e06fig09_bMin", 0.)
   set_default_att(diag_script_info, "e06fig09_bMax", 20.)
   set_default_att(diag_script_info, "e06fig09_multimean", False)
   
   set_default_att(diag_script_info, "e06fig09_latrange", "10S-10N")
   set_default_att(diag_script_info, "e06fig09_start_year", max(y1))
   set_default_att(diag_script_info, "e06fig09_end_year",min(y2))

   flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
   index_mod = ind(flag_mod.gt.0)
   index_obs = ind(flag_mod.eq.0)
   dim_MOD = dimsizes(index_mod)
   dim_OBS = 0

   if (.not.all(ismissing(index_obs))) then
     dim_OBS = dimsizes(index_obs)
     list_obs_mod = names(index_obs)
   end if

         
   aMax = diag_script_info@e06fig09_aMax
   aMin = diag_script_info@e06fig09_aMin
   bMax = diag_script_info@e06fig09_bMax
   bMin = diag_script_info@e06fig09_bMin
   latrange=diag_script_info@e06fig09_latrange
   multimean = diag_script_info@e06fig09_multimean
   if (multimean)
       dimTot=dim_MOD+dim_OBS+1
   else
       dimTot=dim_MOD+dim_OBS
   end if
   year_start=toint(diag_script_info@e06fig09_start_year)
   year_end=toint(diag_script_info@e06fig09_end_year)
  ; make sure path for (mandatory) netcdf output exists

   work_dir = config_user_info@work_dir + "/"
  ; Create work dir
   system("mkdir -p " + work_dir)
   colors = project_style(info0, diag_script_info, "colors")
   dashes = project_style(info0, diag_script_info, "dashes")
   thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
   markers = project_style(info0, diag_script_info, "markers")
end
begin
  aux_title_info = ""
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f09_" + var0)
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
  plotMod=new(dimTot,typeof(names))
  plotColor=new(dimTot,typeof(names))
  plotDash=new(dimTot,typeof(names))
  plotThicks=new(dimTot,"float")
  plotScatters=new(dimTot,typeof(names))
  LineScatters=new(dimTot,typeof(names))
  
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

     dummy1 = time_operations(A0,year_start,year_end,"average","monthlyclim",True)
     dummy2_b=time_operations(A0,year_start,year_end,"extract","monthlyclim",True)
     date= yyyymm_time(year_start,year_end,"integer") 
     yyyy  = date/100 
     mm  = date-yyyy*100                ; mmdd = yyyymmdd%10000
     weights = days_in_month(yyyy, mm)
     dummy2 = dummy2_b(0:11,: )  ; Copy metadata
     do ilev=0,dimsizes(lev0)-1 
        do im = 0, 11
            dummy2(im,ilev ) = dim_stddev_Wrap(dummy2_b(im::12,ilev))
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
     plotMod(imod)=names(imod)
     plotColor(imod)=colors(imod)
     plotDash(imod)=dashes(imod)
     plotThicks(imod)=thicks(imod)
     plotScatters(imod)=markers(imod)
     LineScatters(imod)="Lines" 
end do
 Xdata(:dim_MOD1-1,:,:,:)=(/tmpData/)
LineScatters(index_obs)="Markers"

  ;°°2°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
  if (multimean) then
             Xdata(dimTot-1,:,:,0)=(/dim_avg_n(Xdata(index_mod,:,:,0),0)/)
             plotMod(dimTot-1)="MMM"
             X_MMM_std=(/dim_stddev_n(Xdata(index_mod,:,:,:),0)/)
             plotMod(dimTot-1)="MMM"
             plotColor(dimTot-1)="Black"
             plotDash(dimTot-1)=0.
             plotScatters(dimTot-1)=0.
             LineScatters(dimTot-1)="Lines"
  end if
  Ydata= lev0/100.
  Ydata!0="lev"
  Ydata&lev=lev0/100.
  Ydata@units="hPa"  
;**********************FOURIER ANALYSIS*****************************************
newXXData=Xdata
newXXData=where(ismissing(Xdata),0.,Xdata)
printVarSummary(newXXData)
newXXData!0="line"
newXXData!1="month"
newXXData!2="plev"
newXXData!3="stat"
printVarSummary(newXXData)
;*********************set dimension*******************************************
 ffa=fourier_info(newXXData(line|0,plev|:,month|0:11,stat|0),1,1.0)
 ffa!1="plev"
 ffa&plev=Ydata
 ind_pl=ind((ffa&plev.gt.96).and.(ffa&plev.le.100))
 fase=ffa(1,:,0)
 zlevel =(16*log10(10^3/Ydata))
 max_ind=ind(ffa(0,:,0).eq.(max(ffa(0,ind_pl,0))))
 ffa_m=fourier_info(-newXXData(line|0,plev|:,month|0:11,stat|0),1,1.0)
 ffa_m!1="plev"
 ffa_m&plev=Ydata
 fase_m=ffa_m(1,:,0)
 zlev=  (16*log10(10^3/Ydata))
 max_ind_m=ind(ffa(0,:,0).eq.(max(ffa_m(0,:,0))))
 delete(ind_pl)
 amp_rel=ffa(0,max_ind:,0)/(ffa(0,max_ind,0))
 z_rel=zlev(max_ind:)-zlev(max_ind)

 ph=fase(max_ind:)-fase(max_ind)   ;phase max
 ph_m=fase_m(max_ind_m:)-fase_m(max_ind_m);phase min
 a=(/dimsizes(ph),dimsizes(ph_m)/)
 dime=min(a)
 do i=1,dime-1
       if ((ph(i).le.ph(i-1))).and.(abs(ph(i)-ph(i-1)).ge.abs(ph(i)-ph(i-1)+12)) then
               ph(i)=12+ph(i)
           end if

           if ((ph_m(i).le.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).ge.abs(ph_m(i)-ph_m(i-1)+12)) then
               ph_m(i)=12+ph_m(i)
           end if
           if ((ph_m(i).lt.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).gt.abs(ph_m(i)-ph_m(i-1)+12)) then
               ph_m(i)=12+ph_m(i)
           end if
           if ((ph_m(i).lt.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).gt.abs(ph_m(i)-ph_m(i-1)+12)) then
               ph_m(i)=12+ph_m(i)
           end if

 end do

 phase=(ph(:dime-1)+ph_m(:dime-1))/2


 phase_XX=new((/dimTot,dimsizes(phase)/),typeof(newXXData))
 amp_rel_XX=new((/dimTot,dimsizes(amp_rel)/),typeof(newXXData))
 
 delete(ffa)
 delete(fase)
 delete(max_ind)
 delete(ffa_m)
 delete(fase_m)
 delete(zlev)
 delete(max_ind_m)
 delete(ph)
 delete(ph_m)
 delete(a)
 delete(phase)

;FOURIER ANALYSIS LOOP on MODEL
do imod=0,dimTot-1
 ffa=fourier_info(newXXData(line|imod,plev|:,month|0:11,stat|0),1,1.0)
 ffa!1="plev"
 ffa&plev=Ydata
 ind_pl=ind((ffa&plev.gt.96).and.(ffa&plev.le.100)) 
 fase=ffa(1,:,0)
 zlevel =(16*log10(10^3/Ydata))
 max_ind=ind(ffa(0,:,0).eq.(max(ffa(0,ind_pl,0))))
 printVarSummary(newXXData)
 ffa_m=fourier_info(-newXXData(line|imod,plev|:,month|0:11,stat|0),1,1.0)
 ffa_m!1="plev"
 ffa_m&plev=Ydata
 fase_m=ffa_m(1,:,0)
 zlev=  (16*log10(10^3/Ydata))
 max_ind_m=ind(ffa(0,:,0).eq.(max(ffa_m(0,:,0))))
 delete(ind_pl)
 amp_rel=ffa(0,max_ind:,0)/(ffa(0,max_ind,0))
 z_rel=zlev(max_ind:)-zlev(max_ind)
 ph=fase(max_ind:)-fase(max_ind)   ;phase max
 ph_m=fase_m(max_ind_m:)-fase_m(max_ind_m);phase min
 a=(/dimsizes(ph),dimsizes(ph_m)/)
 dime=min(a)
 
 do i=1,dime-1
           if ((ph(i).le.ph(i-1))).and.(abs(ph(i)-ph(i-1)).ge.abs(ph(i)-ph(i-1)+12)) then
               ph(i)=12+ph(i)
           end if

           if ((ph_m(i).le.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).ge.abs(ph_m(i)-ph_m(i-1)+12)) then
               ph_m(i)=12+ph_m(i)
           end if
           if ((ph_m(i).lt.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).gt.abs(ph_m(i)-ph_m(i-1)+12)) then
               ph_m(i)=12+ph_m(i)
           end if
           if ((ph_m(i).lt.ph_m(i-1))).and.(abs(ph_m(i)-ph_m(i-1)).gt.abs(ph_m(i)-ph_m(i-1)+12)) then
               ph_m(i)=12+ph_m(i)
           end if

  end do

  phase=(ph(:dime-1)+ph_m(:dime-1))/2
  phase_XX(imod,:)=phase+12  ; it seeems the only way to have the sdame results as ey06!!!
  amp_rel_XX(imod,:)=amp_rel
  delete(ffa)
  delete(fase)
  delete(max_ind)
  delete(ffa_m)
  delete(fase_m)
  delete(zlev)
  delete(max_ind_m)
  delete(ph)
  delete(ph_m)
  delete(a)
  delete(phase)

end do
;******************error bar****************************************************
;********************************************************************************
 
  plot_d=new(2,"graphic")
;*************************plotting*******************************************************
     res                      = True                  ; plot mods desired
     ;   res@tiMainString         = models@name(imod)     ; title
     res@gsnDraw              = False
     res@gsnFrame             = False
     res@trXMaxF           = aMax
     res@trXMinF           = aMin
    ;  res@trYMaxF           = aYMax(ipanel)
    ;  res@trYMinF           = aYMin(ipanel)
     res@xyMarkLineModes    = LineScatters

     res@xyLineThicknessF= 4.0
     res@xyLineDashSegLenF=0.13
     res@xyLineColors      = plotColor(:)                 ; line color
     res@xyDashPatterns    = plotDash(:)
    print(plotScatters) 
      res@xyMarkers         =plotScatters(:)
     res@xyMarkerSizeF     =0.01
     res@xyMarkerThicknessF =1.8

     res@trYMaxF              = 12
     res@trYMinF              = 0

     res@vpHeightF            = 0.7                    ; change aspect ratio of plot
     res@vpWidthF             = 0.4

     res@tiYAxisString        = "z-z0 [Km]"    ; y-axis label
     ;print(clima(0))
     ;print(plotMod)
  
      res@tiXAxisString        = "Relative Amplitude"
      plot_d(0) = gsn_csm_xy (wks,amp_rel_XX,z_rel,res)
      res@trXMaxF              = 18
      res@trXMinF              = -2
      res@trXMaxF           = bMax
      res@trXMinF           = bMin
      res@tiXAxisString        = "Phase Lag (months)"
      plot_d(1) =gsn_csm_xy (wks,phase_XX,z_rel,res)



;; Call panel routine
   Pres=True
   gsn_panel(wks,plot_d,(/1,2/),Pres) 
   plotDash(index_obs)=plotScatters(index_obs)
   add_legend(wks,plotMod,plotColor,tofloat(plotDash),LineScatters,plotThicks+2.,"0.2") 

   frame(wks)
;**********************common output file amplitude**************************
     nc_filename = work_dir + "eyring06f09_amplitude" + var0 + ".nc"
     amp_rel_XX@var = var0
     amp_rel_XX@diag_script="eyring06f09"
     amp_rel_XX!0="line"
     amp_rel_XX!1="z_rel"
     amp_rel_XX&line=plotMod
     amp_rel_XX&z_rel=z_rel 
     nc_outfile = ncdf_write(amp_rel_XX, nc_filename)
;**********************common output file phase**************************
     nc_filename1 = work_dir + "eyring06f09_phase" + var0 + ".nc"
     phase_XX@var = var0
     phase_XX@diag_script="eyring06f09"
     phase_XX!0="line"
     phase_XX!1="z_rel"
     phase_XX&line=plotMod
     phase_XX&z_rel=z_rel
     nc_outfile1 = ncdf_write(phase_XX, nc_filename1)
end

