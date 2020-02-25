;EYRING06JGR_FIG14
;; Author: Irene Cionni (ENEA, Italy), Veronika Eyring (DLR, Germany)
;; PROJECT-NAME CCMI
;;#############################################################################
;; Description
;;   Latitudinal distribution of the total column of ozone climatological mean annual cycle ,
;;
;;
;; Optional diag_script_info attributes (diagnostic specific)
;;    @e06fig14_Max:   max limit contour interval
;;    @e06fig14_Min:   min limit contour interval
;;    @e06fig14_Spacing:  contour interval
;;    @e06fig14_multimean:  multimodel mean 
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
   set_default_att(diag_script_info, "e06fig14_Min", 185.)
   set_default_att(diag_script_info, "e06fig14_Max", 500.)
   set_default_att(diag_script_info, "e06fig14_Spacing", 25.)
   set_default_att(diag_script_info, "e06fig14_multimean", False)
   

   flag_mod = where(projects.ne."OBS" .and. projects.ne."obs4mips", 1, 0)
   index_mod = ind(flag_mod.gt.0)
   index_obs = ind(flag_mod.eq.0)
   dim_MOD = dimsizes(index_mod)
   dim_OBS = 0

   if (.not.all(ismissing(index_obs))) then
      dim_OBS = dimsizes(index_obs)
      list_obs_mod = names(index_obs)
   end if

         
   Max = diag_script_info@e06fig14_Max
   Min = diag_script_info@e06fig14_Min
   Spacing= diag_script_info@e06fig14_Spacing
   multimean = diag_script_info@e06fig14_multimean
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
  wks = get_wks("dummy_for_wks", DIAG_SCRIPT, "eyring06f14_" + var0)
  output_filename = wks@fullname
  A0 = read_data(info0[iref])
  printVarSummary(A0)
  lat0=  A0&lat
  delete(A0)

   
  tmpData = new((/dim_MOD+dim_OBS,12,dimsizes(lat0),2/), "double")
  tmpData!0 = "line"
  tmpData!1 = "month"
  tmpData!2 = "lat"
  tmpData!3 = "stat"

  XXdata = new((/dimTot,12,dimsizes(lat0),2/), "double")
  XXdata!0 = "line"
  XXdata!1= "month"
  XXdata!2 = "lat"
  XXdata!3 = "stat"

  XXref= new((/ 12,dimsizes(lat0),2/), "double")
  XXref!0= "month"
  XXref!1 = "lat"
  XXdata!2 = "stat"

  years1=new(dimTot,"integer")
  years2=new(dimTot,"integer")
  plot = new(dimTot, "graphic")
  plot_num = integertochar(ispan(97, 122, 1))

  plotMod=new(dimTot,typeof(names))
  
  ik=0


;******************************************************************
  do imod = 0, dim_MOD1 - 1
    years1(imod)=y1(imod)
    years2(imod)=y2(imod)
    print(imod)
    plotMod(imod)=names(imod)
      ; Extract data
      A0 = read_data(info0[imod])
      dummy1 = time_operations(A0,y1(imod),y2(imod),"average","monthlyclim",True)
      dummy2_b=tofloat(A0)
      date= yyyymm_time(y1(imod),y2(imod),"integer") 
      yyyy  = date/100 
      mm  = date-yyyy*100                ; mmdd = yyyymmdd%10000
      weights = days_in_month(yyyy, mm)
      dummy2 = dummy2_b(0:11,: )  ; Copy metadata
      do ilat=0,dimsizes(lat0)-1
         do im = 0, 11
             dummy2(im,ilat ) = \
               dim_stddev_wgt_Wrap(dummy2_b(im::12,ilat), weights(im::12), 1)
         end do
      end do
      delete(A0)
      delete(weights)
      delete(date)
      delete(yyyy)
      delete(mm)
      delete(dummy2_b)
      dummy2!0=dummy1!0 
       
      tmpData(imod,:,:,0) = (/dummy1/)
      tmpData(imod,:,:,1) = (/dummy2/)
 
      delete(dummy1)
      delete(dummy2)
   end do
   XXdata(:dim_MOD1-1,:,:,:)=(/tmpData/)

   XXref= (/tmpData(iref,:,:,:)/)
         
   XXdata!0 = "line"
   XXdata!1 = "month"
   XXdata!2 = "lat"
   XXdata!3 = "stat"
   XXdata&lat=lat0
   XXdata&month=ispan(1,12,1)

  ;°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°
   if (multimean) then
             XXdata(dimTot-1,:,:,0)=(/dim_avg_n(XXdata(index_mod,:,:,0),0)/)
             plotMod(dimTot-1)="MMM"
             years1(dimTot-1)=max(y1)
             years2(dimTot-1)=min(y2)
   end if

   X_MMM_std=(/dim_stddev_n(XXdata(index_mod,:,:,:),0)/)
   plot_d=new(dimTot,"graphic")
;*************************plotting*******************************************************

      cmap = (/"White","Gray1","Gray75","maroon4","DarkOrchid4","RoyalBlue4", "RoyalBlue","CadetBlue3", "CadetBlue1","DarkOliveGreen3", "Chartreuse3", \
           "GreenYellow", "Yellow","Orange","OrangeRed","OrangeRed4",\
           "Black"/)
   gsn_define_colormap(wks,cmap)
   do ipanel=0, dimTot-1

      res                      = True                  ; plot mods desired
      res@gsnDraw              = False
      res@gsnFrame             = False
      res@lbLabelBarOn         = False

      res@tiMainString         = plotMod(ipanel)+" "+years1(ipanel)+"-"+years2(ipanel)     ; title
      res@tiYAxisFont          = "helvetica"
      res@tiMainFont           = "helvetica"
      res@tmXBLabelFont        = "helvetica"
      res@tmYLLabelFont        = "helvetica"
      res@lbLabelFont          = "helvetica"
      res@cnLevelSelectionMode = "ManualLevels"        ; manual contour levels
      res@cnLevelSpacingF      = Spacing                    ; contour interval
      res@cnMinLevelValF       = Min                   ; min level
      res@cnMaxLevelValF       = Max                 ; max level
      res@cnLineLabelsOn       = False                 ; turn off line labels
      res@cnFillOn             = True                  ; turn on color fill
      res@cnInfoLabelOn        = False                 ; turn off contour label
      res@gsnSpreadColors      = True                  ; use full range of colors
      res@lbLabelAutoStride    = True                  ; optimal labels
      res@gsnLeftString        = ""
      res@gsnRightString       = ""
      res@gsnLeftStringFontHeightF = 0.015
      res@tmXBMode             = "Explicit"            ; explicit labels
      ; location of labels
      res@gsnAddCyclic=True
      res@tmXBValues           = (/0,1,2,3,4,5,6,7,8,9,10,11,12/)
      res@tmXBLabels           = (/"J","F","M","A","M","J","J","A","S","O","N","D","J"/)    ; labels themselves
      res@tmYLValues           = (/-90,-60,-30,0,30,60,90/)
      res@tmYLLabels           = (/"-90","-60","-30","0","30","60","90"/)
      res@vpHeightF            = 0.4                   ; change aspect ratio of plot
      res@vpWidthF             = 0.75
      res@lbLabelFontHeightF   = 0.011
      res@tiYAxisString        = "Latitude"            ; y-axis label
      plot_d(ipanel) = gsn_csm_lat_time (wks,XXdata(line|ipanel,lat|:,month|:,stat|0),res)
      delete(res)


   end do
   ;; Call legend-plotting routine
   Pres = True
   Pres@gsnFrame       = False
   Pres@gsnPanelBottom = 0.2
   Pres@amJust="TopLeft"
;; Call panel routine
   if (dimTot.gt.3) then
      gsn_panel(wks,plot_d,(/1+dimTot/3,3/),Pres)
   else

      gsn_panel(wks,plot_d,(/1,3/),Pres)
   end if
   frame(wks)
 ;*******************OUTPUT FILE**************************************************
    month_names=(/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/) 
    nc_filename = work_dir + "eyring06f14_" + var0 + ".nc"
    Xdata=XXdata(line|:,month|:,lat|:,stat|:)
    Xdata@var = var0
    Xdata@diag_script="eyring06f14"
    Xdata&stat=(/"mean","avg"/)
    Xdata&line=plotMod
    delete(Xdata&month)
    Xdata&month=month_names
    Xdata&lat=lat0
    nc_outfile = ncdf_write(Xdata, nc_filename)

end

