;##############################################################################
;# ANNUAL CYCLE CONTOUR PLOTS OF REGIONAL MEANS                                #
;##############################################################################
;# Author: Irene Cionni(ENEA, Italy)                                          #
;# CRESCENDO project
;##############################################################################
;#                                                                           
;# Description:                                                               
;#    * Creates annual cycle regional plots for individual models and mumltimodle
;#      over the selected project&experiemnt&mip                                    
;#       - Read regions
;#       - select the regions
;#       - calculate annual cycle for each model
;#       - evaluate the MMM over selected project&experiment&mip
;#       - plot the individual models annual cycle 
;#       - plot the MMMs annual cycle and standard deviations
;#                                                                      
;# Required diag_script_info attributes (diagnostics specific):
;#    
;#     * @fig938_region_label:
;#       (/"WNA","ENA","CAM","TSA","SSA","EUM","NAF","CAF","SAF",
;#         "NAS","CAS","EAS","SAS","SEA","AUS"/)
;#        these reagions are defined in function select_region (./diag_scripts/lib/ncl/latlon.ncl
;#        according to a personal comunication of Sin Chan Chou       
;#     * @fig938_project_MMM: projects to average
;#     * @fig938_experiment_MMM: experiments to average
;#     * @fig938_mip_MMM: mip to average
;#
;#        @fig938_mip_MMM must be the same i.e.
;#        @fig938_project_MMM=(/"CMIP5", "CMIP3"/)
;#        @fig938_experiment_MMM=(/"historical", "historical"/)
;#        @fig938_mip_MMM=(/"Amon", "Amon"/)
;#
;# Optional diag_script_info attributes (diagnostics specific):
;#
;#        @fig938_names_MMM   :names in legend  i.e. (/"CMIP5","CMIP3"/) 
;#        @fig938_colors_MMM  : colors (/"red","blue"/)   
;#        @legend_position        : legend position in ncd coordinates     
;#        @fig938_YMin        : minimum Y Axis
;#        @fig938_YMax        : maximum Y Axis 
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;# original file anncycplot.ncl
;#                                                                            #
;#############################################################################
load "$diag_scripts/../interface_scripts/interface.ncl"

load "$diag_scripts/shared/statistics.ncl"
load "$diag_scripts/shared/plot/aux_plotting.ncl"
load "$diag_scripts/shared/plot/style.ncl"
load "/home/irene/ESMValTool/esmvaltool/diag_scripts/regional_downscaling/regional_function.ncl"
;load "$diag_scripts/shared/latlon.ncl"

begin
   enter_msg(DIAG_SCRIPT, "")
   var0 = variable_info[0]@short_name
   units =variable_info[0]@units
   info0 = select_metadata_by_name(input_file_info, var0)
   print(info0)
   dim_MOD1 = ListCount(info0)
   names = metadata_att_as_array(info0, "dataset")
   projects = metadata_att_as_array(info0, "project")
   exps = metadata_att_as_array(info0, "exp")
   ensembles = metadata_att_as_array(info0, "ensemble")
   mips=metadata_att_as_array(info0, "mip")
   y1 = metadata_att_as_array(info0, "start_year")
   y2 = metadata_att_as_array(info0, "end_year") 
   long_name= variable_info[0]@long_name
   if (isatt(variable_info[0], "reference_dataset")) then
     refname = variable_info[0]@reference_dataset
     do imod = 0, dimsizes(names) - 1
         if (names(imod) .eq. refname) then
            iref = imod
          end if
      end do       
   else
     error_msg("f", DIAG_SCRIPT, "", "no reference dataset  " + \
              "(variable_info[0]@reference_dataset) needs to be defined.")
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
   print(LONG_NAME)

   log_info("++++++++++++++++++++++++++++++++++++++++++")
   log_info(DIAG_SCRIPT + " (var: " + var0 + ")")
   log_info("++++++++++++++++++++++++++++++++++++++++++")
  
   ; Set default values for optional diag_script_info attributes

   set_default_att(diag_script_info, "fig938_region_label", "ANN")
   set_default_att(diag_script_info, "styleset", "CMIP5")
   set_default_att(diag_script_info, "fig938_MMM", True)
   set_default_att(diag_script_info, "fig938_project_MMM", "CMIP5")
   set_default_att(diag_script_info, "fig938_experiment_MMM", "historical")
   set_default_att(diag_script_info, "fig938_mip_MMM", "Amon")
   set_default_att(diag_script_info, "fig938_names_MMM", "CMIP5")
   set_default_att(diag_script_info, "fig938_colors_MMM", "red")
   set_default_att(diag_script_info, "legend_position","0.2")
   ;set_default_att(diag_script_info, "Figure_938_YMin", "0") 
   ;set_default_att(diag_script_info, "Figure_938_YMax", "0") 
    colors = project_style(info0, diag_script_info, "colors")
   dashes = project_style(info0, diag_script_info, "dashes")
   thicks = tofloat(project_style(info0, diag_script_info, "thicks"))
   markers = project_style(info0, diag_script_info, "markers") 
   lin_mar_mod=new((/dimsizes(colors)/),"string")
   lin_mar_mod=(/"Lines"/)
 
   if (isatt(diag_script_info, "refModel")) then
        idx_ref = ind(models@name .eq. diag_script_info@refModel)
        if (all(ismissing(idx_ref))) then
            error_msg("f", diag_script, "", "refModel is not included in " + \
                      "model list (namelist)")
        end if
        p_ref = "ref_"
    else
        idx_ref = -1
        p_ref = ""
    end if


   region_label = tostring(diag_script_info@fig938_region_label)
   print(region_label)
   dim_reg=dimsizes(region_label)
   latrange=new((/dimsizes(region_label),2/),"double")
   lonrange=new((/dimsizes(region_label),2/),"double")
   do ir=0,dimsizes(region_label)-1
       print(region_label(ir))
       ;  if (region_label(ir).eq."WNA") then
       ;    region_array = (/30., 67., 190., 250./)
       ;    region_array@name = "WNA"
           ;leave _msg(scriptname, funcname, 10)
        ;   print(region_array)
        ; end if 
       region=select_region1(region_label(ir))
       latrange(ir,:)=region(0:1)
       lonrange(ir,:)=region(2:3)
    end do 



   flag_mod = where(projects.ne."OBS" .and.projects.ne."OBS6".and. projects.ne."obs4mips", 1, 0)
   index_mod = ind(flag_mod.gt.0)
   index_obs = ind(flag_mod.eq.0)
   dim_MOD = dimsizes(index_mod)
   dim_OBS = 0

   if (.not.all(ismissing(index_obs))) then
     dim_OBS = dimsizes(index_obs)
     list_obs_mod = names(index_obs)
   end if


         
   ;   YMax = diag_script_info@fig938_YMax
   ;   XMin = diag_script_info@fig938_YMin
  ; year_start=toint(diag_script_info@fig938_start_year)
  ; year_end=toint(diag_script_info@fig38_end_year)

  ; make sure path for (mandatory) netcdf output exists

  work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)

end
begin
work_dir = config_user_info@work_dir + "/"
  ; Create work dir
  system("mkdir -p " + work_dir)
 wks = get_wks("dummy_for_wks", DIAG_SCRIPT,"Figure9.38_" + var0)
;******************************************************************
; common resources for  plots
    pan=True
    res                        = True
    res@gsnDraw                = False
    res@gsnFrame               = False
    res@vpHeightF= 0.6
    res@vpWidthF= 0.8   
    ;res@txFontHeightF          = 0.018
    months               = ispan(0, 12, 1)
    res@tiMainFontHeightF      = 0.022
    res@tmXBMode               = "Explicit"
    res@tmXBValues             = ispan(0, 12, 1)
    res@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
                       "O", "N", "D", "J"/)
    res@tmXBLabelFontHeightF = 0.035
    res@tmYLLabelFontHeightF = 0.035
    res@tiXAxisString=" "
    res@tiYAxisString=" "
    res@xyLineDashSegLenF=0.5
;;------------res1--------------------------------------------------

    pan1=True
    res1                        = True
    res1@gsnDraw                = False
    res1@gsnFrame               = False
    ;res1@txFontHeightF          = 0.02
    res1@tiMainFontHeightF      = 0.02
    res1@tmXBMode               = "Explicit"
    res1@tmXBValues             = ispan(0, 12, 1)
    res1@tmXBLabels = (/"J", "F", "M", "A", "M", "J", "J", "A", "S", \
                       "O", "N", "D", "J"/)
    res1@vpHeightF= 0.6
    res1@vpWidthF= 0.8   
    res1@tmXBLabelFontHeightF = 0.035
    res1@tmYLLabelFontHeightF = 0.035  
    res1@tiXAxisString=" "
    res1@tiYAxisString=" "
    res1@xyLineDashSegLenF=0.5

;************************select MMM*********************
    dim_MMM=1
    if(isatt(diag_script_info, "fig938_MMM")) then
      if (diag_script_info@fig938_MMM.eq."True")
        if ( (isatt(diag_script_info, "fig938_project_MMM")).and.(isatt(diag_script_info, "fig938_experiment_MMM")).and.(isatt(diag_script_info, "fig938_mip_MMM")))
           if (dimsizes(diag_script_info@fig938_project_MMM).eq.dimsizes(diag_script_info@fig938_experiment_MMM).eq.dimsizes(diag_script_info@fig938_mip_MMM)) then
             dim_MMM=dimsizes(diag_script_info@fig938_project_MMM)
               indic_s=new((/dimsizes(diag_script_info@fig938_project_MMM),dimsizes(names)/),"integer")
               indic_not_s=new((/dimsizes(diag_script_info@fig938_project_MMM),dimsizes(names)/),"integer")
               do ik=0,dimsizes(diag_script_info@fig938_project_MMM)-1
                 ii_s=ind((projects.eq.diag_script_info@fig938_project_MMM(ik)).and.(exps.eq.diag_script_info@fig938_experiment_MMM(ik)).and.(mips.eq.diag_script_info@fig938_mip_MMM(ik)))

                 indic_s(ik,:dimsizes(ii_s)-1)=ii_s
                 delete(ii_s)

                 ii_not_s=ind((projects.ne.diag_script_info@fig938_project_MMM(ik)).or.(exps.ne.diag_script_info@fig938_experiment_MMM(ik)).or.(mips.ne.diag_script_info@fig938_mip_MMM(ik)))
                 indic_not_s(ik,:dimsizes(ii_not_s)-1)=ii_not_s
                 delete(ii_not_s)
              end do
           else
                  error_msg("f", diag_script, "", "diag_script_info@fig938_project_MMM, diag_script_info@fig938_experiment_MMM and diag_script_info@fig938_experiment_MMM must have the same dimension")
           end if
        else
           error_msg("f", diag_script, "", "diag_script_info@fig938_project_MMM, diag_script_info@fig938_experiment_MMM and diag_script_info@fig938_experiment_MMM are not specified")
        end if
       end if
    end if
print(dimsizes(diag_script_info@fig938_project_MMM))
if (dimsizes(diag_script_info@fig938_project_MMM).gt.1) then
   indic_not_sel=get_unique_difference(indic_s,indic_not_s)
else
  ndim=ndtooned(indic_not_s) 
  print(indic_not_s)
  print(ndim)
  ise=ind(.not.ismissing(ndim))
  print(ise)
  indic_not_sel=ndim(ise)
  delete(ise)
  delete(ndim) 
end if
print(indic_not_sel)
;******************************************landsea mask**************************
   a    = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")
   lsdata = a->LSMASK
;********************************************************************************
    plot = new (dim_reg,"graphic")
    plot_MMM = new (dim_reg,"graphic")
    plot_STD = new ((/dim_reg,dim_MMM/),"graphic")   
    do ir=0,dim_reg-1 
       case1=region_label(ir)
        print(case1) 
       if (idx_ref .ne. -1) then
          if(isatt(diag_script_info, "fig938_MMM")) then
             if (diag_script_info@fig938_MMM.eq."True")
                indic_s=where(indic_s.eq.idx_ref,indic_s@_FillValue,indic_s)
                indic_not_sel=where(indic_not_sel.eq.idx_ref,indic_not_sel@_FillValue,indic_not_sel) 
                indic_not_sel1=indic_not_sel(ind(.not.ismissing(indic_not_sel)))
                delete(indic_not_sel)
                indic_not_sel=indic_not_sel1
                delete(indic_not_sel1)
             end if
         end if
         pan@txString =variable_info@long_name+" bias vs. "+models@name(idx_ref)
         D=new((/dim_MOD,13/),"double")
         C=new((/dim_MOD-1,13/),"double")
         models1= new((/dim_MOD-1/),"string")
         colors1 = new((/dim_MOD-1/),typeof(colors)) 
         dashes1 = new((/dim_MOD-1/),typeof(dashes))
         thicks1 = new((/dim_MOD-1/),typeof(thicks))
         markers1 = new((/dim_MOD-1/),typeof(markers))   
         lin_mar_mod1=new((/dim_MOD-1/),"string")
      ;; (b) difference of each model from 'reference'
  ;; CASE 1
      ;case1 = models@name(idx_ref) + " " + models@case_name(idx_ref)
        A0 = read_data(info0[idx_ref])
        ;uniform_coordinates(A0, field_type0)
        printVarSummary(A0)        
        ;; convert to T2Ms, C2Ms
        ;if ((field_type0 .eq. "T2Ds") .or. (field_type0 .eq. "C2Ds"))then
        ;    tmp = day_to_month(A0, field_type0)
        ;    delete(A0)
        ;    A0 = tmp
        ;    delete(tmp)
        ;end if

        ;; Select a longitude sub-set to plot
        ;if (iscoord(A0, "lon") .and. isvar("lonrange") \
        ;    .and. inlist("lonrange", supported_info_transforms)) then
            ;Atmp = A0(time|:, {lat|latrange(ir,0):latrange(ir,1)}, {lon|lonrange(ir,0):lonrange(ir,1)})
             ;Atmp1=calc_cycle_lon(A0)
             Atmp1=A0
 
             lsm  = landsea_mask(lsdata,Atmp1&lat,Atmp&lon)
             Atmp1 = mask(Atmp1,lsm.eq.0,False)
             Atmp1 = mask(Atmp1,lsm.eq.2,False)
             Atmp1 = mask(Atmp1,lsm.eq.4,False)
             delete(lsm)   
            Atmp=area_operations(Atmp1,latrange(ir,0),latrange(ir,1),lonrange(ir,0),lonrange(ir,1),"average",True)
            delete(A0)
            A = Atmp
            delete(Atmp)
            delete(Atmp1) 

        jmod = 0
 do imod = 0, dim_MOD - 1
            if (imod .ne. idx_ref) then
                ;; CASE 2
        ;        case2 = models@name(imod) + " " + models@case_name(imod)
                B0 = read_data(info0[imod])
                ;uniform_coordinates(B0, field_type0)

                ;; convert to T2Ms, C2Ms
                if ((field_type0 .eq. "T2Ds") .or. \
                    (field_type0 .eq. "C2Ds")) then
                    tmp = day_to_month(B0, field_type0)
                    delete(B0)
                    B0 = tmp
                    delete(tmp)
                end if

                   Btmp1=calc_cycle_lon(B0) 
                   lsm  = landsea_mask(lsdata,Btmp1&lat,Btmp1&lon)

                   Btmp1 = mask(Btmp1,lsm.eq.0,False)
                   Btmp1 = mask(Btmp1,lsm.eq.2,False)
                   Btmp1 = mask(Btmp1,lsm.eq.4,False)
                   delete(lsm)
                   Btmp =area_operations(Btmp1,latrange(ir,0),latrange(ir,1),lonrange(ir,0),lonrange(ir,1),"average",True)
                   delete(B0)
                   B = Btmp
                   delete(Btmp)
                   delete(Btmp1)

                ;;---------------------------------------------------------
                ;; difference plot of time(x) vs. latitude(y)

                    ;C = A  ; copy dims, coords
                     C(jmod,:) = (/A - B/)  ; get diff values
                     D(imod,:)=(/A - B/) 
                  if (variable_info@units.eq."K") then
                       UNITS="~S~o~N~C"
                    else
                       UNITS=variable_info@units
                    end if

                delete(B)
                colors1(jmod)=colors(imod) 
                dashes1(jmod)=dashes(imod)
                thicks1(jmod)=thicks(imod) 
                markers1(jmod)=markers(imod)
                models1(jmod)=models@name(imod) 
                lin_mar_mod1(jmod)="Lines"   
             jmod=jmod+1
            end if
        end do  ; imod
        delete(A)
else
       lin_mar_mod1=lin_mar_mod
       C=new((/dim_MOD,13/),"double") 
       do imod = 0, dim_MOD - 1
            ;; CASE 1
            pan@txString=variable_info@long_name
            ;case1 = models@name(imod) + " " + models@case_name(imod)
            ;if ((stringtoint(models@end_year(imod)) \
            ;    - stringtoint(models@start_year(imod))) .gt. 0) then
            ;    yrs_ave1 = models@start_year(imod) + "-" + \
            ;        models@end_year(imod)
            ;end if
            A0 = read_data(info0[imod])
            ;uniform_coordinates(A0, field_type0)

            ;; convert to T2Ms, C2Ms
           ; if ((field_type0 .eq. "T2Ds") .or. (field_type0 .eq. "C2Ds"))then
           ;     tmp = day_to_month(A0, field_type0)
           ;     delete(A0)
           ;     A0 = tmp
           ;     delete(tmp)
           ; end if
           ; Atmp1 = calc_cycle_lon(A0)
            Atmp1 =A0
            lsm  = landsea_mask(lsdata,Atmp1&lat,Atmp1&lon)   
            Atmp1 = mask(Atmp1,lsm.eq.0,False)
            Atmp1 = mask(Atmp1,lsm.eq.2,False)
            Atmp1 = mask(Atmp1,lsm.eq.4,False)
            delete(lsm)
            Atmp=area_operations(Atmp1,latrange(ir,0),latrange(ir,1),lonrange(ir,0),lonrange(ir,1),"average",True)
            delete(A0)
            A = Atmp
            delete(Atmp)
            delete(Atmp1)
            if (units.eq."K") then
                UNITS="~S~o~N~C"
                A=(/A-273.15/)
            else
                UNITS=units
            end if
            printVarSummary(A)
            printVarSummary(C)
                        print(imod)    
            C(imod,0:11)=A
            C(imod,12)=A(0)
            lin_mar_mod1(imod)="Lines"
         end do ;(imod)
         
         colors1=colors
         dashes1=dashes
         thicks1=thicks
         markers1=markers
         models1i=names
      end if
      ;; case contour plots of time(x) vs. latitude(y)
      

      ;res@tiMainString =variable_info@long_name 
      
      min1 = min(C)
      max1 = max(C)
      if(isatt(diag_script_info, "fig938_YMin")) then
         res@trYMinF = diag_script_info@Figure_938_YMin
      else
         res@trYMinF = min1
      end if

      if(isatt(diag_script_info, "fig938_YMax")) then
        res@trYMaxF = diag_script_info@Figure_938_YMax
      else
        res@trYMaxF = max1
      end if

      res@xyLineColors      = colors1  ; line colors
      res@xyLineThicknesses = thicks1        ; line thicknesses
      res@xyDashPatterns    = dashes1       ; line patterns
      res@gsnLeftString = case1
      res@gsnLeftStringFontHeightF=0.04
      res@gsnRightString = UNITS
      res@gsnRightStringFontHeightF=0.04
      res@tfPolyDrawOrder = "Draw"
      plot(ir)  = gsn_csm_xy (wks,months,C,res)
      if(isatt(diag_script_info, "fig938_MMM")) then
          if (diag_script_info@fig938_MMM.eq."True")
              MMM=new((/dimsizes(diag_script_info@fig938_project_MMM)+dimsizes(indic_not_sel),13/),typeof(C))
              STD=new((/dimsizes(diag_script_info@fig938_project_MMM)+dimsizes(indic_not_sel),13/),typeof(C))
              names_MMM=new((/dimsizes(diag_script_info@fig938_project_MMM)+dimsizes(indic_not_sel)/),"string")
              colors_MMM=new((/dimsizes(diag_script_info@fig938_project_MMM)+dimsizes(indic_not_sel)/),typeof(colors))
              thicks_MMM=new((/dimsizes(diag_script_info@fig938_project_MMM)+dimsizes(indic_not_sel)/),typeof(thicks))
              dashes_MMM=new((/dimsizes(diag_script_info@fig938_project_MMM)+dimsizes(indic_not_sel)/),typeof(dashes))
              lin_mar_mod_MMM=new((/dimsizes(diag_script_info@fig938_project_MMM)+dimsizes(indic_not_sel)/),typeof(lin_mar_mod1))
              colori_def=(/"red","blue","green","magenta","yellow"/) 
              do iin=0,dimsizes(diag_script_info@fig938_project_MMM)-1
                 aa=indic_s(iin,ind(.not.ismissing(indic_s(iin,:))))
                 if (idx_ref .ne. -1)
                   D!0="case"
                   MMM(iin,:)=dim_avg_n_Wrap(D(aa,:),0) 
                   STD(iin,:)=dim_stddev_n_Wrap(D(aa,:),0) 
                 else
                   C!0="case"
                   MMM(iin,:)=dim_avg_n_Wrap(C(aa,:),0)
                   STD(iin,:)=dim_stddev_n_Wrap(C(aa,:),0)
                 end if
                 delete(aa)
                 if(isatt(diag_script_info, "fig938_colors_MMM")) then
                     colors_MMM(iin)=diag_script_info@fig938_colors_MMM(iin)
                 else
                     colors_MMM(iin)=colori_def(iin) 
                 end if 
                 thicks_MMM(iin)=3.
                 dashes_MMM(iin)=0
                 lin_mar_mod_MMM(iin)="Lines"
                 if(isatt(diag_script_info, "fig938_names_MMM")) then
                     names_MMM(iin)=diag_script_info@fig938_names_MMM(iin)
                 else
                     names_MMM(iin)=diag_script_info@fig938_project_MMM(iin)+" "+diag_script_info@fig938_experiment_MMM(iin)+" "+diag_script_info@fig938_mip_MMM(iin)
                 end if 
              end do
              if (idx_ref .ne. -1)
                MMM(dimsizes(diag_script_info@fig938_project_MMM):,:)=(/D(indic_not_sel,:)/)
                 colors_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=colors(indic_not_sel)
                 thicks_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=3.
                 dashes_MMM(dimsizes(diag_script_info@fig38_project_MMM):)=dashes(indic_not_sel)
                 lin_mar_mod_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=lin_mar_mod(indic_not_sel)
                 names_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=modelli(indic_not_sel)
              else
                 printVarSummary(MMM(dimsizes(diag_script_info@fig938_project_MMM):,:))
                 printVarSummary(C) 
                 print(indic_not_sel)
                 printVarSummary(C(indic_not_sel,:))
                 MMM(dimsizes(diag_script_info@fig938_project_MMM):,:)=(/C(indic_not_sel,:)/)                 
                 colors_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=colors1(indic_not_sel)
                 thicks_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=3. 
                 dashes_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=dashes1(indic_not_sel)
                 lin_mar_mod_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=lin_mar_mod1(indic_not_sel) 
                 names_MMM(dimsizes(diag_script_info@fig938_project_MMM):)=models1(indic_not_sel)
             end if
          end if
      end if

      if(isatt(diag_script_info, "fig938_YMin")) then
        res1@trYMinF = diag_script_info@fig938_YMin
    else
        res1@trYMinF = min1
    end if

      if(isatt(diag_script_info, "fig938_YMax")) then
        res1@trYMaxF = diag_script_info@fig938_YMax
    else
        res1@trYMaxF = max1
    end if

      res1@xyLineColors      = colors_MMM  ; line colors
      res1@xyLineThicknesses = thicks_MMM        ; line thicknesses
      res1@xyDashPatterns    = dashes_MMM       ; line patterns
      res1@gsnLeftString = case1
      res1@gsnRightString = UNITS
      res1@gsnLeftStringFontHeightF=0.04
      res1@gsnRightStringFontHeightF=0.04
      res1@tfPolyDrawOrder = "Draw"

      plot_MMM(ir)  = gsn_csm_xy (wks,months,MMM,res1)
      xp    = new( (/2*13/), typeof(MMM) )
      yp    = new( (/2*13,dimsizes(diag_script_info@fig938_project_MMM)/), typeof(MMM) )
       do k=0,12
         yp(k,:) = MMM(:dimsizes(diag_script_info@fig938_project_MMM)-1,k) + STD(:dimsizes(diag_script_info@fig938_project_MMM)-1,k)
         xp(k) =months(k) 
         xp(2*13-1-k) = months(k)
         yp(2*13-1-k,:) = MMM(:dimsizes(diag_script_info@fig938_project_MMM)-1,k) - STD(:dimsizes(diag_script_info@fig938_project_MMM)-1,k)
      end do

 
      do iin=0, dimsizes(diag_script_info@fig938_project_MMM)-1
         gsres                   = True                        ; poly res
         gsres@gsFillColor       = colors_MMM(iin) 
         gsres@gsFillOpacityF =0.3
       
         plot_STD(ir,iin)  = gsn_add_polygon (wks,plot_MMM(ir),xp,yp(:,iin),gsres)
      end do

      delete(min1)
      delete(max1)
      delete(C)
      delete(MMM)
      delete(STD)
      delete(xp)
      delete(yp)
    end do  ;region 


    ;; Panel plots
    pan = True
    pan@gsnMaximize = True
    pan@gsnFrame = False
    pan@gsnPaperOrientation = "portrait"
    n_p = dimsizes(plot) /5 
    if (mod(dimsizes(plot), 5) .ne. 0) then
        n_p = n_p + 1
    end if
    gsn_panel(wks, plot, (/n_p, 5/), pan)
    frame(wks)
   ;    add_legenda(wks, models1, colors1, dashes1,lin_mar_mod1, thicks1, diag_script_info@legend_position)
  ;  if(isatt(diag_script_info, "legend_position")) then
   ; else
        print(models1)
        print(colors1)
        print(dashes1)
        print(lin_mar_mod1)
       print(thicks1)
        add_legenda(wks, models1, colors1, dashes1,lin_mar_mod1, thicks1, "0.8") 
    ;end if
    frame(wks)
    if(isatt(diag_script_info, "Figure_938_MMM")) then
       if (diag_script_info@Figure_938_MMM.eq."True")
          pan1 = True
          pan1@gsnMaximize = True
          pan1@gsnFrame = False
          pan1@gsnPaperOrientation = "portrait"
          n_p = dimsizes(plot_MMM) / 5
          if (mod(dimsizes(plot_MMM), 5) .ne. 0) then
             n_p = n_p + 1
          end if
          gsn_panel(wks, plot_MMM, (/n_p, 5/), pan1)
          if(isatt(diag_script_info, "legend_position")) then
             add_legenda(wks, names_MMM, colors_MMM, dashes_MMM,lin_mar_mod_MMM, thicks_MMM, diag_script_info@legend_position)
          else
             add_legenda(wks, names_MMM, colors_MMM, dashes_MMM,lin_mar_mod_MMM, thicks_MMM, "0.3")
          end if
       end if
    end if
    frame(wks)
    delete(plot_MMM)
    delete(res1)
    delete(pan1)


;**********************common output file**************************
     nc_filename = work_dir + "fig938_" + var0 + ".nc"
     Xdata@var = var0
     Xdata@diag_script="fig938"

     printVarSummary(Xdata)
     Xdata&stat=(/"mean","stddev"/)
     Xdata&season=seasons
     Xdata&line=plotMod
     printVarSummary(Xdata) 
     nc_outfile = ncdf_write(Xdata, nc_filename)


end

