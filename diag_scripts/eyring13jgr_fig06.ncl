;##############################################################################
;# eyring13jgr_fig06.ncl                                                      #
;##############################################################################
;# Author: Irene Cionni (ENEA, Italy)                                         #
;# ESMVal project PI: Veronika Eyring (DLR, Germany)                          #
;##############################################################################
;#                                                                            #
;# Description:                                                               #
;#     Calculates 1980 baseline-adjusted stratospheric column ozone time      #
;#     series from 1960-2100.                                                 #
;#                                                                            #
;# Required diag_script_info attributes (diagnostics specific):               #
;#     e13fig06_latrange: min. and max. latitude of the regions (n-element    #
;#         array of 2-element pairs, e.g., (/(/-90,90/), (/-90,-60/)/));      #
;#         one pair of latitudes is required for each season (see below)      #
;#     styleset: style set ("DEFAULT", "CMIP5")                               #
;#     e13fig06_season: seasons (n-element array of strings, "ANN", "JAN",    #
;#         "FEB", "MAR", "DJF", "SON", etc.)                                  #
;#     e13fig06_baseline_adj: do baseline adjustment (string: "True", "False")#
;#     e13fig06_baseline: year for baseline adjustment (e.g., 1980)           #
;#     e13fig06_mod_plot: "MMT" = plot of the MultiModel mean of each         #
;#         scenario and selection "list_chem_mod" and "list_nochem_mod";      #
;#         "IMT" = plot of each single model trend;                           #
;#         "RAW" = plot of each model as raw data                             #
;#     e13fig06_mod_plot_CI: plot confidence interval (string: "True",        #
;#         "False"); for e13fig06_mod_plot = "MMT" only!                      #
;#     e13fig06_mod_plot_PI: plot prediction interval (string: "True",        #
;#         "False"); for e13fig06_mod_plot = "MMT" only!                      #
;#     e13fig06_XMin: min. x-values (start years) for plotting (n-element     #
;#         array, e.g., (/1960., 1960./)); array is required to have the same #
;#         number of elements as "seasons" (see above)                        #
;#     e13fig06_XMax: max. x-values (end years) for plotting (n-element array,#
;#         e.g., (/2010., 2010./)); array is required to have the same number #
;#         of elements as "seasons" (see above)                               #
;#     e13fig06_YMin: min. y-values for plotting (n-element array, e.g.,      #
;#         (/260., 150./)); array is required to have the same number of      #
;#         elements as "seasons" (see above)                                  #
;#     e13fig06_YMax: max. y-values for plotting (n-element array, e.g.,      #
;#         (/330., 500./)); array is required to have the same number of      #
;#         elements as "seasons" (see above)                                  #
;#                                                                            #
;# Optional diag_script_info attributes (diagnostics specific):               #
;#     list_chem_mod: list of models in the group “chem” (array of strings,   #
;#         default = (/”All”))                                                #
;#     list_nochem_mod: list of models in the group “nochem” (array of        #
;#         strings, default = (/"None"/))                                     #
;#     e13fig06_labels_exp_esp: specify experiment name (string: "True",      #
;#         "False"); only if e13fig06_mod_plot = “IMT” or “RAW”!              #
;#                                                                            #
;# Caveats:                                                                   #
;#                                                                            #
;# Modification history:                                                      #
;#    * 20151208-A-laue_ax: code clean-up, moved call to write_references to  #
;#                          beginning of code                                 #
;#    * 20140528-A_GoKl: written based on ancyc_lat.ncl                       #
;#                                                                            #
;##############################################################################
;# The TSAM method is described in: Scinocca, J. F., D. B. Stephenson,        #
;# T. C. Bailey, and J. Austin (2010), Estimates of past and future ozone     #
;# trends frommultimodel simulations using a flexible smoothing spline        #
;# methodology, J. Geophys. Res. Atmos., 115.                                 #
;##############################################################################

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

load "./diag_scripts/lib/ncl/misc_function.ncl"
load "./interface_data/ncl.interface"
load "./interface_scripts/auxiliary.ncl"
load "./interface_scripts/data_handling.ncl"
load "./interface_scripts/messaging.ncl"
load "./diag_scripts/lib/ncl/set_operators.ncl"
load "./diag_scripts/lib/ncl/statistics.ncl"
load "./diag_scripts/lib/ncl/latlon.ncl"
load "./diag_scripts/lib/ncl/regridding.ncl"
load "./diag_scripts/lib/ncl/style.ncl"
load "./diag_scripts/lib/ncl/ensemble.ncl"

;************************************************
begin
    verbosity  = stringtointeger(getenv("ESMValTool_verbosity"))
    info_output("<<<<<<<< Entering " + diag_script, verbosity, 4)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    info_output("plot 06 - " + diag_script + ".ncl (var: " + variables(0) + ")", verbosity, 1)
    info_output("++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
    dim_MOD = dimsizes(models@name)
    var0 = variables(0)
    field_type0 = field_types(0)
    colors=project_style(diag_script_info, "colors")
    dashes=project_style(diag_script_info, "dashes")
    thicks=tofloat(project_style(diag_script_info, "thicks"))

    do idx_mod = 0, dim_MOD - 1
        data0 = read_data(idx_mod,var0, field_type0)
    end do
end

begin
    vardeffile = "interface_data/" + variables + "_info.tmp"
    loadscript(vardeffile)
end

begin
    ;references Information
    write_references(diag_script,      \  ; script name
                     "A_cion_ir",      \  ; authors
                     "",               \  ; contributors
                     "D_eyring13jgr",  \  ; diag_references
                     "",               \  ; obs_references
                     "P_esmval")          ; proj_references

    plot_fields = (/"T3M", "T2Ms"/)
    ind_f = inlist(field_type0, plot_fields)
    if (all(ismissing(ind_f))) then
        print("fatal:can't use field, " + field_type0 + " for diag_script " + diag_script)
        status_exit(1)
    end if
    delete(plot_fields)
    delete(ind_f)

    plot_dir = getenv("ESMValTool_plot_dir")
    diag_script_base = basename(diag_script)
    output_dir = get_output_dir(plot_dir, diag_script_base)
    file_type = getenv("ESMValTool_output_file_type")
    if(ismissing(file_type)) then
        file_type = "ps"
    end if
    file_type = str_lower(file_type)
    color_type = "COLOR"
    if(isatt(diag_script_info, "e13fig06_mod_plot")) then
      mod_plot = diag_script_info@e13fig06_mod_plot
    else
      mod_plot="MMT"
    end if

    wks  = gsn_open_wks (file_type,output_dir+diag_script_base+"_"+mod_plot)
    plot = new(dimsizes(diag_script_info@e13fig06_season),"graphic")
    plot_num = integertochar(ispan(97,122,1))

;-----------------------------------------------------------------
    dim_EXP = 1 ;dimsizes of different experiments
    label_exp = models@experiment(0)
    do n = 1,dimsizes(models@experiment)-1
      if (models@experiment(n).ne.models@experiment(n-1)) then
        dim_EXP = dim_EXP+1
        ltmp = array_append_record(label_exp,models@experiment(n),0)
        delete(label_exp)
        label_exp = ltmp
        delete(ltmp)
      end if
    end do
    enssize =num(ind(models@experiment.eq.models@experiment(0)))
    yearsize = max(toint(models@end_year))-min(toint(models@start_year))+1
    error_bar_PI = new((/dimsizes(diag_script_info@e13fig06_season)*(dim_EXP)*2/),graphic)
    error_bar_CI = new((/dimsizes(diag_script_info@e13fig06_season)*(dim_EXP)*2/),graphic)
    
    if ((isatt(diag_script_info,"list_chem_mod"))) then
      list_chem_mod = diag_script_info@list_chem_mod
      attrib_chem="_CHEM"
      if ((dimsizes(list_chem_mod).eq.1).and.(list_chem_mod.eq."All")) then
          delete(list_chem_mod)
          list_chem_mod = models@name(:enssize)
          attrib_chem=" "
      end if
    else 
      list_chem_mod = models@name 
      attrib_chem=" " 
    end if
    if ((isatt(diag_script_info,"list_nochem_mod"))) then
      list_nochem_mod = diag_script_info@list_nochem_mod
      attrib_nochem="_NOCHEM"
    end if
    if(isatt(diag_script_info, "e13fig06_baseline")) then
      baseline = diag_script_info@e13fig06_baseline
    else
      print("fatal:plot requires diag_script_info attribute 'e13fig06_baseline'")
      status_exit(1)
    end if
    if(isatt(diag_script_info, "e13fig06_labels_exp_esp")) then
          labels_exp_esp = diag_script_info@e13fig06_labels_exp_esp
    else
          labels_exp_esp ="False"
    end if

    if(isatt(diag_script_info, "e13fig06_baseline_adj")) then
      baseline_adj = diag_script_info@e13fig06_baseline_adj
    else
      print("fatal:plot requires diag_script_info attribute 'e13fig06_baseline_adj'")
      status_exit(1)
    end if

    if(isatt(diag_script_info, "e13fig06_mod_plot_PI")) then
      mod_plot_PI = diag_script_info@e13fig06_mod_plot_PI
    else
      mod_plot_PI="False"
    end if

    if(isatt(diag_script_info, "e13fig06_mod_plot_CI")) then
      mod_plot_CI = diag_script_info@e13fig06_mod_plot_CI
    else
      mod_plot_CI="False"
    end if

    if(isvar("latrange")) then
      delete(latrange)
    end if
    if(isatt(diag_script_info, "e13fig06_latrange")) then
      latrange = diag_script_info@e13fig06_latrange
    else
      print("fatal:plot requires diag_script_info attribute 'e13fig06_latrange'")
      status_exit(1)
    end if
    if(isatt(diag_script_info,"e13fig06_season")) then
      season = diag_script_info@e13fig06_season
    else
      if (isStrSubset(field_type0, "z") .eq. False) then
        print("fatal:plot requires diag_script_info attribute 'e13fig06_season'")
        status_exit(1)
      end if
    end if
    if(isatt(diag_script_info, "e13fig06_lev")) then
        lev = toint(diag_script_info@e13fig06_lev)
    end if
    if ((isatt(diag_script_info,"e13fig06_XMin"))) then
       XMin = diag_script_info@e13fig06_XMin
    else
       XMin= new(dimsizes(season),"float")
    end if
    
    if ((isatt(diag_script_info,"e13fig06_XMax"))) then
       XMax = diag_script_info@e13fig06_XMax
    else
       XMax=new(dimsizes(season),"float")
    end if
    
    if ((isatt(diag_script_info,"e13fig06_YMin"))) then
       YMin = diag_script_info@e13fig06_YMin
    else
       YMin=new(dimsizes(season),"float")
    end if
    
    if ((isatt(diag_script_info,"e13fig06_YMax"))) then
       YMax = diag_script_info@e13fig06_YMax
    else
       YMax=new(dimsizes(season),"float")
    end if
   
    if ((isatt(diag_script_info,"e13fig06_Obs"))) then
       Obs=diag_script_info@e13fig06_Obs 
        
    end if 

   if ((isatt(diag_script_info,"e13fig06_Obs_dir"))) then
        Obs_dir=diag_script_info@e13fig06_Obs_dir
        
   end if

  color_fill_exp=(/(/0.90,0.90,0.90/),\
                    (/0.105,0.67,0.8/),\
                    (/0.88,0.88,0.88/),\
                    (/0.66,0.66,0.66/),\
                    (/0.55,0.55,0.55/),\
                    (/0.33,0.33,0.33/),\
                    (/0.22,0.22,0.22/)/)
  
   plot_poly = new(dimsizes(diag_script_info@e13fig06_season)*dim_EXP,"graphic")
   if ((isatt(diag_script_info,"list_nochem_mod"))) then 
      plot_line=  new(dimsizes(diag_script_info@e13fig06_season)*dim_EXP*(dimsizes(list_nochem_mod)+1),"graphic")
   end if
;----------------------------------------------------------------------------
ip=0
in=0
iCI=0
iPI=0
do ipan = 0,dimsizes(season)-1 ;panel loop

    month=(/"JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP","OCT","NOV","DEC"/)
    ind_month=season(ipan)  
    do im=0,11
       if (season(ipan).eq.month(im))
          delete(ind_month) 
          ind_month=im+1
       end if
    end do  
    colors_exp=new(dim_EXP,"string")
    colors_exp_C=new(dim_EXP,"string")  
    colors_exp_P=new(dim_EXP,"string")
    thicks_exp=new(dim_EXP,"float")
    dashs_exp=new(dim_EXP,"integer")
  dashs_exp_nochem=new(dim_EXP,"integer")
;**********************chem var*********************************************
    chem_t       = new((/dim_EXP,dimsizes(list_chem_mod),yearsize/),"float")
    chem_t_IMT       = new((/dim_EXP,dimsizes(list_chem_mod),yearsize/),"float")
    chem_MMT=new((/dim_EXP,5,yearsize/),"float")
    colors_chem=new((/dim_EXP,dimsizes(list_chem_mod)/),typeof(colors)) 
    dashes_chem=new((/dim_EXP,dimsizes(list_chem_mod)/),typeof(dashes))
    thicks_chem=new((/dim_EXP,dimsizes(list_chem_mod)/),typeof(thicks))
    lines_chem=new((/dim_EXP,dimsizes(list_chem_mod)/),typeof(colors))
    exp_chem=new((/dim_EXP,dimsizes(list_chem_mod)/),typeof(models@experiment))
;**********************no chem var*******************************************
    if ((isatt(diag_script_info,"list_nochem_mod"))) then 
       nochem_t       = new((/dim_EXP,dimsizes(list_nochem_mod),yearsize/),"float")
       nochem_t_IMT       = new((/dim_EXP,dimsizes(list_nochem_mod),yearsize/),"float")
       nochem_MMT=new((/dim_EXP,5,yearsize/),"float")
       colors_nochem=new((/dim_EXP,dimsizes(list_nochem_mod)/),typeof(colors))
       dashes_nochem=new((/dim_EXP,dimsizes(list_nochem_mod)/),typeof(dashes))
       thicks_nochem=new((/dim_EXP,dimsizes(list_nochem_mod)/),typeof(thicks))
       lines_nochem=new((/dim_EXP,dimsizes(list_nochem_mod)/),typeof(colors))
       exp_nochem=new((/dim_EXP,dimsizes(list_nochem_mod)/),typeof(models@experiment))
    end if
;**************************************************************************
    time_c=ispan(min(toint(models@start_year)),max(toint(models@end_year)),1)
    iexp = 0
;----------------------------chem------------------------------------
    chem = 0
    time_chem  = new((/dim_EXP,dimsizes(list_chem_mod),yearsize/),float)
    time_chem!0="exp"
    time_chem!1="model"
    time_chem!2="time"
    time_chem&time=(/time_c/)
 ;--------------------------nochem----------------------------------------
    if ((isatt(diag_script_info,"list_nochem_mod"))) then
       nochem = 0
       time_nochem  = new((/dim_EXP,dimsizes(list_nochem_mod),yearsize/),float)
       time_nochem!0="exp"
       time_nochem!1="model"
       time_nochem!2="time"
       time_nochem&time=(/time_c/)
    end if
;----------------------------chem--------------------------------------------------
   chem_t_IMT!2="time"
   chem_t_IMT&time=time_c
   chem_MMT!0="exp"
   chem_MMT!1="case"  
   chem_MMT!2="time"
   chem_MMT&time=time_c

   chem_t!2="time"
   chem_t&time=time_c
    
 ;--------------------------nochem----------------------------------------
    if ((isatt(diag_script_info,"list_nochem_mod"))) then
       nochem_t_IMT!2="time"
       nochem_t_IMT&time=time_c
       nochem_MMT!0="exp"
       nochem_MMT!1="case"
       nochem_MMT!2="time"
       nochem_MMT&time=time_c

       nochem_t!2="time"
       nochem_t&time=time_c
    end if
;******************************************************************************
    do imod = 0,dim_MOD-1
       A0 = extract_data(imod, data0, -1, 0, 0)
       ;uniform_coordinates(A0,field_type0)
       if(isatt(diag_script_info, "e13fig06_lev")) then

          if (max(A0&plev).ge.100000) then
             A0&plev = A0&plev/100
          end if

          A1=A0(time|:,{plev|lev(ipan)},lat|:,lon|:)
          delete(A0)
          A0=A1
          delete(A1) 
       end if 


       ;ttmp= ispan(toint(models@start_year(imod)),toint(models@end_year(imod)),1)
       ;time_chem(iexp,{toint(models@start_year(imod)):toint(models@end_year(imod))}) =(/ttmp/)
       dummy1 = area_operations(A0(time|:,lat|:,lon|:),latrange(ipan,0),latrange(ipan,1),0.,360.,"average",True)
       if (var0.eq."tro3") then
          dummy1 = dummy1*10^-3
       end if
       delete(A0)
       if (season(ipan).eq."ANN") then
          dummy2 = time_operations(dummy1,toint(models@start_year(imod)),toint(models@end_year(imod)),"average","yearly",True)
          ;tmp = time_operations(tmp,toint(models@start_year(imod)),toint(models@end_year(imod)),"average","yearly",True)
          ;delete(tmp)
       else
          dummy2 = time_operations(dummy1,toint(models@start_year(imod)),toint(models@end_year(imod)),"extract",ind_month,True)
          len_of_str = sizeof(stringtochar(ind_month)) - 1
          dummy22=dim_avg_n_Wrap(onedtond(dummy2,(/len_of_str,toint(models@end_year(imod))-toint(models@start_year(imod))+1/)),0)
         delete(dummy2)
         dummy2=dummy22
         delete(dummy22) 
         delete(len_of_str) 
      end if
       if (imod.ge.1) then
          if (models@experiment(imod).ne.models@experiment(imod-1)) then
            iexp = iexp+1
            nochem = 0
            chem = 0
          end if
        end if  
 ;-----------------------------chem----------------------------------------------- 
        if (any(models@name(imod).eq.list_chem_mod)) then
           chem_t(iexp,chem,{toint(models@start_year(imod)):toint(models@end_year(imod))}) = tofloat(dummy2)
           time_chem(iexp,chem,{toint(models@start_year(imod)):toint(models@end_year(imod))})=ispan(toint(models@start_year(imod)),toint(models@end_year(imod)),1)
           dim_m=dimsizes(ispan(toint(models@start_year(imod)),toint(models@end_year(imod)),1)) 
           colors_chem(iexp,chem)=colors(imod)
           dashes_chem(iexp,chem)=dashes(imod)
           thicks_chem(iexp,chem)=2. 
           lines_chem(iexp,chem)="Lines"
           exp_chem(iexp,chem)=models@experiment(imod) 
           chem = chem+1

        end if
;-----------------------------nochem-----------------------------------------------
       if ((isatt(diag_script_info,"list_nochem_mod"))) then
          if (any(models@name(imod).eq.list_nochem_mod)) then
             nochem_t(iexp,nochem,{toint(models@start_year(imod)):toint(models@end_year(imod))}) = tofloat(dummy2)
             time_nochem(iexp,nochem,{toint(models@start_year(imod)):toint(models@end_year(imod))})=ispan(toint(models@start_year(imod)),toint(models@end_year(imod)),1)
             dim_m=dimsizes(ispan(toint(models@start_year(imod)),toint(models@end_year(imod)),1))
             colors_nochem(iexp,nochem)=colors(imod)
             dashes_nochem(iexp,nochem)=dashes(imod)
             thicks_nochem(iexp,nochem)=1.
             lines_nochem(iexp,nochem)="Lines"
             exp_nochem(iexp,nochem)=models@experiment(imod)
             nochem = nochem+1 
          end if
       end if
;****************************************************************************************
       delete(dummy1)
       delete(dummy2)

    end do ;imod
;-------------prepare the input for the r script-----------
 ;;;;;;;;;dat=ozone year model run;;;;;;;
 
    chem_run_t=chem_t
    chem_run_t=where(.not.ismissing(chem_t),1.,0.)
    chem_mod_t=chem_t
    ;chem_mod_t=where(
    
    if ((isatt(diag_script_info,"list_nochem_mod"))) then
       nochem_run_t=nochem_t
       nochem_run_t=where(.not.ismissing(nochem_t),1.,0.)
       nochem_mod_t=nochem_t
    end if


   if ((label_exp(0).eq."historical").and.(dimsizes(label_exp.ge.1)))
       nexp=1
       do iexp=nexp,dim_EXP-1
          do imod=0,dimsizes(list_chem_mod)-1
             do it=0,yearsize-1
                chem_mod_t(iexp,imod,it)=imod+1  
             end do
             chem_t(iexp,imod,:)=(/dim_avg_n_Wrap((/chem_t(0,imod,:),chem_t(iexp,imod,:)/),0)/)
             time_chem(iexp,imod,:)=(/dim_avg_n_Wrap((/time_chem(0,imod,:),time_chem(iexp,imod,:)/),0)/)  
             chem_run_t(iexp,imod,:)=(/dim_sum_n_Wrap((/chem_run_t(0,imod,:),chem_run_t(iexp,imod,:)/),0)/)     
          end do
          if ((isatt(diag_script_info,"list_nochem_mod"))) then
             do imod=0,dimsizes(list_nochem_mod)-1
                do it=0,yearsize-1
                   nochem_mod_t(iexp,imod,it)=imod+1
                end do
                nochem_t(iexp,imod,:)=(/dim_avg_n_Wrap((/nochem_t(0,imod,:),nochem_t(iexp,imod,:)/),0)/)
                time_nochem(iexp,imod,:)=(/dim_avg_n_Wrap((/time_nochem(0,imod,:),time_nochem(iexp,imod,:)/),0)/)
                nochem_run_t(iexp,imod,:)=(/dim_sum_n_Wrap((/nochem_run_t(0,imod,:),nochem_run_t(iexp,imod,:)/),0)/)
             end do
         end if

       end do
    else
       nexp=0
       do iexp=nexp,dim_EXP-1
          do imod=0,dimsizes(list_chem_mod)-1
             do it=0,yearsize-1
                chem_mod_t(iexp,imod,it)=imod+1
             end do
          end do
          if ((isatt(diag_script_info,"list_nochem_mod"))) then
             do imod=0,dimsizes(list_nochem_mod)-1
                do it=0,yearsize-1
                   nochem_mod_t(iexp,imod,it)=imod+1
                end do
             end do
         end if

       end do
    end if 
    do iexp=nexp,dim_EXP-1
        dashs_exp(iexp)=0
        dashs_exp_nochem(iexp)=1
        thicks_exp(iexp)=3.  
        if (((label_exp(iexp)).eq."historical").or.((label_exp(iexp)).eq."refC2")) then
           colors_exp(iexp)="black"
        end if
        if (((label_exp(iexp)).eq."rcp26").or.((label_exp(iexp)).eq."senC2rcp26").or.((label_exp(iexp)).eq."senC2fGHG")) then
           colors_exp(iexp)="(/0.,      0.,      1./)"
           colors_exp_C(iexp)="(/0.6,0.6,0.6/)"
            colors_exp_P(iexp)="(/0.8,0.8,0.8/)"
        end if
        if (((label_exp(iexp)).eq."rcp45").or.((label_exp(iexp)).eq."senC2rcp45")) then
           colors_exp(iexp)="(/0,       0.4757,    0.7373,     1./)"
           ;colors_exp_C(iexp)="(/0,       0.4757,    0.7373,     1./)"
;colors_exp_C(iexp)="(/0,       0.4757,    0.7373,     1./)"

        end if
        if ((label_exp(iexp)).eq."rcp60") then
           colors_exp(iexp)="(/1.,      0.5098,  0.1765/)"
        end if
        if (((label_exp(iexp)).eq."rcp85").or.((label_exp(iexp)).eq."senC2rcp85").or.((label_exp(iexp)).eq."senC2fODS")) then
           colors_exp(iexp)="(/1.,      0.,      0./)"
        end if
    
        dat_chem=ndtooned(chem_t(iexp,:,:))

        dat_t_chem=ndtooned(time_chem(iexp,:,:))
        dat_m_chem=ndtooned(chem_mod_t(iexp,:,:)) 
        dat_r_chem=ndtooned(chem_run_t(iexp,:,:))

        dat_chem_f=dat_chem(ind(.not.ismissing(dat_t_chem)))
        dat_t_chem_f=dat_t_chem(ind(.not.ismissing(dat_t_chem)))
        dat_m_chem_f=dat_m_chem(ind(.not.ismissing(dat_t_chem)))
        dat_r_chem_f=dat_r_chem(ind(.not.ismissing(dat_t_chem)))

        fmtf   = "f8.1,2x,f8.3,2x,f8.1,2x,f8.1,2x,f8.1"                           ; same as "7f7.2"

        dat_input_chem=new((/dimsizes(dat_chem_f),5/),"float")
        dat_input_chem(:,0)=ispan(1,dimsizes(dat_chem_f),1)
        dat_input_chem(:,1)=dat_chem_f
        dat_input_chem(:,2)=dat_t_chem_f
        dat_input_chem(:,3)=dat_m_chem_f
        dat_input_chem(:,4)= dat_r_chem_f
         opt        = True
        opt@title  ="   ozone year model run"
        fou_chem= label_exp(iexp)+"_chem_"+season(ipan)+"_"+latrange(ipan,0)+"_"+latrange(ipan,1)
        opt@fout = output_dir+"RAW_"+fou_chem+".ascii"
        write_matrix (dat_input_chem, fmtf, opt)
;-------------------nochem-----------------------------------------------
       if ((isatt(diag_script_info,"list_nochem_mod"))) then
           dat_nochem=ndtooned(nochem_t(iexp,:,:))
           dat_t_nochem=ndtooned(time_nochem(iexp,:,:))
           dat_m_nochem=ndtooned(nochem_mod_t(iexp,:,:))
           dat_r_nochem=ndtooned(nochem_run_t(iexp,:,:))

           dat_nochem_f=dat_nochem(ind(.not.ismissing(dat_t_nochem)))
           dat_t_nochem_f=dat_t_nochem(ind(.not.ismissing(dat_t_nochem)))
           dat_m_nochem_f=dat_m_nochem(ind(.not.ismissing(dat_t_nochem)))
           dat_r_nochem_f=dat_r_nochem(ind(.not.ismissing(dat_t_nochem)))
    
           fmtf   = "f8.1,2x,f8.3,2x,f8.1,2x,f8.1,2x,f8.1"                           ; same as "7f7.2"

           dat_input_nochem=new((/dimsizes(dat_nochem_f),5/),"float")
           dat_input_nochem(:,0)=ispan(1,dimsizes(dat_nochem_f),1)
           dat_input_nochem(:,1)=dat_nochem_f
           dat_input_nochem(:,2)=dat_t_nochem_f
           dat_input_nochem(:,3)=dat_m_nochem_f
           dat_input_nochem(:,4)= dat_r_nochem_f
           
           opt        = True
           opt@title  ="   ozone year model run"
           fou_nochem= label_exp(iexp)+"_nochem_"+season(ipan)+"_"+latrange(ipan,0)+"_"+latrange(ipan,1)
           opt@fout = output_dir+"RAW_"+fou_nochem+".ascii"
           write_matrix (dat_input_nochem, fmtf, opt)
        end if 
;------------------------------------------------------------------------
;*********      launch r     chem      ***********************
       list_run_chem_mod=list_chem_mod
       list_run_chem_mod=where(.not.ismissing(list_chem_mod),1,0)
       system ("Rscript  ./diag_scripts/lib/R/TSAM.r  "+output_dir+" "+baseline+" "+dimsizes(list_chem_mod)+" "+ min(toint(models@start_year))+" "+ max(toint(models@end_year))+" "+yearsize+" "+fou_chem+" "+str_join(list_chem_mod," ")+" "+str_join(list_run_chem_mod," "))
print("Rscript  ./diag_scripts/lib/R/TSAM.r  "+output_dir+" "+baseline+" "+dimsizes(list_chem_mod)+" "+ min(toint(models@start_year))+" "+ max(toint(models@end_year))+" "+yearsize+" "+fou_chem+" "+str_join(list_chem_mod," ")+" "+str_join(list_run_chem_mod," "))

;**************************************************************

       finp_MMT=output_dir+"TSAM_"+fou_chem+"_"+baseline+"_MMT_estimate.ascii"
       print(finp_MMT)
       ff1=readAsciiTable(finp_MMT,6, "float",0)
       ff1=where(ff1.eq.0,ff1@_FillValue,ff1)
       tt=ff1(:,0)
       if (baseline_adj.eq."True") then
          i_baseline=ind(tt.eq.tofloat(baseline))
          chem_MMT(exp|iexp,{time|toint(tt)},case|:)=(/ff1(:,1:5)-ff1(i_baseline,1)/)
       else
          chem_MMT(exp|iexp,{time|toint(tt)},case|:)=(/ff1(:,1:5)/)
       end if
       delete(ff1)
       delete(tt)
       delete(finp_MMT)
       do imod=0,dimsizes(list_chem_mod)-1
          finp_m_IMT=output_dir+"TSAM_"+fou_chem+"_"+baseline+"_IMT_estimate_"+list_chem_mod(imod)+".ascii"  
          ff2 = readAsciiTable(finp_m_IMT,4, "float",0)
          ff2=where(ff2.eq.0,ff2@_FillValue,ff2)
          tt=ff2(:,0)
          chem_t_IMT(iexp,imod,{toint(tt)})=(/ff2(:,1)/)
          if (baseline_adj.eq."True") then
             i_baseline=ind(tt.eq.tofloat(baseline))   
             chem_t_IMT(iexp,imod,{toint(tt)})=(/ff2(:,1)-ff2(i_baseline,1)/)
          else
             chem_t_IMT(iexp,imod,{toint(tt)})=(/ff2(:,1)/)  
          end if
          delete(ff2)
          delete(tt)
          delete(finp_m_IMT) 
       end do
   
;------------------------------nochem---------------------------------------------
;*********      launch r     nochem      ***********************
       if ((isatt(diag_script_info,"list_nochem_mod"))) then  
          list_run_nochem_mod=list_nochem_mod
          list_run_nochem_mod=where(.not.ismissing(list_nochem_mod),1,0)

          system ("Rscript  ./diag_scripts/lib/R/TSAM.r  "+output_dir+" "+baseline+" "+dimsizes(list_nochem_mod)+" "+ min(toint(models@start_year))+" "+ max(toint(models@end_year))+" "+yearsize+" "+fou_nochem+" "+str_join(list_nochem_mod," ")+" "+str_join(list_run_nochem_mod," "))
print("Rscript  ./diag_scripts/lib/R/TSAM.r  "+output_dir+" "+baseline+" "+dimsizes(list_nochem_mod)+" "+ min(toint(models@start_year))+" "+ max(toint(models@end_year))+" "+yearsize+" "+fou_nochem+" "+str_join(list_nochem_mod," ")+" "+str_join(list_run_nochem_mod," "))
;**************************************************************

          finp_MMT=output_dir+"TSAM_"+fou_nochem+"_"+baseline+"_MMT_estimate.ascii"
          ff1=readAsciiTable(finp_MMT,6, "float",0)
          tt=ff1(:,0)
          ff1=where(ff1.eq.0,ff1@_FillValue,ff1)
          if (baseline_adj.eq."True") then
             i_baseline=ind(tt.eq.tofloat(baseline))
             nochem_MMT(exp|iexp,{time|toint(tt)},case|:)=(/ff1(:,1:5)-ff1(i_baseline,1)/)
          else 
             nochem_MMT(exp|iexp,{time|toint(tt)},case|:)=(/ff1(:,1:5)/)
          end if
          delete(ff1)
          delete(tt)
          delete(finp_MMT)
          do imod=0,dimsizes(list_nochem_mod)-1
             finp_m_IMT=output_dir+"TSAM_"+fou_nochem+"_"+baseline+"_IMT_estimate_"+list_nochem_mod(imod)+".ascii"
             ff2 = readAsciiTable(finp_m_IMT,4, "float",0)
             ff2=where(ff2.eq.0,ff2@_FillValue,ff2)
             tt=ff2(:,0)
             if (baseline_adj.eq."True") then
               i_baseline=ind(tt.eq.tofloat(baseline))
               nochem_t_IMT(iexp,imod,{toint(tt)})=(/ff2(:,1)-ff2(i_baseline,1)/)
             else   
               nochem_t_IMT(iexp,imod,{toint(tt)})=(/ff2(:,1)/)
             end if
             delete(ff2)
             delete(tt)
             delete(finp_m_IMT)
          end do
          delete(dat_nochem_f)
          delete(dat_t_nochem_f)
          delete(dat_m_nochem_f)
          delete(dat_r_nochem_f)
          delete(dat_input_nochem)
        end if ;nochem
      delete(dat_chem_f)
      delete(dat_t_chem_f)
      delete(dat_m_chem_f)
      delete(dat_r_chem_f)
      delete(dat_input_chem)
    end do ;experiment
;******************************************
;plot1 raw data*********
;plot2 IMT***********
;plot3 MMT***************

;--------PLOT---------------------------------------

        ti = ispan(1,(dim_EXP)*3,1)
        if (latrange(ipan,0).ge.0.) then
          latstring1 = latrange(ipan,0)+"N"
        else
          latstring1 = -latrange(ipan,0)+"S"
        end if
       
        if (latrange(ipan,1).ge.0.) then
          latstring2 = latrange(ipan,1)+"N"
        else 
          latstring2 = -latrange(ipan,1)+"S"
        end if

        res                   = True                      ; plot mods desired
        res@tiYAxisString     = "("+variable_info@units+")" ; add units title
        res@tiXAxisString = "Year"     ; y-axis label      
        res@tiMainString      = season(ipan)+" "+latstring1+" "+latstring2+" "+var0 ;
        res@gsnFrame          = False                     ; don't draw yet
        res@gsnDraw           = False  
        ;res@tfPolyDrawOrder = "Predraw"

        if(.not.ismissing(XMax(ipan))) then
           res@trXMaxF         = XMax(ipan)
        end if

        if(.not.ismissing(XMin(ipan))) then
          res@trXMinF         = XMin(ipan)
        end if
        res@vpHeightF       = 0.4                    ; change aspect ratio of plot
        res@vpWidthF        = 0.5   
        if (mod_plot.eq."RAW") then
           if(.not.ismissing(YMax(ipan))) then
                res@trYMaxF         = YMax(ipan)
           else
              if ((isatt(diag_script_info,"list_nochem_mod"))) then
                 vv=max((/max(chem_t),max(nochem_t)/))
              else
                 vv=max(chem_t)
              end if
              res@trYMaxF =vv
           end if

           if(.not.ismissing(YMin(ipan))) then
              res@trYMinF         = YMin(ipan)
           else
              if ((isatt(diag_script_info,"list_nochem_mod"))) then
                 vv=min((/min(chem_t),min(nochem_t)/))
              else
                 vv=min(chem_t)
              end if
              res@trYMinF =vv
           end if
           if ((isatt(diag_script_info,"list_nochem_mod"))) then
              thicks_plot=array_append_record(ndtooned(thicks_chem),ndtooned(thicks_nochem),0)  
              colors_plot=array_append_record(ndtooned(colors_chem),ndtooned(colors_nochem),0) 

              dashes_plot=array_append_record(ndtooned(dashes_chem),ndtooned(dashes_nochem),0)
              models_plot=array_append_record(list_chem_mod,list_nochem_mod,0)
              exp_plot=array_append_record(exp_chem,exp_nochem,0)
              a=onedtond(ndtooned(chem_t),(/dim_EXP*(dimsizes(list_chem_mod)),yearsize/))
          
              b=onedtond(ndtooned(nochem_t),(/dim_EXP*(dimsizes(list_nochem_mod)),yearsize/))
              values_plot=table_attach_rows(a,b,0)
              delete(a)
              delete(b)
              a=onedtond(ndtooned(time_chem),(/dim_EXP*(dimsizes(list_chem_mod)),yearsize/))

              b=onedtond(ndtooned(time_nochem),(/dim_EXP*(dimsizes(list_nochem_mod)),yearsize/))
              times_plot=table_attach_rows(a,b,0)
              delete(a)
              delete(b)
  
           else
              thicks_plot=ndtooned(thicks_chem)
              colors_plot=ndtooned(colors_chem)
              dashes_plot=ndtooned(dashes_chem)
              models_plot=list_chem_mod
              exp_plot=ndtooned(exp_chem)
              a=onedtond(ndtooned(chem_t),(/dim_EXP*(dimsizes(list_chem_mod)),yearsize/))
              values_plot=a
              delete(a)
              a=onedtond(ndtooned(time_chem),(/dim_EXP*(dimsizes(list_chem_mod)),yearsize/))
              times_plot=a
              delete(a) 
           end if
           if (labels_exp_esp.eq."True") then
              res@xyLabelMode            = "Custom"           ; label a line
              res@xyExplicitLabels       = exp_plot            ;text to use
              res@xyLineLabelFontHeightF = 0.010              ; font height
               res@xyLineLabelConstantSpacingF=0.1
              res@xyLineLabelFontColors   = colors_plot              ; label color 
           end if  
              res@xyLineThicknesses=thicks_plot
               res@xyLineColors  = colors_plot                   ; line color
               res@xyDashPatterns = dashes_plot
               plot(ipan)=gsn_csm_xy(wks,times_plot,values_plot,res)

        end if
        if (mod_plot.eq."IMT") then
           if(.not.ismissing(YMax(ipan))) then
              res@trYMaxF         = YMax(ipan)
           else
              if ((isatt(diag_script_info,"list_nochem_mod"))) then
                 vv=max((/max(chem_t_IMT),max(nochem_t_IMT)/))
              else
                 vv=max(chem_t_IMT)
              end if
              res@trYMaxF =vv
           end if

           if(.not.ismissing(YMin(ipan))) then
              res@trYMinF         = YMin(ipan)
           else
              if ((isatt(diag_script_info,"list_nochem_mod"))) then
                 vv=min((/min(chem_t_IMT),min(nochem_t_IMT)/))
              else
                 vv=min(chem_t_IMT)
              end if
              res@trYMinF =vv
          end if   
          if ((isatt(diag_script_info,"list_nochem_mod"))) then
              thicks_plot=array_append_record(ndtooned(thicks_chem),ndtooned(thicks_nochem),0)
              colors_plot=array_append_record(ndtooned(colors_chem),ndtooned(colors_nochem),0)
              exp_plot=array_append_record(ndtooned(exp_chem),ndtooned(exp_nochem),0)
              dashes_plot=array_append_record(ndtooned(dashes_chem),ndtooned(dashes_nochem),0)
              a=onedtond(ndtooned(chem_t_IMT),(/dim_EXP*(dimsizes(list_chem_mod)),yearsize/))

              b=onedtond(ndtooned(nochem_t_IMT),(/dim_EXP*(dimsizes(list_nochem_mod)),yearsize/))
              values_plot=table_attach_rows(a,b,0)
              models_plot=array_append_record(list_chem_mod,list_nochem_mod,0)
              delete(a)
              delete(b)
              a=onedtond(ndtooned(time_chem),(/dim_EXP*(dimsizes(list_chem_mod)),yearsize/))

              b=onedtond(ndtooned(time_nochem),(/dim_EXP*(dimsizes(list_nochem_mod)),yearsize/))
              times_plot=table_attach_rows(a,b,0)
              delete(a)
              delete(b)

           else
              thicks_plot=ndtooned(thicks_chem)
              colors_plot=ndtooned(colors_chem)
              dashes_plot=ndtooned(dashes_chem)
              exp_plot=ndtooned(exp_chem)  
              models_plot=list_chem_mod

              a=onedtond(ndtooned(chem_t_IMT),(/dim_EXP*(dimsizes(list_chem_mod)),yearsize/))
              values_plot=a
              delete(a)
               a=onedtond(ndtooned(time_chem),(/dim_EXP*(dimsizes(list_chem_mod)),yearsize/))
              times_plot=a
              delete(a)
           end if
           if (labels_exp_esp.eq."True") then
              res@xyLabelMode            = "Custom"           ; label a line
              res@xyExplicitLabels       = exp_plot            ;text to use
              res@xyLineLabelFontHeightF = 0.010              ; font height
              res@xyLineLabelConstantSpacingF=0.1
              res@xyLineLabelFontColors   = colors_plot              ; label color
           end if 
           res@xyLineThicknesses=thicks_plot
           res@xyLineColors  = colors_plot                   ; line color
           res@xyDashPatterns = dashes_plot
           plot(ipan)=gsn_csm_xy(wks,times_plot,values_plot,res) 
        end if
        if (mod_plot.eq."MMT") then
           if ((isatt(diag_script_info,"list_nochem_mod"))) then
              thicks_plot=array_append_record(thicks_exp(nexp:),thicks_exp(nexp:),0)
              colors_plot=array_append_record(colors_exp(nexp:),colors_exp(nexp:),0)
              
              dashes_plot=array_append_record(dashs_exp(nexp:),dashs_exp_nochem(nexp:),0)
              values_plot=table_attach_rows(chem_MMT(nexp:,0,:),nochem_MMT(nexp:,0,:),0)
              values_plot_Min_C=table_attach_rows(chem_MMT(nexp:,1,:),nochem_MMT(nexp:,1,:),0)
              values_plot_Max_C=table_attach_rows(chem_MMT(nexp:,2,:),nochem_MMT(nexp:,2,:),0)

              values_plot_Min_P=table_attach_rows(chem_MMT(nexp:,3,:),nochem_MMT(nexp:,3,:),0)
              values_plot_Max_P=table_attach_rows(chem_MMT(nexp:,4,:),nochem_MMT(nexp:,4,:),0)
              models_plot=array_append_record(label_exp(nexp:)+" Chem",label_exp(nexp:)+" NoChem",0)
              times_plot=chem_MMT&time 
              ntimes=dimsizes(times_plot)
           else
              thicks_plot=thicks_exp(nexp:)
              colors_plot=colors_exp(nexp:)
              dashes_plot=dashs_exp(nexp:)
              models_plot=label_exp(nexp:)
              values_plot=chem_MMT(nexp:,0,:)
              values_plot_Min_C=chem_MMT(nexp:,1,:)
              values_plot_Max_C=chem_MMT(nexp:,2,:)

              values_plot_Min_P=chem_MMT(nexp:,3,:)
              values_plot_Max_P=chem_MMT(nexp:,4,:)   
              times_plot=chem_MMT&time
              ntimes=dimsizes(times_plot)   
           end if
           yp_C    = new( (/dimsizes(values_plot(:,0)),2*ntimes/), float )
           xp_C    =  new( (/2*ntimes/), float )   
           do kt=0,ntimes-1
              xp_C(kt) =times_plot(kt)
              yp_C(:,kt) =values_plot_Max_C(:,kt)
              xp_C(2*ntimes-1-kt) =times_plot(kt)
              yp_C(:,2*ntimes-1-kt) = values_plot_Min_C(:,kt)
           end do
           yp_P    = new( (/dimsizes(values_plot(:,0)),2*ntimes/), float )
           xp_P    =  new( (/2*ntimes/), float )
           do kt=0,ntimes-1
              xp_P(kt) =times_plot(kt)
              yp_P(:,kt) =values_plot_Max_P(:,kt)
              xp_P(2*ntimes-1-kt) =times_plot(kt)
              yp_P(:,2*ntimes-1-kt) = values_plot_Min_P(:,kt)
           end do
           res@xyLineThicknesses=thicks_plot
           res@xyLineColors  =colors_plot                    ; line color
           res@xyDashPatterns =dashes_plot
           res@tfPolyDrawOrder    = "Predraw"
            if(.not.ismissing(YMin(ipan))) then
              res@trYMinF         = YMin(ipan)
            end if
            if(.not.ismissing(YMax(ipan))) then
              res@trYMaxF         = YMax(ipan)
            end if
           plot(ipan) =gsn_csm_xy (wks,times_plot,values_plot,res)
  
          do iexp=0,dimsizes(values_plot(:,0))-1 

             gsres                   = True                        ; poly res
             gsres@gsFillColor       = colors_plot(iexp)                  ; color chosen
             
             gsres@gsFillIndex     = dashes_plot(iexp)*17   
             gsres@gsFillScaleF = 0.5
             if (mod_plot_CI.eq."True") then
                  gsres@gsFillOpacityF = 0.3 
                  error_bar_CI(iCI) = gsn_add_polygon (wks,plot(ipan),xp_C,yp_C(iexp,:),gsres)
                iCI=iCI+1
             end if
             if (mod_plot_PI.eq."True") then
                gsres@gsFillOpacityF = 0.1
                error_bar_PI(iPI) = gsn_add_polygon (wks,plot(ipan),xp_P,yp_P(iexp,:),gsres)
                iPI=iPI+1
             end if 

          end do

      end if
     
 end do  ;season

 Pres          =True
 Pres@gsnFrame        = False
 Pres@gsnPanelBottom  = 0.1

  gsn_panel(wks,plot,(/tointeger(dimsizes(season)/2),2/),Pres)
     scatline=new(dimsizes(colors_plot),"string")
     scatline=where(ismissing(scatline),"Lines","Lines")
    if (mod_plot.eq."IMT").or.(mod_plot.eq."RAW") then
       
        add_legenda(wks,models_plot,colors_plot(:dim_EXP-1),dashes_plot(:dim_EXP-1),scatline(:dim_EXP-1),thicks_plot(:dim_EXP-1),"0.1")
    else
        add_legenda(wks,models_plot,colors_plot,dashes_plot,scatline,thicks_plot,"0.1")
    end if

frame(wks)

    info_output(">>>>>>>> Leaving " + diag_script, verbosity, 4)
end 

