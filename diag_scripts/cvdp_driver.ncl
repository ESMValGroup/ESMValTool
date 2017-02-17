;
; CVDP driver script. To run the CVDP at the command line type: ncl driver.ncl
; To run the CVDP at the command line, put it in background mode, and write the terminal output
; to a file named file.out, type: ncl driver.ncl >&! file.out &
; 
;============================================================================================

; A-laue_ax+
  load "./interface_scripts/messaging.ncl"

  diag_script = "cvdp_driver.ncl"
;  outdir            = "./work/plots_CVDP/"   ; location of output files   (must end in a "/")
  outdir = getenv("ESMValTool_plot_dir") + "/"
; A-laue_ax-
                                             ; It is recommended that a new or empty directory be pointed to here 
                                             ; as existing files in outdir can get removed.

  namelists_only       = "False"             ; Set to True to only create the variable namelists. Useful
                                             ; upon running the package for the first time to verify that the correct
                                             ; files are being selected by the package. (See files in namelist_byvar/ directory)
                                             ; Set to False to run the entire package. 
                                             
  obs                  = "True"              ; True = plot OBS (specified in namelist_obs), False = do not
  scale_timeseries     = "True"              ; True = scale timeseries so that x-axis length is comparable across timeseries, False = do not
  output_data          = "True"              ; True = output selected calculated data to a netCDF file. Make sure .nc files from previous CVDP
                                             ;        runs are in outdir or they will get added to/modified.  
  compute_modes_mon    = "False"             ; True = compute DJF, MAM, JJA, SON, Annual and Monthly Atmospheric Modes of Variability 
                                             ; False = do not compute the Monthly Atmospheric Modes of Variability  (saves computation time)                               
;- - - - - - - - - - - - - - - - - -   
  opt_climo         = "Full"                 ; "Full"  = remove climatology based on full record of each simulation,
                                             ; "Custom" = use start and end year of climatological period specified by climo_syear and climo_eyear 
  if (opt_climo.eq."Custom") then
     climo_syear    = 700                    ; remove the climatology/annual cycle based on these years when opt_climo = "Custom"
     climo_eyear    = 719                    ; note: climo_syear and climo_eyear should be within the range of years of all analyzed model runs 
  end if                                     ;       and specified observational datasets.
;- - - - - - - - - - - - - - - - - -   
  colormap          = 0               ; 0 = default colormaps, 1 = colormaps better for color blindness

; A-laue_ax+
;  output_type       = "ps"           ; "ps" or "png". Note that for output_type = "ps", png's are created from the .ps files for web use.
  output_type = getenv("ESMValTool_output_file_type")  ; "ps" or "png". Note that for output_type = "ps",
                                                       ; png's are created from the .ps files for web use.
  if (output_type .ne. "ps" .and. output_type .ne. "png") then
      error_msg("f", diag_script, "", "output_file_type ps and png supported only. Abort.")
  end if
; A-laue_ax-
  
  png_scale         = 1.5             ; Value between .1->5.  Any value > 1 (< 1) increases (decreases) png page size. 
                                      ; When output_type = "png" a value of 1 will result in a png sized 1500 (H) x 1500 (W) before automatic cropping of white space
                                      ; When output_type = "ps"  a value of 1 will result in a png density setting of 144 before automatic cropping of white space
  webpage_title     = "Title goes here"  ; Webpage title   
  
  tar_output        = "False"       ; tar up all output in outdir and remove individual files?                                      
                                    ; Note: ALL files in outdir will be tarred up and then removed fromthe outdir directory. 
                                             
  zp = "./diag_scripts/aux/CVDP/ncl_scripts/"    ; directory path of CVDP NCL scripts. (must end in a "/")
                                    ; Examples: "ncl_scripts/" if all code is local, or on CGD or CISL systems: "~asphilli/CESM-diagnostics/CVDP/ncl_scripts/"
                                    ; Regardless of this setting the following files should be in one directory: namelist, driver.ncl, and namelist_obs. 
                                    ; If pointing to code in ~asphilli make sure the driver script version #s match between the local version of driver.ncl and the 
                                    ; version # in ~asphilli/CESM-diagnostics/CVDP/driver.ncl

; A-laue_ax+  
  ncl_exec = "ncl"       ; this can be changed to a different path if a different version of NCL needs to be used, 
; A-laue_ax-
                         ; such as "/different/path/to/bin/ncl"                                             
;========END USER MODIFICATIONS===========================================================
  version = "3.7.0"
  
  loadscript(zp+"functions.ncl")
  variables = (/"ts","trefht","psl","prect","snowdp", "moc"/)
  do gg = 0,dimsizes(variables)-1
     if (isfilepresent2("obs_"+variables(gg))) then
        system("rm -f obs_"+variables(gg))
     end if 
  end do

  if (isfilepresent2(outdir+"metrics_orig.txt")) then    ; remove metrics_orig.txt file if present
     system("rm -f "+outdir+"metrics_orig.txt")
  end if
  if (isfilepresent2(outdir+"metrics.txt")) then    ; remove metrics.txt file if present
     system("rm -f "+outdir+"metrics.txt")
  end if

  if (.not.isfilepresent2(outdir)) then
     system("mkdir -p "+outdir)
  end if
  if (opt_climo.ne."Custom") then
     climo_syear = -999
     climo_eyear = -999
  end if

  webpage_title = str_sub_str(webpage_title," ","%")
    
  if (namelists_only.eq."True") then
     system("   export OUTDIR="+outdir+";   export OBS="+obs+";   export SCALE_TIMESERIES="+scale_timeseries+";   "+\
         "export OUTPUT_DATA="+output_data+";   export VERSION="+version+";   export PNG_SCALE="+png_scale+";   "+\
         "export OPT_CLIMO="+opt_climo+";   export CLIMO_SYEAR="+climo_syear+";   export CLIMO_EYEAR="+climo_eyear+";   "+\
         "export COMPUTE_MODES_MON="+compute_modes_mon+";   export CVDP_SCRIPTS="+zp+";   export OUTPUT_TYPE="+output_type+";"+\         
         " "+ncl_exec+" "+zp+"namelist.ncl")
     print("Variable namelists have been created. Examine files in namelist_byvar/ directory to verify CVDP file selection.")
     print("Finished: Climate Variability Diagnostics Package")
     exit
  end if             

; If one alters the order that scripts are called in the following system calls one will have to modify the metrics.ncl 
; script to account for this. (metrics.ncl assumes that the metric rows are created in the same order every time.) 
; Adding scripts to the end of this system call is completely fine as that will not alter the order of the 
; calculated metrics.

  ncl_exec = ncl_exec+" -Q"            
  system("   export OUTDIR="+outdir+";   export OBS="+obs+";   export SCALE_TIMESERIES="+scale_timeseries+";   "+\
         "export OUTPUT_DATA="+output_data+";   export VERSION="+version+";   export PNG_SCALE="+png_scale+";   "+\
         "export OPT_CLIMO="+opt_climo+";   export CLIMO_SYEAR="+climo_syear+";   export CLIMO_EYEAR="+climo_eyear+";   "+\
         "export COMPUTE_MODES_MON="+compute_modes_mon+";   export OUTPUT_TYPE="+output_type+";   "+\
         "export COLORMAP="+colormap+";   export CVDP_SCRIPTS="+zp+";"+\     
         " "+str_sub_str(ncl_exec," -Q","")+" "+zp+"namelist.ncl; "+ncl_exec+" "+zp+"sst.indices.ncl; "+ncl_exec+" "+zp+"amo.ncl;"+\         
         " "+ncl_exec+" "+zp+"pdo.ncl; "+ncl_exec+" "+zp+"psl.modes_indices.ncl; "+ncl_exec+" "+zp+"tas.trends_timeseries.ncl;"+\
         " "+ncl_exec+" "+zp+"psl.trends.ncl; "+ncl_exec+" "+zp+"pr.trends_timeseries.ncl; "+ncl_exec+" "+zp+"sst.mean_stddev.ncl;"+\ 
         " "+ncl_exec+" "+zp+"snd.trends.ncl; "+ncl_exec+" "+zp+"psl.mean_stddev.ncl; "+ncl_exec+" "+zp+"tas.mean_stddev.ncl;"+\
         " "+ncl_exec+" "+zp+"pr.mean_stddev.ncl; "+ncl_exec+" "+zp+"sst.trends_timeseries.ncl; "+ncl_exec+" "+zp+"amoc.ncl")                
         
  if (output_type.eq."png") then
     ofiles = systemfunc("ls "+outdir+"*.png")
     do gg = 0,dimsizes(ofiles)-1
        system("convert -trim +repage -border 8 -bordercolor white "+ofiles(gg)+" "+ofiles(gg))
     end do
  else
     ofilesS = systemfunc("ls "+outdir+"*."+output_type)
     ofilesT = systemfunc("ls -l "+outdir+"*."+output_type)
     do gg = 0,dimsizes(ofilesT)-1    ; check for empty .ps files, remove
        filesize =  tofloat(str_get_field(ofilesT(gg),5," "))
        if (filesize.lt.10000) then
           print("Removing: "+ofilesT(gg))
           system("rm -f "+ofilesS(gg))
        end if
     end do
     delete([/filesize,ofilesS,ofilesT/])
     
     ofiles = systemfunc("ls "+outdir+"*."+output_type)
     ofiles_png = str_sub_str(ofiles,"."+output_type,".png")
     d_opt = 144*png_scale
     print("Converting "+output_type+" files to .png")
     do gg = 0,dimsizes(ofiles)-1
        system("convert -density "+d_opt+" -trim +repage -border 8 -bordercolor white -background white -flatten "+ofiles(gg)+" "+ofiles_png(gg))
     end do
     print("Done with "+output_type+"->png conversion")
  end if
  
  system("cp "+zp+"cas-cvdp.png "+outdir)
  system("cp namelist_byvar/* "+outdir)
  system("cp namelist "+outdir)
  if (obs.eq."True") then
     system("cp namelist_obs "+outdir)
  end if
  if (isfilepresent2(outdir+"metrics_orig.txt")) then
     system("   export OUTDIR="+outdir+"; "+ncl_exec+" "+zp+"metrics.ncl")
  end if
  system("   export OUTDIR="+outdir+";   export VERSION="+version+";   export OUTPUT_DATA="+output_data+";   "+\
         "export OPT_CLIMO="+opt_climo+";   export CLIMO_SYEAR="+climo_syear+";   export CLIMO_EYEAR="+climo_eyear+";   "+\
         "export WEBTITLE="+webpage_title+";   export OBS="+obs+";   export CVDP_SCRIPTS="+zp+"; "+ncl_exec+" "+zp+"webpage.ncl")
;-------------------------------
  if (tar_output.eq."True") then
     if (isfilepresent2(outdir+"cvdp.tar")) then
        system("rm -f "+outdir+"cvdp.tar")
     end if
;     print("cd "+outdir+"; tar -cf cvdp.tar *")
     system("cd "+outdir+"; tar -cf cvdp.tar *")
     system("cd "+outdir+"; rm *.png *.ps *.txt *.html *.nc namelist*")
  end if
;-------------------------------
; Cleanup
  do gg = 0,dimsizes(variables)-1
     if (isfilepresent2("obs_"+variables(gg))) then
        system("rm -f obs_"+variables(gg))
     end if 
  end do
  delete([/ofiles,variables,outdir,obs,scale_timeseries,output_data,opt_climo,climo_syear,climo_eyear,png_scale,webpage_title,compute_modes_mon/])
  
  print("Finished: Climate Variability Diagnostics Package")
  
