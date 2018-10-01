######################################################
#-------------ETCCDI preprocessing for HyInt---------#
#-------------E. Arnone (Oct 2017)-------------------#
######################################################
# ABOUT: This function pre-process ETCCDI files obtained with the CRESCENDO_extremeEvents namelist
#        remapping the data from gaussian to lonlat, changing longitude range from 0/360 to -180/180
#        and merging all indices into the HyInt indices file. 

hyint.etccdi.preproc<-function(work_dir,etccdi_dir,etccdi_list_import,cdo_grid,model_idx,season,yrmon="yr") {
# load settings
source('esmvaltool/diag_scripts/aux/hyint/hyint_parameters.r')
for (myname in names(settings)) { temp=get(myname,settings); assign(myname,temp) }

    year1  <- toString(models_start_year[model_idx])
    year2  <- toString(models_end_year[model_idx])
    print(str(c(year1,year2)))
    hyint_file<-getfilename.indices(work_dir,diag_base,model_idx,season)
    etccdi_files<-getfilename.etccdi(etccdi_dir,etccdi_list_import,model_idx,yrmon="yr")
    for (sfile in etccdi_files) {      
      cdo_command<-paste0("cdo -sellonlatbox,-180,180,-90,90  -delvar,time_bnds ",sfile," ",sfile,"_tmp") 
      if (rgrid != F) { 
        cdo_command<-paste0("cdo setgrid,",cdo_grid," -delvar,time_bnds ",sfile," ",sfile,"_tmp")  
      } 
      system(cdo_command)
    }
    
    mv_command<-paste("mv -n ",hyint_file,paste0(hyint_file,"_tmp"))
    etccdi_files_tmp<-paste(etccdi_files,"_tmp",sep="",collapse=" ")
    print(paste0("HyInt: merging ",length(etccdi_files)," ETCCDI files"))
    cdo_command<-paste0("cdo -O merge -sellonlatbox,-180,180,-90,90 ",
			 paste0(hyint_file,"_tmp "),etccdi_files_tmp," ",hyint_file)
    rm_command<-paste("rm ",etccdi_files_tmp)
    print(mv_command)
    print(cdo_command)
    print(rm_command)
    system(mv_command)
    system(cdo_command)
    system(rm_command)
    
return(0)
}


