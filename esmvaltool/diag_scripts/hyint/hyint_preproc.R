######################################################
#---------Regridding preprocessing for HyInt---------#
#-------------E. Arnone (Oct 2017)-------------------#
######################################################

hyint.preproc<-function(work_dir,model_idx,ref_idx,climofile,regfile,rgrid) {

#  absolute axis, remove leap year days, select lonlat box convert to NetCDF4

  # Manage gridding options
  sgrid<-""
  if (rgrid !=F) {sgrid <- paste0("-remapcon2,", rgrid)}

  # regrid/reformat climo file
  print(paste0(diag_base,": pre-processing file: ", climofile))

  cdo_command<-paste(paste0("cdo -L -f nc -a -delete,month=2,day=29"), sgrid, climofile, paste0(regfile,"regtmp"))
  print(cdo_command)
  system(cdo_command)

  cdo_command2<-paste("mv ",paste0(regfile,"regtmp")," ",paste0(regfile,"regtmp1"))
  if (rlonlatdata[1]!=F) {
    cdo_command2<-paste0("cdo sellonlatbox,",paste(rlonlatdata,sep="",collapse=",")," ",paste0(regfile,"regtmp")," ",paste0(regfile,"regtmp1"))
  }
  print(cdo_command2)
  system(cdo_command2)

  cdo_command3<-paste("cdo -f nc4 -copy ", paste0(regfile,"regtmp1"), regfile)
  print(cdo_command3)
  system(cdo_command3)

  rm_command<-paste("rm ", paste0(regfile,"regtmp*"))
  system(rm_command)

  # generate grid file
  gridfile<-getfilename.indices(work_dir,diag_base,model_idx,grid=T)
  grid_command<-paste("cdo griddes ",regfile," > ",gridfile)  
  system(grid_command)

  print(paste0(diag_base,": pre-processed file: ", regfile))

return(0)
}
