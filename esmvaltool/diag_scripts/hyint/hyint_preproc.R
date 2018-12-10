######################################################
#---------Regridding preprocessing for HyInt---------#
#-------------E. Arnone (Oct 2017)-------------------#
######################################################

hyint_preproc <- function(work_dir, model_idx, climofile, regfile) {
  print(paste0(diag_base, ": pre-processing file: ", climofile))

  #  add absolute axis, remove leap year days
  cdo_command <- paste(paste0("cdo -L -f nc -a -delete,month=2,day=29"),
                       climofile, regfile)
  system(cdo_command)

  # generate grid file
  gridfile <- getfilename_indices(work_dir, diag_base, model_idx, grid = T)
  grid_command <- paste("cdo griddes ", regfile, " > ", gridfile)
  system(grid_command)

  print(paste0(diag_base, ": pre-processed file: ", regfile))

  return(0)
}
