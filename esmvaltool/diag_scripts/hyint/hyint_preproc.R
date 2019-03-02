######################################################
#---------Regridding preprocessing for HyInt---------#
#-------------E. Arnone (Oct 2017)-------------------#
######################################################

hyint_preproc <- function(work_dir, model_idx, climofile, regfile) {
  print(paste0(diag_base, ": pre-processing file: ", climofile))

  #  add absolute axis, remove leap year days
  cdo("delete", options = "-L -f nc -a", args = "month=2,day=29",
      input = climofile, output = regfile)

  # generate grid file
  gridfile <- getfilename_indices(work_dir, diag_base, model_idx, grid = T)
  cdo("griddes", input = paste(regfile, ">", gridfile))

  print(paste0(diag_base, ": pre-processed file: ", regfile))

  return(0)
}
