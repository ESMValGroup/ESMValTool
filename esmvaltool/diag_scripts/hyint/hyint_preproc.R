######################################################
#---------Regridding preprocessing for HyInt---------#
#-------------E. Arnone (Oct 2017)-------------------#
######################################################

hyint_preproc <- function(work_dir, model_idx, climofile, regfile) {
  print(paste0(diag_base, ": pre-processing file: ", climofile))

  #  add absolute axis, remove leap year days
  # cdo delete and copy do not like files with whitespace
  tempf <- cdo("addc", args="0", input = climofile)
  cdo("delete", options = "-L -f nc -a", args = "month=2,day=29",
      input = tempf, output = regfile)
  unlink(tempf)

  # generate grid file
  gridfile <- getfilename_indices(work_dir, diag_base, model_idx, grid = T)
  cdo("griddes", input = regfile, stdout = gridfile, noout = T)

  print(paste0(diag_base, ": pre-processed file: ", regfile))

  return(0)
}
