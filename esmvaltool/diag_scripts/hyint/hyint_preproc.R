######################################################
#---------Regridding preprocessing for HyInt---------#
#-------------E. Arnone (Oct 2017)-------------------#
######################################################

hyint_preproc <- function(work_dir,
                          model_idx,
                          ref_idx,
                          climofile,
                          regfile,
                          rgrid) {
  print(paste0(diag_base, ": pre-processing file: ", climofile))

  # add absolute axis, remove leap year days, regrid if needed
  # cdo delete and copy do not like files with whitespace

  if (rgrid != F) {
    if (rgrid == "REF") {
      rgrid <- climofiles[ref_idx]
      gridf <- tempfile()
      cdo("griddes", input = rgrid, stdout = gridf)
    } else {
      gridf <- rgrid
    }
    tempf <- cdo("remapcon2", args = gridf, input = climofile)
    unlink(gridf)
  } else {
    tempf <- cdo("addc", args = "0", input = climofile)
  }

  cdo("-copy",
    options = "-L -f nc -a",
    input = tempf,
    output = regfile
  )

  unlink(tempf)

  # generate grid file
  gridfile <-
    getfilename_indices(work_dir, diag_base, model_idx, grid = T)
  cdo("griddes", input = regfile, stdout = gridfile)

  print(paste0(diag_base, ": pre-processed file: ", regfile))

  return(0)
}
