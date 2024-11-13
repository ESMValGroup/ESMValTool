# #############################################################################
# common_climdex_preprocessing.R
#
# Author: Marit Sandstad (CICERO, Norway)
#       : Christian Wilhelm Mohr (CICERO, Norway)
#
#
# #############################################################################
# Description
#    Common codes to preprocsess climdex files from multiple
#    sources for plotting. This includes creating a common grid
#    Cropping files to the same time span. Regridding, landseamasking
#    and producing timemeans.
#
# Modification history
#    20190506-vonhardenberg_jost:    conversion to ESMValTool2
#    20180725-mohr_christianwilhelm: modification of setTimeForFilesEqual()
#                                    function
#    2017 0920-sandstad_marit:       creation
#
# #############################################################################

##
##
## Method to create an ascii grid file for
## to use to regrid on
## @param idx_dir path of directory containing
## files from which to create the grid
##
create_grid <- function(path = idx_dir, loc = "./gridDef") {
  ## Picking the grid found in the first file to regrid over
  first_file <- list.files(path,
    pattern = paste0(".*", regrid_dataset, ".*\\.nc"),
    full.names = TRUE
  )[1]
  cdo(
    "griddes -delvar,time_bnds",
    input = first_file,
    stdout = loc,
    options = "-s -O"
  )
}

#
# Method to create a landSeaMask on a suitable grid
# @param regrid name w/path of gridfile to use
# to put the landdseamask on
#
create_land_sea_mask <- function(regrid = "./gridDef",
                                 loc = "./",
                                 landmask = "./landSeaMask.nc") {
  # Test if gridfile exists
  # otherwise call function to generate one
  if (!file.exists(regrid)) {
    create_grid(path = loc, loc = regrid)
  }

  ## Making topographic map
  topof <- cdo("topo", options = "-O -f nc")

  ## Regridding the topographic map to chosen grid
  rtopof <-
    cdo("remapcon",
      args = regrid,
      input = topof,
      options = "-O"
    )

  # Set above sea-level gridpoints to missing
  rtopomissf <- cdo("setrtomiss",
    args = "0,9000",
    input = rtopof,
    options = "-O"
  )

  # Set above sea-level gridpoints to 1
  rtopo1posf <- cdo("setmisstoc",
    args = "1",
    input = rtopomissf,
    options = "-O"
  )

  # Set below sea-level gridpoints to missing
  cdo(
    "setrtomiss",
    args = "-9000,0",
    input = rtopo1posf,
    output = landmask,
    options = "-O"
  )

  unlink(c(topof, rtopof, rtopomissf, rtopo1posf))
}

##
## Method crop all index files for a single index
## to the same time period.
## The smallest common time period is chosen
## @param path gives path to location of index files
## @param idx lists the index under consideration
## @param model_list provides the list of selected models for time cropping
## @param time_cropped is the directory to put the time cropped files
## @param max_start is an optional crop start
## @param min_end is an optional crop end
##
set_time_for_files_equal <- function(path, # nolint
                                     idx,
                                     model_list,
                                     time_cropped = "./timeCropped",
                                     max_start = 0,
                                     min_end = 2500) {
  ## Getting a list of all the files for the index
  models_avail <- basename(Sys.glob(file.path(
    path,
    paste(idx, "*.nc", sep = "")
  )))

  ## Selecting only the files from the model list
  models <- vector(mode = "character", length = length(model_list))
  for (i in seq_along(model_list)) {
    models[i] <- models_avail[grep(
      pattern = model_list[i],
      x = models_avail
    )]
  }

  print(models)

  ## Checking if the folder exists and making it if not
  print(time_cropped)
  if (!file.exists(time_cropped)) {
    dir.create(time_cropped)
  }

  ## Arrays to record orginal start and end years
  start <- integer(length(models))
  end <- integer(length(models))

  i <- 1
  # For-loop to find the minimum time interval
  # so we can crop all files to this time interval
  m <- models[1]
  for (m in models) {
    start[i] <- strtoi(substr(m, nchar(m) - 11, nchar(m) - 8))
    end[i] <- strtoi(substr(m, nchar(m) - 6, nchar(m) - 3))

    if (start[i] > max_start) {
      max_start <- start[i]
    }

    if (end[i] < min_end) {
      min_end <- end[i]
    }
    i <- i + 1
  }
  if (max_start >= min_end) {
    print("No time overlap for files")
    print(c(max_start, min_end))
    for (m in models) {
      file.copy(paste0(path, "/", m), paste0(time_cropped, "/", m))
    }
    return(c(max_start, min_end))
  }

  i <- 1
  # For-loop to crop the files
  for (m in models) {
    ## If file is already of appropriate length
    ## Then just copy it over
    if (start[i] == max_start && end[i] == min_end) {
      file.copy(paste0(path, "/", m), paste0(time_cropped, "/", m))
      ## Otherwise do the time cropping
    } else {
      beg <- max_start - start[i]
      sto <- min_end - max_start + beg
      newname <- paste(substr(m, 1, nchar(m) - 12),
        max_start,
        "-",
        min_end,
        ".nc",
        sep = ""
      )
      nco(
        "ncks",
        paste0(
          "-d time,",
          beg,
          ",",
          sto,
          " ",
          path,
          "/",
          m,
          " ",
          time_cropped,
          "/",
          newname
        )
      )
    }
    i <- i + 1
  }
  return(c(max_start, min_end))
}

##
##
## Method that regrids and landseamasks a file
## Timemeaned versions are also produced
## @param idx_raw gives the full path name of the file
## @param regrid gives the file of the grid to regrid on
## @param landmask gives the file that defines the landseamask to be used
##
##
regrid_and_land_sea_mask <- function(idx_raw,
                                     regrid = "./gridDef",
                                     landmask = "./landSeaMask.nc",
                                     regridded = "./Regridded",
                                     land = "./Land",
                                     loc = "./") {
  ## Getting just the raw name of the file
  idx_name <- basename(idx_raw)

  ## If the landmask does not exist, we create one.
  if (!file.exists(landmask)) {
    create_land_sea_mask(
      regrid = regrid,
      loc = loc,
      landmask = landmask
    )
  }

  ## Checking if directories are present and creating them if not:
  if (!dir.exists(regridded)) {
    dir.create(regridded)
  }
  if (!dir.exists(land)) {
    dir.create(land)
  }

  ## Regridding file:
  varname <- strsplit(idx_name, "_")[[1]][1]
  tmpsel <-
    cdo("selvar",
      args = varname,
      input = idx_raw,
      options = "-O"
    )
  cdo(
    "remapcon",
    args = regrid,
    input = tmpsel,
    output = paste0(regridded, "/", idx_name),
    options = "-O"
  )
  unlink(tmpsel)

  ## Applying landseamask:
  cdo(
    "div",
    input = c(paste0(regridded, "/", idx_name), landmask),
    output = paste0(land, "/", idx_name),
    options = "-O"
  )

  ## Also produce timemean:
  ## !! check is this should be subject to some reference period or
  ## time change
  cdo(
    "timmean",
    input = paste0(land, "/", idx_name),
    output = paste0(land, "/tm_", idx_name),
    options = "-O" # nolint
  )
}
