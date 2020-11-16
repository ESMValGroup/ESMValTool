# #############################################################################
# extreme_events.R
#
# Authors: Björn Brötz (DLR, Germany)
#          Marit Sandstad (CICERO, Norway)
#          Christian W. Mohr (CICERO, Norway)
# #############################################################################
# Description
#    Calculate extreme events with plotting functionality
#
# Modification history
#    20190506-vonhardenberg_jost: conversion to ESMValTool2
#    20181006-mohr_christianwilhelm: observation read and sorting fixes
#    20181003-mohr_christianwilhelm: correcting r.interface output for
#                                    observation data.
#    20180725-mohr_christianwilhelm: modification of timeseries_main() and
#                                    climdex selection
#    20180615-mohr_christianwilhelm: more clean up of code
#    20180131-lauer_axel: clean-up of code, adaptation to ESMValTool standards
#                         added tagging support
#    20170920-sandstad_marit: modification to include plotting
#    2016 0414-A_broetz_bjoern: written
# #############################################################################

library(tools)
library(yaml)
library(ncdf4)
library(ncdf4.helpers)
library(scales)
library(RColorBrewer) # nolint

# function to flatten nested lists
flatten_lists <- function(x) {
  if (!inherits(x, "list")) return(x)
  else return(unlist(c(lapply(x, flatten_lists)), recursive = FALSE))
}

provenance_record <- function(infile) {
  xprov <- list(
    ancestors = flatten_lists(as.list(infile)),
    authors = list(
      "broetz_bjoern",
      "sandstad_marit",
      "mohr_christianwilhelm",
      "vonhardenberg_jost"
    ),
    references = list("zhang11wcc"),
    projects = list("crescendo", "c3s-magic"),
    caption = "Extreme events indices",
    statistics = list("other"),
    realms = list("atmos"),
    themes = list("phys"),
    domains = list("global")
  )
  return(xprov)
}

diag_scripts_dir <- Sys.getenv("diag_scripts")
climdex_src <-
  paste0(
    diag_scripts_dir,
    "/extreme_events/climdex.pcic.ncdf/R/ncdf.R"
  ) # nolint
source(paste0(
  diag_scripts_dir,
  "/extreme_events/climdex.pcic.ncdf/R/ncdf.R"
)) # nolint
source(paste0(diag_scripts_dir, "/shared/external.R")) # nolint
source(paste0(diag_scripts_dir, "/extreme_events/cfg_climdex.R")) # nolint
source(paste0(diag_scripts_dir, "/extreme_events/cfg_extreme.R")) # nolint
source(
  paste0(
    diag_scripts_dir,
    "/extreme_events/common_climdex_preprocessing_for_plots.R"
  )
) # nolint
source(paste0(
  diag_scripts_dir,
  "/extreme_events/make_timeseries_plot.R"
)) # nolint
source(paste0(
  diag_scripts_dir,
  "/extreme_events/make_glecker_plot.R"
)) # nolint

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
for (myname in names(settings)) {
  temp <- get(myname, settings)
  assign(myname, temp)
}

list0 <- yaml::read_yaml(settings$input_files[1])
# extract metadata
models_name <- unname(sapply(list0, "[[", "dataset"))
models_ensemble <- unname(sapply(list0, "[[", "ensemble"))
models_start_year <- unname(sapply(list0, "[[", "start_year"))
models_end_year <- unname(sapply(list0, "[[", "end_year"))
models_experiment <- unname(sapply(list0, "[[", "exp"))
models_project <- unname(sapply(list0, "[[", "project"))
diag_base <- unname(sapply(list0, "[[", "diagnostic"))[1]
#### Correction r.interface output correction ####
models_experiment[models_experiment == "No_value"] <- "No-values"

variables <- c()
climofiles <- c()
models <- c()
metadata <- c()

# loop over variables
for (i in seq_along(settings$input_files)) {
  metadata <- yaml::read_yaml(settings$input_files[i])
  models_name <- unname(sapply(metadata, "[[", "dataset"))
  short_name <- unname(sapply(metadata, "[[", "short_name"))
  variables <- c(variables, short_name)
  models <- c(models, models_name)
  climofiles <- c(climofiles, names(metadata))
}

# associated to first climofile
print(paste(diag_base, ": starting routine"))

# create working dirs if they do not exist
work_dir <- settings$work_dir
regridding_dir <- settings$run_dir
plot_dir <- settings$plot_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir,
  recursive = T,
  showWarnings = F
)
dir.create(plot_dir, recursive = T, showWarnings = F)

# setup provenance file and list
provenance_file <-
  paste0(regridding_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

if (anyNA(base_range)) {
  stop("Please choose a base_range!")
}
model_range <- c(
  max(strtoi(models_start_year)),
  min(strtoi(models_end_year))
)
if ((base_range[1] < max(strtoi(models_start_year))) |
  (base_range[2] > min(strtoi(models_end_year)))) {
  stop(
    paste(
      "Base range",
      base_range[1],
      "-",
      base_range[2],
      "outside available model data period",
      model_range[1],
      "-",
      model_range[2]
    )
  )
}
print(paste("Base range:", base_range[1], "-", base_range[2]))

if (anyNA(regrid_dataset)) {
  regrid_dataset <- reference_datasets[1]
  print(paste(
    "Regrid dataset not set, choosing first reference dataset:",
    regrid_dataset
  ))
}

## Find earlier climdex indices in work folder
climdex_files <- list.files(path = work_dir, pattern = "ETCCDI")

# Fix input files removing bounds
print("Removing bounds from preprocessed files")
for (i in seq_along(climofiles)) {
  tmp <- tempfile()
  nco(
    "ncks",
    paste(
      "-C -O -x -v lat_bnds,lon_bnds,time_bnds",
      climofiles[i], tmp
    )
  )
  nco("ncatted", paste("-O -a bounds,time,d,,", tmp))
  nco("ncatted", paste("-O -a bounds,lat,d,,", tmp))
  nco("ncatted", paste("-O -a bounds,lon,d,,", tmp))
  nco(
    "ncatted",
    paste0("-O -a coordinates,", variables[i], ",d,, ", tmp)
  )
  file.copy(tmp, climofiles[i], overwrite = TRUE)
  unlink(tmp)
}

##
## At this stage climdex indices are calculated. This process is extremely
## tedious and check points are in place to check whether the indicese are
## already produced. If the climdex files are there, then this process is
## skipped. Delete the climdex files from the work folder if you wish to have
## the climdex indices recalculated.
##
for (model_idx in c(1:length(models_name))) { # nolint
  author.data <- list(institution = "None", institution_id = "None")
  template <- paste(
    "var_timeres_",
    models_name[model_idx],
    "_",
    models_experiment[model_idx],
    "_",
    models_ensemble[model_idx],
    "_",
    models_start_year[model_idx],
    "01-",
    models_end_year[model_idx],
    "12.nc",
    sep = "",
    collapse = ""
  )
  print("")
  print(paste0(">>>>>>>> Template name: ", template))
  print("")

  idx_select <- unique(c(timeseries_idx, gleckler_idx))

  ## Check point for existing files
  climdex_file_check <- c()
  for (idx in idx_select) {
    if (grepl("mon", idx)) {
      climdex_file_check <- c(
        climdex_file_check,
        paste0(
          idx,
          "_",
          models_name[model_idx],
          "_",
          models_experiment[model_idx],
          "_",
          models_ensemble[model_idx],
          "_",
          models_start_year[model_idx],
          "01-",
          models_end_year[model_idx],
          "12.nc"
        )
      )
    } else {
      climdex_file_check <- c(
        climdex_file_check,
        paste0(
          idx,
          "_",
          models_name[model_idx],
          "_",
          models_experiment[model_idx],
          "_",
          models_ensemble[model_idx],
          "_",
          models_start_year[model_idx],
          "-",
          models_end_year[model_idx],
          ".nc"
        )
      )
    }
  }
  check_control <- vector("logical", length(climdex_file_check))
  n <- 0
  for (chck in climdex_file_check) {
    n <- n + 1
    tmp <- length(grep(chck, climdex_files))
    check_control[n] <- (tmp > 0)
  }

  if (!any(grepl("yr", idx_select))) {
    timeres <- "mon"
    write_plots <- FALSE
  } else if (!any(grepl("mon", idx_select))) {
    timeres <- "annual"
  } else {
    timeres <- "all"
    write_plots <- FALSE
  }

  if (!all(check_control)) {
    print("")
    print(paste0(">>>>>>> Producing Indices for ", models_name[model_idx]))
    print(climofiles[models == models_name[model_idx]])
    print("")
    infiles <- climofiles[models == models_name[model_idx]]
    indices <- sub("ETCCDI.*", "", idx_select)
    # Find best chunk size
    chunk <- 10
    if (!(is.logical(climdex_parallel))) {
      nc <- nc_open(infiles[1])
      chunk <-
        floor(
          (nc$dim$time$len * nc$dim$lon$len * nc$dim$lat$len +
            1000.0) / (climdex_parallel * 1000000)
        )
      chunk <- max(min(100, chunk), 1)
      nc_close(nc)
      print(paste("Chunk size:", chunk))
    }
    create.indices.from.files(
      infiles,
      # nolint
      work_dir,
      template,
      author.data,
      base.range = base_range,
      parallel = climdex_parallel,
      verbose = TRUE,
      climdex.vars.subset = indices,
      climdex.time.resolution = timeres,
      max.vals.millions = chunk,
      src = climdex_src
    )

    # Set provenance for output files
    # Get new list of files after computation
    infiles <- climofiles[models == models_name[model_idx]]
    print("Computing xprov")
    xprov <- provenance_record(infiles)
    climdex_files <- list.files(
      path = work_dir,
      pattern = paste0("ETCCDI.*", models_name[model_idx], ".*\\.nc"),
      full.names = TRUE
    )
    for (fname in climdex_files) {
      print(paste("Provenance for ", fname))
      provenance[[fname]] <- xprov
    }
  }
}

if (write_plots) { # nolint
  #############################
  # A climdex processing section is needed here for observation data.
  # CMORized observation data found in the obs directory,
  # has it's climdex indices calculated,
  # which are then placed in the work/extreme_events directory
  #############################

  ## Splitting models from observations

  ###################################
  #### Produce time series plots ####
  ###################################

  if (anyNA(analysis_range)) {
    analysis_range[1] <- max(strtoi(models_start_year))
    analysis_range[2] <- min(strtoi(models_end_year))
    print(
      paste(
        "Analysis range not defined, assigning model range:",
        analysis_range[1],
        "-",
        analysis_range[2]
      )
    )
  }
  if ((analysis_range[1] < max(strtoi(models_start_year))) |
    (analysis_range[2] > min(strtoi(models_end_year)))) {
    stop(
      paste(
        "Analysis range",
        analysis_range[1],
        "-",
        analysis_range[2],
        "outside available model data period",
        model_range[1],
        "-",
        model_range[2]
      )
    )
  }
  print(paste("Analysis range:", analysis_range[1], "-", analysis_range[2]))

  # These are forced here for testing

  print("------ Model datasets ------")
  print(setdiff(models_name, reference_datasets))
  print("---- Reference datasets ----")
  print(reference_datasets)
  print("----------------------------")
  if (ts_plt) {
    print("")
    print(paste0(">>>>>>>> TIME SERIES PROCESSING INITIATION"))
    plotfiles <- timeseries_main(
      path = work_dir,
      idx_list = timeseries_idx,
      model_list = setdiff(models_name, reference_datasets),
      obs_list = reference_datasets,
      plot_dir = plot_dir,
      normalize = normalize,
      start_yr = analysis_range[1],
      end_yr = analysis_range[2]
    )
    xprov <- provenance_record(list(climofiles))
    for (fname in plotfiles) {
      provenance[[fname]] <- xprov
    }
    # Each timeseries file gets provenance from its reference dataset
    for (model in reference_datasets) {
      ncfiles <- list.files(
        file.path(work_dir, "timeseries"),
        pattern = model,
        full.names = TRUE
      )
      xprov <- provenance_record(list(climofiles[models == model]))
      for (fname in ncfiles) {
        provenance[[fname]] <- xprov
      }
    }
    # The ensemble timeseries get provenance from all model datasets
    ncfiles <- list.files(file.path(work_dir, "timeseries"),
      pattern = "ETCCDI.*ens",
      full.names = TRUE
    )

    ancestors <- sapply(setdiff(models_name, reference_datasets),
      grep,
      climofiles,
      value = TRUE
    )
    xprov <- provenance_record(ancestors)
    for (fname in ncfiles) {
      provenance[[fname]] <- xprov
    }
  }

  ###############################
  #### Produce Gleckler plot ####
  ###############################
  if (glc_plt) {
    print("")
    print(paste0(">>>>>>>> GLECKLER PROCESSING INITIATION"))

    ## Check if Gleckler Array already exists
    nidx <- length(gleckler_idx) # number of indices
    nmodel <- length(models_name) # number of models
    nobs <- length(reference_datasets) # number of observations
    arrayname <- paste0(
      "Gleckler-Array_",
      nidx,
      "-idx_",
      nmodel,
      "-models_",
      nobs,
      "-obs",
      ".RDS"
    )
    arraydirname <- paste0(plot_dir, "/", diag_base, "/", arrayname)
    if (glc_arr) {
      if (file.exists(arraydirname)) {
        file.remove(arraydirname)
      }
      promptinput <- "y"
    }

    if (file.exists(arraydirname)) {
      promptinput <- "n"
    } else {
      promptinput <- "y"
    }

    #### Running gleckler_main ####
    plotfiles <- gleckler_main(
      path = work_dir,
      idx_list = gleckler_idx,
      model_list = setdiff(models_name, reference_datasets),
      obs_list = reference_datasets,
      plot_dir = plot_dir,
      promptinput = promptinput,
      start_yr = analysis_range[1],
      end_yr = analysis_range[2]
    )

    xprov <- provenance_record(list(climofiles))
    for (fname in plotfiles) {
      provenance[[fname]] <- xprov
    }
    ncfiles <- list.files(file.path(work_dir, "gleckler/Gleck*"))
    xprov <- provenance_record(list(climofiles))
    for (fname in ncfiles) {
      provenance[[fname]] <- xprov
    }
  }
}

# Write provenance to file
write_yaml(provenance, provenance_file)
