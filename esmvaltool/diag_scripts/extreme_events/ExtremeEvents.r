# #############################################################################
# ExtremeEvents.r
#
# Authors: Björn Brötz (DLR, Germany)
#          Marit Sandstad (CICERO, Norway)
#          Christian W. Mohr (CICERO, Norway)
# #############################################################################
# Description
#    Calculate extreme events with plotting functionality
#
# Modification history
#    20181006-A_cwmohr : observation read and sorting fixes
#    20181003-A_cwmohr : correcting r.interface output for observation data. 
#    20180725-A_cwmohr : modification of timeseries_main() and climdex selection
#    20180615-A_cwmohr : more clean up of code
#    20180131-A_laue_ax: clean-up of code, adaptation to ESMValTool standards
#                        added tagging support
#    20170920-A_sand_ma: modification to include plotting
#    20160414-A_broe_bj: written
# ############################################################################
library(climdex.pcic.ncdf)
library(tools)
library(yaml)
library(ncdf4)
library(ncdf4.helpers)

cdo <- function(command, args = "", input = "", options = "", output = "",
                stdout = "", noout = F) {
  if (args != "") args <- paste0(",", args)
  if (stdout != "") {
      stdout <- paste0(" > '", stdout, "'")
      noout <- T
  }
  if (input[1] != "") {
    for (i in 1:length(input)) {
      input[i] <- paste0("'", input[i], "'")
    }
    input <- paste(input, collapse = " ")
  }
  output0 <- output
  if (output != "") {
    output <- paste0("'", output, "'")
  } else if ( !noout ) {
    output <- tempfile()
    output0 <- output
  }
  argstr <- paste0(options, " ", command, args, " ", input, " ", output,
                   " ", stdout)
  print(paste("cdo", argstr))
  ret <- system2("cdo", args = argstr)
  if (ret != 0) {
    stop(paste("Failed (", ret, "): cdo", argstr))
  }
  return(output0)
}

nco <- function(cmd, argstr) {
  ret <- system2(cmd, args = argstr)
  if (ret != 0) {
    stop(paste("Failed (", ret, "): ", cmd, " ", argstr))
  }
}

diag_scripts_dir <- Sys.getenv("diag_scripts")
source(paste0(diag_scripts_dir, "/extreme_events/cfg_climdex.r"))
source(paste0(diag_scripts_dir, "/extreme_events/cfg_extreme.r"))
source(paste0(diag_scripts_dir, "/extreme_events/common_climdex_preprocessing_for_plots.r"))
source(paste0(diag_scripts_dir, "/extreme_events/make_timeseries_plot.r"))
source(paste0(diag_scripts_dir, "/extreme_events/make_Glecker_plot2.r"))

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
reference_model <- unname(sapply(list0, "[[", "reference_dataset"))[,1]
diag_base <- unname(sapply(list0, "[[", "diagnostic"))[1]
#### Correction r.interface output correction ####
models_experiment[models_experiment == "No_value"] <- "No-values"

variables <- c()
climofiles <- c()
models <- c()
metadata <- c()

# loop over variables
for (i in 1:length(settings$input_files)) {
    metadata <- yaml::read_yaml(settings$input_files[i])
    models_name <- unname(sapply(metadata, "[[", "dataset"))
    short_name  <- unname(sapply(metadata, "[[", "short_name"))
    variables <- c(variables, short_name)
    models <- c(models, models_name)
    climofiles <- c(climofiles, names(metadata))
}

field_type0 <- "T2Ds"
# associated to first climofile
print(paste(diag_base, ": starting routine"))

# create working dirs if they do not exist
work_dir <- settings$work_dir
regridding_dir <- settings$run_dir
plot_dir <- settings$plot_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir, recursive = T, showWarnings = F)
dir.create(plot_dir, recursive = T, showWarnings = F)

# setup provenance file and list
provenance_file <- paste0(regridding_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

## Find earlier climdex indices in work folder
climdex_files <- list.files(path = work_dir, pattern = "ETCCDI")

# Fix input files removing bounds
print("Removing bounds from preprocessed files")
for (i in 1:length(climofiles)) {
   tmp <- tempfile()
   nco("ncks", paste("-C -O -x -v lat_bnds,lon_bnds,time_bnds",
       climofiles[i], tmp))
   nco("ncatted", paste("-O -a bounds,time,d,,", tmp))
   nco("ncatted", paste("-O -a bounds,lat,d,,", tmp))
   nco("ncatted", paste("-O -a bounds,lon,d,,", tmp))
   nco("ncatted", paste0("-O -a coordinates,", variables[i], ",d,, ", tmp))
   file.copy(tmp, climofiles[i], overwrite = TRUE)
   unlink(tmp)
}

##
## At this stage climdex indices are calculated. This process is extremely tedious and check points are in place to check whether the indicese are already produced. 
## If the climdex files are there, then this process is skipped. Delete the climdex files from the work folder if you wish to have the climdex indices recalculated.
##
for (model_idx in c(1:length(models_name))) {

    author.data <- list(institution = "None", institution_id = "None")
    template <- paste("var_timeres_", models_name[model_idx], "_",
                      models_experiment[model_idx], "_",
                      models_ensemble[model_idx], "_",
                      models_start_year[model_idx],
                      "01-", models_end_year[model_idx],
                      "12.nc", sep = "", collapse = "")
    print("")
    print(paste0(">>>>>>>> Template name: ", template))
    print("")

    idx_select <- unique(c(timeseries_idx, gleckler_idx))
#    climdex.idx.subset  <- unique(idx_df$idxETCCDI[which(idx_df$idxETCCDI_time %in% idx_select)])
#    print("")
#    print(paste0(">>>>>>>> Indices required: ", paste(climdex.idx.subset, collapse = ", ")))
#    print("")

     base.period <- c(max(strtoi(models_start_year)),
                      min(strtoi(models_end_year)))

    ## Check point for existing files
    climdex_file_check <- paste0(idx_select, "_",
                                 models_name[model_idx], "_",
                                 models_experiment[model_idx], "_",
                                 models_ensemble[model_idx], "_",
                                 models_start_year[model_idx], "-",
                                 models_end_year[model_idx])
    check_control <- vector("logical", length(climdex_file_check))
    n <- 0
    for (chck in climdex_file_check) {
      n  <- n + 1
      #print(grep(chck, climdex_files))
      tmp <- length(grep(chck, climdex_files))
      check_control[n] <- (tmp > 0)
    }
    print(check_control)

    if (!all(check_control)) {
      #print(fullpath_filenames)
      print("")
      print(paste0(">>>>>>>> Producing Indices for ", models_name[model_idx]))
      print(climofiles[models == models_name[model_idx]])
      print("")
      create.indices.from.files(climofiles[models == models_name[model_idx]],
                                work_dir, template, author.data,
                                base.range = base.period, parallel = 25,
                                verbose = TRUE, max.vals.millions = 20)
    }
}

if (write_plots) {
############################## 
## A climdex prcessing section is needed here for observation data.
## CMORized observation data found in the obs directory, has it's climdex indices calculated,
## which are then placed in the work/ExtremeEvents directory
############################## 

## Splitting models from observations

#obs_name <- models_name[models_project == "OBS"]
#models_name <- models_name[models_project != "OBS"]

################################### 
#### Produce time series plots ####
################################### 

# These are forced here for testing

  print("------ Model datasets ------")
  print(setdiff(models_name, reference_model))
  print("---- Reference datasets ----")
  print(reference_model)
  print("----------------------------")
  if (chk.ts_plt){
    print("")
    print(paste0(">>>>>>>> TIME SERIE PROCESSING INITIATION"))
    timeseries_main(path = work_dir, idx_list =  timeseries_idx,
                    model_list = setdiff(models_name, reference_model),
                    obs_list = reference_model, plot_dir = plot_dir,
                    normalize = normalize)
  }

############################### 
#### Produce Gleckler plot ####
############################### 
  if (chk.glc_plt) {
    print("")
    print(paste0(">>>>>>>> GLECKLER PROCESSING INITIATION"))

    ## Check if Gleckler Array already exists
    nidx <- length(gleckler_idx) # number of indices
    nmodel <- length(models_name) # number of models
    nobs <- length(reference_model) #number of observations
    ArrayName <- paste0("Gleclker-Array_", nidx, "-idx_",
                        nmodel, "-models_", nobs, "-obs", ".RDS")
    ArrayDirName <- paste0(plot_dir, "/", diag_base, "/", ArrayName)
    if (chk.glc_arr) {
      if (file.exists(ArrayDirName)) {
        file.remove(ArrayDirName)
      }
      promptInput <- "y"
    }

    if (file.exists(ArrayDirName)) {
      promptInput <- "n"
    } else {
      promptInput <- "y"
    }

    #### Running gleckler_main ####
    gleckler_main(path = work_dir, idx_list = gleckler_idx,
                  model_list = setdiff(models_name, reference_model),
                  obs_list = reference_model,
                  plot_dir = plot_dir, promptInput = promptInput)
  }
}
