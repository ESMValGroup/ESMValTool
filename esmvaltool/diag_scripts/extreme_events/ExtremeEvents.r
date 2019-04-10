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
library(ncdf4.helpers)


#source('interface_data/r.interface')

#source('../backup_scripts/r.interface_test')

source(diag_script_cfg)
source('diag_scripts/lib/R/info_output.r')
source('diag_scripts/lib/R/meta_data.r')

## Load in climdex index dataframe
source('nml/cfg_ExtremeEvents/cfg_climdex.r')

## Script for pre-processing climdex data
source('diag_scripts/aux/ExtremeEvents/common_climdex_preprocessing_for_plots.r')

## Scripts for plotting
source('diag_scripts/aux/ExtremeEvents/make_timeseries_plot.r')
source('diag_scripts/aux/ExtremeEvents/make_Glecker_plot2.r')




info_output(paste0("<<<<<<<< Entering ", diag_script), verbosity, 4)
info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
info_output(diag_script, verbosity, 1)
info_output(paste0("var: ", variables), verbosity, 1)
info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

library(tools)
diag_base = file_path_sans_ext(diag_script)

## Create working dirs if they do not exist
out_dir <- paste(work_dir, "/", diag_base, sep = "")
dir.create(out_dir, showWarnings = FALSE)
dir.create(plot_dir, showWarnings = FALSE)
dir.create(file.path(plot_dir, diag_base), showWarnings = FALSE)
dir.create(climo_dir, showWarnings = FALSE)


#### Correction r.interface output correction ####
models_experiment[models_experiment == "No_value"] <- "No-values"

## Find earlier climdex indices in work folder
climdex_files <- list.files(path = out_dir, pattern = "ETCCDI")

##
## At this stage climdex indices are calculated. This process is extremely tedious and check points are in place to check whether the indicese are already produced. 
## If the climdex files are there, then this process is skipped. Delete the climdex files from the work folder if you wish to have the climdex indices recalculated.
##
for (model_idx in c(1:length(models_name))) {
    fullpath_filenames <- character(0)
    for (var_idx in c(1:length(variables))) {
        fullpath_filenames[var_idx] <- interface_get_fullpath(variables[var_idx], field_types[var_idx], model_idx)
        print(fullpath_filenames[var_idx])
    }
    author.data <- list(institution = "None", institution_id = "None")
    template <- paste("var_timeres_", models_name[model_idx], "_", models_experiment[model_idx], "_",
                models_ensemble[model_idx], "_", models_start_year[model_idx],
                "01-", models_end_year[model_idx], "12.nc", sep = "", collapse = "")
    print("")
    info_output(paste0(">>>>>>>> Template name: ", template), verbosity, 4)
    print("")
    idx_select <- unique(c(timeseries_idx, gleckler_idx))
    climdex.idx.subset  <- unique(idx_df$idxETCCDI[which(idx_df$idxETCCDI_time %in% idx_select)])
    print("")
    info_output(paste0(">>>>>>>> Indices required: ", paste(climdex.idx.subset, collapse = ", ")), verbosity, 4)
    print("")
    
    base.period <- c(max(strtoi(models_start_year)), min(strtoi(models_end_year)))
    
    ## Check point for existing files
    climdex_file_check <- paste(idx_select, "_",models_name[model_idx], "_", models_experiment[model_idx], "_",
                                models_ensemble[model_idx], "_", models_start_year[model_idx],"-", models_end_year[model_idx], sep="")
    #print(climdex_file_check)
    check_control <- vector("logical", length(climdex_file_check))
    n = 0
    for(chck in climdex_file_check){
      n  <- n + 1
	    #print(grep(chck, climdex_files))
	    tmp <- length(grep(chck, climdex_files))
      check_control[n] <- (tmp > 0)
    }
    #print(check_control)
    
    if(!all(check_control)){
      #print(fullpath_filenames)
      print("")
      info_output(paste0(">>>>>>>> Producing Indices for ", models_name[model_idx]), verbosity, 4)
      print("")
      create.indices.from.files(fullpath_filenames, out_dir, template, author.data, 
                                   base.range=base.period, parallel = 25, verbose = TRUE, max.vals.millions = 20) # Procuce selected indicesd
    }
}

############################## 
## A climdex prcessing section is needed here for observation data.
## CMORized observation data found in the obs directory, has it's climdex indices calculated,
## which are then placed in the work/ExtremeEvents directory
############################## 

## Splitting models from observations

obs_name <- models_name[models_project == "OBS"]
#obs_name <- c(obs_name, obs_name)
#obs_name <- c(obs_name, obs_name[1])
#obs_name <- obs_name[1]
models_name <- models_name[models_project != "OBS"]



################################### 
#### Produce time series plots ####
################################### 
if(chk.ts_plt){
  print("")
  print("")
  info_output(paste0(">>>>>>>> TIME SERIE PROCESSING INITIATION"), verbosity, 4)
  timeseries_main(path = out_dir, idx_list =  timeseries_idx, model_list = models_name , obs_list = obs_name, plot_dir = paste(plot_dir, "/", diag_base, sep = ""), normalize=normalize)
}


############################### 
#### Produce Gleckler plot ####
############################### 
if(chk.glc_plt){
  print("")
  print("")
  info_output(paste0(">>>>>>>> GLECKLER PROCESSING INITIATION"), verbosity, 4)
    
  ## Check if Gleckler Array already exists
  nidx <- length(gleckler_idx) # number of indices
  nmodel <- length(models_name) # number of models
  nobs <- length(obs_name) #number of observations
  ArrayName <- paste0("Gleclker-Array_", nidx, "-idx_", nmodel,"-models_", nobs, "-obs",  ".RDS")
  ArrayDirName <- paste0(plot_dir, "/", diag_base, "/", ArrayName)
  if(chk.glc_arr){
    if(file.exists(ArrayDirName)){file.remove(ArrayDirName)}
    promptInput <- "y"
  }
  
  if(file.exists(ArrayDirName)){
      promptInput <- "n"
  }else{promptInput <- "y"}
  
  #### Running gleckler_main ####
  gleckler_main(path = out_dir, idx_list = gleckler_idx,
                model_list = models_name, obs_list = obs_name, plot_dir = paste(plot_dir, "/", diag_base, sep = ""), promptInput=promptInput)
  
  print(cfgpar)
  info_output(paste0(">>>>>>>> Leaving ", diag_script), verbosity, 4)
}
