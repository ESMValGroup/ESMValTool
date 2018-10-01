# #############################################################################
# hyint.r
# Authors:       E. Arnone (ISAC-CNR, Italy)
#                J. von Hardenberg (ISAC-CNR, Italy) 
# #############################################################################
# Description
# HyInt is a tool for calculation of the HY-INT index (Giorgi et al. 2011) 
# and additional hydroclimatic indices (Giorgi et al. 2014)
# which allow an estimate of the overall behaviour of the hydroclimatic cycle. 
# The tool calculates also timeseries and trends over selected regions and 
# produces a variety of types of plots including maps and timeseries. The timeseries/trend
# and plotting modules handle also ETCCDI indices data calculated with the climdex library through
# an ad hoc pre-processing.
#  
# Details
# The following indices are calculated based on input daily precipitation data:
# PRY = mean annual precipitation
# INT = mean annual precipitation intensity (intensity during wet days, or simple precipitation intensity index SDII)
# WSL = mean annual wet spell length (number of consecutive days during each wet spell)
# DSL = mean annual dry spell lenght (number of consecutive days during each dry spell)
# PA  = precipitation area (area over which of any given day precipitation occurs) 
# R95 = heavy precipitation index (percent of total precipitation above the 95% percentile of the reference distribution)
# HY-INT = hydroclimatic intensity. HY-INT = normalized(INT) x normalized(DSL).
#
# For EC-Earth data and then extended to any model and observational data, producing plots 
# of data vs. a reference dataset (e.g. ERA-INTERIM). Indices are normalized over a reference 
# period. Both absolute and normalized values are made available: users can select the indices 
# to be stored and plotted. The tool makes extensives use of the cfg_hyint configuration file
# for user selectable options and ease feeding needed inputs (e.g. region boundaries for timeseries 
# or value ranges and labels for figures).     
# 
# Required
# It reads daily precipitation data through ESMValTool. If requested, input precipitation data are pre-processed 
# interpolating on a common grid set by the user in the hyint_parameters file.  
# R libraries:"tools","PCICt","ncdf4","maps"
#
# Optional 
# Several options can be selected via the configuration file, e.g. the provision of an
# external normalization functions for the indices; a reference climatology for the R95 index; type of plots; etc. 
#
# Caveats
# Spatial data selection based on elevation works only with regridding at 320x160 (or by producing by hand grid files at needed resolution) 
#
# Modification history
#    20181001-A_arno_en: converted to latest v2.0
#    20180302-A_arno_en: converted to ESMValTool2  
#    20171206-A_arno_en: modularized version accepting climdex indices  
#    20171010-A_arno_en: modularized version  
#    20170901-A_arno_en: 1st github version  
#
# ############################################################################

library(tools)
library(yaml)
library(climdex.pcic.ncdf)

spath='esmvaltool/diag_scripts/hyint/'
source(paste0(spath,'hyint_functions.R'))
source(paste0(spath,'hyint_metadata.R'))
source(paste0(spath,'hyint_preproc.R'))
source(paste0(spath,'hyint_diagnostic.R'))
source(paste0(spath,'hyint_etccdi_preproc.R'))
source(paste0(spath,'hyint_trends.R'))
source(paste0(spath,'hyint_plot_maps.R'))
source(paste0(spath,'hyint_plot_trends.R'))
source(paste0(spath,'hyint_parameters.r'))

## Do not print warnings
#options(warn=0)var0 <- names(metadata)[1]

# Read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings_file <- args[1]
# settings_file="/work/users/arnone/esmvaltool_output/namelist_hyint_egu2018_1_20180322_222547/run/hyint/main/settings.yml"
# settings_file="/work/users/arnone/esmvaltool_output/namelist_hyint_egu2018_20180416_162130/run/hyint/main/settings.yml"
settings <- yaml::read_yaml(settings_file)
metadata <- yaml::read_yaml(settings$input_files)

# load data from settings 
for (myname in names(settings)) { temp=get(myname,settings); assign(myname,temp) }

# get first variable and list associated to pr variable
variables <- names(metadata)
var0 <- names(metadata)[1]
list0 <- get(var0,metadata) 

# get name of climofile for first variable and list associated to first climofile
climofiles <- names(list0)
climolist0 <- get(climofiles[1],list0)

diag_base = climolist0$diagnostic
print(paste(diag_base,": starting routine"))

if (length(etccdi_dir) != 1) {etcddi_dir = work_dir}
regridding_dir=run_dir
dir.create(plot_dir, recursive=T,  showWarnings = F)
dir.create(work_dir, recursive=T, showWarnings = F)
dir.create(regridding_dir, recursive=T, showWarnings = F)

models_name=unname(sapply(list0, '[[', 'dataset'))
reference_model=unname(sapply(list0, '[[', 'reference_dataset'))[1]
models_start_year=unname(sapply(list0, '[[', 'start_year'))
models_end_year=unname(sapply(list0, '[[', 'end_year'))
models_experiment=unname(sapply(list0, '[[', 'exp'))
models_ensemble=unname(sapply(list0, '[[', 'ensemble'))

# select reference model
ref_idx=which(models_name == reference_model) # select reference dataset; if not available, use last of list
if(length(ref_idx)==0) { ref_idx=length(models_name)}

# setup reference grid file if needed
if (rgrid == "REF") {
  climofile_ref <- climofiles[ref_idx]
  grid_file_ref = paste0(grid_file,ref_idx)
  # generate grid file
  grid_command<-paste("cdo griddes ",climofile_ref," > ",grid_file_ref)
  system(grid_command)
  print(paste(diag_base,": reference grid file ",grid_file_ref))
  # update rgrid tag with actual grid size from grid file
  rgrid = get.cdo.res(grid_file_ref) 
} 
 
## Run regridding and diagnostic
if (write_netcdf) {
  
  # loop through models 
  for (model_idx in c(1:(length(models_name)))) {

    # Setup filenames 
    climofile <- climofiles[model_idx] 
    sgrid <- "noregrid"; if (rgrid != F) {sgrid <- rgrid}
    regfile <- getfilename.regridded(regridding_dir,sgrid,var0,model_idx)

    # If needed, pre-process file regridding, selecting lon/lat region of interest and adding absolute time axis 
    if (run_regridding) {
      if(!file.exists(regfile) | force_regridding) { 
        dummy=hyint.preproc(work_dir,model_idx,ref_idx,climofile,regfile,rgrid)
      } else {   
        gridfile=getfilename.indices(work_dir,diag_base,model_idx,grid=T)
        grid_command<-paste("cdo griddes ",regfile," > ",gridfile)
        system(grid_command)
        print(paste0(diag_base,": data file exists: ", regfile))  
        print(paste0(diag_base,": corresponding grid: ", gridfile))  
      }
    }

    if (run_diagnostic) {    
      # Loop through seasons and call diagnostic
      for (seas in seasons) {
        
        hyint.diagnostic(work_dir,regfile,model_idx,seas,rewrite=force_diagnostic)
      }
    }
  }
}

## Preprocess ETCCDI input files and merge them with HyInt indices
if (write_netcdf & etccdi_preproc) {
  for (model_idx in c(1:(length(models_name)))) {
    gridfile=getfilename.indices(work_dir,diag_base,model_idx,grid=T)
    #grid_file_idx=paste0(grid_file,model_idx)
    dummy<-hyint.etccdi.preproc(work_dir,etccdi_dir,etccdi_list_import,gridfile,model_idx,"ALL",yrmon="yr")
  }
}

## Calculate timeseries and trends
if (write_netcdf & run_timeseries) { 
  for (model_idx in c(1:(length(models_name)))) {
    for (seas in seasons) {
      hyint.trends(work_dir,model_idx,seas) 
    }
  }
}

## Create figures
if (write_plots) { 
  ref_idx=which(models_name == reference_model) # select reference dataset; if not available, use last of list
  if(length(ref_idx)==0) { ref_idx=length(models_name)}
  for (seas in seasons) {
    if (plot_type <= 10) { # Plot maps
      hyint.plot.maps(work_dir,plot_dir,work_dir,ref_idx,seas) 
  } else { # Plot timeseries and trends
      print("HyInt: calling plot.trend")  
      hyint.plot.trends(work_dir,plot_dir,work_dir,ref_idx,seas)  
    }
  }
}
