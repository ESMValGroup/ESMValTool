# #############################################################################
# rainfarm.r
# Authors:       E. Arnone (ISAC-CNR, Italy)
#	         J. von Hardenberg (ISAC-CNR, Italy) 
# #############################################################################
# Description
# ESMValTool diagnostic calling the RainFARM library written in Julia (by von Hardenberg, ISAC-CNR, Italy).
# RainFARM is a stochastic precipitation downscaling method, further adapted for climate downscaling.
#  
# Required
# CDO
# Julia language: https://julialang.org
# RainFARM Julia library: https://github.com/jhardenberg/RainFARM.jl
#
# Optional 
#
# Caveats
#
# Modification history
#    20180508-A_arnone_e: Conversion to v2.0  
#    20170908-A_arnone_e: 1st github version  
#
# ############################################################################

library(tools)
library(yaml)

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
for (myname in names(settings)) { temp=get(myname,settings); assign(myname,temp)}
metadata <- yaml::read_yaml(settings$input_files)

# store needed arguments in lists (outfile is generated based on the input file)
rainfarm_args=list(slope=slope,nens=nens,nf=nf,weights_climo=weights_climo,varname=varname,
		                      conserv_glob=conserv_glob,conserv_smooth=conserv_smooth)
rainfarm_options=list("-s","-e","-n","-w","-v","-g","-c")

# set variable
var0 <- "pr"

# get name of climofile for selected variable and list associated to first climofile
climofiles <- names(metadata)
climolist <- get(climofiles[1],metadata)

diag_base = climolist$diagnostic
print(paste(diag_base,": starting routine"))

# create working dirs if they do not exist
work_dir=settings$work_dir
regridding_dir=settings$run_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir, recursive = T, showWarnings = F)

# extract metadata
models_name=unname(sapply(metadata, '[[', 'dataset'))
reference_model=unname(sapply(metadata, '[[', 'reference_dataset'))[1]
models_start_year=unname(sapply(metadata, '[[', 'start_year'))
models_end_year=unname(sapply(metadata, '[[', 'end_year'))
models_experiment=unname(sapply(metadata, '[[', 'exp'))
models_ensemble=unname(sapply(metadata, '[[', 'ensemble'))

## Loop through input models, apply pre-processing and call RainFARM
for (model_idx in c(1:(length(models_name)))) {
  # Setup parameters and path 
  exp    <- models_name[model_idx]
  year1  <- models_start_year[model_idx]
  year2  <- models_end_year[model_idx]
  infile <- climofiles[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]
  inregname <- file_path_sans_ext(basename(infile))
  inregfile <- paste0(regridding_dir,"/",inregname,".nc")

  ## Pre-process file converting it to NetCDF4 format 
  if(!file.exists(inregfile)) { 
    cdo_command<-paste("cdo -f nc4 -copy ", infile, inregfile)
    print(paste0(diag_base,": pre-processing file: ", infile))
    system(cdo_command)
    print(paste0(diag_base,": pre-processed file: ", inregfile))
  } else {
    print(paste0(diag_base,": data file exists: ", inregfile))  
  }

  ## Call diagnostic
  print(paste0(diag_base,": calling rainfarm"))
  filename <- paste0(work_dir,"/",inregname,"_downscaled") 

  # reformat arguments from parameter_file
  if (rainfarm_args$conserv_glob == T) {rainfarm_args$conserv_glob <- ""}
  if (rainfarm_args$conserv_smooth == T) {rainfarm_args$conserv_smooth <- ""}
 
  # generate weights file if requested
  # (for more information use 'rfweights -h')
  if (rainfarm_args$weights_climo != F) {
    fileweights <- paste0(work_dir,"/",inregname,"_w.nc")
    snf <- ""
    if (rainfarm_args$nf != F) {snf <- paste("-n ",rainfarm_args$nf)}
    command_w<-paste("rfweights -w ",fileweights,snf," -c ",rainfarm_args$weights_climo,inregfile)  
    print(paste0(diag_base,": generating weights file"))
    print(command_w)
    system(command_w)
    rainfarm_args$weights_climo<-fileweights
  }  
  # assign user defined options
  ret <- which(as.logical(rainfarm_args)!=F|is.na(as.logical(rainfarm_args)))
  rargs <- paste(rainfarm_options[ret],rainfarm_args[ret],collapse=" ")
  # call rfarm
  # (for more information use 'rfarm -h')
  command<-paste0("rfarm -o '", filename,"' ",rargs," ",inregfile) 
  print(command)
  system(command)
  print(paste0(diag_base,": downscaled data written to ",paste0(filename,"_*.nc")))
}
