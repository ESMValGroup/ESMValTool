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

# read parameters
source('esmvaltool/diag_scripts/rainfarm/rainfarm_parameters.r')

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
metadata <- yaml::read_yaml(settings$input_files)
for (myname in names(settings)) { temp=get(myname,settings); assign(myname,temp)}

# store needed arguments in lists (outfile is generated based on the input file)
rainfarm_args=list(slope=slope,nens=nens,nf=nf,weights_climo=weights_climo,varname=varname,
		                      conserv_glob=conserv_glob,conserv_smooth=conserv_smooth)
rainfarm_options=list("-s","-e","-n","-w","-v","-g","-c")

# get first variable and list associated to pr variable
var0 <- "pr"
list0 <- metadata

# get name of climofile for first variable and list associated to first climofile
climofiles <- names(list0)
climolist0 <- get(climofiles[1],list0)

diag_base = climolist0$diagnostic
print(paste(diag_base,": starting routine"))

# create working dirs if they do not exist
work_dir=settings$work_dir
regridding_dir=settings$run_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir, recursive = T, showWarnings = F)

# extract metadata
models_name=unname(sapply(list0, '[[', 'model'))
reference_model=unname(sapply(list0, '[[', 'reference_model'))[1]
models_start_year=unname(sapply(list0, '[[', 'start_year'))
models_end_year=unname(sapply(list0, '[[', 'end_year'))
models_experiment=unname(sapply(list0, '[[', 'exp'))
models_ensemble=unname(sapply(list0, '[[', 'ensemble'))

## Loop through input models, apply pre-processing and call RainFARM
for (model_idx in c(1:(length(models_name)))) {
  # Setup parameters and path 
  exp    <- models_name[model_idx]
  year1  <- models_start_year[model_idx]
  year2  <- models_end_year[model_idx]
  infile <- climofiles[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]
  sgrid <- "noregrid"
  if (rgrid != F) {sgrid <- rgrid} 
  inregname <- paste0(exp,"_",model_exp,"_",model_ens,"_",toString(year1),"-",toString(year2),"_",var0,"_",sgrid,
                      "_",paste(rlonlatdata,collapse="-"))
  inregfile <- paste0(regridding_dir,"/",inregname,".nc")

  ## If needed, pre-process file selecting a limited lon/lat region of interest 
  if((!file.exists(inregfile) | force_processing)) { 
    cdo_command<-paste(paste0("cdo -sellonlatbox,",paste(rlonlatdata,sep="",collapse=",")), 
                       infile, paste0(inregfile,"regtmp"))
    cdo_command2<-paste("cdo -f nc4 -copy ", paste0(inregfile,"regtmp"), inregfile)
    rm_command<-paste("rm ", paste0(inregfile,"regtmp"))
    print(paste0(diag_base,": pre-processing file: ", infile))
    system(cdo_command)
    system(cdo_command2)
    system(rm_command)
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
 
  # generate weights file if needed
  # (for more information use 'rfweights -h')
  if (rainfarm_args$weights_climo != F) {
    fileweights <- paste0(work_dir,"/",inregname,"_w.nc")
    snf <- ""
    if (rainfarm_args$nf != F) {snf <- paste("-n ",rainfarm_args$nf)}
    command_w<-paste("rfweights -w ",fileweights,snf," -c ",rainfarm_args$weights_climo,inregfile)  
    print(paste0(diag_base,": generating weights file"))
    print(fileweights)
    system(command_w)
    #print(command_w)
    rainfarm_args$weights_climo<-fileweights
  }  
  ret <- which(as.logical(rainfarm_args)!=F|is.na(as.logical(rainfarm_args)))
  rargs <- paste(rainfarm_options[ret],rainfarm_args[ret],collapse=" ")
  # call rfarm
  # (for more information use 'rfarm -h')
  command<-paste0("rfarm -o '", filename,"' ",rargs," ",inregfile) 
  system(command)
  print(paste0(diag_base,": downscaled data written to ",paste0(filename,"_*.nc")))
}

#info_output(paste0(">>>>>>>> Leaving ", diag_script), verbosity, 4)
