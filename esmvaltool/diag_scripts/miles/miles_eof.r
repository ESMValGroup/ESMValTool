# #############################################################################
# miles_eof.r
# Authors:       P. Davini (ISAC-CNR, Italy) (author of MiLES)
#                J. von Hardenberg (ISAC-CNR, Italy) (ESMValTool adaptation)
# #############################################################################
# Description
# MiLES is a tool for estimating properties of mid-latitude climate originally thought
# for EC-Earth output and then extended to any model data.
# It works on daily 500hPa geopotential height data and it produces climatological figures 
# for the chosen time period. Data are interpolated on a common 2.5x2.5 grid.  
# Model data are compared against ECMWF ERA-INTERIM reanalysis for a standard period (1989-2010).
# It supports analysis for the 4 standard seasons.#
# Required
#
# Optional 
#
# Caveats
#
# Modification history
#
# ############################################################################

library(tools)
library(yaml)

spath='./esmvaltool/diag_scripts/miles/'

source(paste0(spath,'basis_functions.R'))
source(paste0(spath,'eof_figures.R'))
source(paste0(spath,'eof_fast.R'))
source(paste0(spath,'miles_parameters.R'))

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
metadata <- yaml::read_yaml(settings$input_files)
for (myname in names(settings)) { temp=get(myname,settings); assign(myname,temp)}

field_type0 <- "T2Ds"

# get first variable and list associated to pr variable
var0 <- "zg"
list0 <- metadata

# get name of climofile for first variable and list associated to first climofile
climofiles <- names(list0)
climolist0 <- get(climofiles[1],list0)

diag_base = climolist0$diagnostic
print(paste(diag_base,": starting routine"))

# create working dirs if they do not exist
work_dir=settings$work_dir
regridding_dir=settings$run_dir
plot_dir=settings$plot_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir, recursive = T, showWarnings = F)
dir.create(plot_dir, recursive = T, showWarnings = F)

# extract metadata
models_dataset=unname(sapply(list0, '[[', 'dataset'))
models_ensemble=unname(sapply(list0, '[[', 'ensemble'))
models_exp=unname(sapply(list0, '[[', 'exp'))
reference_model=unname(sapply(list0, '[[', 'reference_dataset'))[1]
models_start_year=unname(sapply(list0, '[[', 'start_year'))
models_end_year=unname(sapply(list0, '[[', 'end_year'))
models_experiment=unname(sapply(list0, '[[', 'exp'))
models_ensemble=unname(sapply(list0, '[[', 'ensemble'))

##
## Run it all
##

for (model_idx in c(1:(length(models_dataset)))) {
    exp <- models_exp[model_idx]
    dataset <- models_dataset[model_idx]
    ensemble <- models_ensemble[model_idx]
    year1=models_start_year[model_idx]
    year2=models_end_year[model_idx]
    infile <- climofiles[model_idx]
    for (tele in teles) {
      for (seas in seasons) {
        miles.eofs.fast(dataset=dataset,expid=exp,ens=ensemble,year1=year1,year2=year2,season=seas,tele=tele,z500filename=infile,FILESDIR=work_dir,doforce=TRUE)
      }
    }
}

##
## Make the plots
##
if (write_plots) {
   ref_idx=which(models_dataset == reference_model)
   if(length(ref_idx)==0) {
      ref_idx=length(models_dataset);
   }
   dataset_ref=models_dataset[ref_idx]
   exp_ref=models_exp[ref_idx]
   ensemble_ref=models_ensemble[ref_idx]
   year1_ref=models_start_year[ref_idx]
   year2_ref=models_end_year[ref_idx]

   for (model_idx in c(1:(length(models_dataset)))) {
      if(model_idx != ref_idx) {
	 exp <- models_exp[model_idx]
         dataset <- models_dataset[model_idx]
         ensemble <- models_ensemble[model_idx]
         year1=models_start_year[model_idx]
         year2=models_end_year[model_idx]
         for (tele in teles) {
            for (seas in seasons) {
               miles.eof.figures(expid=exp,year1=year1,year2=year2,dataset=dataset,ens=ensemble,
		  dataset_ref=dataset_ref,expid_ref=exp_ref,year1_ref=year1_ref,ens_ref=ensemble_ref,
                  year2_ref=year2_ref,season=seas,tele=tele,FIGDIR=plot_dir,
                  FILESDIR=work_dir,REFDIR=work_dir)
            }
         }
      }
   }
}

