# #############################################################################
# quantilebias.r
# Authors:       E. Arnone (ISAC-CNR, Italy)
#	         S. Terzago (ISAC-CNR, Italy) 
#	         J. von Hardenberg (ISAC-CNR, Italy) 
# #############################################################################
# Description
# ESMValTool diagnostic for calculation of the precipitation quantile bias
# following Mehran et al. (2014)
#  
# Required
# CDO
#
# Optional 
#
# Caveats
#
# Modification history
#    20180518-A_arnone_e: Written for v2.0  
#
# ############################################################################

library(tools)
library(yaml)

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
metadata <- yaml::read_yaml(settings$input_files)
for (myname in names(settings)) { temp=get(myname,settings); assign(myname,temp)}

# set pr variable and get metadata list associated to it
var0 <- "pr"
list0 <- metadata

# get name of climofile for first variable and list associated to first climofile
climofiles <- names(list0)
climolist0 <- get(climofiles[1],list0)

diag_base = climolist0$diagnostic
print(paste0(diag_base,": starting routine"))

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
  print("======================")
  print(exp)
  print(model_exp)
  print(model_ens)

  inregname <- paste0(exp,"_",model_exp,"_",model_ens,"_",toString(year1),"-",toString(year2),"_",var0)
  inregfile <- infile
  outfile <- paste0(work_dir,"/",inregname,"_",perc_lev,"qb.nc")
  print(paste0(diag_base,": pre-processing file: ", infile))

  ref="/work/datasets/climate/GPCP/pr/data/mon/gpcp_22.nc"                         # reference precipitation data
  model=infile
  orog_tmp="/work/models/cmip5/fx/orog/orog_fx_EC-EARTH_historical_r0i0p0.nc"       # modificare DEM con uno osservato

  perc=perc_lev #"75"
  print(paste0(diag_base,": ",perc," percent quantile"))

  cdo_command=paste("cdo -mulc,86400",model,"tmp_model.nc")
   print(cdo_command)
   system(cdo_command)
  cdo_command=paste("cdo griddes tmp_model.nc > tmp_grid")
   print(cdo_command)
   system(cdo_command)
  cdo_command=paste("cdo remapcon,tmp_grid ",paste0("-selyear,",year1,"/",year2),ref,"tmp_ref.nc")
   print(cdo_command)
   system(cdo_command)

  # Get (75)th percentile of reference dataaset 
  cdo_command=paste(paste0("cdo timpctl,",perc)," tmp_ref.nc  -timmin tmp_ref.nc  -timmax tmp_ref.nc  tmp_ref_perc_p.nc")
   print(cdo_command)
   system(cdo_command)

  # Select points with monthly precipitation greater than 75th perc
  cdo_command=paste("cdo ge tmp_ref.nc  tmp_ref_perc_p.nc tmp_mask_ref.nc")
   print(cdo_command)
   system(cdo_command)
  cdo_command=paste("cdo ge tmp_model.nc   tmp_ref_perc_p.nc tmp_mask_model.nc")
   print(cdo_command)
   system(cdo_command)

  # Precipitation sums
  cdo_command=paste("cdo timsum -mul tmp_mask_ref.nc tmp_ref.nc tmp_ref_sum.nc")
   print(cdo_command)
   system(cdo_command)
  cdo_command=paste("cdo timsum -mul tmp_mask_model.nc  tmp_model.nc tmp_model_sum.nc")
   print(cdo_command)
   system(cdo_command)

  # Quantile bias
  cdo_command=paste("cdo div tmp_model_sum.nc tmp_ref_sum.nc tmp_qb.nc")
   print(cdo_command)
   system(cdo_command)

  # Check with Mehran et al. 2014 
  cdo_command=paste("cdo remapnn,tmp_grid -gtc,5 ",orog_tmp, "tmp_mask_orog.nc")   # modificare DEM con uno osservato  &  la soglia dei 5m
   print(cdo_command)
   system(cdo_command)
  cdo_command=paste("cdo mul tmp_qb.nc tmp_mask_orog.nc tmp_qb_landonly.nc")
   print(cdo_command)
   system(cdo_command)

  # Copy file to output destination and remove temporary files
  mv_command=paste("mv tmp_qb.nc ",outfile)
   print(mv_command)
   system(mv_command)
  rm_command=paste("rm tmp*")
   print(mv_command)
   system(mv_command)
}











