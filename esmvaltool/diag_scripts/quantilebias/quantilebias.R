# #############################################################################
# quantilebias.r
# Authors:       E. Arnone (ISAC-CNR, Italy)
# 	         S. Terzago (ISAC-CNR, Italy)
# 	         J. von Hardenberg (ISAC-CNR, Italy)
# #############################################################################
# Description
# ESMValTool diagnostic for calculation of the precipitation quantile bias
# following Mehran et al. (2014)
#
# Required
# - #CDO
# - observational monthly mean global precipitation climatology
# - (e.g. from the GPCP project: http://gpcp.umd.edu)
#
# Optional
#
# Caveats
#
# Modification history
#    20180926-A_arno_en: Refined for usage as recipe
#    20180518-A_arno_en: Written for v2.0
#
# #############################################################################

library(tools)
library(yaml)

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
metadata <- yaml::read_yaml(settings$input_files)
for (myname in names(settings)) {
  temp <- get(myname, settings)
  assign(myname, temp)
}

# get name of climofile and list associated to first climofile
climofiles <- names(metadata)
climolist <- get(climofiles[1], metadata)

# get variable name
varname <- climolist$short_name

# say hi
diag_base <- climolist$diagnostic
print(paste0(diag_base, ": starting routine"))

# create working dirs if they do not exist
work_dir <- settings$work_dir
regridding_dir <- settings$run_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir, recursive = T, showWarnings = F)
setwd(work_dir)

# extract metadata
models_name <- unname(sapply(metadata, "[[", "dataset"))
reference_model <- unname(sapply(metadata, "[[", "reference_dataset"))[1]
models_start_year <- unname(sapply(metadata, "[[", "start_year"))
models_end_year <- unname(sapply(metadata, "[[", "end_year"))
models_experiment <- unname(sapply(metadata, "[[", "exp"))
models_ensemble <- unname(sapply(metadata, "[[", "ensemble"))

ref_idx <- which(models_name == reference_model)
ref_data_file <- climofiles[ref_idx]

## Loop through input models
for (model_idx in c(1:(length(models_name)))) {
  if (model_idx == ref_idx) {
     next
  }
  # Setup parameters and path
  exp <- models_name[model_idx]
  year1 <- models_start_year[model_idx]
  year2 <- models_end_year[model_idx]
  infile <- climofiles[model_idx]
  model_exp <- models_experiment[model_idx]
  model_ens <- models_ensemble[model_idx]

  inregname <- paste0(exp, "_", model_exp, "_", model_ens, "_",
                      toString(year1), "-", toString(year2), "_", varname)
  outfile <- paste0(work_dir, "/", inregname, "_", perc_lev, "qb.nc")
  outfile_landonly <- paste0(work_dir, "/", inregname, "_",
                             perc_lev, "qb_landonly.nc")
  print(paste0(diag_base, ": pre-processing file: ", infile))

  print(paste0(diag_base, ": ", perc_lev, " percent quantile"))
  system(paste0("cdo selvar,", varname, " ", infile, " tmp_model.nc"))

  # Remap reference onto model grid
  system(paste0("cdo remapcon,tmp_model.nc -selvar,", varname, " ",
                ref_data_file, " tmp_ref.nc"))

  # Get (75)th percentile of reference dataset
  system(paste0("cdo timpctl,", perc_lev, " tmp_ref.nc  -timmin ",
                "tmp_ref.nc -timmax tmp_ref.nc tmp_ref_perc_p.nc"))

  # Select points with monthly precipitation greater than (75)th perc
  system("cdo ge tmp_ref.nc  tmp_ref_perc_p.nc tmp_mask_ref.nc")
  system("cdo ge tmp_model.nc tmp_ref_perc_p.nc tmp_mask_model.nc")

  # Precipitation sums
  system("cdo timsum -mul tmp_mask_ref.nc tmp_ref.nc tmp_ref_sum.nc")
  system("cdo timsum -mul tmp_mask_model.nc tmp_model.nc tmp_model_sum.nc")

  # Calculate quantile bias, set name and attributes
  system("cdo div tmp_model_sum.nc tmp_ref_sum.nc tmp_qb1.nc")
  cdo_command <- paste(
    "cdo ",
    " -setattribute,qb@standard_name='precipitation_quantile_bias'",
    " -setattribute,qb@long_name='Precipitation quantile bias'",
    " -setattribute,qb@units=' '",
    paste0(" -chname,", varname, ",qb"),
    " tmp_qb1.nc tmp_qb.nc"
  )
  print(cdo_command)
  system(cdo_command)

  # Select land only using > 5 meter (Mehran et al. 2014)
  system("cdo -f nc topo tmp_orog.nc")
  system("cdo remapnn,tmp_model.nc -gtc,5 tmp_orog.nc tmp_mask_orog.nc")
  system("cdo mul tmp_qb.nc tmp_mask_orog.nc tmp_qb_landonly.nc")

  # Copy file to output destination and remove temporary files
  mv_command <- paste("mv tmp_qb.nc", outfile)
  print(mv_command)
  system(mv_command)
  mv_command <- paste("mv tmp_qb_landonly.nc", outfile_landonly)
  print(mv_command)
  system(mv_command)
  rm_command <- paste("rm tmp_*")
  print(rm_command)
  system(rm_command)
}

print(paste0(diag_base, ": done."))
