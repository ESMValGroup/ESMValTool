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

cdo <- function(command, args="", input="", options="", output="",
                stdout="") {
  if (args != "") args <- paste0(",", args)
  if (stdout != "") stdout <- paste0(" > ", stdout)
  argstr <- paste0(options, " ", command, args, " ", input, " ", output,
                   " ", stdout)
  print(paste("cdo", argstr))
  ret <- system2("cdo", args = argstr)
  if (ret != 0) {
    stop(paste("Failed (", ret, "): cdo", argstr))
  }
  return(output)
}


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

# say hi
diag_base <- climolist$diagnostic
print(paste0(diag_base, ": starting routine"))

# get variable name
varname <- climolist$short_name

# create working dirs if they do not exist
dir.create(work_dir, recursive = T, showWarnings = F)
setwd(work_dir)

# setup provenance file and list
provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

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
  print(paste0(diag_base, ": pre-processing file: ", infile))

  print(paste0(diag_base, ": ", perc_lev, " percent quantile"))

  # Select variable of interest
  modf <- cdo("selvar", args = varname, input = infile, output = tempfile())

  # Remap reference onto model grid
  reff <- cdo("remapcon", args = modf,
              input = paste0("-selvar,", varname, " ", ref_data_file),
              output = tempfile())

  # Get (75)th percentile of reference dataset
  ref_perc_pf <- cdo("timpctl", args = perc_lev,
           input = paste(reff, " -timmin ", reff, " -timmax ", reff),
           output = tempfile())

  # Select points with monthly precipitation greater than (75)th perc
  mask_reff <- cdo("ge", input = paste(reff, ref_perc_pf), output = tempfile())
  mask_modf <- cdo("ge", input = paste(modf, ref_perc_pf), output = tempfile())

  # Precipitation sums
  ref_sumf <- cdo("timsum", input = paste("-mul", mask_reff, reff),
                  output = tempfile())
  mod_sumf <- cdo("timsum", input = paste("-mul", mask_modf, modf),
                  output = tempfile())

  # Calculate quantile bias, set name and attributes
  qb1f <- cdo("div", input = paste(mod_sumf, ref_sumf), output = tempfile())
  cdo_command <- paste(
    " -setattribute,qb@standard_name='precipitation_quantile_bias'",
    " -setattribute,qb@long_name='Precipitation quantile bias'",
    " -setattribute,qb@units=' '",
    paste0(" -chname,", varname, ",qb")
  )
  cdo(cdo_command, input = qb1f, output = outfile)

  # Remove temporary files
  unlink(c(modf, reff, ref_perc_pf, mask_reff, mask_modf,
           ref_sumf, mod_sumf, qb1f))

  # Set provenance for this output file
  caption <- paste0("Precipitation quantile bias ", perc_lev, "% for years ",
                    year1, " to ", year2, " according to ", exp)
  xbase <- list(list(infile), list("arno_en", "hard_jo"), list("c3s-magic"),
                list("mehran14jgr"), caption, list("perc"), list("atmos"),
                list("phys"), list("global"), ref_data_file)
  names(xbase) <- c("ancestors", "authors", "projects", "references",
                    "caption", "statistics", "realms",
                    "themes", "domains", "reference_dataset")

  # Store provenance in main provenance list
  provenance[[outfile]] <- xbase
}

# Write provenance to file
write_yaml(provenance, provenance_file)

# End of diagnostic
print(paste0(diag_base, ": done."))
