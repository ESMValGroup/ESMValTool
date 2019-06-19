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
#    20180926-arnone_enrico: Refined for usage as recipe
#    20180518-arnone_enrico: Written for v2.0
#
# #############################################################################

library(tools)
library(yaml)

cdo <- function(command, args = "", input = "", options = "", output = "",
                stdout = "", noout = F) {
  if (args != "") args <- paste0(",", args)
  if (stdout != "") stdout <- paste0(" > '", stdout, "'")
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
  modf <- cdo("selvar", args = varname, input = infile)

  # Remap reference onto model grid
  selectf <- cdo("selvar", args = varname, input = ref_data_file)
  reff <- cdo("remapcon", args = modf, input = selectf)

  # Get (X)th percentile of reference dataset
  refminf <- cdo("timmin", input = reff)
  refmaxf <- cdo("timmax", input = reff)
  ref_perc_pf <- cdo("timpctl", args = perc_lev,
                     input = c(reff, refminf, refmaxf))

  # Select points with monthly precipitation greater than (75)th perc
  mask_reff <- cdo("ge", input = c(reff, ref_perc_pf))
  mask_modf <- cdo("ge", input = c(modf, ref_perc_pf))

  # Precipitation sums
  mask_ref2f <- cdo("mul", input = c(mask_reff, reff))
  mask_mod2f <- cdo("mul", input = c(mask_modf, modf))
  ref_sumf <- cdo("timsum", input = mask_ref2f)
  mod_sumf <- cdo("timsum", input = mask_mod2f)

  # Calculate quantile bias, set name and attributes
  qb1f <- cdo("div", input = c(mod_sumf, ref_sumf))
  tempfile <- tempfile()

  temp1f <- cdo("chname", args = paste0(varname, ",qb"), input = qb1f)
  temp2f <- cdo("setattribute", args = "qb@units=' '", input = temp1f)
  temp1f <- cdo("setattribute",
                args = "qb@long_name='Precipitation quantile bias'",
                input = temp2f, output = temp1f)
  cdo("setattribute", args = "qb@standard_name='precipitation_quantile_bias'",
      input = temp1f, output = outfile)

  # Remove temporary files
  unlink(c(modf, reff, ref_perc_pf, mask_reff, mask_modf,
           ref_sumf, mod_sumf, qb1f, refminf, refmaxf, selectf,
           mask_mod2f, mask_ref2f, temp1f, temp2f))

  # Set provenance for this output file
  caption <- paste0("Precipitation quantile bias ", perc_lev, "% for years ",
                    year1, " to ", year2, " according to ", exp)
  xbase <- list(ancestors = list(infile, ref_data_file),
                authors = list("arnone_enrico", "vonhardenberg_jost"),
                projects = list("c3s-magic"), references = list("mehran14jgr"),
                caption = caption, statistics = list("perc"),
                realms = list("atmos"), themes = list("phys"),
                domains = list("global"), reference_dataset = ref_data_file)

  # Store provenance in main provenance list
  provenance[[outfile]] <- xbase
}

# Write provenance to file
write_yaml(provenance, provenance_file)

# End of diagnostic
print(paste0(diag_base, ": done."))
