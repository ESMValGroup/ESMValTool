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

# get path to script and source subroutines (if needed)
args <- commandArgs(trailingOnly = FALSE)
spath <- paste0(dirname(unlist(strsplit(
  grep("--file", args,
    value = TRUE
  ), "="
))[2]), "/")

source(paste0(spath, "quantilebias_functions.R"))
source(paste0(spath, "../shared/external.R")) # nolint

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
dir.create(plot_dir, recursive = T, showWarnings = F)
setwd(work_dir)

# setup provenance file and list
provenance_file <- paste0(run_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

# extract metadata
models_name <- unname(sapply(metadata, "[[", "dataset"))
reference_model <-
  unname(sapply(metadata, "[[", "reference_dataset"))[1]
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

  inregname <- paste0(
    exp,
    "_",
    model_exp,
    "_",
    model_ens,
    "_",
    toString(year1),
    "-",
    toString(year2),
    "_",
    varname
  )
  outfile <-
    paste0(work_dir, "/", inregname, "_", perc_lev, "qb.nc")

  # Select variable of interest
  modf <- cdo("selvar", args = varname, input = infile)

  reff <- ref_data_file

  # Get (X)th percentile of reference dataset
  refminf <- cdo("timmin", input = reff)
  refmaxf <- cdo("timmax", input = reff)
  ref_perc_pf <- cdo("timpctl",
    args = perc_lev,
    input = c(reff, refminf, refmaxf)
  )

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

  temp1f <-
    cdo("chname", args = paste0(varname, ",qb"), input = qb1f)
  temp2f <-
    cdo("setattribute", args = "qb@units=' '", input = temp1f)
  temp1f <- cdo("setattribute",
    args = "qb@long_name='Precipitation quantile bias'",
    input = temp2f,
    output = temp1f
  )
  cdo("setattribute",
    args = "qb@standard_name='precipitation_quantile_bias'",
    input = temp1f,
    output = outfile
  )

  # Remove temporary files
  unlink(
    c(
      modf,
      ref_perc_pf,
      mask_reff,
      mask_modf,
      ref_sumf,
      mod_sumf,
      qb1f,
      refminf,
      refmaxf,
      mask_mod2f,
      mask_ref2f,
      temp1f,
      temp2f
    )
  )


  # Produce figure
  field <- ncdf_opener(outfile, "qb", "lon", "lat", rotate = "full")
  ics_ref <- ics
  ipsilon_ref <- ipsilon

  tmp_figname <- sub(".nc", paste0(".", output_file_type), outfile)
  figname <- sub(work_dir, plot_dir, tmp_figname)

  figure_size <- c(600, 400)
  if (tolower(output_file_type) != "png") {
    figure_size <- c(10, 6)
  }
  graphics_startup(figname, output_file_type, figure_size)

  tmp_levels <- c(0:20) * 0.1
  tmp_colors <- rev(rainbow(30)[1:20])

  # contours
  par(
    cex.main = 1.8,
    cex.axis = 1.4,
    cex.lab = 1.4,
    mar = c(5, 5, 4, 8)
  )
  filled_contour3(
    ics,
    ipsilon,
    field,
    xlab = "Longitude",
    ylab = "Latitude",
    main = paste0(exp),
    levels = tmp_levels,
    col = tmp_colors,
    axes = F,
    asp = 1
  )
  # continents
  map(
    "world",
    regions = ".",
    interior = F,
    exact = F,
    boundary = T,
    add = T,
    col = "black",
    lwd = 2
  )
  axis(1, col = "grey40", at = seq(-180, 180, 45))
  axis(2, col = "grey40", at = seq(-90, 90, 30))

  colorbar_scale <- c(-0.15, -0.08, 0.1, -0.1)
  if (tolower(output_file_type) != "png") {
    colorbar_scale <- c(-0.13, -0.06, 0.1, -0.1)
  }
  image_scale3(
    volcano,
    levels = tmp_levels,
    new_fig_scale = colorbar_scale,
    col = tmp_colors,
    colorbar.label = paste0("QB", perc_lev),
    cex.colorbar = 1.3,
    cex.label = 1.4,
    colorbar.width = 1,
    line.label = 2.9,
    line.colorbar = 1.0,
    extend = F
  )
  graphics_close(figname)

  # Set provenance for this output file
  caption <-
    paste0(
      "Precipitation quantile bias ",
      perc_lev,
      "% for years ",
      year1,
      " to ",
      year2,
      " according to ",
      exp
    )
  xbase <- list(
    ancestors = list(infile, ref_data_file),
    authors = list("arnone_enrico", "vonhardenberg_jost"),
    projects = list("c3s-magic"),
    references = list("mehran14jgr"),
    caption = caption,
    statistics = list("perc"),
    realms = list("atmos"),
    themes = list("phys"),
    domains = list("global"),
    reference_dataset = ref_data_file
  )

  # Store provenance in main provenance list
  provenance[[outfile]] <- xbase
  provenance[[figname]] <- xbase
}

# Write provenance to file
write_yaml(provenance, provenance_file)

# End of diagnostic
print(paste0(diag_base, ": done."))
