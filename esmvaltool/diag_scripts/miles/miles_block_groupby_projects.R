# #############################################################################
# miles_block_groupby_projects.r
# Authors:       P. Davini (ISAC-CNR, Italy) (author of MiLES)
# 	         J. von Hardenberg (ISAC-CNR, Italy) (ESMValTool adaptation)
#                E. Arnone (ISAC-CNR, Italy) (ESMValTool v2.0 adaptation)
#                R. Kazeroni (DLR) (project grouping for TM90/DA98)
# #############################################################################
# Description
# MiLES is a tool for estimating properties of mid-latitude climate.
# It works on daily 500hPa geopotential height data and it produces
# climatological figures for the chosen time period. Data are interpolated
# on a common 2.5x2.5 grid.
# Model data are compared against a reference field such as the
# ECMWF ERA-Interim reanalysis.
#
# Modification history
#   20210122-kazeroni_remi: modified to average over projects (CMIP5, CMIP6)
#   20181203-vonhardenberg_jost: Completed conversion, rlint compliant
#
# ############################################################################

library(tools)
library(yaml)

provenance_record <- function(infile) {
  xprov <- list(
    ancestors = infile,
    authors = list(
      "vonhardenberg_jost", "davini_paolo",
      "arnone_enrico"
    ),
    references = list(
      "davini18", "davini12jclim",
      "tibaldi90tel"
    ),
    projects = list("c3s-magic"),
    caption = "MiLES blocking statistics",
    statistics = list("other"),
    realms = list("atmos"),
    themes = list("phys"),
    domains = list("nh")
  )
  return(xprov)
}

# create Multimodel figure name and folder (no ens and dataset)
fig_multimodel_builder <- function(FIGDIR,
                        dir_name,
                        file_name,
                        dataset,
                        expid,
                        year1,
                        year2,
                        season,
                        output_file_type) {
  # loop on descriptors that are concatenated to create dir and file name
  descriptors <-
    c(dataset, expid, paste0(year1, "-", year2), season)
  for (dcode in descriptors) {
    if (dcode != "NO") {
      FIGDIR <- file.path(FIGDIR, dcode)
      file_name <- paste0(file_name, "_", dcode)
    }
  }

  # add directory name descriptor
  FIGDIR <- file.path(FIGDIR, dir_name)

  # actually dir.exists is in devtools only for R < 3.2,
  # then is included in base package
  if (exists("dir.exists")) {
    if (!dir.exists(FIGDIR)) {
      dir.create(FIGDIR, recursive = T)
    }
  } else {
    dir.create(FIGDIR, recursive = T, showWarnings = F)
  }

  return(file.path(FIGDIR, paste0(file_name, ".", output_file_type)))
}

##
## Set default parameters
##
color_lines <- c("dodgerblue", "darkred", "black", "green", "orange")
indice <- "DA98"
legend_loc <- c(-90, 35)
linewidth <- 4
obs_legend <- c("OBS", "OBS2", "OBS3", "OBS4", "OBS5")
plot_title <- "Instantaneous Blocking: DJFM 1979 - 2000"
transparency <- 0.15 # in [0, 1]
xlabel <- "Longitude"
ylabel <- "Blocked Days (%)"
yrange <- c(0, 35)

diag_scripts_dir <- Sys.getenv("diag_scripts")

source(paste0(diag_scripts_dir, "/miles/basis_functions.R"))
source(paste0(diag_scripts_dir, "/miles/block_fast.R"))
source(paste0(diag_scripts_dir, "/miles/miles_parameters.R"))
source(paste0(diag_scripts_dir, "/shared/external.R")) # nolint

# read settings and metadata files
args <- commandArgs(trailingOnly = TRUE)
settings <- yaml::read_yaml(args[1])
metadata <- yaml::read_yaml(settings$input_files)
for (myname in names(settings)) {
  temp <- get(myname, settings)
  assign(myname, temp)
}

field_type0 <- "T2Ds"

# get first variable and list associated to pr variable
var0 <- "zg"
list0 <- metadata

# get name of climofile for first variable and list
# associated to first climofile
climofiles <- names(list0)
climolist0 <- get(climofiles[1], list0)

diag_base <- climolist0$diagnostic
print(paste(diag_base, ": starting routine"))

# create working dirs if they do not exist
work_dir <- settings$work_dir
regridding_dir <- settings$run_dir
plot_dir <- settings$plot_dir
dir.create(work_dir, recursive = T, showWarnings = F)
dir.create(regridding_dir,
  recursive = T,
  showWarnings = F
)
dir.create(plot_dir, recursive = T, showWarnings = F)

# setup provenance file and list
provenance_file <-
  paste0(regridding_dir, "/", "diagnostic_provenance.yml")
provenance <- list()

# extract metadata
models_dataset <- unname(sapply(list0, "[[", "dataset"))
models_ensemble <- unname(sapply(list0, "[[", "ensemble"))
models_exp <- unname(sapply(list0, "[[", "exp"))
models_projects <- unname(sapply(list0, "[[", "project"))
models_start_year <- unname(sapply(list0, "[[", "start_year"))
models_end_year <- unname(sapply(list0, "[[", "end_year"))
models_experiment <- unname(sapply(list0, "[[", "exp"))
models_ensemble <- unname(sapply(list0, "[[", "ensemble"))

##
## Run it all
##

project_list <- c("OBS", "CMIP5", "CMIP6")
project_obs_list <- c("OBS", "OBS6", "RAWOBS", "native6")
for (model_idx in c(1:(length(models_dataset)))) {
  exp <- models_exp[model_idx]
  dataset <- models_dataset[model_idx]
  ensemble <- models_ensemble[model_idx]
  year1 <- models_start_year[model_idx]
  year2 <- models_end_year[model_idx]
  infile <- climofiles[model_idx]
  project <- models_projects[model_idx]
  if (is.element(project, project_obs_list)) {
    models_projects[model_idx] <- "OBS" # gather all observation data
  }
  for (seas in seasons) {
    filenames <- miles_block_fast(
      year1 = year1,
      year2 = year2,
      expid = exp,
      ens = ensemble,
      dataset = dataset,
      season = seas,
      z500filename = infile,
      FILESDIR = work_dir,
      doforce = TRUE
    )
    # Set provenance for output files
    xprov <- provenance_record(list(infile))
    for (fname in filenames) {
      provenance[[fname]] <- xprov
    }
  }
}
obs_idx <- which(models_projects == "OBS") # isolate OBS datasets

##
## Plotting parameters
##
if (write_plots) {
  color_field <- color_lines
  color_diff <- NULL
  lev_field <- yrange
  lev_diff <- NULL
  lev_hist <- NULL
  legend_unit <- ylabel
  legend_distance <- 3
  title_name <- plot_title
  x_label <- xlabel
  fp <- list(
    color_field = color_field,
    color_diff = color_diff,
    lev_field = lev_field,
    lev_diff = lev_diff,
    lev_hist = lev_hist,
    legend_unit = legend_unit,
    legend_distance = legend_distance,
    title_name = title_name,
    x_label = x_label
  )
  alpha <- transparency*255 # normalized transparency coefficient
  # panels option
  par(
    cex.main = 2,
    cex.axis = 1.5,
    cex.lab = 1.5,
    mar = c(5, 5, 4, 3),
    oma = c(0, 0, 0, 0)
  )
  lwdline <- linewidth
  tm90cols <- fp$color_field

##
## Compute mean and std for datasets grouped by project
##
  field <- indice
  for (season in seasons) {

    # create filenames to handle the provenance
    filenames <- c()
    # create figure names with ad-hoc function
    dataset <- 'Multimodel'
    expid <- 'historical'
    figname <- fig_multimodel_builder(
      plot_dir,
      "Block",
      field,
      dataset,
      expid,
      year1,
      year2,
      season,
      output_file_type
    )
    filenames <- c(filenames, figname)
    field_exp_all <- c() # store field for all datasets
    for (model_idx in c(1:(length(models_dataset)))) {
      expid <- models_exp[model_idx]
      dataset <- models_dataset[model_idx]
      ens <- models_ensemble[model_idx]
      year1 <- models_start_year[model_idx]
      year2 <- models_end_year[model_idx]
      project <- models_projects[model_idx]
      # use file.builder function
      nomefile <- file_builder(
      work_dir,
      "Block",
      "BlockClim",
      dataset,
      expid,
      ens,
      year1,
      year2,
      season
      )
      field_exp <- ncdf_opener(nomefile, namevar = field, rotate = "no")
      assign(paste(field, "_exp", sep = ""), field_exp)
      field_exp_all <- c(field_exp_all, list(field_exp))

      # Set provenance for output files (same as diagnostic files)
      #xprov <- provenance_record(climofiles[model_idx])
      #provenance[[filenames]] <- xprov
    }

    open_plot_device(figname, output_file_type, special = TRUE)
    # rotation to simplify the view (90 deg to the west)
    n <- (-length(ics) / 4)
    ics2 <- c(tail(ics, n), head(ics, -n) + 360)
    text_legend <- c() # text of the legendens <- ''
    i_project <- 1

    for (project in project_list){
      field_exp_mean <- field_exp_all[models_projects == project]
      print(paste(project, "this project"))
      print(paste(length(field_exp_mean), "length exp"))
      if (length(field_exp_mean) > 0) { # plot mean by CMIP5-6-project or every obs dataset available

        if (project != "OBS") { # mean computed only for CMIP5 and 6 datasets
          field_mean <- apply(X=as.data.frame(field_exp_mean), MARGIN = 1, FUN = mean) # mean over datasets of the project
          field_std <- apply(X=as.data.frame(field_exp_mean), MARGIN = 1, FUN = sd) # std
          n_datasets <- length(field_exp_mean) # number of datasets for the project
          mycol <- rgb(aperm(col2rgb(tm90cols[i_project])), max = 255, alpha = alpha) # adjust color for the shaded area (std)
          text_legend <- c(text_legend, paste(paste(project, ":", sep=""), n_datasets, "datasets"))
          field_mean2 <- c(tail(field_mean, n), head(field_mean, -n))
          field_std2 <- c(tail(field_std, n), head(field_std, -n))
          field_mean_up2 <- field_mean2 + field_std2
          field_mean_do2 <- field_mean2 - field_std2
          print(paste(length(field_mean), "length mean"))
          print(paste(length(field_mean2), "length mean2"))

        } else { #start with the first obs data
          field_mean <- unlist(field_exp_mean[1])
          field_mean2 <- c(tail(field_mean, n), head(field_mean, -n))
          text_legend <- c(text_legend, obs_legend[1])
          print(paste(length(field_mean), "length mean"))
          print(paste(length(field_mean2), "length mean2"))
        }

        if (i_project == 1){
          plot(
            ics2,
            field_mean2,
            type = "l",
            lwd = lwdline,
            xlim = c(-90, 270),
            ylim = fp$lev_field,
            main = fp$title_name,
            cex.main = 1.3,
            xlab = fp$x_label,
            ylab = fp$legend_unit,
            frame.plot = FALSE,
            col = tm90cols[1],
            xaxt = "none",
            xaxs = "i",
            yaxs = "i"
          )
          axis(
            1,
            seq(-90, 270, 30),
            labels = c(" ", "60°W", "30°W", "0°", "30°E", "60°E", "90°E", "120°E", "150°E", "180°E", "150°W", "120°W", " ")
          )
          grid(
            nx = 12,
            ny = NULL
          )
        }

        else if (project == "OBS" && length(field_exp_mean) > 1) {
          for (idx in c(2, length(field_exp_mean))) {
            field_mean <- unlist(field_exp_mean[idx])
            field_mean2 <- c(tail(field_mean, n), head(field_mean, -n))
            text_legend <- c(text_legend, obs_legend[idx])
            i_project <- i_project +1
            points(
              ics2,
              field_mean2,
              type = "l",
              lwd = lwdline,
              col = tm90cols[i_project]
            )
          }
        }
        else {
          points(
            ics2,
            field_mean2,
            type = "l",
            lwd = lwdline,
            col = tm90cols[i_project]
          )
        }
        # shaded area
        if (project != "OBS") {
          polygon(c(ics2, rev(ics2)), c(field_mean_do2, rev(field_mean_up2)), col = mycol, border = NA) # shaded area [mean-std; mean+std]
        }
        i_project <- i_project + 1
      }

    }
    legend(
     legend_loc[1],
     legend_loc[2],
      legend = text_legend,
      lwd = lwdline,
      lty = 1,
      col = tm90cols[1:length(project_list)+length(obs_idx)-1],
      bg = "white",
      cex = 1.
    )
    dev.off()
    xprov <- provenance_record(climofiles)
    provenance[[filenames]] <- xprov
  }
}

# Write provenance to file
write_yaml(provenance, provenance_file)
