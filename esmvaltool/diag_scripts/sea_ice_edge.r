interface_data <- file.path(Sys.getenv('ESMValTool_interface_data'), 'r.interface')
source(interface_data)
source('diag_scripts/lib/R/info_output.r')

## Do not print warnings
options(warn=-1)

var0 <- variables[1]
field_type0 <- field_types[1]

info_output(paste0("<<<<<<<< Entering ", diag_script), verbosity, 4)
info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)
info_output(paste0("plot - ", diag_script, " (var: ", variables[1], ")"), verbosity, 1)
info_output("+++++++++++++++++++++++++++++++++++++++++++++++++", verbosity, 1)

library(tools)
diag_base = file_path_sans_ext(diag_script)
library(s2dverification)
library(startR)

## Create working dirs if they do not exist
dir.create(plot_dir, showWarnings = FALSE)
dir.create(file.path(plot_dir, diag_base), showWarnings = FALSE)
dir.create(climo_dir, showWarnings = FALSE)

##
## Run it all
##
print(models_name)
for (model_idx in c(1:length(models_name))) {
    fullpath_filename <- interface_get_fullpath(var0, field_type0, model_idx)
    n_diag_script <- nchar(diag_script)
    diag_script_base <- substr(diag_script, 0, n_diag_script - 2)
    figure_filename <- interface_get_figure_filename(diag_script_base,
                                                 'sea_ice_edge',
                                                 field_type0,
                                                 '',
                                                 model_idx)
    figure_path <- file.path(plot_dir, paste0(figure_filename, '.', output_file_type))
    data <- Start(dat = fullpath_filename,
                  var = var0,
                  time = 'all',
                  lat = 'all',
                  lon = 'all',
                  var_var = 'var_names',
                  return_vars = list(time = NULL,
                                     lat = NULL,
                                     lon = NULL),
                  retrieve = TRUE)

    PlotLayout(PlotStereoMap, c('lat', 'lon'),
               Subset(data, 'time', 1:4),
                  attr(data, 'Variables')$common$lon,
                  attr(data, 'Variables')$common$lat,
                  fileout = figure_path)
}
info_output(paste0(">>>>>>>> Leaving ", diag_script), verbosity, 4)
