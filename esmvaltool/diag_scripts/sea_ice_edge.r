interface_data <- file.path(Sys.getenv('ESMValTool_interface_data'), 'r.interface')
source(interface_data)
source(diag_script_cfg)
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

## Plot function


##
## Run it all
##
print(models_name)
for (model_idx in c(1:length(models_name))) {
    fullpath_filename <- interface_get_fullpath(var0, field_type0, model_idx)
    n_diag_script <- nchar(diag_script)
    diag_script_base <- substr(diag_script, 0, n_diag_script - 2)


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

    data[which(data < concentration_treshold)] <- 0
    data[which(data >= concentration_treshold)] <- 1

    lat = attr(data, 'Variables')$common$lon
    lon = attr(data, 'Variables')$common$lat



    time_dim <- which(names(dim(data)) == 'time')
    data_dim = dim(data)
    data_dim[time_dim] <- 12
    season_mean <- array(0, data_dim)
    for( i in 1:12){
        season_mean[,,i,,] <- Mean1Dim(Season(data, posdim=time_dim, 1, i, i), 3)
    }

    figure_filename <- interface_get_figure_filename(diag_script_base,
                                                 '',
                                                 field_type0,
                                                 '',
                                                 model_idx)
    figure_path <- file.path(plot_dir, diag_base, paste0(figure_filename, '.', output_file_type))
    info_output(figure_path, verbosity, 1)
    title <- paste0('Probability of being inside sea-ice edge (', models_start_year[model_idx],'-', models_end_year[model_idx],')')
    library(RColorBrewer)
    colors<-colorRampPalette(rev(brewer.pal(9, 'Blues')))(10)
    PlotLayout(PlotStereoMap, c('lat', 'lon'),
               Subset(season_mean, 'time', 1:12),
               lat,
               lon,
               bar_limits = c(0, 1),
               cols=colors[1:10],
               col_sup=colors[10],
               col_inf=colors[1],
               filled.continents = TRUE,
               coast_color = 'black',
               titles = as.character(months(attr(data, 'Variables')$common$time[1:12], abbreviate=FALSE)),
               toptitle = title,
               fileout = figure_path)
                   # cols=c('blue', 'white'),
}
info_output(paste0(">>>>>>>> Leaving ", diag_script), verbosity, 4)
