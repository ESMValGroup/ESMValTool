library(yaml)
library(s2dverification)
library(startR)
library(tools)
library(RColorBrewer)
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])

plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
input_files <- read_yaml(params$input_files)

## Do not print warnings
options(warn=-1)

## Create working dirs if they do not exist
dir.create(plot_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(work_dir, recursive=TRUE, showWarnings=FALSE)

## Plot function
for (dataset_id in c(1:length(input_files))) {
    dataset <- input_files[[dataset_id]]
    print(paste0('Computing ', dataset$project, '_', dataset$dataset, '_', dataset$exp, '_', dataset$ensemble))
    data <- Start(dat = dataset$filename,
                  var = dataset$short_name,
                  var_var = 'var_names',
                  time = 'all',
                  lat = 'all',
                  lon = 'all',
                  return_vars = list(time = NULL,
                                     lat = NULL,
                                     lon = NULL),
                  retrieve = TRUE)
    print(str(data))
    print(data[,,1,,1])
    data[which(data < params$concentration_treshold)] <- 0
    data[which(data >= params$concentration_treshold)] <- 1

    lat = attr(data, 'Variables')$common$lon1
    lon = attr(data, 'Variables')$common$lat

    time_dim <- which(names(dim(data)) == 'time')
    data_dim = dim(data)
    data_dim[time_dim] <- 12
    season_mean <- array(0, data_dim)
    print(str(data))
    print(data[,,1,,1])
    for( i in 1:12){
        season_mean[,,i,,] <- Mean1Dim(Season(data, posdim=time_dim, 1, i, i), 3)
    }
    print(str(season_mean))
    print(season_mean[,,6,,1])
    figure_filename <- paste0(dataset$project, '_', dataset$dataset, '_', dataset$exp, '_', dataset$ensemble, '_seaiceedge_', dataset$start_year, '-', dataset$end_year)
    figure_path <- file.path(plot_dir, paste0(figure_filename, '.', params$output_file_type))
    title <- paste0('Probability of being inside sea-ice edge (', dataset$start_year,'-', dataset$end_year,')')


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
}
