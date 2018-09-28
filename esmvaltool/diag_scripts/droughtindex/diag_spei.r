library(yaml)

#Parsing input file paths and creating output dirs
args <- commandArgs(trailingOnly = TRUE)
params <- read_yaml(args[1])
plot_dir <- params$plot_dir
run_dir <- params$run_dir
work_dir <- params$work_dir
## Create working dirs if they do not exist
#dir.create(plot_dir, recursive = TRUE)
#dir.create(run_dir, recursive = TRUE)
#dir.create(work_dir, recursive = TRUE)
print(plot_dir)
print(run_dir)
print(work_dir)

print(params)
names(params)

input_files_per_var <- yaml::read_yaml(params$input_files)
print(input_files_per_var)

var_names <- names(input_files_per_var)
model_names <- lapply(input_files_per_var, function(x) x$model)
model_names <- unname(model_names)
var0 <- lapply(input_files_per_var, function(x) x$short_name)
fullpath_filenames <- names(var0)
var0 <- unname(var0)[1]
experiment <- lapply(input_files_per_var, function(x) x$exp)
experiment <- unlist(unname(experiment))
rcp_scenario <- experiment

model_names <-  lapply(input_files_per_var, function(x) x$model)
model_names <- unlist(unname(model_names))

start_year <- lapply(input_files_per_var, function(x) x$start_year)
start_year <- c(unlist(unname(start_year)))[1]
end_year <- lapply(input_files_per_var, function(x) x$end_year)
end_year <- c(unlist(unname(end_year)))[1]
