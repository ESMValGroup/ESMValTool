library(yaml)
args <- commandArgs(trailingOnly = TRUE)
print(paste0("INFO    Loading settings from ", args[1]))
settings <- yaml::read_yaml(args[1])

print(paste0("INFO    Writing settings to ", settings$setting_name))
yaml::write_yaml(settings, settings$setting_name)
