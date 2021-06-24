import YAML
@info "Starting diagnostic script with" ARGS
config_file = ARGS[1]
cfg = YAML.load_file(config_file)
out_file = cfg["setting_name"]
@info "Copying file to" out_file
Base.Filesystem.cp(config_file, out_file)
@info "Done"
