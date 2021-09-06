begin
    print("INFO    Loading settings from " + getenv("settings"))
    loadscript("$settings")
end
print("INFO Writing " + diag_script_info@setting_name)
n = str_get_nl()
result = "run_dir: " + config_user_info@run_dir + n +\
            "work_dir: " + config_user_info@work_dir + n +\
            "plot_dir: " + config_user_info@plot_dir + n +\
            "log_level: " + config_user_info@log_level + n +\
            "input_files: []" + n

system("echo '" + result + "' > " + diag_script_info@setting_name)
