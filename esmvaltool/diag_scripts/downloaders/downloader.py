from esmvaltool.diag_scripts.shared import  run_diagnostic

def run_my_diagnostic(cfg):

    return

if __name__ == '__main__':
    # always use run_diagnostic() to get the config (the preprocessor
    # nested dictionary holding all the needed information)
    with run_diagnostic() as config:
        # list here the functions that need to run
        run_my_diagnostic(config)
