#!/usr/bin/env python
"""Generate the required user configuration file for ESMValTool."""
import os
import warnings

with warnings.catch_warnings():
    warnings.filterwarnings(
        "ignore",
        message=(
            "\n  Thank you for trying out the new ESMValCore API.\n  Note "
            "that this API is experimental and may be subject to change.\n  "
            "More info: https://github.com/ESMValGroup/ESMValCore/issues/498"),
        category=UserWarning,
        module="esmvalcore.experimental._warnings",
    )
    from esmvalcore.experimental.config._config_object import Config
import yaml

ESMVALTOOL_DIR = os.environ["ESMVALTOOL_DIR"]
LATEST_USER_CONFIG_FILE = os.path.join(ESMVALTOOL_DIR,
                                       "config-user-example.yml")
USER_CONFIG_OUTPUT_PATH = os.environ["USER_CONFIG_PATH"]


def main():
    """Write the updated configuration values to the file defined by
    ``USER_CONFIG_PATH``."""
    # Get the latest configuration values from ESMValTool
    # Write to a dictionary to avoid PosixPath type issues in the YAML file
    # Update with the path to the latest checked out version
    # of config-user-example.yml from the ESMValTool directory.

    config_values = dict(
        Config._load_user_config(LATEST_USER_CONFIG_FILE,
                                 raise_exception=True))

    # Get the configuration values defined in the environment for the
    # ``configure`` task.
    config_values_from_task_env = get_config_values_from_task_env()

    # Update the default configuration values.
    config_values.update(config_values_from_task_env)

    # Write the updated configuration values to the file defined by
    # 'USER_CONFIG_PATH'.
    write_yaml(USER_CONFIG_OUTPUT_PATH, config_values)


def get_config_values_from_task_env():
    """Return the configuration values defined in the environment for the
    ``configure`` task."""
    # Note that 'auxiliary_data_dir', 'download_dir' and
    # 'extra_facets_dir' are set to empty values and cannot currently be
    # configured. However, 'download_dir' is used only when using the
    # automatic download feature via ESMValTool (which we do not intend
    # to use here) and 'extra_facets_dir' is not available in the
    # default configuration file provided by ESMValTool v2.6.0.
    # 'auxiliary_data_dir' is used by some recipes to look for
    # additional datasets, so may need to be configured in the future.
    config_values_from_task_env = {
        "auxiliary_data_dir": "",
        "config_file": USER_CONFIG_OUTPUT_PATH,
        "download_dir": "",
        "drs": {
            "ana4mips": os.environ["DRS_ANA4MIPS"],
            "CMIP3": os.environ["DRS_CMIP3"],
            "CMIP5": os.environ["DRS_CMIP5"],
            "CMIP6": os.environ["DRS_CMIP6"],
            "CORDEX": os.environ["DRS_CORDEX"],
            "native6": os.environ["DRS_NATIVE6"],
            "OBS": os.environ["DRS_OBS"],
            "obs4MIPs": os.environ["DRS_OBS4MIPS"],
            "OBS6": os.environ["DRS_OBS6"],
        },
        "extra_facets_dir": [],
        "max_parallel_tasks": int(os.environ["MAX_PARALLEL_TASKS"]),
        "output_dir": os.environ["OUTPUT_DIR"],
        "remove_preproc_dir": False,
        "rootpath": {
            "ana4mips": os.environ["ROOTPATH_ANA4MIPS"],
            "CMIP3": os.environ["ROOTPATH_CMIP3"],
            "CMIP5": os.environ["ROOTPATH_CMIP5"],
            "CMIP6": os.environ["ROOTPATH_CMIP6"],
            "CORDEX": os.environ["ROOTPATH_CORDEX"],
            "native6": os.environ["ROOTPATH_NATIVE6"],
            "OBS": os.environ["ROOTPATH_OBS"],
            "obs4MIPs": os.environ["ROOTPATH_OBS4MIPS"],
            "OBS6": os.environ["ROOTPATH_OBS6"],
            "RAWOBS": os.environ["ROOTPATH_RAWOBS"],
        },
    }
    return config_values_from_task_env


def write_yaml(file_path, contents):
    """Write the contents specified by ``contents`` to the YAML file specified
    by ``file_path``.

    Parameters
    ----------
    file_path: string
        The full path to the YAML file to write the contents to.
    contents: dictionary
        The contents to write to the YAML file.
    """
    with open(file_path, "w") as file_handle:
        yaml.dump(contents, file_handle, default_flow_style=False)


if __name__ == "__main__":
    main()
