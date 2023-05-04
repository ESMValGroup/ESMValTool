#!/usr/bin/env python
"""Generate the required user configuration file for ESMValTool."""
import os
import pprint

import yaml
from esmvalcore.config._config_object import Config

import esmvaltool


def main():
    """Write the required user configuration file for ESMValTool.

    The default configuration values from the latest version of
    ESMValTool are updated with the configuration values defined in the
    environment for the ``configure`` task and the result is written to
    the file defined by ``USER_CONFIG_PATH`` in the ``flow.cylc`` file.
    """
    # Get the default configuration values from the latest version of
    # ESMValTool.
    config_values = get_default_config_values_from_latest_ESMValTool()

    # Get the configuration values defined in the environment for the
    # 'configure' task.
    config_values_from_task_env = get_config_values_from_task_env()

    # Update the default configuration values.
    config_values.update(config_values_from_task_env)
    user_config_path = os.environ["USER_CONFIG_PATH"]
    config_values["config_file"] = user_config_path

    # Write the updated configuration values to the file defined by
    # 'user_config_path'.
    print(f"Writing the user configuration file to '{user_config_path}' with "
          "values: ")
    pprint.PrettyPrinter().pprint(config_values)
    write_yaml(user_config_path, config_values)


def get_default_config_values_from_latest_ESMValTool():
    """Return the default configuration values from the latest version of
    ESMValTool."""
    latest_esmvaltool_dir = os.environ["ESMVALTOOL_DIR"]
    default_config_file = os.path.join(latest_esmvaltool_dir,
                                       "config-user-example.yml")
    # A dictionary is needed here to avoid the config object from
    # converting paths to 'PosixPath' objects, which causes issues when
    # writing the YAML file.
    config_values = dict(Config._load_user_config(default_config_file))

    # Update the type of the value of 'config_developer_file' from a
    # 'PosixPath' object to a string, to avoid issues when writing the
    # YAML file (all other paths will be overwritten).
    config_developer_file = str(config_values["config_developer_file"])
    config_values["config_developer_file"] = config_developer_file

    print("Default configuration values from the latest version of ESMValTool "
          f"{esmvaltool.__version__} (from '{default_config_file}'): ")
    pprint.PrettyPrinter().pprint(config_values)

    return config_values


def get_config_values_from_task_env():
    """Return the configuration values defined in the environment for the
    ``configure`` task."""
    # Note that 'auxiliary_data_dir' and 'download_dir' are set to empty
    # values and cannot currently be configured. 'auxiliary_data_dir' is
    # used by some recipes to look for additional datasets, so may need
    # to be configured in the future. 'download_dir' is used only when
    # using the automatic download feature via ESMValTool, which is not
    # the intention here, so should not be configurable.
    #
    # In addition, 'check_level' and 'extra_facets_dir' are added to the
    # config object in ESMValTool after loading the user configuration
    # file, but they need updating to avoid issues when writing the YAML
    # file.
    config_values_from_task_env = {
        "auxiliary_data_dir": "",
        "check_level": "DEFAULT",
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
    print("The configuration values defined in the environment for the "
          "'configure' task: ")
    pprint.PrettyPrinter().pprint(config_values_from_task_env)
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
