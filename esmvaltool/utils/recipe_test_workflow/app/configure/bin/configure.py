#!/usr/bin/env python
"""Generate the required user configuration file for ESMValTool."""

import os
import pprint

import yaml
from esmvalcore.config._config_validators import ValidationError, _validators


def main():
    """Write the required user configuration file for ESMValTool.

    The configuration values defined in the environment for the
    ``configure`` task are written to the configuration file defined
    by ``USER_CONFIG_PATH`` in the ``flow.cylc`` file.
    """
    # Get the configuration values defined in the environment for the
    # 'configure' task.
    config_values = get_config_values_from_task_env()

    # Validate the user config file content.
    validate_user_config_file(config_values)

    # Update the configuration from OS environment.
    user_config_path = os.environ["USER_CONFIG_PATH"]
    config_values["config_file"] = user_config_path

    # Write the updated configuration values to the file defined by
    # 'user_config_path'.
    print(
        f"Writing the user configuration file to '{user_config_path}' with "
        "values: "
    )
    pprint.PrettyPrinter().pprint(config_values)
    write_yaml(user_config_path, config_values)


def get_config_values_from_task_env():
    """Get configuration values from environment for configure task."""
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
    print(
        "The configuration values defined in the environment for the "
        "'configure' task: "
    )
    pprint.PrettyPrinter().pprint(config_values_from_task_env)
    return config_values_from_task_env


def validate_user_config_file(user_config_file_content):
    """Validate a user config with ``ESMValCore.config._validators`` functions.

    Parameters
    ----------
    user_config_file_content: dict
        An ESMValTool user configuration file loaded in memory as a Python
        dictionary.

    Raises
    ------
    KeyError
        If ``user_config_file_content`` includes a key for which there is no
        validator listed in ``_validators``,
    ValidationError
        If any of the called validation functions raise a ValidationError.
    """
    errors = [
        "There were validation errors in your user configuration file. See "
        "details below.\n"
    ]
    for user_config_key, usr_config_value in user_config_file_content.items():
        try:
            validation_function = _validators[user_config_key]
        except KeyError as err:
            errors.append(
                f"Key Error for {user_config_key.upper()}. May not be a valid "
                f"ESMValTool user configuration key\nERROR: {err}\n"
            )
        else:
            try:
                print(
                    f"Validating {user_config_key.upper()} with value "
                    f'"{usr_config_value}" using function '
                    f"{validation_function.__name__.upper()}."
                )
                validation_function(usr_config_value)
            except ValidationError as err:
                errors.append(
                    f"Validation error for {user_config_key.upper()} with "
                    f'value "{usr_config_value}"\nERROR: {err}\n'
                )
    if len(errors) > 1:
        raise ValidationError("\n".join(errors))
    print("All validation checks passed.")


def write_yaml(file_path, contents):
    """Write ``contents`` to the YAML file ``file_path``.

    Parameters
    ----------
    file_path: string
        The full path to the YAML file to write the contents to.
    contents: dictionary
        The contents to write to the YAML file.
    """
    with open(file_path, "w", encoding="utf-8") as file_handle:
        yaml.dump(contents, file_handle, default_flow_style=False)


if __name__ == "__main__":
    main()
