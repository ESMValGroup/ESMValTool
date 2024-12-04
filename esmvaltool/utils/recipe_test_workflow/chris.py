# Used esmvt-rtw env with my branch of ESMValCore installed per:
# https://docs.esmvaltool.org/en/latest/quickstart/installation.html#using-the-development-version-of-the-esmvalcore-package
import esmvalcore
from esmvalcore.config import CFG, Config, Session

# Show I'm working and to be friendly
print("Hello Chris!")

print("ESMValCore Version:", esmvalcore.__version__)

# Generate an InvalidConfigParameter Error
# my_list = [1, 2, 3]
# CFG["config_file"] = my_list

# Generate a Deprecation Warning
CFG["config_file"] = "not_a_dir_but_i_still_work"

# Config file not present
print(CFG['config_file'])

print(CFG.__dict__['_mapping'].keys())
