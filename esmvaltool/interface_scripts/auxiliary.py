import logging
import subprocess
import sys

logger = logging.getLogger(__name__)


def ncl_version_check():
    """ @brief Check the NCL version"""
    try:
        cmd = ['ncl', '-V']
        version = subprocess.check_output(cmd)
    except subprocess.CalledProcessError:
        logger.error("Failed to execute '%s'", ' '.join(cmd))
        raise

    version = version.decode(sys.stdout.encoding)

    if version == "6.3.0":
        logger.error("NCL version " + version + " not supported due to a bug "
                     + "(see Known Issues in the ESMValTool user guide)")

    if int(version.split(".")[0]) < 6:
        logger.error("NCL version " + version +
                     " not supported, need version 6.2.0 or higher")

    if int(version.split(".")[0]) == 6 and int(version.split(".")[1]) < 2:
        logger.error("NCL version " + version +
                     " not supported, need version 6.2.0 or higher")
