"""
"""

# Script to install r packages


from __future__ import print_function

from rpy2 import robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import glob, os

def main():
    """Start r package installation."""

    print("Installing the required R libraries.\n")
    with open(os.path.abspath(os.path.dirname(__file__)) + '/r_requirements.txt') as f:
       packageNames = f.read().splitlines()

    if all(rpackages.isinstalled(x) for x in packageNames):
        have_packages = True
    else:
        have_packages = False

    if not have_packages:
        utils = rpackages.importr('utils')
        utils.chooseCRANmirror(ind=1)
        packnames_to_install = [x for x in packageNames if not rpackages.isinstalled(x)]
        if len(packnames_to_install) > 0:
            utils.install_packages(StrVector(packnames_to_install))
    print("Installation has finished.")

if __name__ == "__main__":
    main()
