"""
"""

# Script to install r packages


from __future__ import print_function

from rpy2 import robjects
import rpy2.robjects.packages as rpackages
from rpy2.robjects.vectors import StrVector
import glob, os

import argparse as ap


def main():
    """Start r package installation."""

    parser = ap.ArgumentParser(description="R package installation script")
    parser.add_argument('-f', '--file', type=ap.FileType('r'), help='Requirements file with list of R libraries.')
    args, leftovers = parser.parse_known_args()

    if args.file is None:
        args.file = 'r_requirements.txt'
        print ("Requirements file is not set. Using the default file {}".format(args.file))

    print("Installing the required R libraries.\n")
    with open(os.path.abspath(os.path.dirname(__file__)) + '/' + args.file) as f:
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
