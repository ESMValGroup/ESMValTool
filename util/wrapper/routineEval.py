import yaml
import os
import sys
import tempfile
import shutil
import subprocess
import argparse

from jinja2 import Template

sys.path.insert(0,
                os.path.join(
                    os.path.dirname(os.path.abspath(__file__)),
                    '../nml-utils/generateNML'))

from generateNML import get_namelist

evtRoot = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../..')

runscript = """#!/bin/bash
#SBATCH -J {{ jobid }}
#SBATCH -p prepost
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu 1280
#SBATCH -t 100
#SBATCH -A bd0854
#SBATCH -o {{ jobid }}.o
#SBATCH -e {{ jobid }}.e

module load bb/evt

python main.py namelist.xml

"""

t_runscript = Template(runscript)


class CMIP6Dataset(yaml.YAMLObject):
    yaml_tag = u'!CMIP6Dataset'

    def __init__(self, project, name, product, institute, model, experiment,
                 mip, ensemble, start_year, end_year, variable):
        self.project = project
        self.name = name
        self.product = product
        self.institute = institute
        self.model = model
        self.experiment = experiment
        self.mip = mip
        self.ensemble = ensemble
        self.start_year = start_year
        self.end_year = end_year
        self.variable = variable

    def get_dict(self):
        return self.__dict__

    def __repr__(self):
        return "%s( project=%r, name=%r, product=%r, institute=%r, model=%r, experiment=%r, mip=%r, ensemble=%r, start_year=%r, end_year=%r, variable=%r )" % (
            self.__class__.__name__, self.project, self.name, self.product,
            self.institute, self.model, self.experiment, self.mip,
            self.ensemble, self.start_year, self.end_year, self.variable)


def examplefile():
    e = """--- !CMIP6Dataset
project: ESGF_CMIP6
name: HadGEM3-GC31-LL
product: ScenarioMIP
institute: MOHC
model: HadGEM3-GC31-LL
experiment: ssp245
mip: Amon
ensemble: r1i1p1f1
grid: gn
start_year: 1990
end_year: 1992
variable: prw """
    return e


def process(yfile=examplefile(), verbose=True):
    """According to input file yfile build the temprorary director(y/ies) and submit the jobs
    """
    h = yaml.load_all(yfile)
    lfd = 0
    user = os.getenv('USER')
    for d in h:
        lfd += 1
        tmpbase = tempfile.mkdtemp(dir='/scratch/b/{0}/rEval'.format(
            user))  # Needs to be mounted on compute nodes
        tmpdir = os.path.join(tmpbase, 'ESMValTool')
        shutil.copytree(
            evtRoot,
            tmpdir,
            symlinks=False,
            ignore=shutil.ignore_patterns('.git', '.git*'))
        #shutil.copytree(evtRoot, tmpdir, symlinks=False, ignore=shutil.ignore_patterns('.git', '.git*', '*'))
        with open(os.path.join(tmpdir, 'namelist.xml'), 'w') as f:
            f.write(get_namelist(**d.get_dict()))
        s_runscript = os.path.join(tmpdir, 'runscript')
        with open(s_runscript, 'w') as f:
            f.write(
                t_runscript.render(jobid="rEval{0}".format(str(lfd).zfill(3))))
        if verbose:
            print(tmpdir)
        subprocess.Popen('sbatch runscript', cwd=tmpdir, shell=True)


def main():
    """Executes the ESMValTool for routine evaluation on CMIP6-Data described in the input file. The ESMValTool must be configured beforehand."""
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument(
        '-f',
        '--infile',
        dest='infile',
        help='Path to file in yaml-format describing the dataset.',
        type=argparse.FileType('r'))
    parser.add_argument(
        '-e',
        '--example',
        dest='example',
        action='store_true',
        help='Output example input file.')

    args = parser.parse_args()

    if args.example:
        print(examplefile())
        return

    process(yfile=args.infile)


if __name__ == "__main__":
    main()
