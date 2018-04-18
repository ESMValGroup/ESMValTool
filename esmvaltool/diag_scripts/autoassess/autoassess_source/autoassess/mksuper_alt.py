#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-
"""
(C) Crown Copyright 2017, the Met Office

Python script to download seasonal and annual means from MOOSE and
mean them together to form seasonal and annual supermeans (using ppmnvar)
"""

# What if meaning Decembers? What happens to year label then?

import os
import shutil


class UMFileError(Exception):
    pass


class UMFile(object):
    '''Class to handle UM files'''

    _moose = 'moose:/crum'
    _years = ['ann']
    _seasons = 'ndjfmamjjasond'
    _months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun',
               'jul', 'aug', 'sep', 'oct', 'nov', 'dec']

    def __init__(self, jobid, startyear, nyears, period):
        self.id = jobid
        self.startyear = startyear
        self.nyears = nyears
        self.period = period
        self.moose = '{0._moose}/{0.id}'.format(self)

    @property
    def endyear(self):
        '''Define last year to mean from first year and number of years'''
        return self.startyear + self.nyears - 1

    @property
    def id5(self):
        '''Set str to be last 5 characters of jobid'''
        return self.id[-5:]

    @property
    def season(self):
        '''Return file name used by UM to represent meaned period'''
        if self.period in self._years:
            seasonstr = '1201'
        else:
            seasonstr = self.period
        return seasonstr

    @property
    def streamcode(self):
        '''Return MASS stream used for source data'''
        if self.period in self._years:
            stream = 'y'
        elif self.period in self._seasons:
            stream = 's'
        elif self.period in self._months:
            stream = 'm'
        return stream

    @property
    def yearcode(self):
        '''Return MASS stream used for meaned data'''
        if self.nyears <= 9:
            stream = chr(self.nyears+48)
        elif self.nyears <= 10:
            stream = 'a'
        elif self.nyears <= 19:
            stream = 'b'
        elif self.nyears <= 24:
            stream = 'k'
        elif self.nyears <= 29:
            stream = 's'
        elif self.nyears <= 39:
            stream = 't'
        elif self.nyears <= 49:
            stream = 'q'
        elif self.nyears <= 99:
            stream = 'l'
        elif self.nyears <= 249:
            stream = 'u'
        elif self.nyears <= 499:
            stream = 'w'
        else:
            stream = 'd'
        return stream

    @property
    def url_file(self):
        template = '{0.moose}/ap{0.streamcode}.pp'
        return template.format(self)

    @property
    def url_mean(self):
        template = '{0.moose}/am{0.yearcode}.pp'
        return template.format(self)

    def get_file(self, year):
        """Generate a um file name"""
        template = '{0.id5}a.p{0.streamcode}{1!s}{0.season}.pp'
        return template.format(self, year)

    def get_mean(self, filestr):
        """
        Generate a mean file name
        filestr: m = mean, v = variance, s = standard deviation
        """
        template = '{0.id5}a.{1}{0.yearcode}{0.endyear!s}{0.period}.pp'
        return template.format(self, filestr)

    def generate_files(self):
        (syear, nyear) = (self.startyear, self.nyears)
        return [self.get_file(syear+i) for i in range(0, nyear)]

    def generate_urls(self):
        return [os.path.join(self.url_file, x) for x in self.generate_files()]


class ArchiveError(Exception):
    pass


class Archive(object):

    def __init__(self, listfiles, target):
        self._files = listfiles
        self._url = target

    @property
    def _strfiles(self):
        return ' '.join(self._files)

    @property
    def _prtfiles(self):
        return '\n'.join(self._files)

    def arch_cmd(self, moocmd):
        template = 'moo {1} -f {0._strfiles} {0._url}'
        return template.format(self, moocmd)

    def print_msg(self, moocmd):
        template = '{1} the following files to {0._url}\n{0._prtfiles}'
        return template.format(self, moocmd)

    def put_files(self):
        print(self.print_msg('Archiving'))
        ret = os.system(self.arch_cmd('put'))
        if ret != 0:
            print('Error archiving files to MASS-R')
            print('Return code = {0}'.format(ret))
            raise ArchiveError('Error archiving file to MASS')

    def get_files(self):
        print(self.print_msg('Retrieving'))
        ret = os.system(self.arch_cmd('get'))
        if ret != 0:
            print('Error retrieving files from MASS-R.')
            print('Return code = {0}'.format(ret))
            raise ArchiveError('Error retrieving files from MASS')


class PPMeanVarError(Exception):
    pass


class PPMeanVar(object):

    _exec = '/opt/ukmo/utils/bin/run_ppmnvar'
    _listfile = 'ppmnvar.list'

    def __init__(self, infiles, meanfile, varfile):
        self._files = infiles
        self._mean = meanfile
        self._var = varfile

    @property
    def _prtfiles(self):
        return '\n'.join(self._files)

    def mean(self):
        print('Making supermean {0._mean} from:\n{0._prtfiles}'.format(self))

        # Populate input file for ppmnvar
        with open(self._listfile, 'w') as outf:
            lastfile = self._files[-1]
            for infile in self._files:
                outf.write(infile)
                if infile == lastfile:
                    outf.write(' {0._mean} {0._var}'.format(self))
                outf.write('\n')

        # Call run_ppmnvar
        ret = os.system('{0._exec} -F {0._listfile} -V'.format(self))
        if ret != 0:
            raise PPMeanVarError('Error meaning files ({0})'.format(ret))


class STDevError(Exception):
    pass


class STDev(object):

    # Can this be written in Python using raw PP Iris interface instead?
    _pro = '/home/h02/hadvg/cma/general/idl/procedures/var_to_stdev.pro'
    _idl = 'stdev_idl.pro'
    _log = 'mksuper_idl.log'

    def __init__(self, varfile, stdevfile):
        self._var = varfile
        self._std = stdevfile

    def make_stdev(self):
        print('Making standard deviation {0._std}'.format(self))
        idlcode = 'status=1\n'
        idlcode += '.compile {0._pro}\n'
        idlcode += '!PP_LARGEFILE_EXTEND=0\n'
        idlcode += 'var_to_stdev, "./{0._var}", "./{0._std}", status=status\n'
        idlcode += 'exit, status=status\n'
        with open(self._idl, 'w') as idlf:
            idlf.write(idlcode.format(self))
        # Make sure (T)IDL environment is not contaminated
        os.environ['IDL_STARTUP'] = ' '
        os.environ['IDL_PP_PATH'] = '.'
        ret = os.system('tidl -quiet {0._idl} 2>{0._log}'.format(self))
        if ret != 0:
            print('Error converting variances into standard deviations:')
            print('From variance file {0._var}'.format(self))
            print('To standard deviation file {0._std}'.format(self))
            print('Using temporary idl code {0._idl}'.format(self))
            print('Check {0._log}'.format(self))
            print('Return code = {0}'.format(ret))
            raise STDevError('Error converting variances to standard devs')


class MakeSuperError(Exception):
    pass


class MakeSuper(object):

    _worktemp = 'mksuper_{0}_{1!s}_{2!s}'

    def __init__(self, workdir, jobid, nyears, startyear, seasons):
        tempdir = self._worktemp.format(jobid, startyear, nyears)
        self._workdir = os.path.join(workdir, tempdir)
        self._jobid = jobid
        self._nyears = nyears
        self._startyear = startyear
        self._seasons = seasons

    @property
    def _cwd(self):
        return os.path.join(self._workdir, self._season)

    def _mkcwd(self):
        # Remove directory if it exists
        if os.path.exists(self._cwd):
            shutil.rmtree(self._cwd)
        # Create directory
        os.makedirs(self._cwd)
        # Enter directory
        os.chdir(self._cwd)

    def _mksuper(self, stdev=True):
        """The main part of the program"""

        # Create object to generate UM file names
        try:
            umfiles = UMFile(self._jobid, self._startyear,
                             self._nyears, self._season)
        except UMFileError:
            raise

        # Retrieve files to mean from MASS
        try:
            urlarray = umfiles.generate_urls()
            Archive(urlarray, '.').get_files()
        except ArchiveError:
            raise

        # Make a supermean using ppmnvar
        try:
            filearray = umfiles.generate_files()
            meanfile = umfiles.get_mean('m')
            varfile = umfiles.get_mean('v')
            PPMeanVar(filearray, meanfile, varfile).mean()
            arch_files = [meanfile]
        except PPMeanVarError:
            raise

        # Convert variances to standard deviations if required
        if stdev:
            try:
                stdevfile = umfiles.get_mean('s')
                STDev(varfile, stdevfile).make_stdev()
                arch_files.append(stdevfile)
            except STDevError:
                raise

        # Archive the supermeans to MASS
        try:
            Archive(arch_files, umfiles.url_mean).put_files()
        except ArchiveError:
            raise

    def means(self, stdev=True):
        # Loop over seasons
        for self._season in self._seasons:
            # Create and enter working directory
            self._mkcwd()
            # Create supermean for given season
            self._mksuper(stdev=stdev)
            # Remove working directory and contents
            shutil.rmtree(self._cwd)
        # Tidy up - fails if not empty
        os.rmdir(self._workdir)


class DecDate(object):

    _dec_temp = '{0!s:4s}1201T0000Z'

    def __init__(self, datec):
        self._datec = str(datec)

    @property
    def year(self):
        return int(self._datec[0:4])

    @property
    def dec_year(self):
        return self._dec_temp.format(self.year)

    @property
    def dec_after(self):
        datec_yr = self.year
        if self.dec_year < self._datec:
            datec_yr += 1
        return datec_yr

    @property
    def dec_before(self):
        datec_yr = self.year
        if self._datec < self.dec_year:
            datec_yr -= 1
        return datec_yr


def main():
    workdir = os.environ['SCRATCH']
    jobid = os.environ['JOBID']
    basis = os.environ['BASIS']
    tcycle = os.environ['TCYCLE']
    ncycle = os.environ['NCYCLE']
    seasons = os.environ['SEASONS'].split()
    years = os.environ['YEARS'].split()

    # Calculate year of first 0Z 1st December after basis
    basis_yr = DecDate(basis).dec_after

    # Calculate year of last 0Z 1st December before current model time
    cycle_yr = DecDate(ncycle).dec_before

    # Only run if 0Z 1st December within bounds of this cycle and next cycle
    if tcycle < DecDate(cycle_yr).dec_year <= ncycle:
        for year in years:
            (x, y) = year.split(':')
            (fyear, nyears) = (int(x), int(y))
            if basis_yr + fyear + nyears == cycle_yr:
                startyear = basis_yr + fyear + 1
                mksuper = MakeSuper(workdir, jobid, nyears, startyear, seasons)
                mksuper.means()

if __name__ == '__main__':
    main()
