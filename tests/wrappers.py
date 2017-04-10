"""
module that contains wrapper classes that enable mainly to construct filenames
for the random data generator dependent on the model classes specified in the namelist files
"""

import os

class Wrapper(object):
    def __init__(self, s):
        self._h = s.lstrip().split()
        self._check_version()
        self._append_coordinates = False  # specifies if lat/lon fields should be appended as 2D variables when random data is generated
        self._append_cellsize = False # same for cellsize

    def _check_version(self):
        assert self._h[0].strip().upper() == self._ref.upper()

    def get_dimension_tag(self, dim):
        if dim == 2:
            return 'T2Ms'
        elif dim == 3:
            return 'T3M'
        else:
            assert False

    def replace_dimension(self, s, dim):
        return s.replace('@{DIM}', self.get_dimension_tag(dim))


class CMIP5_base(Wrapper):
    """
    wrapper class to construct CMIP5 type directory structures
    """
    def __init__(self, s):
        """
        Parameters
        ----------
        s : str
            string for CMIP5 class like specified in namelist
            CMIP5  MPI-ESM-LR Amon historical r1i1p1 2000 2004  @{MODEL_PATH}
        """
        super(CMIP5_base, self).__init__(s)

        self.name = self._h[1].strip()
        self.mip = self._h[2].strip()
        self.experiment = self._h[3].strip()
        self.ens = self._h[4].strip()
        self.start = int(self._h[5].strip())  # start year
        self.stop = int(self._h[6].strip()) # stop year
        self.dir = self._h[7].strip()

    def get_oname(self, s, varname, dim=None):
        """
        Parameters
        ----------
        s : str
            string as specified in namelist
        varname : str
            name of variable
        """
        assert dim is not None
        #@{VAR_FILE}
        sep ='_'
        filename = varname + sep + self.mip + sep + self.name + sep + self.experiment + sep + self.ens + '.nc'
        o = s.replace('@{VAR_FILE}', filename)
        return self.replace_dimension(o, dim)




class CMIP5(CMIP5_base):
    def __init__(self, s):
        self._ref = 'CMIP5'
        super(CMIP5, self).__init__(s)

class CMIP5fx_w(CMIP5_base):
    def __init__(self, s):
        self._ref = 'CMIP5_fx'
        super(CMIP5fx_w, self).__init__(s)

class CMIP5_ETHZ(CMIP5_base):
    def __init__(self, s):
        self._ref = 'CMIP5_ETHZ'
        super(CMIP5_ETHZ, self).__init__(s)

    def get_oname(self, s, varname, dim=None):
        assert dim is not None
        #@{VAR_FILE}
        sep ='_'

        dirname = self.experiment + os.sep + self.mip + os.sep + varname + os.sep + self.name + os.sep + self.ens + os.sep
        filename = varname + sep + self.mip + sep + self.name + sep + self.experiment + sep + self.ens + '.nc'
        o = s.replace('@{VAR_FILE}', dirname + filename)
        return self.replace_dimension(o, dim)





class OBS_base(Wrapper):

    def __init__(self, s):

        super(OBS_base, self).__init__(s)

        self.name = self._h[1].strip()
        self.mip = self._h[2].strip()
        self.experiment = self._h[3].strip()

        self.start = int(self._h[4].strip())  # start year
        self.stop = int(self._h[5].strip()) # stop year
        self.dir = self._h[6].strip()

    def get_oname(self, s, varname, dim=None):
        """
        Parameters
        ----------
        s : str
            string as specified in namelist
        varname : str
            name of variable
        """
        assert dim is not None
        #@{VAR_FILE}
        sep ='_'

        filename = 'OBS' + sep + self.name + sep + self.mip + sep + self.experiment + sep + '@{DIM}' + sep + varname + sep + '123456-123456' + '.nc'

        #~ varname + sep + self.mip + sep + self.name + sep + self.experiment

        #~ filename = varname + sep + self.mip + sep + self.name + sep + self.experiment
        o = s.replace('@{VAR_FILE}', filename)
        return self.replace_dimension(o, dim)

#~ sic_reanaly_HadISST_20130524.nc


class OBS(OBS_base):
    def __init__(self, s):
        self._ref = 'OBS'
        super(OBS, self).__init__(s)

class OBS_gridfile(OBS_base):
    def __init__(self, s):
        self._ref = 'OBS_gridfile'
        super(OBS_gridfile, self).__init__(s)
        self._append_coordinates = True
        self._append_cellsize = True


class OBS4MIPS(OBS_base):
    def __init__(self, s):
        self._ref = 'obs4mips'
        super(OBS4MIPS, self).__init__(s)

    def get_oname(self, s, varname, dim=None):
        assert dim is not None
        sep ='_'

        filename = self.name + os.sep + varname + sep + self.name + sep + self.mip + sep + self.experiment + sep + 'abcd' + '.nc'
        o = s.replace('@{VAR_FILE}', filename)
        return self.replace_dimension(o, dim)


class ANA4MIPS(OBS_base):
    def __init__(self, s):
        self._ref = 'ana4mips'
        super(ANA4MIPS, self).__init__(s)

    def get_oname(self, s, varname, dim=None):
        assert dim is not None
        sep ='_'

        filename = self.name + os.sep + varname + sep + self.mip + sep + self.experiment + sep + self.name + sep + 'abcd' + '.nc'
        o = s.replace('@{VAR_FILE}', filename)
        return self.replace_dimension(o, dim)
