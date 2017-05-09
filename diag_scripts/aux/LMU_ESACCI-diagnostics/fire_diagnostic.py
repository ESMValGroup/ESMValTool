import os
import subprocess
from netCDF4 import Dataset
from ESMValMD import ESMValMD
from diagnostic import BasicDiagnostics


class FireDiagnostic(BasicDiagnostics):
    """
    class to implement fire diagnostics,
    like e.g. global means, global differences, RMSD etc.

    TODO implement testing for this diagnostic

    """

    def __init__(self, **kwargs):
        super(FireDiagnostic, self).__init__(**kwargs)

        self._project_info = {}
        self._modtype = None
        self._reftype = None
        self._plot_dir = '.' + os.sep
        self._work_dir = '.' + os.sep

        self._vartype = 'burned area'  # default value as there must be one
        self.output_type = 'png'  # default ouput file type
        self._changed = False

    def run_diagnostic(self):
        """
        running the diagnostics
        """
        super(FireDiagnostic, self).run_diagnostic()

    def _specific_diag(self):
        """
        Diagnostic management
        """

    def write_data(self, plot=True):
        """
        write data
        """
        super(FireDiagnostic, self).write_data()

    def _load_model_data(self):
        """ load model data """

        mod_info = Dataset(self._mod_file)
        try:
            lat = mod_info.dimensions['lat'].size
            lon = mod_info.dimensions['lon'].size
        except:  # regridding required in any case
            lat = -1
            lon = -1
        mod_info.close()

        if not (lat == self.reso["lat"] and lon == self.reso["lon"]):

            grid = self.resoLUT[str(self.reso["lat"]) + "-" +
                                str(self.reso["lon"])]

            newfile = self._mod_file + "." + grid + "built.nc"
            newfile = newfile.split("/")
            newdir = (self._work_dir if self._work_dir[-1] ==
                      os.sep else self._work_dir + os.sep) +\
                "AUX_Files_fire_ESACCI"
            newfile = newdir + os.sep + newfile[-1]

            if not os.path.exists(newfile):
                tempfile = self._aggregate_resolution(
                    self._mod_file, grid, remove=False)
                subprocess.call(["mkdir", newdir])
                subprocess.call(['cp', tempfile, newfile])
                os.remove(tempfile)

            self._mod_file = newfile

        # load data
        proj_var = self._project_info['RUNTIME']['currDiag'].get_variables()[0]
        self._mod_data = self._load_cmip_generic(self._mod_file, proj_var)

    def _load_observation_data(self):
        """ load obs data """

        mod_info = Dataset(self._ref_file)
        try:
            lat = mod_info.dimensions['lat'].size
            lon = mod_info.dimensions['lon'].size
        except:  # regridding required in any case
            lat = -1
            lon = -1
        mod_info.close()

        reso = {"lat": lat, "lon": lon}

        if str(lat) + "-" + str(lon) not in self.resoLUT.keys():

            grid = self.resoLUT[str(self.reso["lat"]) + "-" +
                                str(self.reso["lon"])]

            newfile = self._ref_file + "." + grid + "built.nc"
            newfile = newfile.split("/")
            newdir = (self._work_dir if self._work_dir[-1] ==
                      os.sep else self._work_dir + os.sep) +\
                "AUX_Files_fire_ESACCI"
            newfile = newdir + os.sep + newfile[-1]

            if not os.path.exists(newfile):
                tempfile = self._aggregate_resolution(
                    self._ref_file, grid, remove=False)
                subprocess.call(["mkdir", newdir])
                subprocess.call(['cp', tempfile, newfile])
                os.remove(tempfile)

            self._ref_file = newfile

            reso = self.reso

        self.reso = reso

        proj_var = self._project_info['RUNTIME']['currDiag'].get_variables()[0]
        if proj_var == "burntArea":
            self._ref_data = self._load_cci_generic(self._ref_file, proj_var)
        else:
            assert False, 'Not supported yet'

        # cdo package does not provide correct scaling?!
        # offline cdo provides correct scaling!
        self._ref_data.data = self._ref_data.data * 100
