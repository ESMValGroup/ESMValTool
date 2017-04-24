import os
import subprocess
from netCDF4 import Dataset
from ESMValMD import ESMValMD
from diagnostic import BasicDiagnostics


class fAPARDiagnostic(BasicDiagnostics):
    """
    class to implement albedo diagnostics,
    like e.g. global means, global differences, RMSD etc.

    TODO implement testing for this diagnostic

    """

    def __init__(self, **kwargs):
        super(fAPARDiagnostic, self).__init__(**kwargs)

        self._project_info = {}
        self._modtype = None
        self._reftype = None
        self._plot_dir = '.' + os.sep
        self._work_dir = '.' + os.sep

        self._vartype = 'albedo'  # default value as there must be one
        self.output_type = 'png'  # default ouput file type
        self._changed = False

    def run_diagnostic(self):
        """
        running the diagnostics
        """
        super(fAPARDiagnostic, self).run_diagnostic()

    def _specific_diag(self, percentile=True):
        """
        Diagnostic management
        """

    def write_data(self, plot=True):
        """
        write data
        """
        super(fAPARDiagnostic, self).write_data()

    def _load_model_data(self):
        """ load albedo model data """

        edited = False

        newfile = self._mod_file + ".T85built.nc"
        newfile = newfile.split("/")
        newdir = (self._work_dir if self._work_dir[-1] ==
                  os.sep else self._work_dir + os.sep) + \
            "AUX_Files_fAPAR_QA4ECV"
        newfile = newdir + os.sep + newfile[-1]

        mod_info = Dataset(self._ref_file)
        try:
            lat = mod_info.dimensions['lat'].size
            lon = mod_info.dimensions['lon'].size
        except:  # regridding required in any case
            lat = -1
            lon = -1
        mod_info.close()

        if not ((lat == 128 and lon == 256) or
                (lat == 96 and lon == 192) or
                (lat == 18 and lon == 36) or
                (lat == 6 and lon == 12)):
            # TODO add diffs

            if not os.path.exists(newfile):
                tempfile = self._aggregate_resolution(
                    self._mod_file, "T85", remove=False)
                subprocess.call(["mkdir", newdir])
                subprocess.call(['cp', tempfile, newfile])
                os.remove(tempfile)

            self._mod_file_E = newfile
            edited = True

        # load data
        self._mod_data = self._load_cmip_generic(
                self._mod_file_E if edited else self._mod_file,
                self._project_info['RUNTIME']['currDiag'].get_variables()[0])
        self._mod_data.unit = "-"

    def _load_observation_data(self):
        """ load obs data """
        newfile = self._ref_file + ".T85built.nc"
        newfile = newfile.split("/")
        newdir = (self._work_dir if self._work_dir[-1] ==
                  os.sep else self._work_dir + os.sep) + \
            "AUX_Files_fAPAR_QA4ECV"
        newfile = newdir + os.sep + newfile[-1]

        mod_info = Dataset(self._ref_file)
        try:
            lat = mod_info.dimensions['lat'].size
            lon = mod_info.dimensions['lon'].size
        except:  # regridding required in any case
            lat = -1
            lon = -1
        mod_info.close()

        if not ((lat == 128 and lon == 256) or
                (lat == 96 and lon == 192) or
                (lat == 18 and lon == 36) or
                (lat == 6 and lon == 12)):
            # TODO add diffs

            if not os.path.exists(newfile):
                tempfile = self._aggregate_resolution(
                    self._ref_file, "T85", remove=False)
                subprocess.call(["mkdir", newdir])
                subprocess.call(['cp', tempfile, newfile])
                os.remove(tempfile)

            self._ref_file = newfile

        proj_var = self._project_info['RUNTIME']['currDiag'].get_variables()[0]
        if proj_var == "fAPAR":
            self._ref_data = self._load_cci_generic(self._ref_file, proj_var)
        else:
            assert False, 'Not supported yet'
        self._ref_data.unit = "-"
