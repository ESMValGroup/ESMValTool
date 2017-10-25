"""
This modules provides generic functions needed for ESMValTool
performance metrics testing
"""

import os
from xml.dom import minidom

from easytest import EasyTest


class ESMValToolTest(EasyTest):
    """
    main class for all ESMValTool tests
    """

    def __init__(self, **kwargs):
        """
        output_directory : str
            a default output directory is set as ESMVALROOT/work/plots
            but user can specify custom output directory
        """
        self.nml = kwargs.pop('nml', None)
        assert self.nml is not None, 'Namelist needs to be provided!'

        exe = 'python main.py'  # the command to call ESMValTool
        assert exe is not None, 'Executable needs to be given!'

        self.esmval_dir = kwargs.pop('esmval_dir', None)
        assert self.esmval_dir is not None, 'esmval_dir directory needs to be given'

        xmldoc = minidom.parse(os.path.join(self.esmval_dir, self.nml))
        nml_plot_dir = xmldoc.getElementsByTagName('plot_dir')[0].childNodes[0].nodeValue.strip()

        # output directory as defined in nml
        default_output_dir = os.path.join(self.esmval_dir, nml_plot_dir)
        output_directory = kwargs.pop('output_directory', default_output_dir)
        self.refdir_root = self.esmval_dir + 'testdata' + os.sep
        super(ESMValToolTest, self).__init__(exe,
                                             args=[self.nml],
                                             output_directory=output_directory,
                                             checksum_exclude=['ps', 'epsi', 'eps'],
                                             **kwargs)

    def generate_reference_data(self):
        """
        generate reference data by executing the namelist once and then copy
        results to the output directory
        """
        self._execute(wdir=self.esmval_dir)
        if not os.path.exists(self.refdirectory):
            self._copy_output()

    def run_nml(self):
        self._execute(wdir=self.esmval_dir)

    def _copy_output(self):
        """ copy entire result output to reference data directory """
        shutil.copytree(self.output_directory, self.refdirectory)
