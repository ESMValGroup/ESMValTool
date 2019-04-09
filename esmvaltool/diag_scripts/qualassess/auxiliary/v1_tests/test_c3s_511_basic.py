import pytest
import os
import sys
import iris

# ESMValTool main.py environments

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '../../../../interface_scripts'))
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '../../../lib/python'))
sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '../../../../diag_scripts'))

from auxiliary import info, error, print_header, ncl_version_check
## from climate import climate
from optparse import OptionParser
import datetime
import projects
# import pdb
import reformat
import xml.sax
import xml_parsers

# Define ESMValTool version
version = "1.1.0"
os.environ['0_ESMValTool_version'] = version

# watermarks (None, default)
os.environ['0_ESMValTool_watermark'] = "None"
###

sys.path.insert(0, os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..'))
from c3s_511_basic import __Diagnostic_skeleton__, Basic_Diagnostic
import warnings


class TestDiagnosticSkeleton:
    def setup(self):
        self.S = __Diagnostic_skeleton__()

    def test_init(self):
        self.S.__init__()
        assert isinstance(self.S.__project_info__, dict)
        assert os.path.isdir(self.S.__plot_dir__)
        assert os.path.isdir(self.S.__work_dir__)
        assert self.S.__varname__ == "var"
        assert self.S.__output_type__ == "png"
        assert self.S.__regions__ == {"example": (10, 20, -10, -20)}
        assert self.S.__verbosity_level__ == 0
        assert self.S.__debug_info__ == "No debug info"
        assert isinstance(self.S.__config__, dict)
        assert isinstance(self.S.__basetags__, list)
        assert isinstance(self.S.__infiles__, list)
        assert isinstance(self.S.authors, list)
        assert isinstance(self.S.diagname, str)
        assert self.S.data is None

    def test_set_info(self):
        with pytest.warns(UserWarning):
            self.S.set_info()

    def test_read_data(self):
        with pytest.warns(UserWarning):
            a = self.S.read_data()

    def test_run_diagnostic(self):
        with pytest.warns(UserWarning):
            a = self.S.run_diagnostic()

    def test_do_overview(self):
        with pytest.warns(UserWarning):
            a = self.S.__do_overview__()

    def test_do_mean_var(self):
        with pytest.warns(UserWarning):
            a = self.S.__do_mean_var__()

    def test_do_trends(self):
        with pytest.warns(UserWarning):
            a = self.S.__do_trends__()

    def test_do_extremes(self):
        with pytest.warns(UserWarning):
            a = self.S.__do_extremes__()

    def test_do_maturity_matrix(self):
        with pytest.warns(UserWarning):
            a = self.S.__do_maturity_matrix__()

    def test_do_gcos_requirements(self):
        with pytest.warns(UserWarning):
            a = self.S.__do_gcos_requirements__()

    def test_prepare_report(self):
        with pytest.warns(UserWarning):
            a = self.S.__prepare_report__()


class TestBasicDiagnostic:
    def setup(self):
        self.S = Basic_Diagnostic()

    def test_init(self):
        self.S.__init__()
        assert "data" in dir(self.S)

    def test_read_data(self):
        self.S.read_data()
        assert isinstance(self.S.data, iris.cube.Cube)

    def test_read_data_mock(self):
        assert False

    def test_read_data_mock(self):
        self.S.read_data_mock()
        assert isinstance(self.S.data, iris.cube.Cube)

    def test_do_report(self):
        flist = ["diag_scripts/aux/C3S_511/example_images/albedo_QA4ECV_all_models_regionalized_smean_ts.png"]
        testdict = {"resolution" : {"dx":"1 deg", "dy":"1 deg"}, "time":{"period":{"start":"2000-01-01", "end":"2000-12-31"}}, "attribute1":"none"}

        do_report(flist, "figure list test")
        do_report(testdict, "dictionary test")


