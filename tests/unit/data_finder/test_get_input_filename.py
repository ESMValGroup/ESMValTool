"""Unit tests for :func:`esmvaltool._data_finder.regrid._stock_cube`"""

import unittest
import unittest.mock as mock

from esmvaltool._data_finder import get_input_filename


class TestGetInputDirnameTemplate(unittest.TestCase):
    """Tests for get_start_end_year function"""

    def setUp(self):
        self.cfg = {
            'input_dir': 'input_dir',
            'input_file': 'input_file',
        }
        self.var = {
            'project': 'project',
            'start_year': 2000,
            'end_year': 2001,
        }

    @mock.patch("esmvaltool._data_finder.get_project_config")
    @mock.patch("esmvaltool._data_finder.get_input_dirname_template")
    def test_basic_case(self, get_input_dirname_template, get_project_config):
        """Test case with no drs"""
        get_input_dirname_template.side_effect = \
          lambda var, root_path, drs: [root_path + '/template']
        get_project_config.return_value = self.cfg
        self.assertEqual(
            get_input_filename(
                self.var,
                '/root/path',
                {'project': 'project_path', 'default': 'default_path'}
            ),
            '/root/path/template/input_file'
        )

    @mock.patch("esmvaltool._data_finder.get_project_config")
    @mock.patch("esmvaltool._data_finder.get_input_dirname_template")
    def test_start_end_year(self, get_input_dirname_template,
                                 get_project_config):
        """Test case with no drs"""
        get_input_dirname_template.side_effect = \
          lambda var, root_path, drs: root_path + '/template'
        self.cfg ['input_file'] = 'input_file_*.nc'
        get_project_config.return_value = self.cfg 
        self.assertEqual(
            get_input_filename(
                self.var,
                '/root/path1',
                {'project': 'project_path', 'default': 'default_path'}
            ),
            '/root/path1/template/input_file_200001-200112.nc',
        )
