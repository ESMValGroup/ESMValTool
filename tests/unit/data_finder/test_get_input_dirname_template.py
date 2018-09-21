"""Unit tests for :func:`esmvaltool._data_finder.regrid._stock_cube`"""

import unittest
import unittest.mock as mock

from esmvaltool._data_finder import get_input_dirname_template


class TestGetInputDirnameTemplate(unittest.TestCase):
    """Tests for get_start_end_year function"""

    def setUp(self):
        self.cfg = {
            'input_dir': 'input_dir',
        }
        self.var = {
            'project': 'project',
        }

    @mock.patch("esmvaltool._data_finder.replace_tags")
    @mock.patch("esmvaltool._data_finder.get_project_config")
    def test_basic_case(self, get_project_config, replace_tags):
        """Test case with no drs"""
        replace_tags.side_effect = lambda input_dir, var: input_dir
        get_project_config.return_value = self.cfg
        self.assertListEqual(
            get_input_dirname_template(
                self.var,
                {'project':'/root/path'},
                {'project': 'project_path', 'default': 'default_path'}
            ),
            ['/root/path/input_dir']
        )

    @mock.patch("esmvaltool._data_finder.replace_tags")
    @mock.patch("esmvaltool._data_finder.get_project_config")
    def test_use_default(self, get_project_config, replace_tags):
        """Test case with no drs"""
        replace_tags.side_effect = lambda input_dir, var: input_dir
        get_project_config.return_value = self.cfg
        self.assertListEqual(
            get_input_dirname_template(
                self.var,
                {'default':'/default/path'},
                {'default': 'default_path'}
            ),
            ['/default/path/input_dir']
        )
