"""Unit tests for :func:`esmvaltool._data_finder.regrid._stock_cube`"""

import unittest

from esmvaltool._data_finder import get_input_dirname_template


class TestGetInputDirnameTemplate(unittest.TestCase):
    """Tests for get_start_end_year function"""

    def setUp(self):
        self.var = {
            'project': 'project',
        }

    def test_basic_case(self):
        """Test case with no drs"""
        self.assertListEqual(
            get_input_dirname_template(self.var, '/root/path', None),
            ['/root/path/var']
        )

if __name__ == '__main__':
    unittest.main()