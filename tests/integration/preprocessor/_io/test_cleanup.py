"""Integration tests for :func:`esmvalcore.preprocessor._io.cleanup`"""

import os
import tempfile
import unittest

from esmvalcore.preprocessor import _io


class TestCleanup(unittest.TestCase):
    """Tests for :func:`esmvalcore.preprocessor._io.cleanup`"""

    def setUp(self):
        self.temp_paths = []
        descriptor, temp_file = tempfile.mkstemp('.nc')
        os.close(descriptor)
        self.temp_paths.append(temp_file)
        self.temp_paths.append(tempfile.mkdtemp())

    def tearDown(self):
        for path in self.temp_paths:
            if os.path.isfile(path):
                os.remove(path)
            elif os.path.isdir(path):
                os.rmdir(path)

    def test_cleanup(self):
        """Test cleanup"""
        _io.cleanup([], self.temp_paths)
        for path in self.temp_paths:
            self.assertFalse(os.path.exists(path))

    def test_cleanup_when_files_removed(self):
        """Test cleanup works even with missing files or folders"""
        self.tearDown()
        _io.cleanup([], self.temp_paths)
        for path in self.temp_paths:
            self.assertFalse(os.path.exists(path))
