# -*- coding: utf-8 -*-

# This file is part of ESMValTool


"""
Tests are implemented using *assert* statements
"""

import sys
import os
import glob

import unittest
import tempfile
from nose.tools import assert_raises

class TestLauncher(unittest.TestCase):

    def setUp(self):
        # implement here everything you would like to see happen BEFORE a test is executed

        # to allow that test find the ESMValTool modules, we add here pathes to the system path
        esmval_path = os.path.dirname(os.path.realpath(__file__)) + os.sep + '..' + os.sep
        sys.path.append(esmval_path)

        # temporary directory for output
        self.tmpdir = tempfile.mkdtemp() + os.sep

    def tearDown(self):
        # implement here everything you would like to see happen AFTER a test was executed
        pass

    def test_python_launcher_init(self):
        from interface_scripts.launchers import py_launcher
        L = py_launcher()
        self.assertFalse(L.execute_as_shell)
        self.assertEqual(L.lang, 'PY ')

    def test_python_launcher_exectue_shell(self):
        # execute launcher in SHELL and check if output files existing
        from interface_scripts.launchers import py_launcher

        # generate some dummy python script that should be executed
        # this script will just save some file
        testoutput = self.tmpdir + 'test_result_shell.json'

        script = self.tmpdir + 'testscript.py'
        o = open(script, 'w')
        o.write('import json\n')
        o.write("d = {'A':1,'B':2}\n")
        o.write("json.dump(d, open('" + testoutput + "', 'w'))\n")
        o.close()

        L = py_launcher(execute_as_shell=True)
        project_info = {}
        L.execute(script, project_info, 0, False)
        self.assertTrue(os.path.exists(testoutput))

    def test_python_launcher_execute_script(self):
        # execute launcher in SCRIPT and check if output files existing
        from interface_scripts.launchers import py_launcher
        sys.path.append(self.tmpdir)

        # generate some dummy python script that should be executed
        # this script will just save some file
        testoutput = self.tmpdir + 'test_result_script.json'

        script = self.tmpdir + 'testscript_script.py'
        o = open(script, 'w')
        o.write('import json\n')
        o.write('def main(project_info):\n')
        o.write("    d = {'A':1,'B':2}\n")
        o.write("    json.dump(d, open('" + testoutput + "', 'w'))\n")
        o.close()

        L = py_launcher(execute_as_shell=False)
        project_info = {}
        L.execute(script, project_info, 0, False)
        self.assertTrue(os.path.exists(testoutput))

    def test_python_launcher_missing_script(self):
        # test here that error handling works: script NOT available!
        from interface_scripts.launchers import py_launcher
        L = py_launcher(execute_as_shell=False)
        project_info = {}
        with self.assertRaises(ValueError):
            L.execute('nothing.py', project_info, 0, False)

    def test_python_launcher_execute_script_with_errors(self):
        # execute launcher in SCRIPT and check if output files existing
        from interface_scripts.launchers import py_launcher
        sys.path.append(self.tmpdir)

        # generate some dummy python script that should be executed
        # this script will just save some file

        script = self.tmpdir + 'testscript_script_fail.py'
        o = open(script, 'w')
        o.write('import json\n')
        o.write('def main(project_info):\n')
        o.write("  here we introduce some syntax error this script should then fail")
        o.close()

        L = py_launcher(execute_as_shell=False)
        project_info = {}
        with self.assertRaises(ValueError):
            L.execute(script, project_info, 0, False)




if __name__ == "__main__":
    unittest.main()


