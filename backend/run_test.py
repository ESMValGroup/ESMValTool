# coding=utf-8
"""
Script to run the tests for EarthDiagnostics and generate the code coverage report
"""
import coverage
import unittest
import os
cov = coverage.Coverage()
cov.set_option("run:branch", True)
cov.start()
import test
suite = unittest.TestLoader().discover('test')
unittest.TextTestRunner(verbosity=2).run(suite)
cov.stop()
cov.save()
source_files = list()
for path, dirs, files in os.walk('backend'):
    for filename in files:
        if filename.endswith('.py') and filename.startswith('test'):
            source_files.append(os.path.join(path, filename))
cov.report(source_files)
cov.html_report(source_files)