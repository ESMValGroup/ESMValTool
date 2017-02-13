# coding=utf-8
"""
Script to run the tests and generating the code coverage report
"""
import coverage
import unittest
import os
cov = coverage.Coverage(omit='test/*')
cov.set_option("run:branch", True)

cov.start()
suite = unittest.TestLoader().discover('test')
unittest.TextTestRunner(verbosity=2).run(suite)
cov.stop()

cov.save()
cov.report()
cov.html_report()
