"""
Module to run all unittests
The module also checks which namelists are existing 
and analyzes for which namelists tests still need to be implemented
"""


import glob
import os

class RunTests(object):
    def __init__(self):
        self.diagnostics = self._get_diagnostic_list()
        print self.diagnostics

    def _get_diagnostic_list(self):
        nml_dir = self._get_nml_dir()
        L = glob.glob(nml_dir + '*.xml')
        r = []
        for l in L:
            r.append(os.path.basename(l))
        return r

    def _get_nml_dir(self):
        """ get directory with namelists """
        return os.path.dirname(os.path.abspath(__file__)) + os.sep + '..' + os.sep + 'nml' + os.sep

    def test_core(self):
        # supposed to run the tests for the backend
        print('Backend tests not integrated yet')

    def test_diagnostics(self):
        for d in self.diagnostics:
            self._test_diagnostic(d)

    def _test_diagnostic(self, d):
        """
        run test for a single diagnostic

        Parameters
        ----------
        d : str
            name of diagnostic file (namelist filename)
        """
        print('Testing diagnostic: ' + d)



def main():
    R = RunTests()
    #R.test_diagnostics()
    #R.test_core()

if __name__ == '__main__':
    main()




