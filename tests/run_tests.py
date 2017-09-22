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
        self.check_test_availability()
        self.info()

    def info(self):
        print('*** EXISTING TESTS ***')
        print('For the following diagnostics tests are available:')
        for i, d in enumerate(self.valid_diagnostics):
            print('%i. %s' % (i,d))

        print('')
        print('*** MISSING TESTS ***')
        print('For the following diagnostics tests are missing:')
        for i, d in enumerate(self.missing_diagnostics):
            print('%i. %s' % (i,d))

    def check_test_availability(self):
        """ checks if tests are available and generates two lists """
        self.valid_diagnostics = []
        self.missing_diagnostics = []
        for d in self.diagnostics:
            if self._test_exists(d):
                self.valid_diagnostics.append(d)
            else:
                self.missing_diagnostics.append(d)
        self.valid_diagnostics.append('python_dummy')
        assert len(self.missing_diagnostics) + len(self.valid_diagnostics) == len(self.diagnostics) + 1


    def _get_test_dir(self, d):
        return os.path.dirname(os.path.abspath(__file__)) + os.sep + 'test_diagnostics' + os.sep + 'test_' + os.path.splitext(d)[0]

    def _test_exists(self, d):
        fname = self._get_test_dir(d) + os.sep +  d
        return os.path.exists(fname)

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

    def test_diagnostics(self, ignore=[]):
        self.results = {}
        for d in self.valid_diagnostics:
            if d not in ignore:
                r = self._test_diagnostic(d)
                self.results.update({d : r})

    def _test_diagnostic(self, d):
        """
        run test for a single diagnostic

        Parameters
        ----------
        d : str
            name of diagnostic file (namelist filename)
        """
        print('Testing diagnostic: ' + d)
        
        # at the moment, the command line is used to launch the test. At a later stage this could
        # be done differently, also replacing nosetests with pytest would be a good idea
        cmd = 'nosetests ' + self._get_test_dir(d)
        r = os.system(cmd)
        if r == 0:
            return True
        else:
            return False

    def print_results(self):
        print('')
        print('*********')
        print('Results of diagnostic testing')
        for d in self.results.keys():
            print(d, self.results[d])

    def check_assert(self):
        """
        issue assert in case that tests failed
        """
        for d in self.results.keys():
            if not self.results[d]:
                assert False, 'Some diagnostics failed!'


def main():
    R = RunTests()

    # specify here currently diagnostics that should be irgnore for testing
    # if doing so, always specify a reason why this is done
    ignore=[]
    #ignore.append('namelist_Evapotranspiration.xml')  # NCL errors happen!


    R.test_diagnostics(ignore=ignore)
    #R.test_core()
    R.print_results()
    R.check_assert()

if __name__ == '__main__':
    main()




