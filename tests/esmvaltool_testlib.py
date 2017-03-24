"""
This modules provides generic functions needed for ESMValTool
performance metrics testing
"""

import os
import unittest
from xml.etree import ElementTree

from easytest import EasyTest
from dummydata import Model3, Model2
from wrappers import *
import glob


class ESMValTestDiagnostic(unittest.TestCase):
    def __init__(self, *args, **kwargs):
        super(ESMValTestDiagnostic, self).__init__(*args, **kwargs)
        pass

    def read_reffiles(self, fname, tdir='plot'):
        """
        read reference files from given file
        the routine also ensure that trailing whitespaces
        like e.g. \n are removed.

        Parameters
        ----------
        fname : str
            name of file with required files (one filename per line)
        tdir : str
            specifier for target directory ['plot']
        """
        assert os.path.exists(fname)
        reffiles=[]
        with open(fname, 'rU') as f:
            for line in f:
                reffiles.append((tdir, line.split()[0]))
        return reffiles


class ESMValToolTest(EasyTest):
    """
    main class for all ESMValTool tests
    """

    def __init__(self, **kwargs):
        """
        output_directory : str
            a default output directory is set as ESMVALROOT/work/plots
            but user can specify custom output directory
        """
        self._nml = kwargs.pop('nml', None)
        assert self._nml is not None, 'Namelist needs to be provided!'

        exe = 'python main.py'  # the command to call ESMValTool
        assert exe is not None, 'Executable needs to be given!'

        #~ self.esmval_dir = kwargs.pop('esmval_dir', None)
        #~ assert self.esmval_dir is not None, 'esmval_dir directory needs to be given'

        # directory of current file
        fdir = os.path.dirname(os.path.realpath(__file__))
        self.esmval_dir = fdir + os.sep + '..'



        self._default_dir = os.path.abspath(fdir + os.sep + 'data' + os.sep + os.path.basename(self._nml))
        self._default_output_dir = self._default_dir + os.sep + 'output' + os.sep
        self._default_input_dir = self._default_dir + os.sep + 'input' + os.sep

        self.wrk_dir = self._default_output_dir +  'wrk_dir' + os.sep
        self.plot_dir = self._default_output_dir +  'plot_dir' + os.sep
        self.climo_dir = self._default_output_dir +  'climo_dir' + os.sep

        self._gen_dirs(self.wrk_dir)
        self._gen_dirs(self.plot_dir)
        self._gen_dirs(self.climo_dir)

        # set reference filenames
        self._set_reffiles(kwargs.pop('files', None), subdir=kwargs.get('subdirectory', None))


        #~ xmldoc = minidom.parse(self.nml)
        #~ nml_plot_dir = xmldoc.getElementsByTagName('plot_dir')[0].childNodes[0].nodeValue.strip()

        #~ output_directory = kwargs.pop('output_directory', self._default_output_dir)

        self._modify_nml()
        #~ self.nml = self._nml

        super(ESMValToolTest, self).__init__(exe,
                                             args=[self.nml],
                                             output_directory=self._default_output_dir,
                                             #refdirectory = self._default_input_dir,
                                             files = self.ref_files,
                                             checksum_exclude=['ps', 'epsi', 'eps'],
                                             basedir=self.esmval_dir,
                                             **kwargs)

    def run_tests(self, **kwargs):
        # overwrites routine by Easytest
        # run first all Easytest routines and then run internal tests

        EasyTest.run_tests(self, **kwargs)
        # check that log files are existing and have a size > 0 bytes
        if not self.test_logfile_exists():
            self.sucess = False

    def _size_gt_0(self, f):
        # check that filesize > 0 bytes
        if os.path.getsize(f) > 0.:
            return True
        else:
            print('ERROR: The filesize of the file ' + f + ' is 0 bytes !!!')
            return False

    def test_logfile_exists(self):
        # test if logfile exists and contains content
        files = glob.glob(self.wrk_dir + 'refs-acknows*.log')
        if len(files) == 0:
            return False
        for f in files:
            if not self._size_gt_0(f):
                return False

        return True



    def get_field_definitions(self):
        """
        routine to specify the structure of the sample data to be generated
        needs to be specified in the child class!
        """
        r = {}
        #~ rpath = self._default_input_dir
        #~ # TODO: document variable field
        #~ r.update({'ta' : {'method' : 'uniform', 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        #~ r.update({'pr' : {'method' : 'constant', 'constant' : 20., 'filename' : rpath + os.sep + '@{VAR_FILE}', 'ndim' : 3}})
        return r


    def _set_reffiles(self, x, subdir=None):
        if x is None:
            self.ref_files = None
        else:
            self.ref_files=[]
            for l in x:
                if l[0] == 'plot':
                    if subdir is None:
                        s = self.plot_dir + l[1]
                    else:
                        s = self.plot_dir + subdir + os.sep + l[1]
                elif l[0] == 'work':
                    s = self.wrk_dir + l[1]
                else:
                    assert False, 'Invalid directory specification: ' + l[0]
                self.ref_files.append(s)


    def _gen_dirs(self, d):
        if not os.path.exists(d):
            os.makedirs(d)


    def _get_model_class(self, s):
        # return model class name
        return s.lstrip().split()[0]

    def _set_model_class(self, s):
        """
        store model class as attribute; currently only ONE model class for all models is supported
        """
        mc = self._get_model_class(s).upper()
        if mc == 'CMIP5':
            M = CMIP5(s)
        elif mc == 'CMIP5_fx':
            M = CMIP5fx_w(s)
        elif mc == 'CMIP5_ETHZ':
            M = CMIP5_ETHZ(s)
        elif mc == 'OBS':
            M = OBS(s)
        elif mc == 'OBS_GRIDFILE':
            M = OBS_gridfile(s)
        elif mc == 'OBS4MIPS':
            M = OBS4MIPS(s)
        elif mc == 'ANA4MIPS':
            M = ANA4MIPS(s)
        else:
            assert False, 'Model class not supported yet! ' + mc
        return M




    def _replace_modelpath(self, s):
        o = s.replace('@{MODEL_PATH}', self._default_input_dir)
        o = o.replace('@{MODELPATH}', self._default_input_dir)
        return o


    def _modify_nml(self):
        """
        modify namelist content and store in new list
        """
        #http://stackoverflow.com/questions/6523886/find-and-replace-values-in-xml-using-python

        # read namelist
        #~ xmldoc = minidom.parse(self._nml)
        tree=ElementTree.parse(self._nml)

        # change global attributes
        tree.find('GLOBAL/wrk_dir').text = self.wrk_dir
        tree.find('GLOBAL/plot_dir').text = self.plot_dir
        tree.find('GLOBAL/climo_dir').text = self.climo_dir

        # change attributes in MODEL section
        models = tree.find('MODELS')
        for m in models.getchildren():
            themodel = self._set_model_class(m.text)
            m.text = self._replace_modelpath(m.text)

            # for some diagnostics the direct link to the data file is part of the models tag. This is replaced here (required for sea ice diagnosic)
            m.text = self._replace_datafile_name(m.text, themodel)

            self._generate_input_fields(model=themodel, start_year=themodel.start, stop_year=themodel.stop)

        # change attributes in each of the diagnostic sections
        D=tree.find('DIAGNOSTICS')
        for d in D.getchildren():
            if d.tag == 'diag':
                for x in d.getchildren():
                    if x.tag == 'model':
                        themodel = self._set_model_class(x.text)
                        x.text = self._replace_modelpath(x.text)
                        x.text = self._replace_datafile_name(x.text, themodel)
                        self._generate_input_fields(model=themodel, start_year=themodel.start, stop_year=themodel.stop)


        # write results
        self.nml = self._default_dir + os.sep + os.path.basename(self._nml)[:-4] + '_new.xml'
        tree.write(self.nml)

    def _replace_datafile_name(self, s, model):

        replaced=None
        if True:
            # limitation: currently only name of FIRST variable will be used ???
            var_defs = self.get_field_definitions()  # needs to come from child class
            for k in var_defs.keys():
                oname = var_defs[k].get('filename')
                oname = model.get_oname(oname, k, dim=var_defs[k]['ndim'])
                if 'DATAFILE' in s:
                    s = s.replace('@{DATAFILE}', oname)
                    replaced = True

        #~ if replaced is not None:
            #~ print '****THETEXT: ', s
            #~ assert False
        return s


    def _generate_input_fields(self, force=False, model=None, start_year=None, stop_year=None, size='tiny'):
        """
        routine which generates input fields
        requires a wrapper function from the child class
        """
        var_defs = self.get_field_definitions()  # needs to come from child class
        cnt = 1
        assert start_year is not None
        assert stop_year is not None
        for k in var_defs.keys():
            method = var_defs[k].get('method')
            constant = var_defs[k].get('constant', None)
            oname = var_defs[k].get('filename')

            if var_defs[k]['ndim'] == 2:
                M = Model2
            elif var_defs[k]['ndim'] == 3:
                M = Model3
            else:
                assert False

            # replace tags in output filename according to model class (e.g. CMIP5)
            if model is not None:
                oname = model.get_oname(oname, k, dim=var_defs[k]['ndim'])
            else:
                assert False, 'Should probably not happen'


            # generate ouput directory if not existing yet
            if not os.path.exists(os.path.dirname(oname)):
                os.makedirs(os.path.dirname(oname))

            if force:
                if os.path.exists(oname):
                    os.remove(oname)
            if not os.path.exists(oname):
                X = M(method=method, constant=constant, oname=oname, var=k, start_year=start_year, stop_year=stop_year, size=size, append_coordinates=model._append_coordinates, append_cellsize=model._append_cellsize)
                del X
            cnt += 1


    def generate_reference_data(self):
        """
        generate reference data by executing the namelist once and then copy
        results to the output directory
        """
        self._execute(wdir=self.esmval_dir)
        if not os.path.exists(self.refdirectory):
            self._copy_output()

    def run_nml(self):
        self._execute()

    def _copy_output(self):
        """ copy entire result output to reference data directory """
        shutil.copytree(self.output_directory, self.refdirectory)
