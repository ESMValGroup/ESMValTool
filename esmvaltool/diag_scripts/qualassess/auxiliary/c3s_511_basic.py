"""
Basic implementation for diagnostics into ESMValTool
"""
# used modules
import inspect
import iris
import iris.coord_categorisation
from iris.analysis import stats as cubestats
import iris.coords as icoords
import cf_units as unit
import os
import sys
import shutil
import errno
import numpy as np
import random, string
import collections
import matplotlib
matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
from scipy import stats
import datetime
from .libs import c3s_511_util as utils
from .libs.customErrors import ImplementationError, ConfigurationError, PathError, EmptyContentError
import warnings
from .libs.reporting import do_report as report
from .plots.matrices import do_smm_table
from .plots.matrices import do_gcos_table
from .plots.basicplot import Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup
from .libs.MD_old.ESMValMD import ESMValMD
import logging
from pprint import pprint
from .libs.predef.ecv_lookup_table import ecv_lookup

logger = logging.getLogger(os.path.basename(__file__))

class __Diagnostic_skeleton__(object):
    """
    Basic class to implement any kind of diagnostic
    """

    def __init__(self, **kwargs):    
        super(__Diagnostic_skeleton__, self).__init__(**kwargs)
        """
        Default values to experiment with the diagnostics
        """
        
        self.diagname = "Diagnostic_skeleton.py"
        
        # config
        self.__cfg__ = None
        self.__logger__ = logging.getLogger(os.path.basename(__file__))

        self.__project_info__ = None  # empty project info
        self.__plot_dir__ = '.' + os.sep  # default plot directory
        self.__work_dir__ = '.' + os.sep  # default work dir

        self.__varname__ = 'var'  # default value
        self.__output_type__ = 'png'  # default ouput file type
        self.__regions__ = {}  # default regions

        self.__verbosity_level__ = 0  # default information during runtime
        self.__debug_info__ = "No debug info"  # default debug information
        self.__config__ = dict()  # default configuration input

        self.__basetags__ = []
        self.__infile__ = None
        self.__inpath__= None
        
        self.__time_period__ = None
        self.__dataset_id__ = None
        
        self.authors = ["A_muel_bn", "A_hass_bg", "A_laue_ax",
                        "A_broe_bj", "A_mass_fr", "A_nico_nd",
                        "A_schl_mn", "A_bock_ls"]  # TODO fill in
        self.diagname = "Diagnostic_skeleton"
        self.CDS_ID = "CDS_ID_needed"
        
        self.colormaps=dict({"default":"jet"})#"binary"})
        self.__latex_output__ = False
        self.levels = None

        self.sp_data = None
        self.mp_data = None
        self.var3D = None
        self.__basic_filename__ = "someFile"
        self.__dimensions__ = []
        self.level_dim = None
        self.__time_read__ = []
        self.__avg_timestep__ = None
    
    def set_info(self, **kwargs):
        self.__cfg__ = kwargs.get('cfg', None)
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        raise ImplementationError("set_info","This method has to be implemented.")
        return

    def read_data(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        raise ImplementationError("read_data","This method has to be implemented.")
        return

    def run_diagnostic(self):
        self.__logger__.debug("There is nothing to run! This is the empty skeleton.")
        warnings.warn("Implementation Warning", UserWarning)
        pass
    
    def __do_overview__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={},filename="do_overview_default")
#        raise ImplementationError("__do_overview__","This method has to be implemented.")
        return

    def __do_mean_var__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={},filename="do_mean_var_default")
        raise ImplementationError("__do_mean_var__","This method has to be implemented.")
        return

    def __do_trends__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={},filename="do_trends_default")
        raise ImplementationError("__do_trends__","This method has to be implemented.")
        return

    def __do_extremes__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + "is not implemented!")
        self.__do_report__(content={},filename="do_extremes_default")
        warnings.warn("Implementation Warning", UserWarning)
        return
    
    def __do_sectors__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={},filename="do_sectors_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_maturity_matrix__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={},filename="do_maturity_matrix_default")
        raise ImplementationError("__do_maturity_matrix__","This method has to be implemented.")
        return
    
    def __do_app_perf_matrix__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={},filename="do_app_perf_matrix_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_gcos_requirements__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={},filename="do_gcos_requirements_default")
        raise ImplementationError("__do_gcos_requirements__","This method has to be implemented.")
        return
    
    def __do_esm_evaluation__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={},filename="do_esmevaluation_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

#    def __do_report__(self,*args,**kwargs):
#        this_function = inspect.currentframe().f_code.co_name
#        self.__logger__.error(this_function + " is not implemented!")
#        raise ImplementationError("__do_report__","This method has to be implemented.")
#        return
    
    def __do_report__(self,**kwargs):
        """
        reporting function for use with single diagnostics
        """
        content = kwargs.get('content', [])
        if not isinstance(content, (list, dict)):
            raise TypeError("content", "Element is not a list, nor a dict.")
            
        rand_str = lambda n: ''.join([random.choice(string.ascii_lowercase) for i in range(n)])
            
        filename = kwargs.get('filename', rand_str(10))
        if not isinstance(filename, str):
            raise TypeError("filename", "Element is not a string.")
            
        report(content,
               filename,
               self.__work_dir__,
               ecv = ecv_lookup(self.__varname__),
               dataset = "".join(str(self.__dataset_id__[0])),
               signature = self.CDS_ID,
               latex_opts=self.__latex_output__)
        return
    
    def __file_anouncement__(self,subdir="input",expfile="end.txt",protofile="none.txt",function="no_fun"):
        
        expected_input = self.__work_dir__ + os.sep + subdir + os.sep + ("_".join(self.__basic_filename__) if type(self.__basic_filename__) is list else self.__basic_filename__) + expfile
        if not os.path.isfile(expected_input):
            try:
                os.makedirs(os.path.dirname(expected_input))
            except OSError as exc: # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
            shutil.copy2(os.path.dirname(os.path.realpath(__file__)) + "/libs/predef/" + protofile, expected_input)
            #TODO Improve formatting? Was nicer in V1
            self.__logger__.info("************************************** WARNING **************************************")
            self.__logger__.info("Expected " + function + " input file " + expected_input + " not found!")
            self.__logger__.info("Created dummy file instead. Please fill in appropriate information and rerun!")
            self.__logger__.info("(This won't fail if you do not do so and produce an empty output!!!!)")
            self.__logger__.info("************************************** WARNING **************************************")

            found = False
        else:
            self.__logger__.info("Processing " + function + " input file: " + expected_input) 
            found = True
            
        return expected_input, found

class Basic_Diagnostic_SP(__Diagnostic_skeleton__):
    """
    class to implement basic diagnostics, like e.g. global means,
    global differences, RMSD etc.
    """

    def __init__(self, **kwargs):
        super(Basic_Diagnostic_SP, self).__init__(**kwargs)
        self.diagname = "Basic_Diagnostic_SP.py"
        
        # save GCOS requirements
        self.__gcos_dict__ = dict()
        self.__gcos_dict__.update({"Frequency":{"value":None, "unit":None}})
        self.__gcos_dict__.update({"Resolution":{"value":None, "unit":None}})
        self.__gcos_dict__.update({"Accuracy":{"value":None, "unit":None}})
        self.__gcos_dict__.update({"Stability":{"value":None, "unit":None}})

    def set_info(self, **kwargs):
        """
        gather information for diagnostic
        """

        self.__cfg__ = kwargs.get('cfg', None)
        
        #######################################################################
        #TODO delete after debugging
#        self.__logger__.info("Object content pre:")
        self.li3 = [method_name for method_name in dir(self) if callable(getattr(self, method_name))].copy()
#        self.__logger__.info(np.sort(self.li3))
        self. li1 = self.__dict__.copy().keys()
#        self.__logger__.info(np.sort([*self.li1]))
        #######################################################################        

        if not isinstance(self.__cfg__, dict) or len(self.__cfg__)==0:
            raise EmptyContentError("__cfg__", "Element is empty.")
            
        try:
            self.__latex_output__ = self.__cfg__.show_latex
        except:
            pass
        
        try:
            self.levels = self.__cfg__.levels
        except:
            pass

        self.__plot_dir__ = self.__cfg__['plot_dir']
        self.__work_dir__ = self.__cfg__['work_dir']

        self.__infile__ = list(self.__cfg__['input_data'].keys())
        if not len(self.__infile__) == 1:
            raise ConfigurationError("self.__infile__", "There should be only one infile!")
        self.__infile__ = self.__infile__[0]

        fileinfo = self.__cfg__['input_data'][self.__infile__]
        self.__varname__ = fileinfo['short_name']
#
#        self.__output_type__ = self.__cfg__['output_file_type']  # default ouput file type for the basic diagnostics
#        if not self.__output_type__ == 'png':
#            raise ConfigurationError("self.__output_type__", "Only png is currently supported.")
#        self.__regions__ = {'Germany_2001-2005':{'latitude':(47,56),'longitude':(5,16),'time':(datetime.datetime(2001,1,1),datetime.datetime(2005,12,31))}}  # default regions
#
#        # for metadata
#        self.__basetags__ = [self.__varname__] + []# TODO transport tags from namelist
#
        self.__time_period__ = "-".join([str(fileinfo['start_year']),str(fileinfo['end_year'])])
        
        try:
            # model
            self.__dataset_id__ = [fileinfo["cmor_table"], fileinfo["dataset"], fileinfo["mip"], fileinfo["exp"], fileinfo["ensemble"], fileinfo["short_name"]]
        except:
            # obs
            self.__dataset_id__ = [fileinfo] # TODO adjust to OBS
#
#        self.__basic_filename__ = "_".join(self.__dataset_id__ + [self.__time_period__])
#
        self.__dimensions__ = np.array(["time", "latitude", "longitude"]) # TODO: get from cube
#
#        # TODO: for testing purpose (should come from CDS)
#        self.CDS_ID = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
        self.__output_type__ = self.CDS_ID + "." + self.__output_type__

    def read_data(self):

        """
        reads data
        """ 
        try:
            if os.path.isfile(self.__infile__):
                self.sp_data = iris.load_cube(self.__infile__)
                # get dimensions and correct if needed
                sp_dimensions = [c.name() for c in self.sp_data.coords()]
                aux_dimensions = [c.name() for c in self.sp_data.aux_coords]
                sp_dimensions = [item for item in sp_dimensions if item not in aux_dimensions]
                extra_dimensions = [item for item in sp_dimensions if item not in self.__dimensions__]
                self.__dimensions__ = sp_dimensions
                if not self.sp_data.coord('latitude').has_bounds():
                    self.sp_data.coord('latitude').guess_bounds()
                if not self.sp_data.coord('longitude').has_bounds():
                    self.sp_data.coord('longitude').guess_bounds()
                if self.sp_data.units == "no-unit":
                    self.sp_data.units = '1'
            else:
                self.__logger__.error("self.__infile__ is not set to valid path.")
                raise PathError(self.__infile__,"This file is not accessible.")
        except:
            self.__logger__.error("self.__infile__ is not readable.")
            raise PathError("Basic_Diagnostic_SP.__init__", "self.__infile__ is not readable: " + self.__infile__)
            
        if len(extra_dimensions) == 0:
            self.var3D=False
            self.level_dim = None
        elif len(extra_dimensions) == 1:
            self.var3D=True 
            self.level_dim = extra_dimensions[0]
            self.__logger__.info("available levels are: " + str(self.sp_data.coord(self.level_dim)))
        else:
            assert False, "Error in data dimensions. Too many dimensions!"
        
        # Human readable time   
        time_coord = self.sp_data.coord("time")
        
        time_read = unit.num2date(time_coord.points,time_coord.units.name,time_coord.units.calendar)
        
        self.__time_read__ = time_read
        
        return
    
    def run_diagnostic(self):
        self.__do_overview__()
        self.__do_mean_var__()
        self.__do_trends__()
#        self.__do_extremes__()
#        self.__do_sectors__()
#        self.__do_maturity_matrix__()
#        self.__do_gcos_requirements__()
#        self.__do_esm_evaluation__()
#        self.__do_app_perf_matrix__()
        
        
        #######################################################################
        #TODO delete after debugging
#        self.__logger__.info("Object content post:")
        li4 = [method_name for method_name in dir(self) if callable(getattr(self, method_name))]
#        self.__logger__.info(np.sort(li4))
        li2 =self.__dict__.keys()
#        self.__logger__.info(np.sort([*li2]))

        self.__logger__.info("Object content diff:")        
#        self.__logger__.info(np.sort([x for x in [*li2] if x not in [*self.li1] + ['li1']]))
        self.__logger__.info({y:self.__dict__[y] for y in np.sort([x for x in [*li2] if x not in [*self.li1] + ['li1']])})
        self.__logger__.info(np.sort([x for x in li4 if x not in self.li3]))
        #######################################################################
        pass

    def __do_overview__(self):
        
        this_function = "overview"
        
        if not self.var3D:
            list_of_plots=self.__overview_procedures_2D__(cube=self.sp_data)
        else: 
            list_of_plots=[]
            for lev in self.levels:
                loc_cube=self.sp_data.extract(iris.Constraint(coord_values={str(self.level_dim):lambda cell: cell==lev}))
                if loc_cube is None:
                    raise ValueError("Invalid input: at least one of " + 
                                     "the requested levels seem " + 
                                     "not to be available in the analyzed " + 
                                     "cube.")
                lop=self.__overview_procedures_2D__(loc_cube,level=lev)
                list_of_plots = list_of_plots + lop
        
        # dimension information
        lon_range = self.sp_data.coord("longitude").points
        lat_range = self.sp_data.coord("latitude").points
        tim_range = self.sp_data.coord("time").points
        t_info = str(self.sp_data.coord("time").units)
        
        lon_range_spec = utils.__minmeanmax__(lon_range)
        lat_range_spec = utils.__minmeanmax__(lat_range)
        tim_range_spec = utils.__minmeanmax__(tim_range)
        
        tim_range_spec_read = unit.num2date(tim_range_spec,self.sp_data.coord("time").units.name,self.sp_data.coord("time").units.calendar)

        lon_freq = np.diff(lon_range)
        lat_freq = np.diff(lat_range)
        tim_freq = np.diff(tim_range)
        
        lon_freq_spec = utils.__minmeanmax__(lon_freq)
        lat_freq_spec = utils.__minmeanmax__(lat_freq)
        tim_freq_spec = utils.__minmeanmax__(tim_freq)
        
        overview_dict=collections.OrderedDict()
        
        overview_dict.update({'longitude range [' + str(self.sp_data.coord("longitude").units) + ']': collections.OrderedDict([("min",str(lon_range_spec[0])), ("max",str(lon_range_spec[2]))])})
        overview_dict.update({'longitude frequency [' + str(self.sp_data.coord("longitude").units) + ']': collections.OrderedDict([("min",str(lon_freq_spec[0])), ("average",str(lon_freq_spec[1])), ("max",str(lon_freq_spec[2]))])})
        overview_dict.update({'latitude range [' + str(self.sp_data.coord("latitude").units) + ']': collections.OrderedDict([("min",str(lat_range_spec[0])), ("max",str(lat_range_spec[2]))])})
        overview_dict.update({'latitude frequency [' + str(self.sp_data.coord("latitude").units) + ']': collections.OrderedDict([("min",str(lat_freq_spec[0])), ("average",str(lat_freq_spec[1])), ("max",str(lat_freq_spec[2]))])})
        overview_dict.update({'temporal range': collections.OrderedDict([("min",str(tim_range_spec_read[0])), ("max",str(tim_range_spec_read[2]))])})
        overview_dict.update({'temporal frequency [' + t_info.split(" ")[0] + ']': collections.OrderedDict([("min",str(tim_freq_spec[0])), ("average",str(round(tim_freq_spec[1],2))), ("max",str(tim_freq_spec[2]))])})
        
        try:
            lev_range = self.sp_data.coord(self.level_dim).points
            lev_range_spec = utils.__minmeanmax__(lev_range)
            lev_freq_spec = utils.__minmeanmax__(np.diff(lev_range))
            overview_dict.update({'level range [' + str(self.sp_data.coord(self.level_dim).units) + ']': collections.OrderedDict([("min",str(lev_range_spec[0])), ("max",str(lev_range_spec[2]))])})
            overview_dict.update({'level frequency [' + str(self.sp_data.coord(self.level_dim).units) + ']': collections.OrderedDict([("min",str(lev_freq_spec[0])), ("average",str(lev_freq_spec[1])), ("max",str(lev_freq_spec[2]))])})
        except:
            pass    
        
        # save GCOS requirements
        self.__gcos_dict__.update({"Frequency":{"value":round(tim_freq_spec[1],2), "unit":t_info.split(" ")[0]}})
        self.__gcos_dict__.update({"Resolution":{"value":round(np.mean([lon_freq_spec[1],lat_freq_spec[1]]),2), "unit":str(self.sp_data.coord("longitude").units)}})
        
        # save values for other diagnostics
        self.__avg_timestep__ = tim_freq_spec[1]
        
        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_overview_input",
                                      expfile="_overview.txt",
                                      protofile="empty.txt",
                                      function=this_function)
            
        if found:    
            self.__do_report__(content={"listtext":overview_dict,"plots":list_of_plots,"freetext":expected_input}, filename=this_function.upper())
        else:
            self.__do_report__(content={"listtext":overview_dict,"plots":list_of_plots}, filename=this_function.upper())

        return
    
    def __overview_procedures_2D__(self,cube=None,level=None):
        
        if cube is None:
            cube=self.sp_data
            
        if level is not None:
            basic_filename = self.__basic_filename__ + "_lev" + str(level)
            dataset_id = [self.__dataset_id__[0],"at","level",str(level),str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[0]]
        
        list_of_plots=[]
        
        reg_dimensions =  [item for item in self.__dimensions__ if item not in set([self.level_dim])]
        maxnumtemp = float(len(cube.coord("time").points))

        try:
            sp_masked_vals = cube.collapsed("time",iris.analysis.COUNT,function=lambda values: values.mask)
            sp_masked_vals.data.mask = sp_masked_vals.data==maxnumtemp
            sp_masked_vals.data = sp_masked_vals.data*0.+1.
        except:
            sp_masked_vals = cube.collapsed("time",iris.analysis.MEAN) *0.+1. #errors when data is not correctly masked!

        for d in reg_dimensions:
            
            long_left_over = [rd for rd in reg_dimensions if rd!=d]
            short_left_over = np.array([sl[0:3] for sl in long_left_over])
                        
            num_available_vals = cube.collapsed(d,iris.analysis.COUNT,function=lambda values: values>=cube.collapsed(reg_dimensions,iris.analysis.MIN).data)
            if d in ["time"]:
                frac_available_vals = iris.analysis.maths.divide(num_available_vals,maxnumtemp)
            else:
                sp_agg_array=sp_masked_vals.collapsed(d,iris.analysis.SUM)   
                frac_available_vals = iris.analysis.maths.divide(num_available_vals,sp_agg_array.data)
            
            try:
                # plotting routine
                filename = self.__plot_dir__ + os.sep + basic_filename + "_frac_avail_" + "_".join(short_left_over) + "." + self.__output_type__
                list_of_plots.append(filename)
                frac_available_vals.rename("availability as a fraction")
                x=Plot2D(frac_available_vals)
                
                fig = plt.figure()
                (fig,ax,_) = plot_setup(d=d,fig=fig)
                x.plot(ax=ax, vminmax=[0.,1.], color=self.colormaps, color_type="Sequential", title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
                fig.savefig(filename)
                plt.close(fig.number)
                
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + ['DM_global', 'C3S_overview'],
                         str('Overview on ' + "/".join(long_left_over) + ' availablility of ' + ecv_lookup(self.__varname__) + ' for the data set ' + " ".join(dataset_id) + ' (' + self.__time_period__ + '). NA-values are shown in grey.'),
                         '#C3S' + 'frav' + "".join(short_left_over) + self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors)
                
            except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    self.__logger__.error(exc_type, fname, exc_tb.tb_lineno)
                    self.__logger__.error("Availability in " + d)
                    self.__logger__.error('Warning: blank figure!')
                    
                    x=Plot2D_blank(frac_available_vals)
                        
                    fig = plt.figure()
                    
                    (fig,ax,_) = plot_setup(d=d,fig=fig)
                    x.plot(ax=ax, color=self.colormaps, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
                    fig.savefig(filename)
                    plt.close(fig.number)
                    
                    ESMValMD("meta",
                             filename,
                             self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                             str('Availability plot for the data set "' + "_".join(dataset_id) + '" (' + self.__time_period__ + '); Data can not be displayed due to cartopy error!'),
                             '#C3S' + 'frav' + "".join(short_left_over) + self.__varname__,
                             self.__infile__,
                             self.diagname,
                             self.authors)
            
        del sp_masked_vals
        del num_available_vals
        del frac_available_vals
        
        # histogram plot of available measurements
        all_data = cube.copy()
         
        try:
            # plotting routine
            filename = self.__plot_dir__ + os.sep + basic_filename + "_hist_all_vals" + "." + self.__output_type__
            list_of_plots.append(filename)
            x=PlotHist(all_data)
            fig = x.plot(title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_overview'],
                     str('Full spatio-temporal histogram of ' + self.__varname__ + ' for the data set ' + " ".join(dataset_id) + ' (' + self.__time_period__ + ').'),
                     '#C3S' + 'histall' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
        except:
            self.__logger__.error('Something is probably wrong with the plotting routines/modules!')
        
        del all_data
        
        return list_of_plots    

    def __do_mean_var__(self):
        
        this_function = "mean and variability"
        
        if not self.var3D:
            list_of_plots=self.__mean_var_procedures_2D__(cube=self.sp_data)
        else: 
            list_of_plots=[]
            for lev in self.levels:
                loc_cube=self.sp_data.extract(iris.Constraint(coord_values={str(self.level_dim):lambda cell: cell==lev}))
                lop=self.__mean_var_procedures_2D__(loc_cube,level=lev)
                list_of_plots = list_of_plots + lop
            lop=self.__mean_var_procedures_3D__()
            list_of_plots = list_of_plots + lop

        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_mean_var_input",
                                      expfile="_mean_var.txt",
                                      protofile="empty.txt",
                                      function=this_function)
            
        if found:    
            self.__do_report__(content={"plots":list_of_plots,"freetext":expected_input}, filename=this_function.upper())
        else:
            self.__do_report__(content={"plots":list_of_plots}, filename=this_function.upper())

        return

    def __mean_var_procedures_2D__(self,cube=None,level=None):
                
        if cube is None:
            cube=self.sp_data
            
        if level is not None:
            basic_filename = self.__basic_filename__ + "_lev" + str(level)
            dataset_id = [self.__dataset_id__[0],"at","level",str(level),str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[0]]
        
        list_of_plots = []
        
        reg_dimensions =  [item for item in self.__dimensions__ if item not in set([self.level_dim])]

        
        maths = ["MEAN",
                 "STD_DEV",
#                 "LOG_COV",
                 "PERCENTILE",
                 "CLIMATOLOGY"]
        
        percentiles = [1.,5.,10.,50.,90.,95.,99.]
        
        for d in reg_dimensions:
            
            long_left_over = [rd for rd in reg_dimensions if rd!=d]
            short_left_over = np.array([sl[0:3] for sl in long_left_over])
        
            if d in ["time"]:
                temp_series_1d = cube.collapsed(long_left_over, iris.analysis.MEAN, weights=iris.analysis.cartography.area_weights(cube))
        
                filename = self.__plot_dir__ + os.sep + basic_filename + "_" + "temp_series_1d" + "." + self.__output_type__
                list_of_plots.append(filename)
                
                x=Plot1D(temp_series_1d)
                        
                fig = plt.figure()

                ax = [plt.subplot(1,1,1)]
                fig.set_figheight(1.2*fig.get_figheight())
                x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
                fig.savefig(filename)
                plt.close(fig.number)

                        
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                         str("/".join(long_left_over).title() + ' aggregated ' + d + ' series' + ' values of ' + ecv_lookup(self.__varname__) + ' for the data set ' + " ".join(dataset_id) + ' (' + self.__time_period__ + ').'),
                         '#C3S' + d + 'series' + "".join(short_left_over) + self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors)
                
                del temp_series_1d
        
            mean_std_cov=collections.OrderedDict()
            disp_min_max=collections.OrderedDict()
            disp_min_max.update({"abs_vals":np.array([np.nan])})
            disp_min_max.update({"diff_vals":np.array([np.nan])})

            for m in maths:
                
                if m == "PERCENTILE":
                    
                    try:
                        perc = cube.collapsed(d, iris.analysis.__dict__["W" + m], percent=percentiles,weights=iris.analysis.cartography.area_weights(cube))
                        
                        precentile_list=list()
                        
                        for p in percentiles:
                            
                            loc_data = perc.extract(iris.Constraint(weighted_percentile_over_time=p))
                            
                            precentile_list.append(loc_data)
                               
                            try:
                                disp_min_max.update({"abs_vals":np.nanpercentile(np.concatenate([disp_min_max["abs_vals"],np.percentile(loc_data.data.compressed(),[5,95])]),[0,100])})
                            except:
                                disp_min_max.update({"abs_vals":np.nanpercentile(np.concatenate([disp_min_max["abs_vals"],np.percentile(loc_data.data,[5,95])]),[0,100])})

                        mean_std_cov.update({m:precentile_list})
                            
                        del precentile_list
                        del perc
                    except:
                        pass
                    
                elif m == "CLIMATOLOGY":
                    
                    try:
                        clim = cube.copy()
                        iris.coord_categorisation.add_month_number(clim, d, name='month_num')
                        clim_comp = clim.aggregated_by('month_num',iris.analysis.MEAN)
                        
                        clim_anom = clim.copy()
                        
                        del clim
                        
                        for mn in range(len(clim_anom.coord('month_num').points)):
                        
                            idx = clim_comp.coord('month_num').points.tolist().index(clim_anom.coord('month_num').points[mn])
                            clim_anom.data[mn,:,:] = clim_anom.data[mn,:,:]-clim_comp.data[idx,:,:]
                        
                        iris.coord_categorisation.add_year(clim_anom, d, name='yr')
                        clim_anom = clim_anom.aggregated_by('yr', iris.analysis.MEAN)
                        
                        clim_comp_list=list()
                        
                        for mon in np.sort(clim_comp.coord('month_num').points):
                    
                            loc_data=clim_comp.extract(iris.Constraint(month_num=mon))
                            
                            clim_comp_list.append(loc_data)
                            
                            try:
                                disp_min_max.update({"abs_vals":np.nanpercentile(np.concatenate([disp_min_max["abs_vals"],np.nanpercentile(loc_data.data.compressed(),[5,95])]),[0,100])})
                            except:
                                disp_min_max.update({"abs_vals":np.nanpercentile(np.concatenate([disp_min_max["abs_vals"],np.nanpercentile(loc_data.data,[5,95])]),[0,100])})

                        
                        mean_std_cov.update({m:clim_comp_list})
                        
                        del clim_comp
                        del clim_comp_list
                        
                        clim_anom_list=list()
                        
                        year_list = np.sort(clim_anom.coord('yr').points)
                        
                        for y in year_list:
                            
                            loc_data=clim_anom.extract(iris.Constraint(yr=y))
                            
                            try:
                                pot_min_max = np.concatenate(
                                        [disp_min_max["diff_vals"],
                                        np.nanpercentile(loc_data.data.compressed(),[5,95])
                                        ])
                            except:
                                pot_min_max = np.concatenate(
                                        [disp_min_max["diff_vals"],
                                        np.nanpercentile(loc_data.data,[5,95])
                                        ])
                                
                            disp_min_max.update({"diff_vals":list(set(np.nanpercentile(
                                    np.concatenate([pot_min_max,-1.*pot_min_max])
                                    ,[0,100])))})

                            if len(disp_min_max["diff_vals"])==1:
                                disp_min_max["diff_vals"]*=2
    
                            clim_anom_list.append(loc_data)
                            
                        mean_std_cov.update({"mean anomalies from " + m:clim_anom_list})
                                                            
                        del clim_anom
                        del clim_anom_list
                            
                    except:
                        pass
                
                elif m == "MEAN":
                    
                    loc_data=cube.collapsed(d, iris.analysis.__dict__[m],weights=iris.analysis.cartography.area_weights(cube))
        
                    mean_std_cov.update({m:loc_data})
                    
                    disp_min_max.update({"abs_vals":np.nanpercentile(np.concatenate([disp_min_max["abs_vals"],np.nanpercentile(loc_data.data.compressed(),[5,95])]),[0,100])})
                    
                elif m == "STD_DEV":
        
                    mean_std_cov.update({m:utils.weighted_STD_DEV(cube,d,weights=iris.analysis.cartography.area_weights(cube))})
                                 
                elif m == "LOG_COV": 
                    
                    mean_std_cov.update({m:iris.analysis.maths.log(iris.analysis.maths.divide(mean_std_cov["STD_DEV"],mean_std_cov["MEAN"]))})
         
                else:
                    raise ConfigurationError(m,"This maths functionality has to be implemented in _do_mean_var_.") 
            
            for m in mean_std_cov.keys():
                
                # plotting routine
                vminmax=None
                ctype=None
                
                if np.any([m_typ in m for m_typ in ["MEAN","PERCENTILE", "CLIMATOLOGY"]]):
                    vminmax=disp_min_max["abs_vals"]
                    ctype="Data"
                    
                if np.any([m_typ in m for m_typ in ["STD_DEV"]]):
                    ctype="Sequential"
                
                if np.any([m_typ in m for m_typ in ["anomalies"]]):
                    vminmax=disp_min_max["diff_vals"]
                    ctype="Diverging"
                
                if mean_std_cov[m] is not None:
                # this needs to be done due to an error in cartopy
                    filename = self.__plot_dir__ + os.sep + basic_filename + "_" + "_".join(m.split(" ") )+ "_" + "_".join(short_left_over) + "." + self.__output_type__
                    list_of_plots.append(filename)
                    
                    try:
                        x=Plot2D(mean_std_cov[m])
                        
                        caption = str("/".join(long_left_over).title() + ' ' + m.lower() + ' maps of ' + ecv_lookup(self.__varname__) + ' for the data set ' + " ".join(dataset_id) + ' (' + self.__time_period__ + ').')
                        
                        try:
                            numfigs=len(mean_std_cov[m])
                        except:
                            numfigs=1
                            
                        fig = plt.figure()
                        
                        (fig,ax,caption) = plot_setup(d=d,m=m,numfigs=numfigs,fig=fig,caption=caption)
                        
                        x.plot(ax=ax, color=self.colormaps, color_type=ctype,
                               title=" ".join([self.__dataset_id__[indx] for 
                                               indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")",
                               vminmax=vminmax, ext_cmap="both")
                        fig.savefig(filename)
                        plt.close(fig.number)
                    
                        ESMValMD("meta",
                                 filename,
                                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                                 caption + " NA-values are shown in grey.",
                                 '#C3S' + m + "".join(short_left_over) + self.__varname__,
                                 self.__infile__,
                                 self.diagname,
                                 self.authors)
                        
                    except Exception as e:
                        exc_type, exc_obj, exc_tb = sys.exc_info()
                        fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                        self.__logger__.error(exc_type, fname, exc_tb.tb_lineno)
                        self.__logger__.error(mean_std_cov[m])
                        self.__logger__.error('Warning: blank figure!')
                        
                        x=Plot2D_blank(mean_std_cov[m])
                        
                        try:
                            numfigs=len(mean_std_cov[m])
                        except:
                            numfigs=1
                            
                        fig = plt.figure()
                        
                        (fig,ax,caption) = plot_setup(d=d,m=m,numfigs=numfigs,fig=fig,caption=caption)
                        
                        x.plot(ax=ax, color=self.colormaps, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
                        fig.savefig(filename)
                        plt.close(fig.number)
                    
                        del mean_std_cov[m]
                        
                        ESMValMD("meta",
                                 filename,
                                 self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                                 str("/".join(long_left_over).title() + ' ' + m.lower() + ' values of ' + ecv_lookup(self.__varname__) + ' for the data set ' + " ".join(dataset_id) + ' (' + self.__time_period__ + '); Data can not be displayed due to cartopy error!'),
                                 '#C3S' + m + "".join(short_left_over) + self.__varname__,
                                 self.__infile__,
                                 self.diagname,
                                 self.authors)
                    
            del mean_std_cov
            
        return list_of_plots

    def __do_trends__(self):
        
        this_function = "trends and stability"
        
        if not self.var3D:
            list_of_plots=self.__trends_procedures_2D__(cube=self.sp_data)
        else: 
            list_of_plots=[]
            for lev in self.levels:
                loc_cube=self.sp_data.extract(iris.Constraint(coord_values={str(self.level_dim):lambda cell: cell==lev}))
                lop=self.__trends_procedures_2D__(loc_cube,level=lev)
                list_of_plots = list_of_plots + lop
        
        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_trends_input",
                                      expfile="_trends.txt",
                                      protofile="empty.txt",
                                      function=this_function)
            
        if found:    
            self.__do_report__(content={"plots":list_of_plots,"freetext":expected_input}, filename=this_function.upper())
        else:
            self.__do_report__(content={"plots":list_of_plots}, filename=this_function.upper())
        
        # update gcos
        self.__gcos_dict__.update({"Accuracy":{"value":None, "unit":None}})
        self.__gcos_dict__.update({"Stability":{"value":None, "unit":None}})
        
        return
    
    def __trends_procedures_2D__(self,cube=None,level=None):
        
        if cube is None:
            cube=self.sp_data
            
        if level is not None:
            basic_filename = self.__basic_filename__ + "_lev" + str(level)
            dataset_id = [self.__dataset_id__[0],"at","level",str(level),str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[0]]
        
        list_of_plots = []
        
        # simple linear trend (slope) and p-values
        _,S,_,P = utils.__temporal_trend__(cube, pthres=1.01)      
        
        signtrends=(P.data<=0.05)*np.sign(S.data)
        ST=S.copy()
        ST.data=signtrends
        
        try:
            # plotting routines
            x=Plot2D(S)
            
            vminmax = np.array([-1,1])*np.max(np.abs(np.nanpercentile(S.data.compressed(),[5,95])))
            
            filename = self.__plot_dir__ + os.sep + basic_filename + "_trend." + self.__output_type__
            list_of_plots.append(filename)
            
            fig = plt.figure()
            (fig,ax,_) = plot_setup(fig=fig)
            x.plot(ax=ax, color=self.colormaps, ext_cmap="both", color_type="Diverging", vminmax=vminmax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' slope values of ' + ecv_lookup(self.__varname__) + ' temporal trends per decade for the data set ' + " ".join(dataset_id) + ' (' + self.__time_period__ + '). NA-values are shown in grey.'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
            
            x=Plot2D(ST)
                        
            filename = self.__plot_dir__ + os.sep + basic_filename + "_signtrend." + self.__output_type__
            list_of_plots.append(filename)
            
            fig = plt.figure()
            (fig,ax,_) = plot_setup(fig=fig)
            x.plot(ax=ax, color=self.colormaps, color_type="Diverging", vminmax=[-1.,1.], title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' significant slope signs (p<=0.05) of ' + ecv_lookup(self.__varname__) + ' temporal trends per decade for the data set ' + " ".join(dataset_id) + ' (' + self.__time_period__ + ')'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
            
            x=Plot2D(P)
    
            filename = self.__plot_dir__ + os.sep + basic_filename + "_pvals." + self.__output_type__
            list_of_plots.append(filename)
            
            fig = plt.figure()
            (fig,ax,_) = plot_setup(fig=fig)
            x.plot(ax=ax, color=self.colormaps, ext_cmap="both", color_type="Sequential", color_reverse=True, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' p-values for slopes of ' + ecv_lookup(self.__varname__) + ' temporal trends per decade for the data set ' + " ".join(dataset_id) + ' (' + self.__time_period__ + '). NA-values are shown in grey.'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
            
        except Exception as e:
            exc_type, exc_obj, exc_tb = sys.exc_info()
            fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
            self.__logger__.error(exc_type, fname, exc_tb.tb_lineno)
            self.__logger__.info("Trend")
            self.__logger__.info('Warning: blank figure!')
            
            x=Plot2D_blank(S)
            
            fig = plt.figure()
            (fig,ax,_) = plot_setup(fig=fig)
            x.plot(ax=ax, color=self.colormaps, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)
            
            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                     str('Trend plots for the data set "' + "_".join(dataset_id) + '" (' + self.__time_period__ + '); Data can not be displayed due to cartopy error!'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)
        
        del P
        del ST
        del S
        

        return list_of_plots

#    def __do_trends__(self):
#
#        this_function = "trends & stability"
#
#        list_of_plots = []
#
#        # simple linear trend (slope) and p-values
#        _,S,_,P = utils.__temporal_trend__(self.sp_data, pthres=1.01)
#
#        try:
#            # plotting routines
#            x=Plot2D(S)
#
#            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_trend." + self.__output_type__
#            list_of_plots.append(filename)
#
#            fig = plt.figure()
#            ax = [plt.subplot(1,1,1)]
#            fig.set_figheight(1.2*fig.get_figheight())
#            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
#            fig.savefig(filename)
#            plt.close(fig.number)
#
#            ESMValMD("meta",
#                     filename,
#                     self.__basetags__ + ['DM_global', 'C3S_trend'],
#                     str("Latitude/Longitude" + ' slope values of ' + self.__varname__ + ' temporal trends per decade for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
#                     '#C3S' + 'temptrend' + self.__varname__,
#                     self.__infile__,
#                     self.diagname,
#                     self.authors)
#
#            x=Plot2D(P)
#
#            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_pvals." + self.__output_type__
#            list_of_plots.append(filename)
#
#            fig = plt.figure()
#            ax = [plt.subplot(1,1,1)]
#            fig.set_figheight(1.2*fig.get_figheight())
#            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
#            fig.savefig(filename)
#            plt.close(fig.number)
#
#            ESMValMD("meta",
#                     filename,
#                     self.__basetags__ + ['DM_global', 'C3S_trend'],
#                     str("Latitude/Longitude" + ' p-values for slopes of ' + self.__varname__ + ' temporal trends per decade for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
#                     '#C3S' + 'temptrend' + self.__varname__,
#                     self.__infile__,
#                     self.diagname,
#                     self.authors)
#
#        except:
#            print('no figure done #004')
#
#        del P
#
#        # linear trend (slope),breakpoints, and actual data after homogenization
#        TempStab = utils.__TS_of_cube__(self.sp_data,
#                                        dates=self.__tim_read__,
#                                        breakpoint_method="CUMSUMADJ",
#                                        max_num_periods=3,
#                                        periods_method="autocorr",
#                                        temporal_resolution=self.__avg_timestep__,
#                                        min_avail_pts=2)
#
#        # plotting routines
#        try:
#            # plotting the slope after break point correction
#            x=Plot2D(TempStab["slope"])
#
#            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_min_slope." + self.__output_type__
#            list_of_plots.append(filename)
#
#            fig = plt.figure()
#            ax = [plt.subplot(1,1,1)]
#            fig.set_figheight(1.2*fig.get_figheight())
#            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
#            fig.savefig(filename)
#            plt.close(fig.number)
#
#            ESMValMD("meta",
#                     filename,
#                     self.__basetags__ + ['DM_global', 'C3S_trend'],
#                     str("Latitude/Longitude" + ' slope values of ' + self.__varname__ + ' temporal trends per decade after breakpoint detection for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
#                     '#C3S' + 'temptrend' + self.__varname__,
#                     self.__infile__,
#                     self.diagname,
#                     self.authors)
#
#            # plotting number of breakpoints
#            x=Plot2D(TempStab["number_breakpts"])
#
#            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_num_bps." + self.__output_type__
#            list_of_plots.append(filename)
#
#            fig = plt.figure()
#            ax = [plt.subplot(1,1,1)]
#            fig.set_figheight(1.2*fig.get_figheight())
#            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
#            fig.savefig(filename)
#            plt.close(fig.number)
#
#            ESMValMD("meta",
#                     filename,
#                     self.__basetags__ + ['DM_global', 'C3S_trend'],
#                     str("Latitude/Longitude" + ' number of breakpoints of ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
#                     '#C3S' + 'temptrend' + self.__varname__,
#                     self.__infile__,
#                     self.diagname,
#                     self.authors)
#        except:
#            print('no figure done #005')
#
#        # plotting the version (2=deseason, 1=detrend, 0=neither, -1=not enough data available, -2=something went wrong)
##        x=Plot2D(TempStab["version"])
##        figS = x.plot(summary_plot=True, title=" ".join([self.__dataset_id__[idx] for idx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
##
##        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_ver_trend." + self.__output_type__
##        list_of_plots.append(filename)
##        figS.savefig(filename)
##        plt.close(figS)
##
##        ESMValMD("meta",
##                 filename,
##                 self.__basetags__ + ['DM_global', 'C3S_trend'],
##                 str("Latitude/Longitude" + ' versions of trend estimation of ' + self.__varname__ + ' temporal trends per decade after breakpoint detection for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
##                 '#C3S' + 'temptrend' + self.__varname__,
##                 self.__infile__,
##                 self.diagname,
##                 self.authors)
#        try:
#            x=Plot2D(TempStab["slope"]-S)
#
#            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_hom_mean_diff." + self.__output_type__
#            list_of_plots.append(filename)
#
#            fig = plt.figure()
#            ax = [plt.subplot(1,1,1)]
#            fig.set_figheight(1.2*fig.get_figheight())
#            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
#            fig.savefig(filename)
#            plt.close(fig.number)
#
#            ESMValMD("meta",
#                     filename,
#                     self.__basetags__ + ['DM_global', 'C3S_trend'],
#                     str("Latitude/Longitude" + ' slope difference after breakpoint reduction for ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
#                     '#C3S' + 'temptrend' + self.__varname__,
#                     self.__infile__,
#                     self.diagname,
#                     self.authors)
#
#            x=Plot2D(cubestats.pearsonr(TempStab["homogenized"], self.sp_data, corr_coords="time"))
#
#            filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_hom_corr." + self.__output_type__
#            list_of_plots.append(filename)
#
#            fig = plt.figure()
#            ax = [plt.subplot(1,1,1)]
#            fig.set_figheight(1.2*fig.get_figheight())
#            x.plot(ax=ax, title=" ".join([self.__dataset_id__[indx] for indx in [0,2,1,3]]) + " (" + self.__time_period__ + ")")
#            fig.savefig(filename)
#            plt.close(fig.number)
#
#            ESMValMD("meta",
#                     filename,
#                     self.__basetags__ + ['DM_global', 'C3S_trend'],
#                     str("Latitude/Longitude" + ' correlation after breakpoint reduction for ' + self.__varname__ + ' for the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')'),
#                     '#C3S' + 'temptrend' + self.__varname__,
#                     self.__infile__,
#                     self.diagname,
#                     self.authors)
#        except:
#            print('no figure done #006')
#
#        del TempStab
#
#        # produce report
#        self.__do_report__(content={"plots":list_of_plots}, filename=this_function.upper())
#
#        # update gcos
#        self.__gcos_dict__.update({"Accuracy":{"value":None, "unit":None}})
#        self.__gcos_dict__.update({"Stability":{"value":None, "unit":None}})
#
#        return
#
#
#    def __do_maturity_matrix__(self):
#
#        this_function = "System maturity matrix"
#
#        expected_input = self.__work_dir__ + os.sep + "smm_input" + os.sep + self.__basic_filename__ + "_smm_expert.csv"
#        if not os.path.isfile(expected_input):
#            try:
#                os.makedirs(os.path.dirname(expected_input))
#            except OSError as exc: # Guard against race condition
#                if exc.errno != errno.EEXIST:
#                    raise
#            shutil.copy2(os.path.dirname(os.path.realpath(__file__)) + "/libs3/predef/empty_smm_expert.csv", expected_input)
#            print("************************************** WARNING **************************************")
#            print(("Expected " + this_function + " input file " + expected_input + " not found!"))
#            print("Created dummy file instead. Please fill in appropriate values and rerun!")
#            print("(This won't fail if you do not do so and produce a white matrix!!!!)")
#            print("************************************** WARNING **************************************")
#            captionerror = True
#        else:
#            print(("Processing " + this_function + " input file: " + expected_input))
#            captionerror = False
#
#        # plotting routines
#        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_" + "".join(this_function.split()) + "." + self.__output_type__
#        fig = do_smm_table(expected_input, os.path.dirname(os.path.realpath(__file__)) + "/libs3/predef/smm_definitions.csv")
#        fig.savefig(filename)
#        plt.close(fig)
#
#        if captionerror:
#            caption = "The expected input file: " + expected_input + " was not found and an empty dummy file created, therefore this plot is blank. Please edit the expected file!"
#        else:
#            caption = str(this_function + ' for the variable ' + self.__varname__ + ' in the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')')
#
#
#        ESMValMD("meta",
#                 filename,
#                 self.__basetags__ + ['C3S_SMM'],
#                 caption,
#                 '#C3S' + 'SMM' + self.__varname__,
#                 self.__infile__,
#                 self.diagname,
#                 self.authors)
#
#        # produce report
#        self.__do_report__(content={"plots":[filename]}, filename="".join(this_function.upper().split()))
#
#        return
#
#
#    def __do_gcos_requirements__(self):
#
#        this_function = "GCOS requirements"
#
#        # plotting routines
#        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + "_" + "".join(this_function.split()) + "." + self.__output_type__
#        fig = do_gcos_table(self.__varname__, self.__gcos_dict__, os.path.dirname(os.path.realpath(__file__)) + "/example_csvs/example_gcos_reference.csv")
#        fig.savefig(filename)
#        plt.close(fig)
#
#        caption = str(this_function + ' for the variable ' + self.__varname__ + ' in the data set "' + "_".join(self.__dataset_id__) + '" (' + self.__time_period__ + ')')
#
#        ESMValMD("meta",
#                 filename,
#                 self.__basetags__ + ['C3S_GCOS'],
#                 caption,
#                 '#C3S' + 'GCOS' + self.__varname__,
#                 self.__infile__,
#                 self.diagname,
#                 self.authors)
#
#        # produce report
#        self.__do_report__(content={"plots":[filename]}, filename="".join(this_function.upper().split()))
#
#        return
#
#
#    def __spatiotemp_subsets__(self,dict_of_regions=None):
#        """
#        produces spatial subset data sets for further calculation
#        """
#
#        if dict_of_regions is None:
#            dict_of_regions = self.__regions__
#            print(dict_of_regions)
#
#        subset_cubes={}
#
#        for R in list(dict_of_regions.keys()):
#            if all([k_dim in self.__dimensions__ for k_dim in list(dict_of_regions[R].keys())]):
#                loc_subset=self.sp_data
#                for dim in list(dict_of_regions[R].keys()):
#                    if dim=='latitude':
#                        r_min,r_max=np.sort(dict_of_regions[R]['latitude'])
#                        loc_subset=loc_subset.extract(iris.Constraint(latitude=lambda point: r_min <= point <= r_max))
#                    if dim=='longitude':
#                        r_min,r_max=np.sort(dict_of_regions[R]['longitude'])
#                        if r_min<0:
#                            if r_max<0:
#                                r_min+=360
#                                r_max+=360
#                            else:
#                                r_min+=360
#                                if r_min>r_max:
#                                    raise ImplementationError("__spatial_subsets__","This method is not yet implemented for regions crossing the 0 degrees longitude.")
#                        loc_subset=loc_subset.extract(iris.Constraint(longitude=lambda point: r_min <= point <= r_max))
#                    if dim=='time':
#                        r_min,r_max=np.sort(dict_of_regions[R]['time'])
#                        loc_subset=loc_subset.extract(iris.Constraint(time=lambda cell: r_min <= cell.point  <= r_max))
#                subset_cubes.update({R:loc_subset})
#            else:
#                print("Region " + R + " specifications not specified correctly: " + str(dict_of_regions[R]) + "!")
#        return subset_cubes
