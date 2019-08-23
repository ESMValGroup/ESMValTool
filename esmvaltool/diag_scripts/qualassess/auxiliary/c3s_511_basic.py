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
import random
import string
import collections
import csv
import matplotlib
#matplotlib.use('Agg')
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import dask
#from scipy import stats
import datetime
from json import load,dump
#from memory_profiler import profile

from .libs import c3s_511_util as utils
from .libs.customErrors import \
    ImplementationError, ConfigurationError, PathError, EmptyContentError
import warnings
from .libs.reporting import do_report as report
from .libs.reporting import do_full_report as full_report
from .plots.matrices import do_smm_table
from .plots.matrices import do_gcos_table
from .plots.matrices import do_eval_table
from .plots.basicplot import \
    Plot2D, PlotHist, Plot2D_blank, Plot1D, PlotScales, plot_setup
from .libs.MD_old.ESMValMD import ESMValMD
import logging
from pprint import pprint
from .libs.predef.ecv_lookup_table import ecv_lookup
from .libs.predef.color_lookup_table import color_lookup

logger = logging.getLogger(os.path.basename(__file__))

# For warning suppression. Comment for development.
warnings.simplefilter("ignore")


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
        # TODO explain all attributes
        self.__cfg__ = None
        self.__logger__ = logging.getLogger(os.path.basename(__file__))

        self.__plot_dir__ = '.' + os.sep  # default plot directory
        self.__work_dir__ = '.' + os.sep  # default work dir
        self.__report_dir__ = self.__work_dir__ + os.sep + "c3s_511" + os.sep

        self.__varname__ = 'var'  # default value
        self.__output_type__ = 'png'  # default ouput file type
        self.__regions__ = {}  # default regions

        self.__basetags__ = []
        self.__infile__ = None

        self.__time_period__ = None
        self.__dataset_id__ = None

        self.authors = ["A_muel_bn", "A_hass_bg", "A_laue_ax",
                        "A_broe_bj", "A_mass_fr", "A_nico_nd",
                        "A_schl_mn", "A_bock_ls"]  # TODO fill in
        self.diagname = "Diagnostic_skeleton"
        self.CDS_ID = "CDS_ID_needed"

        self.colormaps = dict({"default": "binary"})
        self.__latex_output__ = False
        self.levels = [None]
        self.log_level = False
        self.log_data = False

        self.sp_data = None
        self.map_area_frac = None
        self.mp_data = None
        self.var3D = None
        self.__basic_filename__ = "someFile"
        self.__dimensions__ = []
        self.level_dim = None
        self.__time_read__ = []
        self.__avg_timestep__ = None
        
        self.__requested_diags__ = []
        self.__report_directories__ = None
        self.__report_order__ = []
        
        self.__extremes__ = {}
        self.__extremes_regions__ = {}
        
        self.reporting_structure = collections.OrderedDict()

    def set_info(self, **kwargs):
        self.__cfg__ = kwargs.get('cfg', None)
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        raise ImplementationError(
            "set_info",
            "This method has to be implemented.")
        return

    def read_data(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        raise ImplementationError(
            "read_data",
            "This method has to be implemented.")
        return

    def run_diagnostic(self):
        self.__logger__.debug(
            "There is nothing to run! This is the empty skeleton.")
        warnings.warn("Implementation Warning", UserWarning)
        pass

    def __do_overview__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={}, filename="do_overview_default")
        raise ImplementationError(
            "__do_overview__",
            "This method has to be implemented.")
        return

    def __do_mean_var__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={}, filename="do_mean_var_default")
        raise ImplementationError(
            "__do_mean_var__",
            "This method has to be implemented.")
        return

    def __do_trends__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={}, filename="do_trends_default")
        raise ImplementationError(
            "__do_trends__",
            "This method has to be implemented.")
        return

    def __do_extremes__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + "is not implemented!")
        self.__do_report__(content={}, filename="do_extremes_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_sectors__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={}, filename="do_sectors_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_maturity_matrix__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={}, filename="do_maturity_matrix_default")
        raise ImplementationError(
            "__do_maturity_matrix__",
            "This method has to be implemented.")
        return

    def __do_app_perf_matrix__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={}, filename="do_app_perf_matrix_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

    def __do_gcos_requirements__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={}, filename="do_gcos_requirements_default")
        raise ImplementationError(
            "__do_gcos_requirements__",
            "This method has to be implemented.")
        return

    def __do_esm_evaluation__(self):
        this_function = inspect.currentframe().f_code.co_name
        self.__logger__.error(this_function + " is not implemented!")
        self.__do_report__(content={}, filename="do_esmevaluation_default")
        warnings.warn("Implementation Warning", UserWarning)
        return

#    def __do_report__(self,*args,**kwargs):
#        this_function = inspect.currentframe().f_code.co_name
#        self.__logger__.error(this_function + " is not implemented!")
#        raise ImplementationError("__do_report__","This method has to be implemented.")
#        return

    def __do_report__(self, **kwargs):
        """
        reporting function for use with single diagnostics
        """
        content = kwargs.get('content', [])
        if not isinstance(content, (list, dict)):
            raise TypeError("content", "Element is not a list, nor a dict.")

        def rand_str(n): return ''.join(
            [random.choice(string.ascii_lowercase) for i in range(n)])

        filename = kwargs.get('filename', rand_str(10))
        if not isinstance(filename, str):
            raise TypeError("filename", "Element is not a string.")

        report(content,
               filename,
               self.__work_dir__,
               ecv=ecv_lookup(self.__varname__),
               dataset="".join(str(self.__dataset_id__[1])),
               signature=self.CDS_ID,
               latex_opts=self.__latex_output__)
        return

    def __do_full_report__(self):
        """
        reporting function for use with all diagnostics
        """
        
        self.__gather_content__()
        
        for content in self.reporting_structure:
            if not isinstance(self.reporting_structure[content], (list, dict)):
                raise TypeError("contents of reporting_structure",
                                "Elements are not a list, nor a dict.")

        filename = "Full Assessment"
        if not isinstance(filename, str):
            raise TypeError("filename", "Element is not a string.")

        full_report(self.reporting_structure,
                    filename,
                    self.__work_dir__,
                    ecv=ecv_lookup(self.__varname__),
                    dataset="".join(str(self.__dataset_id__[1])),
                    signature=self.CDS_ID,
                    latex_opts=self.__latex_output__)
        return
    
    def __gather_content__(self):
        """
        gather all other diagnostic information for reporting
        """
        self.__logger__.info("Reporting structure and order:")
        self.__logger__.info(self.reporting_structure)
        self.__logger__.info(self.__report_order__)
        for theme in self.__report_order__:
            directory = [s for s in self.__cfg__["input_files"] if theme in s]
            if len(directory) == 1: 
                if os.path.isdir(directory[0]):
                    self.__logger__.info("gathering from directory " + 
                                         directory[0])
                    repstruct = load(open(directory[0] + os.sep + 
                                          'c3s_511' + os.sep + 
                                          'custom_reporting.json'))
#                    repstruct = yaml.load(open(directory[0] + os.sep + 
#                                               'c3s_511' + os.sep + 
#                                               'custom_reporting.yml'),
#                                    Loader=yamlordereddictloader.SafeLoader)
                    utils.dict_merge(self.reporting_structure, repstruct)
            else:
                raise ConfigurationError("self.__report_order__", 
                                         "Requested theme " + theme +
                                         " could not be found!")
        
    def __file_anouncement__(
            self,
            subdir="input",
            expfile="end.txt",
            protofile="none.txt",
            function="no_fun"):

        expected_input = self.__work_dir__ + os.sep + subdir + os.sep + (
            "_".join(
                self.__basic_filename__) if isinstance(
                self.__basic_filename__,
                list) else self.__basic_filename__) + expfile
        if not os.path.isfile(expected_input):
            try:
                os.makedirs(os.path.dirname(expected_input))
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
            shutil.copy2(
                os.path.dirname(
                    os.path.realpath(__file__)) +
                "/libs/predef/" +
                protofile,
                expected_input)
            self.__logger__.info(
                "\n************************************** " +
                "WARNING" +
                " **************************************\n" + 
                "Expected " + function + " input file " + expected_input +
                " not found!\n" + 
                "Created dummy file instead. Please fill in appropriate " + 
                "information and rerun!\n" + 
                "(This won't fail if you do not do " + 
                "so and produce an empty output!!!!)" +
                "\n************************************** " + 
                "WARNING" + 
                " **************************************\n")

            found = False
        else:
            self.__logger__.info(
                "Processing " +
                function +
                " input file: " +
                expected_input)
            found = True

        return expected_input, found

    def __spatiotemp_subsets__(self, cube, dict_of_regions=None):
        """
        produces spatial subset data sets for further calculation
        """

        if dict_of_regions is None:
            dict_of_regions = self.__regions__

        subset_cubes = {}

        for R in dict_of_regions.keys():
            if all([k_dim in self.__dimensions__ for 
                    k_dim in dict_of_regions[R].keys()]):
                loc_subset = cube.copy()
                for dim in dict_of_regions[R].keys():
                    if dim == 'latitude':
                        r_min, r_max = np.sort(dict_of_regions[R]['latitude'])
                        try:  # iris v2
                            loc_subset = loc_subset.extract(iris.Constraint(
                                latitude=lambda cell: r_min <= cell.point <= r_max))
                        except BaseException:  # iris v1
                            loc_subset = loc_subset.extract(iris.Constraint(
                                latitude=lambda point: r_min <= point <= r_max))
                    if dim == 'longitude':
                        r_min, r_max = np.sort(dict_of_regions[R]['longitude'])
                        loc_subset = loc_subset.intersection(
                            longitude=(r_min, r_max))
                        loc_subset = loc_subset.intersection(
                            longitude=(-180, 180))
                    if dim == 'time':
                        r_min, r_max = np.sort(dict_of_regions[R]['time'])
                        try:  # iris v2
                            loc_subset = loc_subset.extract(iris.Constraint(
                                time=lambda cell: r_min <= cell.point <= r_max))
                        except BaseException:  # iris v1
                            tu = loc_subset.coord(
                                "time").units.origin.split(" ")
                            tu_origin = datetime.datetime.strptime(
                                tu[-2], "%Y-%m-%d")  # we assume 00:00:00 as a starttime
                            tu_inc = tu[0]
                            if tu_inc in ["day", "days"]:
                                loc_subset = loc_subset.extract(
                                    iris.Constraint(
                                        time=lambda point: (
                                            r_min -
                                            tu_origin).days <= point <= (
                                            r_max -
                                            tu_origin).days))
                            elif tu_inc in ["second", "seconds"]:
                                loc_subset = loc_subset.extract(
                                    iris.Constraint(
                                        time=lambda point: (
                                            r_min -
                                            tu_origin).total_seconds() <= point <= (
                                            r_max -
                                            tu_origin).total_seconds()))
                            else:
                                raise ValueError("Time coordinate has a " +
                                                 "time increment that was " +
                                                 "not expected: " + tu_inc)
                                
                subset_cubes.update({R: loc_subset})
            else:
                logger.error("Region " + R +
                             " specifications not specified correctly: " +
                             str(dict_of_regions[R]) + "!")
        #BAS: TODO the warning below can be made more helpfull, printing the 
        # key of the empty subsets.
        if any([subset_cubes[sc] is None for sc in subset_cubes]):
            raise ValueError(
                "Could not calculate all subsets. " + 
                "Some are non-overlapping! " +
                str(subset_cubes))

        return subset_cubes
    
    def __apply_fun2cube__(self, cube, dims=None, function=None,
                           incl_weights = True, **kwargs):
        """
        applies function to a sliced cube (memory saving)
        """
        self.__logger__.info("====================================================")
        self.__logger__.info("Running apply_fun2_cube")
#        self.__logger__.info(cube)
        self.__logger__.info(dims)
        self.__logger__.info(kwargs)
        self.__logger__.info("still lazy?")
        self.__logger__.info(cube.has_lazy_data())
        
        try:
#        for i in [1]:
            if "time" not in dims:
    #        if "latitude" in dims:
                latlon_list = []
                
                for latlon in cube.slices(["latitude","longitude"]):
                    
    #                self.__logger__.info(latlon)
                    
                    latlon.remove_coord("day_of_month")
                    latlon.remove_coord("day_of_year")
                    latlon.remove_coord("month_number")
                    latlon.remove_coord("year")
                    
                    if incl_weights and "latitude" in dims:
                        iris.analysis.maths.multiply(latlon, self.map_area_frac, in_place=True)
                        latlon.standard_name = cube.standard_name
                        latlon.long_name = cube.long_name
                        
    #                self.__logger__.info(latlon)
    
                    latlon = latlon.collapsed(dims,
                                              function,
                                              **kwargs)
                    latlon_list.append(latlon)
            
                cube_list = iris.cube.CubeList(latlon_list)
    #            self.__logger__.info(cube_list)
    #            self.__logger__.info([c.coords for c in cube_list])
                
                new_cube = cube_list.merge_cube()
                
                self.__logger__.info("going path one")
        
#        except: 
            else:
                new_cube = cube.collapsed(dims, function, **kwargs)
                self.__logger__.info("going path two")
        except:
            new_cube = cube.collapsed(dims, function, **kwargs)
            self.__logger__.info("going path three")
            
#        self.__logger__.info(new_cube)
            
#        self.__logger__.info(function)
        self.__logger__.info("still lazy?")
        self.__logger__.info(cube.has_lazy_data())
        self.__logger__.info("====================================================")
        
        return new_cube

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
        self.__gcos_dict__.update(
                {"Frequency": {"value": None, "unit": None}})
        self.__gcos_dict__.update(
                {"Resolution": {"value": None, "unit": None}})
        self.__gcos_dict__.update(
                {"Accuracy": {"value": None, "unit": None}})
        self.__gcos_dict__.update(
                {"Stability": {"value": None, "unit": None}})

    def set_info(self, **kwargs):
        """
        gather information for diagnostic
        """

        self.__cfg__ = kwargs.get('cfg', None)
        
        #######################################################################
#        # TODO delete after debugging
##        self.__logger__.info("Object content pre:")
#        self.li3 = [method_name for method_name in dir(
#            self) if callable(getattr(self, method_name))].copy()
##        self.__logger__.info(np.sort(self.li3))
#        self. li1 = self.__dict__.copy().keys()
##        self.__logger__.info(np.sort([*self.li1]))
        #######################################################################

        # TODO use self.__cfg__.pop() here!!!
        if not isinstance(self.__cfg__, dict) or len(self.__cfg__) == 0:
            raise EmptyContentError("__cfg__", "Element is empty.")

        try:
            self.__latex_output__ = self.__cfg__['show_latex']
        except BaseException:
            pass

        try:
            self.levels = self.__cfg__['levels']
        except BaseException:
            pass
        
        try:
            self.log_level = self.__cfg__['log_level']
        except BaseException:
            pass
        
        try:
            self.log_data = self.__cfg__['log_data']
        except BaseException:
            pass
        
        try:
            self.__requested_diags__ = self.__cfg__['requests']
        except BaseException:
            pass
        
        try:
            self.__report_order__ = self.__cfg__['order']
        except BaseException:
            self.__report_order__ = self.__cfg__["input_files"]
            
        try:
            # reads extremes setup from recipe
            self.__extremes__ = dict({"min_measurements": self.__cfg__["minimal_number_measurements"],
                                      "which_percentile": self.__cfg__["which_percentile"],
                                      "window_size": self.__cfg__["window_size"],
                                      "extreme_events": self.__cfg__["extreme_events"],
                                      "multiprocessing": self.__cfg__["multiprocessing"],
                                      })
        except BaseException:
            pass
            
        self.__plot_dir__ = self.__cfg__['plot_dir']
        self.__work_dir__ = self.__cfg__['work_dir']
        self.__report_dir__ = self.__work_dir__ + os.sep + "c3s_511" + os.sep

        self.__infile__ = list(self.__cfg__['input_data'].keys())
        if not len(self.__infile__) == 1:
            raise ConfigurationError(
                "self.__infile__",
                "There should be only one infile!")
        self.__infile__ = self.__infile__[0]

        fileinfo = self.__cfg__['input_data'][self.__infile__]
        self.__varname__ = fileinfo['short_name']
        
        try:
            for c in color_lookup(self.__varname__).keys():
                matplotlib.cm.get_cmap(color_lookup(self.__varname__)[c])
            self.colormaps = color_lookup(self.__varname__)
        except:
            logging.warning("There is no usable specification of data " + 
                            "colors. Falling back to default.")

#        self.__output_type__ = self.__cfg__['output_file_type']  # default ouput file type for the basic diagnostics
#        if not self.__output_type__ == 'png':
#            raise ConfigurationError("self.__output_type__", "Only png is currently supported.")
        self.__regions__ = {
            'Europe_2000': {
                'latitude': (30, 75),
                'longitude': (-10, 35),
                'time': (datetime.datetime(2000, 5, 1),
                         datetime.datetime(2000, 9, 30)
                         )}}  # default region
#
#        # for metadata
#        self.__basetags__ = [self.__varname__] + []# TODO transport tags from namelist
#
        self.__time_period__ = "-".join(
            [str(fileinfo['start_year']), str(fileinfo['end_year'])])

        try:
            # model
            self.__dataset_id__ = [
                fileinfo["cmor_table"],
                fileinfo["dataset"],
                fileinfo["mip"],
                fileinfo["exp"],
                fileinfo["ensemble"],
                fileinfo["short_name"]]
        except BaseException:
            # obs/reanalysis
            self.__dataset_id__ = [
                fileinfo["cmor_table"],
                fileinfo["dataset"],
                fileinfo["mip"],
                fileinfo["type"],
                fileinfo["version"],
                fileinfo["short_name"]]

        self.__basic_filename__ = "_".join(
            self.__dataset_id__ + [self.__time_period__])

        self.__dimensions__ = np.array(["time", "latitude", "longitude"])  
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
                sp_dimensions = [item for item in sp_dimensions if 
                                 item not in aux_dimensions]
                extra_dimensions = [item for item in sp_dimensions if 
                                    item not in self.__dimensions__]
                self.__dimensions__ = sp_dimensions
                if not self.sp_data.coord('latitude').has_bounds():
                    self.sp_data.coord('latitude').guess_bounds()
                if not self.sp_data.coord('longitude').has_bounds():
                    self.sp_data.coord('longitude').guess_bounds()
                if self.sp_data.units == "no-unit":
                    self.sp_data.units = '1'
            else:
                self.__logger__.error(
                    "self.__infile__ is not set to valid path.")
                raise PathError(
                    self.__infile__,
                    "This file is not accessible.")
        except BaseException:
            self.__logger__.error("self.__infile__ is not readable.")
            raise PathError(
                "Basic_Diagnostic_SP.__init__",
                "self.__infile__ is not readable: " +
                self.__infile__)

        if len(extra_dimensions) == 0:
            self.var3D = False
            self.level_dim = None
        elif len(extra_dimensions) == 1:
            self.var3D = True
            self.level_dim = extra_dimensions[0]
            self.__logger__.info("available levels are: " +
                                 str(self.sp_data.coord(self.level_dim)))
        else:
            assert False, "Error in data dimensions. Too many dimensions!"

        # Human readable time
        time_coord = self.sp_data.coord("time")

        time_read = unit.num2date(
            time_coord.points,
            time_coord.units.name,
            time_coord.units.calendar)

        self.__time_read__ = time_read
        
        # provide weights
        
        self.map_area_frac = iris.analysis.cartography.area_weights(
                next(self.sp_data.slices(["latitude","longitude"])))
        self.map_area_frac = dask.array.from_array(self.map_area_frac,chunks="auto")
        
        # save timestep for other diagnostics
        tim_range = self.sp_data.coord("time").points
        tim_freq = np.diff(tim_range)
        tim_freq_spec = utils.__minmeanmax__(tim_freq)
        self.__avg_timestep__ = tim_freq_spec

        return


    def run_diagnostic(self, cfg = None):
#        self.sp_data = self.__spatiotemp_subsets__(self.sp_data)['Europe_2000']
        self.__do_overview__() # granular, requests in overview code
        self.__do_mean_var__() # granular, requests in mean_var code
        if "lintrend" in self.__requested_diags__:
            self.__do_trends__()
        self.__do_extremes__()
        self.__do_sectors__()
        if "SMM" in self.__requested_diags__:
            self.__do_maturity_matrix__()
        if "GCOS" in self.__requested_diags__:
            self.__do_gcos_requirements__()
        if "ESM" in self.__requested_diags__:
            self.__do_esm_evaluation__()
        self.__do_app_perf_matrix__()
        
        if "report" in self.__requested_diags__:
            self.__do_full_report__()
        else:    
            dump(self.reporting_structure,
                      open(self.__report_dir__ + 'custom_reporting.json', 'w'))
#            yaml.dump(self.reporting_structure,
#                      open(self.__report_dir__ + 'custom_reporting.yml', 'w'),
#                      Dumper=yamlordereddictloader.SafeDumper,
#                      default_flow_style=False)

        #######################################################################
#        # TODO delete after debugging
##        self.__logger__.info("Object content post:")
#        li4 = [method_name for method_name in dir(
#            self) if callable(getattr(self, method_name))]
##        self.__logger__.info(np.sort(li4))
#        li2 = self.__dict__.keys()
##        self.__logger__.info(np.sort([*li2]))
#
#        self.__logger__.info("Object content diff:")
## self.__logger__.info(np.sort([x for x in [*li2] if x not in [*self.li1]
## + ['li1']]))
#        self.__logger__.info({y: self.__dict__[y] for y in np.sort([x for x in [*li2] if x not in [*self.li1] + ['li1']])})
#        self.__logger__.info(np.sort([x for x in li4 if x nopasst in self.li3])) 
#        self.__logger__.info("attributes")
#        self.__logger__.info(li2)
#        for l in li2:
#            self.__logger__.info(self.__dict__[l])
#            self.__logger__.info(type(self.__dict__[l]))
#            self.__logger__.info("-----------------------------------------------------------")
        #######################################################################
        pass

    def __do_overview__(self):

        this_function = "overview"
        
        list_of_plots = []

        if not self.var3D:
            lop = self.__overview_procedures_2D__(cube=self.sp_data)
            list_of_plots = list_of_plots + lop
            lop = self.__add_overview_procedures_2D__(cube=self.sp_data)
            list_of_plots = list_of_plots + lop
        else:
            list_of_plots = []
            for lev in self.levels:
                loc_cube = self.sp_data.extract(iris.Constraint(
                    coord_values={str(self.level_dim): lambda cell: cell == lev}))
#                if loc_cube is None:
#                    raise ValueError("Invalid input: at least one of " +
#                                     "the requested levels seem " +
#                                     "not to be available in the analyzed " +
#                                     "cube.")
                lop = self.__overview_procedures_2D__(loc_cube, level=lev)
                list_of_plots = list_of_plots + lop
                lop = self.__add_overview_procedures_2D__(loc_cube, level=lev)
                list_of_plots = list_of_plots + lop
            lop = self.__overview_procedures_3D__()
            list_of_plots = list_of_plots + lop
            lop = self.__add_overview_procedures_3D__()
            list_of_plots = list_of_plots + lop

        overview_dict = None

        if "info" in self.__requested_diags__:
    
            # dimension information
            lon_range = self.sp_data.coord("longitude").points
            lat_range = self.sp_data.coord("latitude").points
            tim_range = self.sp_data.coord("time").points
            t_info = str(self.sp_data.coord("time").units)
    
            lon_range_spec = utils.__minmeanmax__(lon_range)
            lat_range_spec = utils.__minmeanmax__(lat_range)
            tim_range_spec = utils.__minmeanmax__(tim_range)
    
            tim_range_spec_read = unit.num2date(
                tim_range_spec,
                self.sp_data.coord("time").units.name,
                self.sp_data.coord("time").units.calendar)
    
            lon_freq = np.diff(lon_range)
            lat_freq = np.diff(lat_range)
    
            lon_freq_spec = utils.__minmeanmax__(lon_freq)
            lat_freq_spec = utils.__minmeanmax__(lat_freq)
    
            overview_dict = collections.OrderedDict()
    
            overview_dict.update({'longitude range [' + str(self.sp_data.coord("longitude").units) + ']': collections.OrderedDict(
                [("min", str(lon_range_spec[0])), ("max", str(lon_range_spec[2]))])})
            overview_dict.update({'longitude frequency [' + str(self.sp_data.coord("longitude").units) + ']': collections.OrderedDict(
                [("min", str(lon_freq_spec[0])), ("average", str(lon_freq_spec[1])), ("max", str(lon_freq_spec[2]))])})
            overview_dict.update({'latitude range [' + str(self.sp_data.coord("latitude").units) + ']': collections.OrderedDict(
                [("min", str(lat_range_spec[0])), ("max", str(lat_range_spec[2]))])})
            overview_dict.update({'latitude frequency [' + str(self.sp_data.coord("latitude").units) + ']': collections.OrderedDict(
                [("min", str(lat_freq_spec[0])), ("average", str(lat_freq_spec[1])), ("max", str(lat_freq_spec[2]))])})
            overview_dict.update({'temporal range': collections.OrderedDict(
                [("min", str(tim_range_spec_read[0])), ("max", str(tim_range_spec_read[2]))])})
            overview_dict.update({'temporal frequency [' + t_info.split(" ")[0] + ']': collections.OrderedDict(
                [("min", str(self.__avg_timestep__[0])), ("average", str(round(self.__avg_timestep__[1], 2))), ("max", str(self.__avg_timestep__[2]))])})
    
            try:
                lev_range = self.sp_data.coord(self.level_dim).points
                lev_range_spec = utils.__minmeanmax__(lev_range)
                lev_freq_spec = utils.__minmeanmax__(np.diff(lev_range))
                overview_dict.update({'level range [' + str(self.sp_data.coord(self.level_dim).units) + ']': collections.OrderedDict(
                    [("min", str(lev_range_spec[0])), ("max", str(lev_range_spec[2]))])})
                overview_dict.update({'level frequency [' + str(self.sp_data.coord(self.level_dim).units) + ']': collections.OrderedDict(
                    [("min", str(lev_freq_spec[0])), ("average", str(lev_freq_spec[1])), ("max", str(lev_freq_spec[2]))])})
            except BaseException:
                pass

            # save GCOS requirements
            self.__gcos_dict__.update(
                    {"Frequency": {"value": round(self.__avg_timestep__[1], 2),
                                   "unit": t_info.split(" ")[0]}})
            self.__gcos_dict__.update(
                    {"Resolution": {"value": round(np.mean([lon_freq_spec[1],
                                                            lat_freq_spec[1]]), 2),
                                    "unit": str(self.sp_data.coord(
                                            "longitude").units)}})

        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_overview_input",
                                      expfile="_overview.txt",
                                      protofile="empty.txt",
                                      function=this_function)

        if found:
#            self.__do_report__(
#                content={
#                    "listtext": overview_dict,
#                    "plots": list_of_plots,
#                    "freetext": expected_input},
#                filename=this_function.upper())
            if not overview_dict is None:
                self.reporting_structure.update(
                        {"Overview": 
                            {"listtext": overview_dict,
                             "plots": list_of_plots,
                             "freetext": expected_input}})
            else:
                self.reporting_structure.update(
                        {"Overview": 
                            {"plots": list_of_plots,
                             "freetext": expected_input}})
        else:
#            self.__do_report__(
#                content={
#                    "listtext": overview_dict,
#                    "plots": list_of_plots},
#                filename=this_function.upper())
            if not overview_dict is None:
                self.reporting_structure.update(
                        {"Overview": 
                            {"listtext": overview_dict,
                             "plots": list_of_plots}})
            else:
                self.reporting_structure.update(
                        {"Overview": 
                            {"plots": list_of_plots}})
    
        return

    def __overview_procedures_2D__(self, cube=None, level=None):

        if cube is None:
            cube = self.sp_data

        if level is not None:
            basic_filename = self.__basic_filename__ + "_lev" + str(level)
            dataset_id = [self.__dataset_id__[1], "at", "level", str(
                level), str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[1]]

        list_of_plots = []

        reg_dimensions = [item for item in self.__dimensions__ if 
                          item not in set([self.level_dim])]
        maxnumtemp = float(len(cube.coord("time").points))

        if "availability" in self.__requested_diags__:
    
            try:
                sp_masked_vals = cube.collapsed(
                    "time", iris.analysis.COUNT,
                    function=lambda values: values.mask)
                sp_masked_vals.data.mask = sp_masked_vals.data == maxnumtemp
                sp_masked_vals.data = sp_masked_vals.data * 0. + 1.
            except BaseException:
                sp_masked_vals = cube.collapsed(
                    "time", iris.analysis.MEAN) * 0. + 1.  # errors when data is not correctly masked!
    
            for d in reg_dimensions:
    
                long_left_over = [rd for rd in reg_dimensions if rd != d]
                short_left_over = np.array([sl[0:3] for sl in long_left_over])
    
                num_available_vals = cube.collapsed(d, iris.analysis.COUNT,
                                                    function=lambda values: values >= cube.collapsed(reg_dimensions, iris.analysis.MIN).data)
                if d in ["time"]:
                    frac_available_vals = iris.analysis.maths.divide(
                        num_available_vals, maxnumtemp)
                else:
                    sp_agg_array = sp_masked_vals.collapsed(d, iris.analysis.SUM)
                    frac_available_vals = iris.analysis.maths.divide(
                        num_available_vals, sp_agg_array.data)
    
                try:
                    # plotting routine
                    filename = self.__plot_dir__ + os.sep + basic_filename + \
                        "_frac_avail_" + "_".join(short_left_over) + "." + \
                        self.__output_type__
                    list_of_plots.append(filename)
                    frac_available_vals.rename("availability as a fraction")
                    x = Plot2D(frac_available_vals)
    
                    fig = plt.figure()
                    (fig, ax, _) = plot_setup(d=d, fig=fig)
                    x.plot(ax=ax,
                           vminmax=[0., 1.],
                           color={"Sequential": "YlGn"},
                           color_type="Sequential",
                           title=" ".join([self.__dataset_id__[idx] for 
                                           idx in [0, 2, 1, 3]]) + " (" + 
                        self.__time_period__ + ")")
                    fig.savefig(filename)
                    plt.close(fig.number)
    
                    ESMValMD("meta",
                             filename,
                             self.__basetags__ + ['DM_global', 'C3S_overview'],
                             str('Overview on ' + "/".join(long_left_over) +
                                 ' availablility of ' + 
                                 ecv_lookup(self.__varname__) +
                                 ' for the data set ' + " ".join(dataset_id) +
                                 ' (' + self.__time_period__ +
                                 '). ' + 
                                 ('NA-values and values of 0 and below are shown ' + 
                                  'in grey.' if self.log_data else
                                  'NA-values are shown in grey.')),
                             '#C3S' + 'frav' + "".join(short_left_over) +
                             self.__varname__,
                             self.__infile__,
                             self.diagname,
                             self.authors)
    
                except Exception as e:
                    exc_type, exc_obj, exc_tb = sys.exc_info()
                    fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
                    self.__logger__.error(exc_type, fname, exc_tb.tb_lineno)
                    self.__logger__.error("Availability in " + d)
                    self.__logger__.error('Warning: blank figure!')
    
                    x = Plot2D_blank(frac_available_vals)
    
                    fig = plt.figure()
    
                    (fig, ax, _) = plot_setup(d=d, fig=fig)
                    x.plot(ax=ax,
                           color=self.colormaps,
                           title=" ".join([self.__dataset_id__[indx] for 
                                           indx in [0, 2, 1, 3]]) + " (" + 
                                 self.__time_period__ + ")")
                    fig.savefig(filename)
                    plt.close(fig.number)
    
                    ESMValMD("meta",
                             filename,
                             self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                             str('Availability plot for the data set "' +
                                 "_".join(dataset_id) + '" (' +
                                 self.__time_period__ +'); Data can not be ' + 
                                 'displayed due to cartopy error!'),
                             '#C3S' + 'frav' + "".join(short_left_over) +
                             self.__varname__,
                             self.__infile__,
                             self.diagname,
                             self.authors)
    
            del sp_masked_vals
            del num_available_vals
            del frac_available_vals

        if "histogram" in self.__requested_diags__:
            # histogram plot of available measurements
    
            try:
                # plotting routine
                filename = self.__plot_dir__ + os.sep + basic_filename + \
                    "_hist_all_vals" + "." + self.__output_type__
                list_of_plots.append(filename)

                x = PlotHist(cube)
                fig = x.plot(dat_log=self.log_data)
                fig.savefig(filename)
                plt.close(fig)
    
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + ['DM_global', 'C3S_overview'],
                         str('Full spatio-temporal histogram of ' +
                             self.__varname__ + ' for the data set ' +
                             " ".join(dataset_id) + ' (' + self.__time_period__ +
                             ').'),
                         '#C3S' + 'histall' + self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors)
                         
            except BaseException:
                self.__logger__.error('Something is probably wrong with the ' + 
                                      'plotting routines/modules!')
    
        return list_of_plots
    
    def __add_overview_procedures_2D__(self, cube=None, level=None):
        return []
    
    def __overview_procedures_3D__(self, cube=None):
        return []
    
    def __add_overview_procedures_3D__(self, cube=None):
        return []

    def __do_mean_var__(self):
        
        this_function = "mean and variability"

        list_of_plots = []
        
        if not self.var3D:
            lop = self.__mean_var_procedures_2D__(cube=self.sp_data)
            list_of_plots = list_of_plots + lop
            if "custom" in self.__requested_diags__:
                lop = self.__add_mean_var_procedures_2D__(cube=self.sp_data)
                list_of_plots = list_of_plots + lop
        else:
            list_of_plots = []
            for lev in self.levels:
                loc_cube = self.sp_data.extract(
                        iris.Constraint(coord_values={str(self.level_dim): lambda cell: cell == lev}))
                lop = self.__mean_var_procedures_2D__(loc_cube, level=lev)
                list_of_plots = list_of_plots + lop
                if "custom" in self.__requested_diags__:
                    lop = self.__add_mean_var_procedures_2D__(loc_cube, level=lev)
                    list_of_plots = list_of_plots + lop
            lop = self.__mean_var_procedures_3D__()
            list_of_plots = list_of_plots + lop
            if "custom" in self.__requested_diags__:
                lop = self.__add_mean_var_procedures_3D__()
                list_of_plots = list_of_plots + lop

        if len(list_of_plots)>0:
            # produce report
            expected_input, found = \
                self.__file_anouncement__(subdir="c3s_511/single_mean_var_input",
                                          expfile="_mean_var.txt",
                                          protofile="empty.txt",
                                          function=this_function)
    
            if found:
                self.reporting_structure.update(
                        {"Mean and Variability": 
                            {"plots": list_of_plots,
                             "freetext": expected_input}})
            else:
                self.reporting_structure.update(
                        {"Mean and Variability": 
                            {"plots": list_of_plots}})
        
        return

#    @profile
    def __mean_var_procedures_2D__(self, cube=None, level=None):

        try:
            if cube is None:
                cube = self.sp_data
    
            if level is not None:
                basic_filename = self.__basic_filename__ + "_ledataset_idv" + str(level)
                dataset_id = [self.__dataset_id__[1], "at", "level", str(
                    level), str(self.sp_data.coord(self.level_dim).units)]
            else:
                basic_filename = self.__basic_filename__
                dataset_id = [self.__dataset_id__[1]]
    
            reg_dimensions = [item for item in self.__dimensions__ if 
                              item not in set([self.level_dim])]
        except:
            return []
        
        data_info = []
        list_of_plots = []


        if "mean" in self.__requested_diags__:
            
            for d in reg_dimensions:
                
                long_left_over = [rd for rd in reg_dimensions if rd != d]
                
                plotcube = utils.dask_weighted_mean_wrapper(cube,self.map_area_frac,dims=d)

                vminmax = utils.minmax_cubelist([plotcube], [5, 95]) 
                
                data_info.append(dict({"name":"mean",
                                       "data": plotcube,
                                       "dim":d,
                                       "llo":long_left_over,
                                       "ctype":"Data",
                                       "vminmax":vminmax,
                                       "mmtype":"abs",
                                       }))
    
            del plotcube
            
        if "std_dev" in self.__requested_diags__:
            
            for d in reg_dimensions:
                
                long_left_over = [rd for rd in reg_dimensions if rd != d]

                plotcube = utils.dask_weighted_stddev_wrapper(cube,self.map_area_frac,dims=d)
                
                vminmax = utils.minmax_cubelist([plotcube], [5, 95]) 
                
                data_info.append(dict({"name":"std_dev",
                                       "data": plotcube,
                                       "dim":d,
                                       "llo":long_left_over,
                                       "ctype":"Sequential",
                                       "vminmax":vminmax,
                                       "mmtype":None,
                                       }))
    
            del plotcube
            
        if "percentiles" in self.__requested_diags__:
            
            percentiles = [1., 5., 10., 50., 90., 95., 99.]
            
            for d in ["time"]:
                
                long_left_over = [rd for rd in reg_dimensions if rd != d]
            
                perc_comp_dict = utils.lazy_percentiles(cube,percentiles)
                
                plotcubes = list(perc_comp_dict.values())
                
                vminmax = utils.minmax_cubelist(plotcubes, [5, 95]) 
                
                data_info.append(dict({"name":"percentiles",
                                       "data": plotcubes,
                                       "dim":d,
                                       "llo":long_left_over,
                                       "ctype":"Data",
                                       "vminmax":vminmax,
                                       "mmtype":"abs",
                                       }))
    
            del plotcubes
            
        if "climatology" in self.__requested_diags__:
            
            for d in ["time"]:
                
                long_left_over = [rd for rd in reg_dimensions if rd != d]
            
                clim_comp_dict = utils.lazy_climatology(cube,'month_number')
                
                plotcubes = list(clim_comp_dict.values())
                
                vminmax = utils.minmax_cubelist(plotcubes, [5, 95])    
                
                data_info.append(dict({"name":"climatology",
                                       "data": plotcubes,
                                       "dim":d,
                                       "llo":long_left_over,
                                       "ctype":"Data",
                                       "vminmax":vminmax,
                                       "mmtype":"abs",
                                       }))
    
            del plotcubes
            
        if "anomalies" in self.__requested_diags__:
            
            for d in ["time"]:
                
                long_left_over = [rd for rd in reg_dimensions if rd != d]
            
                anom_comp_dict = utils.lazy_climatology(cube,'year')
                
                mean = cube.collapsed(d,iris.analysis.MEAN)
                
                plotcubes =  [pc - mean  for pc in list(anom_comp_dict.values())]
                
                vminmax = utils.minmax_cubelist(plotcubes, [1/3 *100, 2/3 *100], symmetric=True) 
                
                data_info.append(dict({"name":"anomalies",
                                       "data": plotcubes,
                                       "dim":d,
                                       "llo":long_left_over,
                                       "ctype":"Diverging",
                                       "vminmax":vminmax,
                                       "mmtype":"diff",
                                       }))
    
            del plotcubes
            
        if "time_series" in self.__requested_diags__:
            
            for d in ["time"]:
                
                long_left_over = [rd for rd in reg_dimensions if rd != d]
            
#                plotcubes_m = utils.dask_weighted_mean_wrapper(
#                        cube,self.map_area_frac, dims=long_left_over)
#                plotcubes_std = utils.dask_weighted_stddev_wrapper(
#                        cube,self.map_area_frac, dims=long_left_over)

#                plotcubes_std_p = plotcubes_m + plotcubes_std
#                plotcubes_std_m = plotcubes_m - plotcubes_std
                percentiles = [1,5,25,50,75,95,99]
                
                pltcb = utils.lazy_percentiles(cube,percentiles,dims = long_left_over)
                
                pltcb = list(pltcb.values())
    
                filename = self.__plot_dir__ + os.sep + basic_filename + \
                    "_" + "temp_series_1d" + "." + self.__output_type__
                list_of_plots.append(filename)
    
                x = Plot1D(iris.cube.CubeList(pltcb))
    
                fig = plt.figure()
    
                ax = [plt.subplot(1, 1, 1)]
                fig.set_figheight(1.2 * fig.get_figheight())
                x.plot(ax=ax,
                       title=["{}%".format(p) for p in percentiles],
                       colors=["skyblue",
                               "blue",
                               "darkorchid",
                               "red",
                               "darkorchid",
                               "blue",
                               "skyblue",
                               ])
                fig.savefig(filename, bbox_inches='tight')
                plt.close(fig.number)
    
                ESMValMD("meta",
                         filename,
                         self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                         str("/".join(long_left_over).title() + 
                             ' aggregated ' + d + ' series' + ' values of ' +
                             ecv_lookup(self.__varname__) + 
                             ' for the data set ' + " ".join(dataset_id) +
                             ' (' + self.__time_period__ + ').'),
                         '#C3S' + d + 'series' + "".join(
                             np.array([sl[0:3] for sl in long_left_over])) +
                         self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors)
            
        # Check mins and maxs
            
        abs_mins = [dask.array.atleast_1d(di["vminmax"].min()) for di in
                    data_info if di["mmtype"]=="abs"]
        abs_maxs = [dask.array.atleast_1d(di["vminmax"].max()) for di in
                    data_info if di["mmtype"]=="abs"]
        [di.update({"vminmax":[dask.array.concatenate(abs_mins).min().compute(),
                               dask.array.concatenate(abs_maxs).max().compute()]}) 
         for di in data_info if di["mmtype"]=="abs"]
        
        diff_mins = [dask.array.atleast_1d(di["vminmax"].min()) for di in
                     data_info if di["mmtype"]=="diff"]
        diff_maxs = [dask.array.atleast_1d(di["vminmax"].max()) for di in
                     data_info if di["mmtype"]=="diff"]
        [di.update({"vminmax":[dask.array.concatenate(diff_mins).min().compute(),
                               dask.array.concatenate(diff_maxs).max().compute()]}) 
         for di in data_info if di["mmtype"]=="diff"]
        
        # plotting (2D maps only)
        for di in data_info:
            
            filename = self.__plot_dir__ + os.sep + basic_filename + \
                    "_" + "_" + di["name"] + "_" + \
                    "_".join(np.array([sl[0:3] for sl in di["llo"]])) + \
                    "." + self.__output_type__
            list_of_plots.append(filename)

            try:
#            for i in [0]:
                x = Plot2D(di["data"])

                caption = str("/".join(di["llo"]).title() +
                              ' ' +
                              di["name"] +
                              ' maps of ' +
                              ecv_lookup(self.__varname__) +
                              ' for the data set ' +
                              " ".join(dataset_id) + ' (' +
                              self.__time_period__ + ').')

                try:
                    numfigs = len(di["data"])
                except BaseException:
                    numfigs = 1

                fig = plt.figure()

                (fig, ax, caption) = plot_setup(d=di["dim"], m=di["name"],
                                                numfigs=numfigs,
                                                fig=fig,
                                                caption=caption)

                x.plot(ax=ax,
                       color=self.colormaps,
                       color_type=di["ctype"],
                       title=" ".join([self.__dataset_id__[indx] for
                                       indx in [0, 2, 1, 3]]) + \
                             " (" + self.__time_period__ + ")",
                       vminmax=di["vminmax"],
                       ext_cmap="both",
                       dat_log = self.log_data)
                fig.savefig(filename)
                plt.close(fig.number)

                ESMValMD("meta",
                         filename,
                         self.__basetags__ + 
                         ['DM_global', 'C3S_mean_var'],
                         caption +  
                         ('NA-values and values of 0 and below are shown ' + 
                          'in grey.' if self.log_data else
                          'NA-values are shown in grey.'),
                         '#C3S' + "mean_var" + "".join(
                             np.array([sl[0:3] for sl in di["llo"]])) + \
                         self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors)

            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(
                    exc_tb.tb_frame.f_code.co_filename)[1]
                self.__logger__.error(
                    exc_type, fname, exc_tb.tb_lineno)
                self.__logger__.error(di["name"])
                self.__logger__.error('Warning: blank figure!')

                x = Plot2D_blank(di["data"])

                fig = plt.figure()

                (fig, ax, caption) = plot_setup(d=d, m=di["name"],
                                                numfigs=numfigs,
                                                fig=fig,
                                                caption=caption)

                x.plot(ax=ax,
                       color=self.colormaps,
                       title=" ".join([self.__dataset_id__[indx] for 
                                       indx in [0, 2, 1, 3]]) + \
                             " (" + self.__time_period__ + ")")
                fig.savefig(filename)
                plt.close(fig.number)

                ESMValMD("meta",
                         filename,
                         self.__basetags__ +
                         ['DM_global', 'C3S_mean_var'],
                         str("/".join(long_left_over).title() +
                             ' ' + di["name"] + ' values of ' +
                             ecv_lookup(self.__varname__) + 
                             ' for the data set ' + 
                             " ".join(dataset_id) +
                             ' (' + self.__time_period__ + 
                             '); Data can ' +
                             'not be displayed due to cartopy error!'),
                         '#C3S' + "mean_var" + "".join(
                             np.array([sl[0:3] for sl in di["llo"]])) +
                         self.__varname__,
                         self.__infile__, 
                         self.diagname, 
                         self.authors)

        return list_of_plots
    
    def __add_mean_var_procedures_2D__(self, cube=None, level=None):
        return []
    
#    @profile
    def __mean_var_procedures_3D__(self, cube=None):
        
        try:
            if cube is None:
                cube=self.sp_data
                
            basic_filename = self.__basic_filename__
            dataset_id = self.__dataset_id__[1]
              
            data_info = []
            list_of_plots = []
            
            reg_dimensions = [item for item in self.__dimensions__ if 
                              item not in set([self.level_dim])]
            
        except:
            return []

        if "mean" in self.__requested_diags__:
            
            for d in reg_dimensions:
                
                long_agg = [rd for rd in reg_dimensions if rd != d]
                long_left_over = [d,self.level_dim]
            
                plotcube = utils.dask_weighted_mean_wrapper(cube,self.map_area_frac,dims=long_agg)
                
                vminmax = utils.minmax_cubelist([plotcube], [5, 95])
                
                data_info.append(dict({"name":"mean",
                                       "data": plotcube,
                                       "dim":d,
                                       "llo":long_left_over,
                                       "ctype":"Data",
                                       "vminmax":vminmax,
                                       "mmtype":"abs",
                                       }))
        
            del plotcube
    
        if "std_dev" in self.__requested_diags__:
            
            for d in reg_dimensions:
                
                long_agg = [rd for rd in reg_dimensions if rd != d]
                long_left_over = [d,self.level_dim]
            
                plotcube = utils.dask_weighted_stddev_wrapper(cube,self.map_area_frac,dims=long_agg)
                
                vminmax = utils.minmax_cubelist([plotcube], [5, 95])
                
                data_info.append(dict({"name":"std_dev",
                                       "data": plotcube,
                                       "dim":d,
                                       "llo":long_left_over,
                                       "ctype":"Sequential",
                                       "vminmax":vminmax,
                                       "mmtype":None,
                                       }))
    
            del plotcube

        for di in data_info:
        
            filename = self.__plot_dir__ + os.sep + basic_filename + \
                    "_" + "_" + di["name"] + "_" + \
                    "_".join(np.array([sl[0:3] for sl in di["llo"]])) + \
                    "." + self.__output_type__
            
            list_of_plots.append(filename)

            try:
                x = Plot2D(di["data"])

                caption = str("/".join(di["llo"]).title() +
                              ' ' +
                              di["name"] +
                              ' maps of ' +
                              ecv_lookup(self.__varname__) +
                              ' for the data set ' +
                              dataset_id + ' (' +
                              self.__time_period__ + ').')

                try:
                    numfigs = len(di["data"])
                except BaseException:
                    numfigs = 1

                fig = plt.figure()

                (fig, ax, caption) = plot_setup(d="levels", m=di["name"],
                                                numfigs=numfigs,
                                                fig=fig,
                                                caption=caption)

                x.plot(ax=ax,
                       color=self.colormaps,
                       color_type=di["ctype"],
                       title=" ".join([self.__dataset_id__[indx] for
                                       indx in [0, 2, 1, 3]]) + \
                             " (" + self.__time_period__ + ")",
                       vminmax=di["vminmax"],
                       y_log=self.log_level,
                       dat_log = self.log_data,
                       ext_cmap="both")
                fig.savefig(filename)
                plt.close(fig.number)

                ESMValMD("meta",
                         filename,
                         self.__basetags__ + 
                         ['DM_global', 'C3S_mean_var'],
                         caption + 
                         ('NA-values and values of 0 and below are shown ' + 
                          'in grey.' if self.log_data else
                          'NA-values are shown in grey.'),
                         '#C3S' + "mymean" + "".join(
                             np.array([sl[0:3] for sl in di["llo"]])) + \
                         self.__varname__,
                         self.__infile__,
                         self.diagname,
                         self.authors)

            except Exception as e:
                exc_type, exc_obj, exc_tb = sys.exc_info()
                fname = os.path.split(
                    exc_tb.tb_frame.f_code.co_filename)[1]
                self.__logger__.error(
                    exc_type, fname, exc_tb.tb_lineno)
                self.__logger__.error(di["name"])
                self.__logger__.error('Warning: blank figure!')

                x = Plot2D_blank(di["data"])

                fig = plt.figure()

                (fig, ax, caption) = plot_setup(d=d, m=di["name"],
                                                numfigs=numfigs,
                                                fig=fig,
                                                caption=caption)

                x.plot(ax=ax,
                       color=self.colormaps,
                       title=" ".join([self.__dataset_id__[indx] for 
                                       indx in [0, 2, 1, 3]]) + \
                             " (" + self.__time_period__ + ")")
                fig.savefig(filename)
                plt.close(fig.number)

                ESMValMD("meta",
                         filename,
                         self.__basetags__ +
                         ['DM_global', 'C3S_mean_var'],
                         str("/".join(long_left_over).title() +
                             ' ' + di["name"] + ' values of ' +
                             ecv_lookup(self.__varname__) + 
                             ' for the data set ' + 
                             " ".join(dataset_id) +
                             ' (' + self.__time_period__ + 
                             '); Data can ' +
                             'not be displayed due to cartopy error!'),
                         '#C3S' + "mymean" + "".join(
                             np.array([sl[0:3] for sl in di["llo"]])) +
                         self.__varname__,
                         self.__infile__, 
                         self.diagname, 
                         self.authors)
                    
        return list_of_plots
    
    def __add_mean_var_procedures_3D__(self, cube=None):
        return []

    def __do_trends__(self):

        this_function = "trends and stability"
        
        list_of_plots = []

        if not self.var3D:
            lop = self.__trends_procedures_2D__(cube=self.sp_data)
            list_of_plots = list_of_plots + lop
            lop = self.__add_trend_procedures_2D__(cube=self.sp_data)
            list_of_plots = list_of_plots + lop
        else:
            list_of_plots = []
            for lev in self.levels:
                loc_cube = self.sp_data.extract(iris.Constraint(
                    coord_values={str(self.level_dim): lambda cell: cell == lev}))
                lop = self.__trends_procedures_2D__(loc_cube, level=lev)
                list_of_plots = list_of_plots + lop
            lop = self.__trend_procedures_3D__()
            list_of_plots = list_of_plots + lop
            lop = self.__add_trend_procedures_3D__()
            list_of_plots = list_of_plots + lop

        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_trends_input",
                                      expfile="_trends.txt",
                                      protofile="empty.txt",
                                      function=this_function)

        if found:
#            self.__do_report__(
#                content={
#                    "plots": list_of_plots,
#                    "freetext": expected_input},
#                filename=this_function.upper())
            self.reporting_structure.update(
                    {"Trends": 
                        {"plots": list_of_plots,
                         "freetext": expected_input}})
        else:
#            self.__do_report__(
#                content={
#                    "plots": list_of_plots},
#                filename=this_function.upper())
            self.reporting_structure.update(
                    {"Trends": 
                        {"plots": list_of_plots}})

        # update gcos
        self.__gcos_dict__.update({"Accuracy": {"value": None, "unit": None}})
        self.__gcos_dict__.update({"Stability": {"value": None, "unit": None}})

        return

    def __trends_procedures_2D__(self, cube=None, level=None):

        import resource
        
        before0 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        
        if cube is None:
            cube = self.sp_data

        if level is not None:
            basic_filename = self.__basic_filename__ + "_lev" + str(level)
            dataset_id = [self.__dataset_id__[1], "at", "level", str(
                level), str(self.sp_data.coord(self.level_dim).units)]
        else:
            basic_filename = self.__basic_filename__
            dataset_id = [self.__dataset_id__[1]]

        list_of_plots = []

        # simple linear trend (slope) and p-values
        _, S, _, P = utils.temporal_trend(cube, pthres=1.01)

        signtrends = (P.data <= 0.05) * np.sign(S.data)
        ST = S.copy()
        ST.data = signtrends
        
        after1 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("full trend without plotting" + ":")
        logger.info(str(round((after1-before0)/1024.,2)) + "MB")
        
        try:
            # plotting routines
            x = Plot2D(S)

            vminmax = np.array([-1, 1]) * np.max(np.abs(
                    np.nanpercentile(S.data.compressed(), [5, 95])))

            filename = self.__plot_dir__ + os.sep + \
                basic_filename + "_trend." + self.__output_type__
            list_of_plots.append(filename)

            fig = plt.figure()
            (fig, ax, _) = plot_setup(fig=fig)
            x.plot(ax=ax,
                   color=self.colormaps,
                   ext_cmap="both",
                   color_type="Diverging",
                   vminmax=vminmax,
                   title=" ".join([self.__dataset_id__[indx] for 
                                   indx in [0, 2, 1, 3]]) + " (" + 
                    self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)

            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' slope values of ' +
                         ecv_lookup(self.__varname__) +
                         ' temporal trends per decade for the data set ' +
                         " ".join(dataset_id) + ' (' + self.__time_period__ +
                         '). NA-values are shown in grey.'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)

            x = Plot2D(ST)

            filename = self.__plot_dir__ + os.sep + \
                basic_filename + "_signtrend." + self.__output_type__
            list_of_plots.append(filename)

            fig = plt.figure()
            (fig, ax, _) = plot_setup(fig=fig)
            x.plot(ax=ax,
                   color=self.colormaps,
                   color_type="Diverging",
                   vminmax=[-1.,1.],
                   title=" ".join([self.__dataset_id__[indx] for 
                                   indx in [0, 2, 1, 3]]) + " (" + 
                     self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)

            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + 
                         ' significant slope signs (p<=0.05) of ' +
                         ecv_lookup(self.__varname__) +
                         ' temporal trends per decade for the data set ' +
                         " ".join(dataset_id) + ' (' + self.__time_period__ +
                         ')'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)

            x = Plot2D(P)

            filename = self.__plot_dir__ + os.sep + \
                basic_filename + "_pvals." + self.__output_type__
            list_of_plots.append(filename)

            fig = plt.figure()
            (fig, ax, _) = plot_setup(fig=fig)
            x.plot(ax=ax,
                   color=self.colormaps,
                   ext_cmap="both",
                   color_type="Sequential",
                   color_reverse=True,
                   title=" ".join([self.__dataset_id__[indx] for 
                                   indx in [0, 2, 1, 3]]) + " (" + 
                     self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)

            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_trend'],
                     str("Latitude/Longitude" + ' p-values for slopes of ' +
                         ecv_lookup(self.__varname__) +
                         ' temporal trends per decade for the data set ' +
                         " ".join(dataset_id) + ' (' + self.__time_period__ +
                         '). NA-values are shown in grey.'),
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

            x = Plot2D_blank(S)

            fig = plt.figure()
            (fig, ax, _) = plot_setup(fig=fig)
            x.plot(ax=ax,
                   color=self.colormaps,
                   title=" ".join([self.__dataset_id__[indx] for 
                                   indx in [0, 2, 1, 3]]) + " (" + 
                    self.__time_period__ + ")")
            fig.savefig(filename)
            plt.close(fig.number)

            ESMValMD("meta",
                     filename,
                     self.__basetags__ + ['DM_global', 'C3S_mean_var'],
                     str('Trend plots for the data set "' +
                         "_".join(dataset_id) + '" (' + self.__time_period__ +
                         '); Data can not be displayed due to cartopy error!'),
                     '#C3S' + 'temptrend' + self.__varname__,
                     self.__infile__,
                     self.diagname,
                     self.authors)

        del P
        del ST
        del S
        
        after2 = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
        logger.info("full trend only plotting" + ":")
        logger.info(str(round((after2-after1)/1024.,2)) + "MB")

        return list_of_plots
    
    def __add_trend_procedures_2D__(self, cube=None, level=None):
        return []
    
    def __trend_procedures_3D__(self, cube=None):
        return []    
    
    def __add_trend_procedures_3D__(self, cube=None):
        return []

    def __do_maturity_matrix__(self):

        this_function = "System maturity matrix"

        self.__file_anouncement__(
            subdir="c3s_511/single_smm_input",
            expfile="_SMM_CORE_CLIMAX_c3s511_Adapted_v5_0.xlsx",
            protofile="SMM_CORE_CLIMAX_c3s511_Adapted_v5_0.xlsx",
            function=this_function)

        self.__file_anouncement__(
            subdir="c3s_511/single_smm_input",
            expfile="_SMM_Guide_for_USERS_C3S_511_v1.pdf",
            protofile="SMM_Guide_for_USERS_C3S_511_v1.pdf",
            function=this_function)

        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_smm_input",
                                      expfile="_smm_expert.csv",
                                      protofile="empty_smm_expert.csv",
                                      function=this_function)

        if not found:
            suggestion = "Please make use of " + \
                "SMM_Guide_for_USERS_C3S_511_v1.pdf and " + \
                "SMM_CORE_CLIMAX_c3s511_Adapted_v5_0.xlsx for producing " + \
                "the requested file!"
            caption = "The expected input file: " + expected_input + \
                " was not found and an empty dummy file created, " + \
                "therefore this plot is blank. Please edit the expected file!"
            self.__logger__.info(suggestion)
        else:
            with open(expected_input, "r") as file_read:
                if not any([lit.isdigit() for lit in list(file_read.read())]):
                    caption = "The expected input file: " + expected_input + \
                        " was found but empty, therefore this plot " + \
                        "is blank. Please edit the requested file!"
                else:
                    caption = str(
                        this_function +
                        ' for the variable ' +
                        ecv_lookup(
                            self.__varname__) +
                        ' in the data set "' +
                        self.__dataset_id__[1] +
                        '" (' +
                        self.__time_period__ +
                        ')')

        # plotting routines
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + \
            "_" + "".join(this_function.split()) + "." + self.__output_type__
        fig = do_smm_table(
            expected_input,
            os.path.dirname(
                os.path.realpath(__file__)) +
            "/libs/predef/smm_definitions.csv")
        fig.savefig(filename)
        plt.close(fig)

        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['C3S_SMM'],
                 caption,
                 '#C3S' + 'SMM' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)

        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_smm_input",
                                      expfile="_smm_add.csv",
                                      protofile="empty.txt",
                                      function=this_function)

        if found:
#            self.__do_report__(
#                content={
#                    "plots": [filename],
#                    "freetext": expected_input},
#                filename=this_function.upper())
            self.reporting_structure.update(
                    {"System Maturity Matrix": 
                        {"plots": [filename],
                         "freetext": expected_input}})
        else:
#            self.__do_report__(
#                content={
#                    "plots": [filename]},
#                filename=this_function.upper())
            self.reporting_structure.update(
                    {"System Maturity Matrix": 
                        {"plots": [filename]}})

        return

    def __do_gcos_requirements__(self):

        this_function = "GCOS requirements"

        # plotting routines
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + \
            "_" + "".join(this_function.split()) + "." + self.__output_type__
        subdir = "c3s_511/single_gcos_input"
        if not os.path.isdir(self.__work_dir__ + os.sep + subdir):
            os.mkdir(self.__work_dir__ + os.sep + subdir)
        fig = do_gcos_table(
            self.__varname__,
            self.__gcos_dict__,
            os.path.dirname(
                os.path.realpath(__file__)) +
            "/libs/predef/GCOS_specific_info.csv",
            self.__work_dir__ +
            os.sep +
            subdir +
            os.sep +
            self.__basic_filename__ +
            "_gcos_values_editable.csv")
        fig.savefig(filename)
        plt.close(fig)

        caption = str(
            this_function +
            ' for the variable ' +
            ecv_lookup(
                self.__varname__) +
            ' in the data set "' +
            self.__dataset_id__[1] +
            '" (' +
            self.__time_period__ +
            ').')

        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['C3S_GCOS'],
                 caption,
                 '#C3S' + 'GCOS' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)

        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir=subdir,
                                      expfile="_gcos.txt",
                                      protofile="empty.txt",
                                      function=this_function)

        if found:
#            self.__do_report__(
#                content={
#                    "plots": [filename],
#                    "freetext": expected_input},
#                filename=this_function.upper())
            self.reporting_structure.update(
                    {"GCOS requirements": 
                        {"plots": [filename],
                         "freetext": expected_input}})
        else:
#            self.__do_report__(
#                content={
#                    "plots": [filename]},
#                filename=this_function.upper())
            self.reporting_structure.update(
                    {"GCOS requirements": 
                        {"plots": [filename]}})

        return

    def __do_esm_evaluation__(self):

        this_function = "ESM evaluation"

        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_esmeval_input",
                                      expfile="_esmeval_expert.csv",
                                      protofile="empty_esmeval_expert.csv",
                                      function=this_function)

        # PART 1 of ESM evaluation: table
        # read in the ESM evaluation grading csv file
        esm_eval_input = os.path.dirname(os.path.realpath(
            __file__)) + "/lib/predef/example_eval_expert.csv"

        # calculate the length of the dataset
        ecv_length = int(
                self.__time_period__[5:10]) - int(self.__time_period__[0:4]
                ) + 1

        # plotting routines
        filename = self.__plot_dir__ + os.sep + self.__basic_filename__ + \
            "_" + "".join(this_function.split()) + "." + self.__output_type__
        fig = do_eval_table(
            self.__varname__,
            expected_input,
            os.path.dirname(
                os.path.realpath(__file__)) +
            "/libs/predef/example_eval_data.csv",
            ecv_length)
        fig.savefig(filename)
        plt.close(fig)

        caption = str(
            this_function +
            ' for the variable ' +
            ecv_lookup(
                self.__varname__) +
            ' in the data set "' +
            self.__dataset_id__[1] +
            '" (' +
            self.__time_period__ +
            ')' +
            ' (Petrol: data set is recommended for this application; ' + 
            'Brown: data set is not recommended for this application; ' + 
            'Grey: no decision about applicability of ' +
            'the data set can be made (e.g. uncertainty too high)).')

        # PART 2 of ESM evaluation: bullet point list
        # read in the ESM evaluation csv file
        esm_eval_input = os.path.dirname(os.path.realpath(
            __file__)) + "/libs/predef/esmeval_expert.csv"

        # Create a list. Each item of the list will be itself a list of strings, corresponding either to the
        # headers or to the ESM evaluation entries for the different ECVs
        contents = list()
        with open(esm_eval_input, 'r') as csvfile:
            s = csv.reader(csvfile, delimiter=";", skipinitialspace=True)
            for row in s:
                contents.append(row)

        # convert the list into an array
        esmeval_data = np.asarray(contents)

        # check the number of entries for the ECV in question, and write all of the available entries in an ordered dictionary
        # for easy output in the reports
        esmeval_dict = collections.OrderedDict()

        for num_entries in range(0, len(np.nonzero(
                esmeval_data == self.__varname__)[0])):
            insert_dict = collections.OrderedDict()
            for column in range(1, len(esmeval_data[0, :])):
                insert_dict.update({esmeval_data[0, column]:
                    esmeval_data[np.nonzero(esmeval_data ==
                                            self.__varname__)[0][num_entries],
                    column]})
            esmeval_dict.update({'R' + str(num_entries + 1): insert_dict})

        ESMValMD("meta",
                 filename,
                 self.__basetags__ + ['C3S_ESMeval'],
                 caption,
                 '#C3S' + 'ESMeval' + self.__varname__,
                 self.__infile__,
                 self.diagname,
                 self.authors)

        # produce report
        expected_input, found = \
            self.__file_anouncement__(subdir="c3s_511/single_esmeval_input",
                                      expfile="_esmeval.txt",
                                      protofile="empty.txt",
                                      function=this_function)

        if found:
#            self.__do_report__(
#                content={
#                    "listtext": esmeval_dict,
#                    "plots": [filename],
#                    "freetext": expected_input},
#                filename=this_function.upper())
            self.reporting_structure.update(
                    {"ESM evaluation": 
                        {"listtext": esmeval_dict,
                         "plots": [filename],
                         "freetext": expected_input}})
        else:
#            self.__do_report__(
#                content={
#                    "listtext": esmeval_dict,
#                    "plots": [filename]},
#                filename=this_function.upper())
            self.reporting_structure.update(
                    {"ESM evaluation": 
                        {"listtext": esmeval_dict,
                         "plots": [filename],
                         "freetext": expected_input}})

        return
