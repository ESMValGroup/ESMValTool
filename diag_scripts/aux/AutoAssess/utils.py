# module of utility functions required by various modules
import numbers
import datetime

class MetricsError(Exception): pass

def combine_metrics( mtrc_list ):
    '''
    Combine multiple metrics dictionaries into a single dict.
    Check for duplicates and raise an exception if there are any.

    Arguments:
        mtrc_list - a sequence of metrics dictionaries.  Each dictionary 
                    contains name-value pairs for each metric

    Returns:
        metrics - a single dictionary containing the combined metrics

    todo?: check for gaps against a control list (likely csv file)

    '''

    metrics = dict()

    for mdct in mtrc_list:
        for name, value in mdct.items():
            # Check this is a valid name-number pair
            # (It's generally recommended not to test if a python object is a 
            # string, but in this case it makes sense as if it is not a string
            # it suggests the area code has done something wrong.)
            assert isinstance(name, basestring), \
                    "Metric name is not a string: %s" % (name,)
            assert isinstance(value, numbers.Number), \
                    "Value for metric %s is not a number: %s" % (name,value)

            # Check is not a duplicate
            if name in metrics:
                msg = 'Metric name %s appears more than once' % (name,)
                raise DuplicateError, msg

            metrics[name] = value

    return metrics


def make_run_dict():
    '''
    Generate a run dictionary for a single run, which can be passed to 
    top.top() for testing metrics code for an individual area.  It can also be 
    used to test lower-level components of an area, since these can take the 
    same run dict as an argument (see land_surface.top.top() for example).

    '''

    runs = make_run_dict_multi()
    return runs[0]


def make_run_dict_multi():
    '''
    Generate a list of runs which can be passed to top.multi() for testing
    plotting code for an individual area.

    todo: find out from Paul why no ss_* dirs
          then add third run to list for more general testing (not just a vs b)
          also write instructions for making it easy to import this module (put
          __init__ file in maverick dirs to make it easier?)

    todo?: write a test function which calls a chunk of metrics code with a
           standard run, and checks the validity of the returned metrics (a sort
           of unit test).  Not sure how useful this would be.
    '''

    run1 = {'runid': 'amzgg',
            'longname': 'amzgg',
            'run_type': 'AMIP',
            'ocean_model': 'NEMO',
            'start': 1982.0,
            'nyear': 25.0,
            'validation_note': '/data/cr1/hadtq/valnote/antia_v_amzgg',
            'summary_file': '/data/cr1/hadtq/valnote/antia_v_amzgg/summary_global.csv',
            'radiation_table': '/data/cr1/hadtq/valnote/antia_v_amzgg/rad.txt',
            'ss_annual': '',
            'ss_seasonal': '',
            'ss_monthly': '',
            'ss_daily': '',
            'ss_spatiot': '',
            'data_root': '/project/hadgem3/data/amzgg',
            'supermean_root': '/project/hadgem3/data/amzgg/supermeans/amzgga.ms2006',
            'from_daily': datetime.date(1982, 1, 1),
            'to_daily': datetime.date(1986, 12, 1),
            'from_monthly': datetime.date(1982, 1, 1),
            'to_monthly': datetime.date(2006, 12, 1),
            'from_annual': datetime.date(1981, 1, 1),
            'to_annual': datetime.date(2005, 1, 1) 
           }

    run2 = {'runid': 'antia',
            'longname': 'antia',
            'run_type': 'AMIP',
            'ocean_model': 'NEMO',
            'start': 1982.0,
            'nyear': 25.0,
            'validation_note': '/data/cr1/hadtq/valnote/antia_v_amzgg',
            'summary_file': '/data/cr1/hadtq/valnote/antia_v_amzgg/summary_global.csv',
            'radiation_table': '/data/cr1/hadtq/valnote/antia_v_amzgg/rad.txt',
            'ss_annual': '',
            'ss_monthly': '',
            'ss_seasonal': '',
            'ss_daily': '',
            'ss_spatiot': '',
            'data_root': '/project/hadgem3/data/antia',
            'supermean_root': '/project/hadgem3/data/antia/supermeans/antiaa.ms2006',
            'from_daily': datetime.date(1982, 1, 1),
            'to_daily': datetime.date(1986, 12, 1),
            'from_monthly': datetime.date(1982, 1, 1),
            'to_monthly': datetime.date(2006, 12, 1),
            'from_annual': datetime.date(1981, 1, 1),
            'to_annual': datetime.date(2005, 1, 1)
           }

    runs = [run1, run2]

    return runs

