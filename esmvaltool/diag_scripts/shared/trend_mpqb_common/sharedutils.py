'''
# Author: Bas Crezee

# This script contains a few functions for utilizing parallelization along one axis of datasets
# that fit fully into memory.
'''

import multiprocessing
import numpy as np
from os import cpu_count
from concurrent.futures import ProcessPoolExecutor, as_completed


def wrap1dfunc(array,*args,**kwargs):
    # This function needs to follow calling signature for func1ds provided to np.apply_along_axis
    myfunc1d = kwargs.pop('diag1d')
    splitaxis = kwargs.pop('splitaxis')
    a,b=np.split(array,2,axis=splitaxis)
    return myfunc1d(a,b,*args,**kwargs)

def parallel_apply_along_axis(func1d, axis, arr, *args, **kwargs):
    """
    Like numpy.apply_along_axis(), but takes advantage of multiple
    cores.
    Taken from: https://stackoverflow.com/questions/45526700/easy-parallelization-of-numpy-apply-along-axis
    
    Caveats:
     - func1d can not be a lambda function, since that is unpicklable
     - only axis=0 works at the moment
     - number of workers is hard-coded
    """
    # Add func1d to the kwargs, since it might have to be passed on
    kwargs.update({'diag1d' : func1d})
    

    
    # Distinguish between case of single arr, or two arrs
    if type(arr)==tuple:
        print("Executing function on two arrays")
        # By convention, concatenate them along the given axis
        arr = np.concatenate(arr,axis=axis)
        # Pass on this axis in the kwargs to the function
        kwargs.update({'splitaxis' : axis})
        # Redefine the function
        func_to_call = wrap1dfunc
    else:
        print("Executing function on one array")
        func_to_call = kwargs.pop('diag1d')
    
    # Now get input dimensions
    indim = arr.ndim
    
    # Effective axis where apply_along_axis() will be applied by each
    # worker (any non-zero axis number would work, so as to allow the use
    # of `np.array_split()`, which is only done on axis 0):
    effective_axis = 1 if axis == 0 else axis
    if effective_axis != axis:
        arr = arr.swapaxes(axis, effective_axis)
        
    # Use 2/3 of available CPUs
    n_workers = int(cpu_count() / 1.5)

    # Chunks for the mapping (only a few chunks):
    chunks = [(func_to_call, effective_axis, sub_arr, args, kwargs)
              for sub_arr in np.array_split(arr, n_workers)]

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        results = executor.map(unpacking_apply_along_axis,chunks)

    result = np.concatenate(list(results),axis=axis)
    
    # Swap back but only if multiple values were returned
    if effective_axis != axis and result.ndim==indim:
        result = result.swapaxes(axis, effective_axis)
    return result

def unpacking_apply_along_axis(all_args):
    """
    Like numpy.apply_along_axis(), but with arguments in a tuple
    instead.

    This function is useful with multiprocessing.Pool().map(): (1)
    map() only handles functions that take a single argument, and (2)
    this function can generally be imported from a module, as required
    by map().
    Taken from: https://stackoverflow.com/questions/45526700/easy-parallelization-of-numpy-apply-along-axis
    """
    (func1d, axis, arr, args, kwargs) = all_args
    return np.apply_along_axis(func1d, axis, arr, *args, **kwargs)
