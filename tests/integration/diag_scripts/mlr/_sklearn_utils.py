"""Testing utilities for custom :mod:`sklearn` functionalities.

Parts of this code have been copied from :mod:`sklearn`.

License: BSD 3-Clause License

Copyright (c) 2007-2020 The scikit-learn developers.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

* Neither the name of the copyright holder nor the names of its
  contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""

import warnings

import numpy as np
import pytest
import scipy as sp
from numpy.testing import assert_allclose, assert_array_equal
from sklearn.base import BaseEstimator


class FailingClassifier(BaseEstimator):
    """Classifier that deliberately fails when using fit()."""

    FAILING_PARAMETER = 42

    def __init__(self, parameter=None):
        """Initialize class instance."""
        self.parameter = parameter

    def fit(self, *_, **__):
        """Fit."""
        if self.parameter == FailingClassifier.FAILING_PARAMETER:
            raise ValueError("Failing classifier failed as required")

    @staticmethod
    def predict(x_data):
        """Predict."""
        return np.zeros(x_data.shape[0])

    @staticmethod
    def score(*_, **__):
        """Score."""
        return 0.0


def assert_warns(warning_class, func, *args, **kw):
    """Test that a certain warning occurs.

    Parameters
    ----------
    warning_class : the warning class
        The class to test for, e.g. UserWarning.

    func : callable
        Callable object to trigger warnings.

    *args : the positional arguments to `func`.

    **kw : the keyword arguments to `func`

    Returns
    -------
    result : the return value of `func`

    """
    with warnings.catch_warnings(record=True) as warn:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        # Trigger a warning.
        result = func(*args, **kw)
        if hasattr(np, 'FutureWarning'):
            # Filter out numpy-specific warnings in numpy >= 1.9
            warn = [e for e in warn
                    if e.category is not np.VisibleDeprecationWarning]

        # Verify some things
        if not len(warn) > 0:
            raise AssertionError("No warning raised when calling %s"
                                 % func.__name__)

        found = any(warning.category is warning_class for warning in warn)
        if not found:
            raise AssertionError("%s did not give warning: %s( is %s)"
                                 % (func.__name__, warning_class, warn))
    return result


def assert_warns_message(warning_class, message, func, *args, **kw):
    # very important to avoid uncontrolled state propagation
    """Test that a certain warning occurs and with a certain message.

    Parameters
    ----------
    warning_class : the warning class
        The class to test for, e.g. UserWarning.

    message : str or callable
        The message or a substring of the message to test for. If callable,
        it takes a string as the argument and will trigger an AssertionError
        if the callable returns `False`.

    func : callable
        Callable object to trigger warnings.

    *args : the positional arguments to `func`.

    **kw : the keyword arguments to `func`.

    Returns
    -------
    result : the return value of `func`

    """
    with warnings.catch_warnings(record=True) as warn:
        # Cause all warnings to always be triggered.
        warnings.simplefilter("always")
        if hasattr(np, 'FutureWarning'):
            # Let's not catch the numpy internal DeprecationWarnings
            warnings.simplefilter('ignore', np.VisibleDeprecationWarning)
        # Trigger a warning.
        result = func(*args, **kw)
        # Verify some things
        if not len(warn) > 0:
            raise AssertionError("No warning raised when calling %s"
                                 % func.__name__)

        found = [issubclass(warning.category, warning_class) for warning in
                 warn]
        if not any(found):
            raise AssertionError("No warning raised for %s with class "
                                 "%s"
                                 % (func.__name__, warning_class))

        message_found = False
        # Checks the message of all warnings belong to warning_class
        for index in [i for i, x in enumerate(found) if x]:
            # substring will match, the entire message with typo won't
            msg = warn[index].message  # For Python 3 compatibility
            msg = str(msg.args[0] if hasattr(msg, 'args') else msg)
            if callable(message):  # add support for certain tests
                check_in_message = message
            else:
                def check_in_message(msg):
                    return message in msg

            if check_in_message(msg):
                message_found = True
                break

        if not message_found:
            raise AssertionError("Did not receive the message you expected "
                                 "('%s') for <%s>, got: '%s'"
                                 % (message, func.__name__, msg))

    return result


def assert_raise_message(exceptions, message, function, *args, **kwargs):
    """Test whether the message is in an exception.

    Given an exception, a callable to raise the exception, and
    a message string, tests that the correct exception is raised and
    that the message is a substring of the error thrown. Used to test
    that the specific message thrown during an exception is correct.

    Parameters
    ----------
    exceptions : exception or tuple of exception
        An Exception object.

    message : str
        The error message or a substring of the error message.

    function : callable
        Callable object to raise error.

    *args : the positional arguments to `function`.

    **kwargs : the keyword arguments to `function`.

    """
    try:
        function(*args, **kwargs)
    except exceptions as exc:
        error_message = str(exc)
        if message not in error_message:
            raise AssertionError("Error message does not include the expected "
                                 "string: %r. Observed error message: %r" %
                                 (message, error_message))
    else:
        # concatenate exception names
        if isinstance(exceptions, tuple):
            names = " or ".join(e.__name__ for e in exceptions)
        else:
            names = exceptions.__name__

        raise AssertionError("%s not raised by %s" %
                             (names, function.__name__))


def assert_allclose_dense_sparse(x_arr, y_arr, rtol=1e-07, atol=1e-9,
                                 err_msg=''):
    """Assert allclose for sparse and dense data.

    Both x_arr and y_arr need to be either sparse or dense, they
    can't be mixed.

    Parameters
    ----------
    x_arr : {array-like, sparse matrix}
        First array to compare.

    y_arr : {array-like, sparse matrix}
        Second array to compare.

    rtol : float, default=1e-07
        relative tolerance; see numpy.allclose.

    atol : float, default=1e-9
        absolute tolerance; see numpy.allclose. Note that the default here is
        more tolerant than the default for numpy.testing.assert_allclose, where
        atol=0.

    err_msg : str, default=''
        Error message to raise.

    """
    if sp.sparse.issparse(x_arr) and sp.sparse.issparse(y_arr):
        x_arr = x_arr.tocsr()
        y_arr = y_arr.tocsr()
        x_arr.sum_duplicates()
        y_arr.sum_duplicates()
        assert_array_equal(x_arr.indices, y_arr.indices, err_msg=err_msg)
        assert_array_equal(x_arr.indptr, y_arr.indptr, err_msg=err_msg)
        assert_allclose(x_arr.data, y_arr.data, rtol=rtol, atol=atol,
                        err_msg=err_msg)
    elif not sp.sparse.issparse(x_arr) and not sp.sparse.issparse(y_arr):
        # both dense
        assert_allclose(x_arr, y_arr, rtol=rtol, atol=atol, err_msg=err_msg)
    else:
        raise ValueError("Can only compare two sparse matrices, not a sparse "
                         "matrix and an array.")


def _convert_container(container, constructor_name, columns_name=None):
    """Convert container."""
    if constructor_name == 'list':
        return list(container)
    if constructor_name == 'tuple':
        return tuple(container)
    if constructor_name == 'array':
        return np.asarray(container)
    if constructor_name == 'sparse':
        return sp.sparse.csr_matrix(container)
    if constructor_name == 'dataframe':
        pandas = pytest.importorskip('pandas')
        return pandas.DataFrame(container, columns=columns_name)
    if constructor_name == 'series':
        pandas = pytest.importorskip('pandas')
        return pandas.Series(container)
    if constructor_name == 'index':
        pandas = pytest.importorskip('pandas')
        return pandas.Index(container)
    if constructor_name == 'slice':
        return slice(container[0], container[1])
    if constructor_name == 'sparse_csr':
        return sp.sparse.csr_matrix(container)
    if constructor_name == 'sparse_csc':
        return sp.sparse.csc_matrix(container)
    raise TypeError(f"Constructor name '{constructor_name}' not supported")
