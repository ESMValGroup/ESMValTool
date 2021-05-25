"""Integration tests for functions of custom :mod:`sklearn` functionalities.

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

# pylint: disable=arguments-differ
# pylint: disable=invalid-name
# pylint: disable=no-self-use
# pylint: disable=protected-access
# pylint: disable=redefined-outer-name
# pylint: disable=too-few-public-methods
# pylint: disable=too-many-arguments

from copy import copy, deepcopy

import numpy as np
import pytest
import scipy.sparse as sp
from sklearn import datasets
from sklearn.base import BaseEstimator
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import KernelPCA
from sklearn.exceptions import FitFailedWarning, NotFittedError
from sklearn.linear_model import LinearRegression
from sklearn.metrics import (
    explained_variance_score,
    make_scorer,
    mean_absolute_error,
    mean_squared_error,
)
from sklearn.model_selection import LeaveOneGroupOut, ShuffleSplit
from sklearn.svm import SVC

from esmvaltool.diag_scripts.mlr.custom_sklearn import (
    _DEFAULT_TAGS,
    AdvancedPipeline,
    AdvancedRFE,
    AdvancedRFECV,
    FeatureSelectionTransformer,
    _check_fit_params,
    _determine_key_type,
    _fit_and_score_weighted,
    _get_fit_parameters,
    _is_pairwise,
    _map_features,
    _num_samples,
    _rfe_single_fit,
    _safe_indexing,
    _safe_split,
    _safe_tags,
    _score_weighted,
    _split_fit_kwargs,
    _update_transformers_param,
    cross_val_score_weighted,
    get_rfecv_transformer,
    perform_efecv,
)

from ._sklearn_utils import (
    FailingClassifier,
    _convert_container,
    assert_allclose_dense_sparse,
    assert_raise_message,
    assert_warns,
    assert_warns_message,
)

# _determine_key_type


TEST_DETERMINE_KEY_TYPE = [
    (0, 'int'),
    ('0', 'str'),
    (True, 'bool'),
    (np.bool_(True), 'bool'),
    ([0, 1, 2], 'int'),
    (['0', '1', '2'], 'str'),
    ((0, 1, 2), 'int'),
    (('0', '1', '2'), 'str'),
    (slice(None, None), None),
    (slice(0, 2), 'int'),
    (np.array([0, 1, 2], dtype=np.int32), 'int'),
    (np.array([0, 1, 2], dtype=np.int64), 'int'),
    (np.array([0, 1, 2], dtype=np.uint8), 'int'),
    ([True, False], 'bool'),
    ((True, False), 'bool'),
    (np.array([True, False]), 'bool'),
    ('col_0', 'str'),
    (['col_0', 'col_1', 'col_2'], 'str'),
    (('col_0', 'col_1', 'col_2'), 'str'),
    (slice('begin', 'end'), 'str'),
    (np.array(['col_0', 'col_1', 'col_2']), 'str'),
    (np.array(['col_0', 'col_1', 'col_2'], dtype=object), 'str'),
]


@pytest.mark.parametrize('key,dtype', TEST_DETERMINE_KEY_TYPE)
def test_determine_key_type(key, dtype):
    """Test working ``_determine_key_type``."""
    assert _determine_key_type(key) == dtype


def test_determine_key_type_error():
    """Test failing ``_determine_key_type``."""
    with pytest.raises(ValueError, match="No valid specification of the"):
        _determine_key_type(1.0)


def test_determine_key_type_slice_error():
    """Test failing ``_determine_key_type``."""
    with pytest.raises(TypeError, match="Only array-like or scalar are"):
        _determine_key_type(slice(0, 2, 1), accept_slice=False)


# _safe_indexing


@pytest.mark.parametrize(
    'array_type', ['list', 'array', 'sparse', 'dataframe'],
)
@pytest.mark.parametrize(
    'indices_type', ['list', 'tuple', 'array', 'series', 'slice'],
)
def test_safe_indexing_2d_container_axis_0(array_type, indices_type):
    """Test ``_safe_indexing`` with 2D container."""
    indices = [1, 2]
    if indices_type == 'slice' and isinstance(indices[1], int):
        indices[1] += 1
    array = _convert_container([[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type)
    indices = _convert_container(indices, indices_type)
    subset = _safe_indexing(array, indices, axis=0)
    assert_allclose_dense_sparse(
        subset, _convert_container([[4, 5, 6], [7, 8, 9]], array_type)
    )


X_DATA_TOY = np.arange(9).reshape((3, 3))


@pytest.mark.parametrize('array_type', ['list', 'array', 'series'])
@pytest.mark.parametrize(
    'indices_type', ['list', 'tuple', 'array', 'series', 'slice'],
)
def test_safe_indexing_1d_container(array_type, indices_type):
    """Test ``_safe_indexing`` with 1D container."""
    indices = [1, 2]
    if indices_type == 'slice' and isinstance(indices[1], int):
        indices[1] += 1
    array = _convert_container([1, 2, 3, 4, 5, 6, 7, 8, 9], array_type)
    indices = _convert_container(indices, indices_type)
    subset = _safe_indexing(array, indices, axis=0)
    assert_allclose_dense_sparse(
        subset, _convert_container([2, 3], array_type)
    )


@pytest.mark.parametrize('array_type', ['array', 'sparse', 'dataframe'])
@pytest.mark.parametrize(
    'indices_type', ['list', 'tuple', 'array', 'series', 'slice'],
)
@pytest.mark.parametrize('indices', [[1, 2], ['col_1', 'col_2']])
def test_safe_indexing_2d_container_axis_1(array_type, indices_type, indices):
    """Test ``_safe_indexing`` with 2D container."""
    # validation of the indices
    # we make a copy because indices is mutable and shared between tests
    indices_converted = copy(indices)
    if indices_type == 'slice' and isinstance(indices[1], int):
        indices_converted[1] += 1

    columns_name = ['col_0', 'col_1', 'col_2']
    array = _convert_container(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type, columns_name
    )
    indices_converted = _convert_container(indices_converted, indices_type)

    if isinstance(indices[0], str) and array_type != 'dataframe':
        err_msg = ("Specifying the columns using strings is only supported "
                   "for pandas DataFrames")
        with pytest.raises(ValueError, match=err_msg):
            _safe_indexing(array, indices_converted, axis=1)
    else:
        subset = _safe_indexing(array, indices_converted, axis=1)
        assert_allclose_dense_sparse(
            subset, _convert_container([[2, 3], [5, 6], [8, 9]], array_type)
        )


@pytest.mark.parametrize('array_read_only', [True, False])
@pytest.mark.parametrize('indices_read_only', [True, False])
@pytest.mark.parametrize('array_type', ['array', 'sparse', 'dataframe'])
@pytest.mark.parametrize('indices_type', ['array', 'series'])
@pytest.mark.parametrize(
    'axis, expected_array',
    [(0, [[4, 5, 6], [7, 8, 9]]), (1, [[2, 3], [5, 6], [8, 9]])],
)
def test_safe_indexing_2d_read_only_axis_1(array_read_only, indices_read_only,
                                           array_type, indices_type, axis,
                                           expected_array):
    """Test ``_safe_indexing`` with 2D container."""
    array = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    if array_read_only:
        array.setflags(write=False)
    array = _convert_container(array, array_type)
    indices = np.array([1, 2])
    if indices_read_only:
        indices.setflags(write=False)
    indices = _convert_container(indices, indices_type)
    subset = _safe_indexing(array, indices, axis=axis)
    assert_allclose_dense_sparse(
        subset, _convert_container(expected_array, array_type)
    )


@pytest.mark.parametrize('array_type', ['list', 'array', 'series'])
@pytest.mark.parametrize('indices_type', ['list', 'tuple', 'array', 'series'])
def test_safe_indexing_1d_container_mask(array_type, indices_type):
    """Test ``_safe_indexing`` with 1D container."""
    indices = [False] + [True] * 2 + [False] * 6
    array = _convert_container([1, 2, 3, 4, 5, 6, 7, 8, 9], array_type)
    indices = _convert_container(indices, indices_type)
    subset = _safe_indexing(array, indices, axis=0)
    assert_allclose_dense_sparse(
        subset, _convert_container([2, 3], array_type)
    )


@pytest.mark.parametrize('array_type', ['array', 'sparse', 'dataframe'])
@pytest.mark.parametrize('indices_type', ['list', 'tuple', 'array', 'series'])
@pytest.mark.parametrize(
    'axis, expected_subset',
    [(0, [[4, 5, 6], [7, 8, 9]]),
     (1, [[2, 3], [5, 6], [8, 9]])],
)
def test_safe_indexing_2d_mask(array_type, indices_type, axis,
                               expected_subset):
    """Test ``_safe_indexing`` with 2D container."""
    columns_name = ['col_0', 'col_1', 'col_2']
    array = _convert_container(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type, columns_name
    )
    indices = [False, True, True]
    indices = _convert_container(indices, indices_type)

    subset = _safe_indexing(array, indices, axis=axis)
    assert_allclose_dense_sparse(
        subset, _convert_container(expected_subset, array_type)
    )


@pytest.mark.parametrize(
    'array_type, expected_output_type',
    [('list', 'list'), ('array', 'array'),
     ('sparse', 'sparse'), ('dataframe', 'series')],
)
def test_safe_indexing_2d_scalar_axis_0(array_type, expected_output_type):
    """Test ``_safe_indexing`` with 2D container."""
    array = _convert_container([[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type)
    indices = 2
    subset = _safe_indexing(array, indices, axis=0)
    expected_array = _convert_container([7, 8, 9], expected_output_type)
    assert_allclose_dense_sparse(subset, expected_array)


@pytest.mark.parametrize('array_type', ['list', 'array', 'series'])
def test_safe_indexing_1d_scalar(array_type):
    """Test ``_safe_indexing`` with 1D container."""
    array = _convert_container([1, 2, 3, 4, 5, 6, 7, 8, 9], array_type)
    indices = 2
    subset = _safe_indexing(array, indices, axis=0)
    assert subset == 3


@pytest.mark.parametrize(
    'array_type, expected_output_type',
    [('array', 'array'), ('sparse', 'sparse'), ('dataframe', 'series')],
)
@pytest.mark.parametrize('indices', [2, 'col_2'])
def test_safe_indexing_2d_scalar_axis_1(array_type, expected_output_type,
                                        indices):
    """Test ``_safe_indexing`` with 1D container."""
    columns_name = ['col_0', 'col_1', 'col_2']
    array = _convert_container(
        [[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type, columns_name
    )

    if isinstance(indices, str) and array_type != 'dataframe':
        err_msg = ("Specifying the columns using strings is only supported "
                   "for pandas DataFrames")
        with pytest.raises(ValueError, match=err_msg):
            _safe_indexing(array, indices, axis=1)
    else:
        subset = _safe_indexing(array, indices, axis=1)
        expected_output = [3, 6, 9]
        if expected_output_type == 'sparse':
            # sparse matrix are keeping the 2D shape
            expected_output = [[3], [6], [9]]
        expected_array = _convert_container(
            expected_output, expected_output_type
        )
        assert_allclose_dense_sparse(subset, expected_array)


@pytest.mark.parametrize('array_type', ['list', 'array', 'sparse'])
def test_safe_indexing_none_axis_0(array_type):
    """Test ``_safe_indexing`` with None."""
    x_data = _convert_container([[1, 2, 3], [4, 5, 6], [7, 8, 9]], array_type)
    x_data_subset = _safe_indexing(x_data, None, axis=0)
    assert_allclose_dense_sparse(x_data_subset, x_data)


def test_safe_indexing_pandas_no_matching_cols_error():
    """Test ``_safe_indexing`` with pandas."""
    pd = pytest.importorskip('pandas')
    err_msg = "No valid specification of the columns."
    x_data = pd.DataFrame(X_DATA_TOY)
    with pytest.raises(ValueError, match=err_msg):
        _safe_indexing(x_data, [1.0], axis=1)


@pytest.mark.parametrize('axis', [None, 3])
def test_safe_indexing_error_axis(axis):
    """Test ``_safe_indexing`` error."""
    with pytest.raises(ValueError, match="'axis' should be either 0"):
        _safe_indexing(X_DATA_TOY, [0, 1], axis=axis)


@pytest.mark.parametrize('x_constructor', ['array', 'series'])
def test_safe_indexing_1d_array_error(x_constructor):
    """Test ``_safe_indexing`` error."""
    # check that we are raising an error if the array-like passed is 1D and
    # we try to index on the 2nd dimension
    x_data = list(range(5))
    if x_constructor == 'array':
        x_constructor = np.asarray(x_data)
    elif x_constructor == 'series':
        pd = pytest.importorskip("pandas")
        x_constructor = pd.Series(x_data)

    err_msg = "'x_data' should be a 2D NumPy array, 2D sparse matrix or pandas"
    with pytest.raises(ValueError, match=err_msg):
        _safe_indexing(x_constructor, [0, 1], axis=1)


def test_safe_indexing_container_axis_0_unsupported_type():
    """Test ``_safe_indexing`` error."""
    indices = ["col_1", "col_2"]
    array = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
    err_msg = "String indexing is not supported with 'axis=0'"
    with pytest.raises(ValueError, match=err_msg):
        _safe_indexing(array, indices, axis=0)


# _num_samples


def test_retrieve_samples_from_non_standard_shape():
    """Test ``_num_samples``."""
    class TestNonNumericShape:
        """Non-numeric shape."""

        def __init__(self):
            """Init."""
            self.shape = ("not numeric",)

        def __len__(self):
            """Length."""
            return len([1, 2, 3])

    x_data = TestNonNumericShape()
    assert _num_samples(x_data) == len(x_data)

    # Check that it gives a good error if there's no __len__
    class TestNoLenWeirdShape:
        """Weird shape with no length."""

        def __init__(self):
            """Init."""
            self.shape = ("not numeric",)

    with pytest.raises(TypeError, match="Expected sequence or array-like"):
        _num_samples(TestNoLenWeirdShape())


# _check_fit_params


@pytest.mark.parametrize('indices', [None, [1, 3]])
def test_check_fit_params(indices):
    """Test ``_check_fit_params``."""
    x_data = np.random.randn(4, 2)
    fit_params = {
        'list': [1, 2, 3, 4],
        'array': np.array([1, 2, 3, 4]),
        'sparse-col': sp.csc_matrix([1, 2, 3, 4]).T,
        'sparse-row': sp.csc_matrix([1, 2, 3, 4]),
        'scalar-int': 1,
        'scalar-str': 'xxx',
        'None': None,
    }
    result = _check_fit_params(x_data, fit_params, indices)
    indices_ = indices if indices is not None else list(range(x_data.shape[0]))

    for key in ['sparse-row', 'scalar-int', 'scalar-str', 'None']:
        assert result[key] is fit_params[key]

    assert result['list'] == _safe_indexing(fit_params['list'], indices_)
    np.testing.assert_array_equal(
        result['array'], _safe_indexing(fit_params['array'], indices_)
    )
    assert_allclose_dense_sparse(
        result['sparse-col'],
        _safe_indexing(fit_params['sparse-col'], indices_)
    )


# _safe_tags


class NoTagsEstimator:
    """Estimators with no tags."""


class MoreTagsEstimator:
    """Estimator with ``_more_tags``."""

    def _more_tags(self):
        """Return more tags."""
        return {"allow_nan": True}


@pytest.mark.parametrize(
    'estimator,err_msg',
    [
        (BaseEstimator(), 'The key xxx is not defined in _get_tags'),
        (NoTagsEstimator(), 'The key xxx is not defined in _DEFAULT_TAGS'),
    ],
)
def test_safe_tags_error(estimator, err_msg):
    """Test ``_safe_tags`` with error."""
    # Check that safe_tags raises error in ambiguous case.
    with pytest.raises(ValueError, match=err_msg):
        _safe_tags(estimator, key="xxx")


@pytest.mark.parametrize(
    'estimator,key,expected_results',
    [
        (NoTagsEstimator(), None, _DEFAULT_TAGS),
        (NoTagsEstimator(), 'allow_nan', _DEFAULT_TAGS['allow_nan']),
        (MoreTagsEstimator(), None, {**_DEFAULT_TAGS, **{'allow_nan': True}}),
        (MoreTagsEstimator(), 'allow_nan', True),
        (BaseEstimator(), None, _DEFAULT_TAGS),
        (BaseEstimator(), 'allow_nan', _DEFAULT_TAGS['allow_nan']),
        (BaseEstimator(), 'allow_nan', _DEFAULT_TAGS['allow_nan']),
    ],
)
def test_safe_tags_no_get_tags(estimator, key, expected_results):
    """Test ``_safe_tags`` without ``_get_tags``."""
    assert _safe_tags(estimator, key=key) == expected_results


# _is_pairwise


def test_is_pairwise():
    """Test ``_is_pairwise``."""
    # Simple checks for _is_pairwise
    pca = KernelPCA(kernel='precomputed')
    with pytest.warns(None) as record:
        assert _is_pairwise(pca)
    assert not record

    # Pairwise attribute that is not consistent with the pairwise tag
    class IncorrectTagPCA(KernelPCA):
        """Class with incorrect _pairwise attribute."""

        _pairwise = False

    pca = IncorrectTagPCA(kernel='precomputed')
    msg = "_pairwise attribute is inconsistent with tags."
    with pytest.warns(FutureWarning, match=msg):
        assert not _is_pairwise(pca)

    # The _pairwise attribute is present and set to True while pairwise tag is
    # not present
    class TruePairwise(BaseEstimator):
        """Class without pairwise tag."""

        _pairwise = True

    true_pairwise = TruePairwise()
    with pytest.warns(FutureWarning, match=msg):
        assert _is_pairwise(true_pairwise)

    # Pairwise attribute is not defined thus tag is used
    est = BaseEstimator()
    with pytest.warns(None) as record:
        assert not _is_pairwise(est)
    assert not record


# _safe_split


def test_safe_split():
    """Test ``_safe_split``."""
    clf = SVC()
    clfp = SVC(kernel="precomputed")

    iris = datasets.load_iris()
    (x_data, y_data) = (iris.data, iris.target)
    kernel = np.dot(x_data, x_data.T)

    cv = ShuffleSplit(test_size=0.25, random_state=0)
    train, test = list(cv.split(x_data))[0]

    x_train, y_train = _safe_split(clf, x_data, y_data, train)
    kernel_train, y_train2 = _safe_split(clfp, kernel, y_data, train)
    np.testing.assert_array_almost_equal(kernel_train, np.dot(x_train,
                                                              x_train.T))
    np.testing.assert_array_almost_equal(y_train, y_train2)

    x_test, y_test = _safe_split(clf, x_data, y_data, test, train)
    kernel_test, y_test2 = _safe_split(clfp, kernel, y_data, test, train)
    np.testing.assert_array_almost_equal(kernel_test, np.dot(x_test,
                                                             x_train.T))
    np.testing.assert_array_almost_equal(y_test, y_test2)


# _fit_and_score_weighted


def test_fit_and_score_weighted_failing():
    """Test if ``_fit_and_score_weighted`` fails as expected."""
    # Create a failing classifier to deliberately fail
    failing_clf = FailingClassifier(FailingClassifier.FAILING_PARAMETER)

    # Dummy X data
    x_data = np.arange(1, 10)
    scorer = make_scorer(mean_squared_error)
    fit_and_score_args = [failing_clf, x_data, None, scorer, None, None, None,
                          None]

    # Passing error score to trigger the warning message
    fit_and_score_kwargs = {'error_score': 42}

    # Check if the warning message type is as expected
    assert_warns(FitFailedWarning, _fit_and_score_weighted,
                 *fit_and_score_args, **fit_and_score_kwargs)

    # Since we're using FailingClassifier, our error will be the following
    error_message = "ValueError: Failing classifier failed as required"

    # The warning message we're expecting to see
    warning_message = ("Estimator fit failed. The score on this train-test "
                       "partition for these parameters will be set to %f. "
                       "Details: \n%s" % (fit_and_score_kwargs['error_score'],
                                          error_message))

    def test_warn_trace(msg):
        """Traceback in warning message."""
        assert 'Traceback (most recent call last):\n' in msg
        split = msg.splitlines()  # note: handles more than '\n'
        mtb = split[0] + '\n' + split[-1]
        return warning_message in mtb

    # Check traceback is included
    assert_warns_message(FitFailedWarning, test_warn_trace,
                         _fit_and_score_weighted, *fit_and_score_args,
                         **fit_and_score_kwargs)

    # Check return of error_score in case of failed fit
    result = _fit_and_score_weighted(*fit_and_score_args,
                                     **fit_and_score_kwargs)
    assert isinstance(result, int)
    assert result == fit_and_score_kwargs['error_score']

    # Check if exception was raised, with default error_score='raise'
    fit_and_score_kwargs = {'error_score': 'raise'}
    assert_raise_message(ValueError, "Failing classifier failed as required",
                         _fit_and_score_weighted, *fit_and_score_args,
                         **fit_and_score_kwargs)

    # Wrong parameter type for error_score
    fit_and_score_kwargs = {'error_score': 'wrong_type'}
    with pytest.raises(ValueError):
        _fit_and_score_weighted(*fit_and_score_args, **fit_and_score_kwargs)

    assert failing_clf.score() == 0.0


X_DATA = np.arange(5).reshape(5, 1)
Y_DATA = np.array([0, 1, 2, 3, -1])
TRAIN = np.array([1, 2, 3, 4])
TEST = np.array([0])


TEST_FIT_AND_SCORE_WEIGHTED_NO_WEIGHTS = [
    (make_scorer(mean_absolute_error), 2.5),
    (make_scorer(mean_squared_error), 6.25),
    (make_scorer(explained_variance_score), 1.0),
]


@pytest.mark.parametrize('scorer,output',
                         TEST_FIT_AND_SCORE_WEIGHTED_NO_WEIGHTS)
def test_fit_and_score_weighted_no_weights(scorer, output):
    """Test ``_fit_and_score_weighted`` without weights."""
    clf = LinearRegression()
    fit_and_score_weighted_args = [clf, X_DATA, Y_DATA, scorer, TRAIN, TEST]
    fit_and_score_weighted_kwargs = {
        'parameters': {'normalize': False},
        'fit_params': None,
    }

    result = _fit_and_score_weighted(*fit_and_score_weighted_args,
                                     **fit_and_score_weighted_kwargs)
    np.testing.assert_allclose(clf.coef_, [-0.5])
    np.testing.assert_allclose(clf.intercept_, 2.5)
    assert isinstance(result, float)
    np.testing.assert_allclose(result, output)


SAMPLE_WEIGHTS = np.array([1.0, 1.0, 1.0, 1.0, 0.0])
TEST_FIT_AND_SCORE_WEIGHTED_WEIGHTS = [
    (make_scorer(mean_absolute_error), 0.0),
    (make_scorer(mean_squared_error), 0.0),
    (make_scorer(explained_variance_score), 1.0),
]


@pytest.mark.parametrize('scorer,output', TEST_FIT_AND_SCORE_WEIGHTED_WEIGHTS)
def test_fit_and_score_weighted_weights(scorer, output):
    """Test ``_fit_and_score_weighted`` with weights."""
    clf = LinearRegression()
    fit_and_score_weighted_args = [clf, X_DATA, Y_DATA, scorer, TRAIN, TEST]
    fit_and_score_weighted_kwargs = {
        'parameters': {'normalize': False},
        'fit_params': {'sample_weight': SAMPLE_WEIGHTS},
        'sample_weights': SAMPLE_WEIGHTS,
    }

    result = _fit_and_score_weighted(*fit_and_score_weighted_args,
                                     **fit_and_score_weighted_kwargs)
    np.testing.assert_allclose(clf.coef_, [1.0])
    np.testing.assert_allclose(clf.intercept_, 0.0, atol=1e-10)
    assert isinstance(result, float)
    np.testing.assert_allclose(result, output, atol=1e-10)


# _get_fit_parameters


STEPS_1 = [('a', 1)]
STEPS_2 = [('a', 1), ('b', 0)]
TEST_GET_FIT_PARAMETERS = [
    ({'a': 1}, STEPS_1, ValueError),
    ({'a': 1, 'a__b': 1}, STEPS_1, ValueError),
    ({'a__x': 1}, [], ValueError),
    ({'a__x': 1}, STEPS_1, {'a': {'x': 1}}),
    ({'a__x': 1, 'a__y': 2}, STEPS_1, {'a': {'x': 1, 'y': 2}}),
    ({'a__x': 1, 'a__y__z': 2}, STEPS_1, {'a': {'x': 1, 'y__z': 2}}),
    ({'a__x': 1, 'b__y': 2}, STEPS_1, ValueError),
    ({'a__x': 1, 'b__y': 2}, STEPS_2, {'a': {'x': 1}, 'b': {'y': 2}}),
]


@pytest.mark.parametrize('kwargs,steps,output', TEST_GET_FIT_PARAMETERS)
def test_get_fit_parameters(kwargs, steps, output):
    """Test retrieving of fit parameters."""
    if isinstance(output, type):
        with pytest.raises(output):
            _get_fit_parameters(kwargs, steps, 'x')
        return
    params = _get_fit_parameters(kwargs, steps, 'x')
    assert params == output


# _score_weighted


def test_score_weighted_failing():
    """Test if ``_score_weighted`` fails as expected."""
    error_message = "Scoring must return a number, got None"

    def two_params_scorer(*_, **__):
        """Scorer function."""
        return None

    score_args = [None, None, None, two_params_scorer]
    assert_raise_message(ValueError, error_message, _score_weighted,
                         *score_args)


TEST_SCORE_WEIGHTED_NO_WEIGHTS = [
    (make_scorer(mean_absolute_error), 1.2),
    (make_scorer(mean_squared_error), 2.0),
    (make_scorer(explained_variance_score), 0.0),
]


@pytest.mark.parametrize('scorer,output', TEST_SCORE_WEIGHTED_NO_WEIGHTS)
def test_score_weighted_no_weights(scorer, output):
    """Test ``_score_weighted`` without weights."""
    clf = LinearRegression()
    clf.fit(X_DATA, Y_DATA)
    np.testing.assert_allclose(clf.coef_, [0.0], atol=1e-10)
    np.testing.assert_allclose(clf.intercept_, 1.0)
    result = _score_weighted(clf, X_DATA, Y_DATA, scorer)
    assert isinstance(result, float)
    np.testing.assert_allclose(result, output, atol=1e-10)


TEST_SCORE_WEIGHTED_WEIGHTS = [
    (make_scorer(mean_absolute_error), 0.0),
    (make_scorer(mean_squared_error), 0.0),
    (make_scorer(explained_variance_score), 1.0),
]


@pytest.mark.parametrize('scorer,output', TEST_SCORE_WEIGHTED_WEIGHTS)
def test_score_weighted_weights(scorer, output):
    """Test ``_score_weighted`` with weights."""
    clf = LinearRegression()
    clf.fit(X_DATA, Y_DATA, sample_weight=SAMPLE_WEIGHTS)
    np.testing.assert_allclose(clf.coef_, [1.0])
    np.testing.assert_allclose(clf.intercept_, 0.0, atol=1e-10)
    result = _score_weighted(clf, X_DATA, Y_DATA, scorer,
                             sample_weights=SAMPLE_WEIGHTS)
    assert isinstance(result, float)
    np.testing.assert_allclose(result, output, atol=1e-10)


# _split_fit_kwargs


def test_split_fit_kwargs():
    """Test ``_split_fit_kwargs``."""
    fit_kwargs = {
        'a': [1, 2],
        'b': [1, 2],
        'sample_weight': [1, 2],
        'a__sample_weight__b': [1, 2],
        'sample_weight_eval_set': [1, 2],
        'a__sample_weight_eval_set__b': [1, 2],
    }
    train_idx = 0
    test_idx = 1

    (train_kwargs, test_kwargs) = _split_fit_kwargs(fit_kwargs, train_idx,
                                                    test_idx)

    assert train_kwargs is not fit_kwargs
    assert test_kwargs is not fit_kwargs
    assert len(train_kwargs) == len(fit_kwargs)
    assert len(test_kwargs) == len(fit_kwargs)
    for key in train_kwargs:
        assert train_kwargs[key] is not fit_kwargs[key]
    for key in test_kwargs:
        assert test_kwargs[key] is not fit_kwargs[key]

    expected_train_kwargs = {
        'a': [1, 2],
        'b': [1, 2],
        'sample_weight': 1,
        'a__sample_weight__b': 1,
        'sample_weight_eval_set': [1, 2],
        'a__sample_weight_eval_set__b': [1, 2],
    }
    assert train_kwargs == expected_train_kwargs

    expected_test_kwargs = {
        'a': [1, 2],
        'b': [1, 2],
        'sample_weight': 2,
        'a__sample_weight__b': 2,
        'sample_weight_eval_set': [1, 2],
        'a__sample_weight_eval_set__b': [1, 2],
    }
    assert test_kwargs == expected_test_kwargs


# _rfe_single_fit


@pytest.fixture
def advanced_rfe():
    """``AdvancedRFE`` object."""
    rfe_kwargs = {
        'estimator': LinearRegression(),
        'n_features_to_select': 1,
    }
    return AdvancedRFE(**rfe_kwargs)


def test_rfe_single_fit(advanced_rfe):
    """Test ``_rfe_single_fit``."""
    x_data = np.array(
        [[0.0, 0.0, 0.0],
         [2.0, 0.0, 1.0],
         [3.0, 0.0, -2.0],
         [4.0, 0.0, -4.0]],
    )
    y_data = np.array([0.0, 1.0, -2.0, -4.0])
    sample_weights = np.array([1.0, 1.0, 0.0, 1.0])
    train = np.array([0, 1, 2])
    test = np.array([3])
    scorer = make_scorer(mean_absolute_error)

    # No weights
    scores = _rfe_single_fit(advanced_rfe, advanced_rfe.estimator, x_data,
                             y_data, train, test, scorer)
    assert advanced_rfe.n_features_ == 1
    np.testing.assert_array_equal(advanced_rfe.ranking_, [2, 3, 1])
    np.testing.assert_array_equal(advanced_rfe.support_,
                                  [False, False, True])
    est = advanced_rfe.estimator_
    assert isinstance(est, LinearRegression)
    np.testing.assert_allclose(est.coef_, [1.0])
    np.testing.assert_allclose(est.intercept_, 0.0, atol=1e-10)
    np.testing.assert_allclose(scores, [0.0, 0.0, 0.0], atol=1e-10)

    # With weights
    scores = _rfe_single_fit(advanced_rfe, advanced_rfe.estimator, x_data,
                             y_data, train, test, scorer,
                             sample_weight=sample_weights)
    assert advanced_rfe.n_features_ == 1
    np.testing.assert_array_equal(advanced_rfe.ranking_, [1, 3, 2])
    np.testing.assert_array_equal(advanced_rfe.support_,
                                  [True, False, False])
    est = advanced_rfe.estimator_
    assert isinstance(est, LinearRegression)
    np.testing.assert_allclose(est.coef_, [0.5])
    np.testing.assert_allclose(est.intercept_, 0.0, atol=1e-10)
    np.testing.assert_allclose(scores, [4.8, 4.8, 6.0])


# _map_features


TEST_MAP_FEATURES = [
    ([True, True, True], [0, 1, 2]),
    ([1, 1, 1], [0, 1, 2]),
    ([True, True, False], [0, 1]),
    ([1, 1, 0], [0, 1]),
    ([True, False, True], [0, 1]),
    ([1, 0, 1], [0, 1]),
    ([False, True, True], [0, 1]),
    ([0, 1, 1], [0, 1]),
    ([True, False, False], [0]),
    ([1, 0, 0], [0]),
    ([False, True, False], [0]),
    ([0, 1, 0], [0]),
    ([False, False, True], [0]),
    ([0, 0, 1], [0]),
    ([False, False, False], []),
    ([0, 0, 0], []),
]

FEATURES = [0, 1, 2]


@pytest.mark.parametrize('support,output', TEST_MAP_FEATURES)
def test_map_features(support, output):
    """Test ``_map_features``."""
    new_features = _map_features(FEATURES, support)
    np.testing.assert_array_equal(new_features, output)


# _update_transformers_param


class NoPipeline(BaseEstimator):
    """No pipeline."""

    def __init__(self, test_transformers=1):
        """Initialize instance."""
        self.test_transformers = test_transformers

    def fit(self, *_):
        """Fit method."""
        return self


def test_update_transformers_fail_no_pipeline():
    """Test ``_update_transformers_param`` expected fail."""
    msg = 'estimator is not a Pipeline or AdvancedPipeline'
    est = NoPipeline()
    with pytest.raises(TypeError, match=msg):
        _update_transformers_param(est, [])


def test_update_transformers_fail_no_column_transformer():
    """Test ``_update_transformers_param`` expected fail."""
    msg = 'pipeline step is not a ColumnTransformer'
    est = AdvancedPipeline([('no_pipeline', NoPipeline())])
    with pytest.raises(TypeError, match=msg):
        _update_transformers_param(est, [])


def test_update_transformers_param_lin():
    """Test ``_update_transformers_param``."""
    est = LinearRegression()
    params = deepcopy(est.get_params())
    _update_transformers_param(est, [])
    assert est.get_params() == params


def test_update_transformers_column_transformer():
    """Test ``_update_transformers_param``."""
    trans1 = [('1', 'drop', [0, 1, 2])]
    trans2 = [('2a', 'drop', [0]), ('2b', 'passthrough', [1, 2])]
    est = AdvancedPipeline([
        ('trans1', ColumnTransformer(trans1)),
        ('trans2', ColumnTransformer(trans2)),
        ('passthrough', 'passthrough'),
    ])
    params = deepcopy(est.get_params())
    _update_transformers_param(est, np.array([True, False, True]))
    new_params = est.get_params()
    assert new_params != params
    assert new_params['trans1__transformers'] == [('1', 'drop', [0, 1])]
    assert new_params['trans2__transformers'] == [
        ('2a', 'drop', [0]),
        ('2b', 'passthrough', [1]),
    ]


# cross_val_score_weighted


def test_cross_val_score_weighted():
    """Test ``cross_val_score_weighted``."""
    sample_weights = np.array([1.0, 1.0, 0.0, 1.0, 1.0, 0.0])
    cv_score_kwargs = {
        'estimator': LinearRegression(),
        'x_data': np.arange(6).reshape(6, 1),
        'y_data': np.array([0, 1, 1000, 0, -1, -1000]),
        'groups': ['A', 'A', 'A', 'B', 'B', 'B'],
        'scoring': 'neg_mean_absolute_error',
        'cv': LeaveOneGroupOut(),
        'fit_params': {'sample_weight': sample_weights},
        'sample_weights': sample_weights,
    }
    scores = cross_val_score_weighted(**cv_score_kwargs)
    np.testing.assert_allclose(scores, [-2.0, -4.0])


# get_rfecv_transformer


def test_get_rfecv_transformer_not_fitted():
    """Test ``get_rfecv_transformer`` expected fail."""
    rfecv = AdvancedRFECV(LinearRegression())
    msg = ('RFECV instance used to initialize FeatureSelectionTransformer '
           'must be fitted')
    with pytest.raises(NotFittedError, match=msg):
        get_rfecv_transformer(rfecv)


def test_get_rfecv_transformer():
    """Test ``get_rfecv_transformer``."""
    rfecv = AdvancedRFECV(LinearRegression(), cv=2)
    x_data = np.arange(30).reshape(10, 3)
    y_data = np.arange(10)
    rfecv.fit(x_data, y_data)
    transformer = get_rfecv_transformer(rfecv)
    assert isinstance(transformer, FeatureSelectionTransformer)
    assert rfecv.n_features_ == transformer.n_features
    np.testing.assert_allclose(rfecv.grid_scores_, transformer.grid_scores)
    np.testing.assert_allclose(rfecv.ranking_, transformer.ranking)
    np.testing.assert_allclose(rfecv.support_, transformer.support)
    assert len(transformer.grid_scores) == 3
    assert len(transformer.ranking) == 3
    assert len(transformer.support) == 3


# perform_efecv


def test_perform_efecv():
    """Test ``perform_efecv``."""
    x_data = np.array([
        [0, 0, 0],
        [1, 1, 0],
        [2, 0, 2],
        [0, 3, 3],
        [4, 4, 4],
        [4, 4, 0],
    ])
    y_data = np.array([1, 0, 3, -5, -3, -3])

    (best_est, transformer) = perform_efecv(LinearRegression(), x_data, y_data,
                                            cv=2)

    assert isinstance(best_est, LinearRegression)
    np.testing.assert_allclose(best_est.coef_, [1.0, -2.0])
    np.testing.assert_allclose(best_est.intercept_, 1.0)

    assert isinstance(transformer, FeatureSelectionTransformer)
    assert transformer.n_features == 2
    assert len(transformer.grid_scores) == 7
    np.testing.assert_array_equal(transformer.ranking, [1, 1, 2])
    np.testing.assert_array_equal(transformer.support, [True, True, False])
