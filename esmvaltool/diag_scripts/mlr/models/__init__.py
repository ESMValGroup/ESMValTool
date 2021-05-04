"""Base class for MLR models.

Example recipe
--------------
The :ref:`MLR main diagnostic script<api.esmvaltool.diag_scripts.mlr.main>`
provides an interface for using MLR models in recipes. The following recipe
shows a typical example on how to setup MLR recipes/diagnostics with the
following properties:

#. Setup an MLR model with target variable ``y`` (using the tag ``Y``) and
   three predictors ``x1``, ``x2`` and ``latitude`` (with tags ``X1``, ``X2``
   and ``latitude``, respectively). The target variable needs the attribute
   ``var_type: label``; the predictors ``x1`` and ``x2`` the attribute
   ``var_type: feature``.  The coordinate feature ``latitude`` is added via the
   option ``coords_as_features: [latitude]``.
#. Suppose ``y`` and ``x1`` are 3D fields (pressure, latitude, longitude);
   ``x2`` is a 2D field (latitude, longitude). Thus, it is necessary to add the
   attribute ``broadcast_from: [1, 2]`` to it (see ``dim_map`` parameter in
   :func:`iris.util.broadcast_to_shape` for details).  In order to consider
   multiple climate models (``A``, ``B`` and ``C``) at once, the option
   ``group_datasets_by_attributes: [dataset]`` is necessary.  Otherwise the
   diagnostic will complain about duplicate data.
#. For the prediction, data from dataset ``D`` is used (with
   ``var_type: prediction_input``). For the feature ``X1`` additional input
   error (with ``var_type: prediction_input_error``) is used.

   .. code-block:: yaml

      diag_feature_x1:
        variables:
          feature:
            ... # specify project, mip, start_year, end_year, etc.
            short_name: x1
            var_type: feature
            tag: X1
            additional_datasets:
              - {dataset: A, ...}
              - {dataset: B, ...}
              - {dataset: C, ...}
          prediction_input:
            ... # specify project, mip, start_year, end_year, etc.
            short_name: x1
            var_type: prediction_input
            tag: X1
            additional_datasets:
              - {dataset: D, ...}
          prediction_input_error:
            ... # specify project, mip, start_year, end_year, etc.
            short_name: x1Stderr
            var_type: prediction_input_error
            tag: X1
            additional_datasets:
              - {dataset: D, ...}
        scripts:
          null

      diag_feature_x2:
        variables:
          feature:
            ... # specify project, mip, start_year, end_year, etc.
            short_name: x2
            var_type: feature
            broadcast_from: [1, 2]
            tag: X2
            additional_datasets:
              - {dataset: A, ...}
              - {dataset: B, ...}
              - {dataset: C, ...}
          prediction_input:
            ... # specify project, mip, start_year, end_year, etc.
            short_name: x2
            var_type: prediction_input
            broadcast_from: [1, 2]
            tag: X2
            additional_datasets:
              - {dataset: D, ...}
        scripts:
          null

      diag_label:
        variables:
          label:
            ... # specify project, mip, start_year, end_year, etc.
            short_name: y
            var_type: label
            tag: Y
            additional_datasets:
              - {dataset: A, ...}
              - {dataset: B, ...}
              - {dataset: C, ...}
        scripts:
          null

#. In this example, a
   `GBRT model
   <https://scikit-learn.org/stable/modules/ensemble.html
   #gradient-tree-boosting>`_ (with ``mlr_model_type: gbr_sklearn``) is used.
   Parameters for this are specified via ``parameters_final_regressor``. Apart
   from the best-estimate prediction, the estimated MLR model error
   (``save_mlr_model_error: test``) and the propagated prediction input error
   (``save_propagated_errors: true``) are returned.
#. With ``postprocess.py``, the global mean of the best estimate prediction and
   the corresponding errors (MLR model + propagated input error) are calculted.

   .. code-block:: yaml

      diag_mlr_gbrt:
        scripts:
          mlr:
            script: mlr/main.py
            ancestors: [
               'diag_label/y',
               'diag_feature_*/*',
            ]
            coords_as_features: [latitude]
            group_datasets_by_attributes: [dataset]
            mlr_model_name: GBRT
            mlr_model_type: gbr_sklearn
            parameters_final_regressor:
              learning_rate: 0.1
              n_estimators: 100
            save_mlr_model_error: test
            save_propagated_errors: true
          postprocess:
            script: mlr/postprocess.py
            ancestors: ['diag_mlr_gbrt/mlr']
            ignore:
              - {var_type: null}
            mean: [pressure, latitude, longitude]

#. Plots of the global distribution (latitude, longitude) are created with
   ``plot.py`` after calculating the mean over the pressure coordinate using
   ``preprocess.py``.

   .. code-block:: yaml

      diag_plot:
        scripts:
          preprocess:
            script: mlr/preprocess.py
            ancestors: ['diag_mlr_gbrt/mlr']
            collapse: [pressure]
            ignore:
              - {var_type: null}
          plot:
            script: mlr/plot.py
            ancestors: ['diag_plot/preprocess']
            plot_map:
               plot_kwargs:
                 cbar_label: 'Y'
                 cbar_ticks: [0, 1, 2, 3]
                 vmin: 0
                 vmax: 3

All datasets must have the attribute ``var_type`` which specifies the type of
the dataset.  Possible values are ``feature`` (independent variables used for
training/testing), ``label`` (dependent variables, y-axis),
``prediction_input`` (independent variables used for prediction of dependent
variables, usually observational data), ``prediction_input_error`` (standard
error of the ``prediction_input`` data, optional) or ``prediction_reference``
(`true` values for the ``prediction_input`` data, optional). In addition, all
datasets must habe the attribute ``tag``, which specifies the name of
variable/diagnostic. All datasets can be converted to new units in the loading
step by specifying the key ``convert_units_to`` in the respective dataset(s).

Training data
-------------
All groups (specified in ``group_datasets_by_attributes``, if desired) given
for ``label`` datasets must also be given for the ``feature`` datasets. Within
these groups, all ``feature`` and ``label`` datasets must have the same shape,
except the attribute ``broadcast_from`` is set to a list of suitable coordinate
indices to map this dataset to regular datasets (see parameter ``dim_map`` in
:func:`iris.util.broadcast_to_shape`).

Prediction data
---------------
All ``tag`` s specified for ``prediction_input`` datasets must also be given
for the ``feature`` datasets (except ``allow_missing_features`` is set to
``True``).  Multiple predictions can be specified by ``prediction_name``.
Within these predictions, all ``prediction_input`` datasets must have the same
shape, except the attribute ``broadcast_from`` is given. Errors in the
prediction input data can be specified by ``prediction_input_error``. If given,
these errors are used to calculate errors in the final prediction using linear
error propagation given by `LIME <https://arxiv.org/abs/1602.04938>`_.
Additionally, `true` values for ``prediction_input`` can be specified with
``prediction_reference`` datasets (together with the respective
``prediction_name``). This allows an evaluation of the performance of the MLR
model by calculating residuals (`true` minus predicted values).

Available MLR models
--------------------
MLR models are subclasses of this base class. A list of all available MLR
models can be found :ref:`here <availableMLRModels>`. To add a new MLR model,
create a new file in ``esmvaltool/diag_scripts/mlr/models/`` with a child class
of :class:`esmvaltool.diag_scripts.mlr.models.MLRModel` decorated with
:meth:`esmvaltool.diag_scripts.mlr.models.MLRModel.register_mlr_model`.

.. _MLRModeloptionalparameters:

Optional parameters for class initialization
--------------------------------------------
accept_only_scalar_data: bool (default: False)
    If set to ``True``, only accept scalar input data. Should be used together
    with the option ``group_datasets_by_attributes``.
allow_missing_features: bool (default: False)
    Allow missing features in the training data.
cache_intermediate_results: bool (default: True)
    Cache the intermediate results of the pipeline's transformers.
categorical_features: list of str
    Names of features which are interpreted as categorical features (in
    contrast to numerical features).
coords_as_features: list of str
    If given, specify a list of coordinates which should be used as features.
dtype: str (default: 'float64')
    Internal data type which is used for all calculations, see
    `<https://docs.scipy.org/doc/numpy/user/basics.types.html>`_ for a list of
    allowed values.
fit_kwargs: dict
    Optional keyword arguments for the pipeline's ``fit()`` function.  These
    arguments have to be given for each step of the pipeline separated by two
    underscores, i.e. ``s__p`` is the parameter ``p`` for step ``s``.
group_datasets_by_attributes: list of str
    List of dataset attributes which are used to group input data for
    ``feature`` s and ``label`` s. For example, this is necessary if the MLR
    model should consider multiple climate models in the training phase. If
    this option is not given, specifying multiple datasets with identical
    ``var_type`` and ``tag`` entries results in an error. If given, all the
    input data is first grouped by the given attributes and then checked for
    uniqueness within this group. After that, all groups are stacked to form a
    single set of training data.
imputation_strategy: str (default: 'remove')
    Strategy for the imputation of missing values in the features. Must be one
    of ``'remove'``, ``'mean'``, ``'median'``, ``'most_frequent'`` or
    ``'constant'``.
log_level: str (default: 'info')
    Verbosity for the logger. Must be one of ``'debug'``, ``'info'``,
    ``'warning'`` or ``'error'``.
mlr_model_name: str
    Human-readable name of the MLR model instance (e.g used for labels).
n_jobs: int (default: 1)
    Maximum number of jobs spawned by this class. Use ``-1`` to use all
    processors. More details are given `here
    <https://scikit-learn.org/stable/glossary.html#term-n-jobs>`_.
output_file_type: str (default: 'png')
    File type for the plots.
parameters: dict
    Parameters used for the whole pipeline. Have to be given for each step of
    the pipeline separated by two underscores, i.e. ``s__p`` is the parameter
    ``p`` for step ``s``.
parameters_final_regressor: dict
    Parameters used for the **final** regressor. If these parameters are
    updated using the function :meth:`update_parameters`, the new names have to
    be given for each step of the pipeline separated by two underscores, i.e.
    ``s__p`` is the parameter ``p`` for step ``s``.
pca: bool (default: False)
    Preprocess numerical input features using PCA. Parameters for this pipeline
    step can be given via the ``parameters`` argument.
plot_dir: str (default: ~/plots)
    Root directory to save plots.
plot_units: dict
    Replace specific units (keys) with other text (values) in plots.
savefig_kwargs: dict
    Keyword arguments for :func:`matplotlib.pyplot.savefig`.
seaborn_settings: dict
    Options for :func:`seaborn.set` (affects all plots).
standardize_data: bool (default: True)
    Linearly standardize numerical input data by removing mean and scaling to
    unit variance.
sub_dir: str
    Create additional subdirectory for output in ``work_dir`` and ``plot_dir``.
test_size: float (default: 0.25)
    If given, randomly exclude the desired fraction of input data from training
    and use it as test data.
weighted_samples: dict
    If specified, use weighted samples whenever possible. The given keyword
    arguments are directly passed to
    :func:`esmvaltool.diag_scripts.mlr.get_all_weights` to calculate the sample
    weights. By default, area weights and time weights are used.
work_dir: str (default: ~/work)
    Root directory to save all other files (mainly ``*.nc`` files).

"""

import importlib
import logging
import os
from copy import deepcopy
from inspect import getfullargspec
from pprint import pformat

import iris
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from cf_units import Unit
from joblib import Parallel, delayed
from lime.lime_tabular import LimeTabularExplainer
from matplotlib.ticker import ScalarFormatter
from scipy.stats import shapiro
from sklearn import metrics
from sklearn.compose import ColumnTransformer
from sklearn.decomposition import PCA
from sklearn.exceptions import NotFittedError
from sklearn.impute import SimpleImputer
from sklearn.inspection import plot_partial_dependence
from sklearn.model_selection import (
    GridSearchCV,
    LeaveOneGroupOut,
    LeaveOneOut,
    train_test_split,
)
from sklearn.preprocessing import StandardScaler

from esmvaltool.diag_scripts import mlr
from esmvaltool.diag_scripts.mlr.custom_sklearn import (
    AdvancedPipeline,
    AdvancedRFECV,
    AdvancedTransformedTargetRegressor,
    cross_val_score_weighted,
    get_rfecv_transformer,
    perform_efecv,
)
from esmvaltool.diag_scripts.shared import (
    ProvenanceLogger,
    group_metadata,
    io,
    select_metadata,
)

logger = logging.getLogger(os.path.basename(__file__))


class MLRModel():
    """Base class for MLR models."""

    _CLF_TYPE = None
    _MODELS = {}
    _MLR_MODEL_TYPE = None

    @staticmethod
    def _load_mlr_models():
        """Load MLR models from :mod:`esmvaltool.diag_scripts.mlr.models`."""
        current_path = os.path.dirname(os.path.realpath(__file__))
        models_path = os.path.join(current_path)
        for (root, _, model_files) in os.walk(models_path):
            for model_file in model_files:
                rel_path = ('' if root == models_path else os.path.relpath(
                    root, models_path))
                module = os.path.join(rel_path,
                                      os.path.splitext(model_file)[0])
                try:
                    importlib.import_module(
                        'esmvaltool.diag_scripts.mlr.models.{}'.format(
                            module.replace(os.sep, '.')))
                except ImportError:
                    pass

    @classmethod
    def register_mlr_model(cls, mlr_model_type):
        """Add MLR model (subclass of this class) (decorator)."""
        logger.debug("Found available MLR model '%s'", mlr_model_type)

        def decorator(subclass):
            """Decorate subclass."""
            subclass._MLR_MODEL_TYPE = mlr_model_type
            cls._MODELS[mlr_model_type] = subclass
            return subclass

        return decorator

    @classmethod
    def create(cls, mlr_model_type, *args, **kwargs):
        """Create desired MLR model subclass (factory method)."""
        cls._load_mlr_models()
        if not cls._MODELS:
            raise NotImplementedError(
                f"Cannot initialize new MLR model with type "
                f"'{mlr_model_type}', no MLR models found. Please add "
                f"subclasses of {cls} in new files under 'esmvaltool/"
                f"diag_scripts/mlr/models/' decorated by 'esmvaltool."
                f"diag_scripts.mlr.models.{cls.__name__}."
                f"register_mlr_model()'")
        if mlr_model_type not in cls._MODELS:
            raise NotImplementedError(
                f"MLR model type '{mlr_model_type}' not found in 'esmvaltool/"
                f"diag_scripts/mlr/models/'")
        subclass = cls._MODELS[mlr_model_type]
        logger.info(
            "Initialized MLR model with type '%s' and final regressor %s",
            mlr_model_type, subclass._CLF_TYPE)
        return subclass(*args, **kwargs)

    def __init__(self, input_datasets, **kwargs):
        """Initialize class members.

        Parameters
        ----------
        input_datasets : list of dict
            List of dataset metadata used as data for the MLR model.
        **kwargs
            Optional keyword arguments, see next sections.

        Raises
        ------
        NotImplementedError
            Class is initialized directly without the use of its factory
            function ``create()``.
        ValueError
            Invalid data given.

        """
        self._check_clf()

        # Private attributes
        self._cfg = deepcopy(kwargs)
        self._clf = None
        self._lime_explainer = None
        self._data = {}
        self._data['pred'] = {}
        self._datasets = {}
        self._classes = {}
        self._parameters = {}

        # Set default settings
        self._set_default_settings()

        # Seaborn
        sns.set(**self._cfg.get('seaborn_settings', {}))

        # Adapt output directories
        self._cfg['mlr_work_dir'] = os.path.join(self._cfg['work_dir'],
                                                 self._cfg['sub_dir'])
        self._cfg['mlr_plot_dir'] = os.path.join(self._cfg['plot_dir'],
                                                 self._cfg['sub_dir'])
        if not os.path.exists(self._cfg['mlr_work_dir']):
            os.makedirs(self._cfg['mlr_work_dir'])
            logger.info("Created %s", self._cfg['mlr_work_dir'])
        if not os.path.exists(self._cfg['mlr_plot_dir']):
            os.makedirs(self._cfg['mlr_plot_dir'])
            logger.info("Created %s", self._cfg['mlr_plot_dir'])

        # Load datasets, classes and training data
        self._load_input_datasets(input_datasets)
        self._load_classes()
        self._load_data()

        # Create pipeline (with all preprocessor steps and final regressor)
        self.reset_pipeline()
        if self._cfg['parameters']:
            logger.debug("Using parameter(s): %s", self._cfg['parameters'])
        self.update_parameters(**self._cfg['parameters'])

        # Log successful initialization
        logger.info("Initialized MLR model (using at most %i processes)",
                    self._cfg['n_jobs'])
        logger.debug("With parameters")
        logger.debug(pformat(self.parameters))

    @property
    def categorical_features(self):
        """numpy.ndarray: Categorical features."""
        return self.features[self._classes['features'].categorical]

    @property
    def data(self):
        """dict: Input data of the MLR model."""
        return self._data

    @property
    def features(self):
        """numpy.ndarray: Features of the input data."""
        return self._classes['features'].index.values

    @property
    def features_after_preprocessing(self):
        """numpy.ndarray: Features of the input data after preprocessing."""
        x_train = self.get_x_array('train')
        y_train = self.get_y_array('train')
        try:
            self._check_fit_status('Calculating features after preprocessing')
        except NotFittedError:
            self._clf.fit_transformers_only(x_train, y_train,
                                            **self.fit_kwargs)
        x_trans = self._clf.transform_only(x_train)
        features = self.features
        n_features_may_drop = False
        if 'feature_selection' in self._clf.named_steps:
            support = self._clf.named_steps['feature_selection'].support
            features = features[support]
            n_features_may_drop = True
        if 'pca' in self._clf.named_steps:
            categorical_features = np.array([
                f for f in features if f in self.categorical_features])
            n_numerical_features = x_trans.shape[1] - categorical_features.size
            features = [
                f'Principal component {idx}'
                for idx in range(n_numerical_features)
            ]
            features.extend(categorical_features)
            n_features_may_drop = True
        if not n_features_may_drop and x_trans.shape[1] != self.features.size:
            logger.warning(
                "Number of features decreased from %i to %i during "
                "preprocessing for unknown reasons (neither feature selection "
                "using recursive feature elimination nor PCA is performed)",
                self.features.size, x_trans.shape[1])
            features = [
                f'Unknown feature {idx}' for idx in range(x_trans.shape[1])
            ]
        return np.array(features, dtype='str')

    @property
    def features_types(self):
        """pandas.Series: Types of the features."""
        return self._classes['features'].types

    @property
    def features_units(self):
        """pandas.Series: Units of the features."""
        return self._classes['features'].units

    @property
    def fit_kwargs(self):
        """dict: Keyword arguments for :meth:`fit`."""
        fit_kwargs = self._cfg['fit_kwargs']
        fit_kwargs = self._update_fit_kwargs(fit_kwargs)
        verbosity_kwargs = self._get_verbosity_parameters(self._clf.fit)
        for (key, val) in verbosity_kwargs.items():
            fit_kwargs.setdefault(key, val)
        return fit_kwargs

    @property
    def group_attributes(self):
        """numpy.ndarray: Group attributes of the input data."""
        return self._classes['group_attributes']

    @property
    def label(self):
        """str: Label of the input data."""
        return self._classes['label'].index.values[0]

    @property
    def label_units(self):
        """str: Units of the label."""
        return self._classes['label'].units.values[0]

    @property
    def mlr_model_type(self):
        """str: MLR model type."""
        return self._MLR_MODEL_TYPE

    @property
    def numerical_features(self):
        """numpy.ndarray: Numerical features."""
        return self.features[~self._classes['features'].categorical]

    @property
    def parameters(self):
        """dict: Parameters of the complete MLR model pipeline."""
        return self._parameters

    def efecv(self, **kwargs):
        """Perform exhaustive feature elimination using cross-validation.

        Parameters
        ----------
        **kwargs : keyword arguments, optional
            Additional options for :func:`esmvaltool.diag_scripts.mlr.
            custom_sklearn.cross_val_score_weighted`.

        """
        logger.info(
            "Performing exhaustive feature elimination using cross-validation "
            "with final regressor %s on %i training points (thiy may take a "
            "while...)", self._CLF_TYPE,
            len(self.data['train'].index))

        # Get fit parameters
        fit_kwargs = deepcopy(self.fit_kwargs)
        keys_to_remove = []
        for key in fit_kwargs:
            if key.endswith('eval_set'):
                keys_to_remove.append(key)
        for key in keys_to_remove:
            logger.warning(
                "Fit parameter '%s' is not supported for efecv()", key)
            fit_kwargs.pop(key)

        # Get other keyword arguments
        kwargs = deepcopy(kwargs)
        verbosity_kwargs = self._get_verbosity_parameters(
            cross_val_score_weighted)
        for (key, val) in verbosity_kwargs.items():
            kwargs.setdefault(key, val)
        kwargs.setdefault('n_jobs', self._cfg['n_jobs'])
        kwargs['fit_params'] = fit_kwargs
        kwargs['sample_weights'] = self._get_sample_weights('train')
        if kwargs.get('cv') == 'logo':
            kwargs.update(self._get_logo_cv_kwargs())

        # Exhaustive feature selection
        (self._clf, transformer) = perform_efecv(
            self._clf, self.data['train'].x, self.get_y_array('train'),
            **kwargs)
        self._clf.steps.insert(0, ('feature_selection', transformer))

        # Log results
        new_features = self.features[transformer.support]
        logger.info(
            "Exhaustive feature elimination was successful, %i of the %i "
            "features remain", new_features.size, self.features.size)
        logger.info("Old features: %s", self.features)
        logger.info("New features: %s", new_features)
        logger.info("Successfully fitted MLR model on %i training point(s)",
                    len(self.data['train'].index))
        logger.debug("Pipeline steps:")
        logger.debug(pformat(list(self._clf.named_steps.keys())))
        logger.debug("Parameters:")
        logger.debug(pformat(self.parameters))

        # LIME
        self._load_lime_explainer()

    def export_prediction_data(self, filename=None):
        """Export all prediction data contained in `self._data`.

        Parameters
        ----------
        filename : str, optional (default: '{data_type}_{pred_name}.csv')
            Name of the exported files.

        """
        for pred_name in self.data['pred']:
            self._save_csv_file('pred', filename, pred_name=pred_name)

    def export_training_data(self, filename=None):
        """Export all training data contained in `self._data`.

        Parameters
        ----------
        filename : str, optional (default: '{data_type}.csv')
            Name of the exported files.

        """
        for data_type in ('all', 'train', 'test'):
            self._save_csv_file(data_type, filename)

    def fit(self):
        """Fit MLR model.

        Note
        ----
        Specifying keyword arguments for this function is not allowed here
        since :attr:`features_after_preprocessing` might be altered by
        that. Use the keyword argument ``fit_kwargs`` during class
        initialization instead.

        """
        logger.info(
            "Fitting MLR model with final regressor %s on %i training "
            "point(s)", self._CLF_TYPE, len(self.data['train'].index))

        # Create MLR model with desired parameters and fit it
        self._clf.fit(self.data['train'].x, self.data['train'].y,
                      **self.fit_kwargs)
        self._parameters = self._get_clf_parameters()
        logger.info("Successfully fitted MLR model on %i training point(s)",
                    len(self.data['train'].index))
        logger.debug("Pipeline steps:")
        logger.debug(pformat(list(self._clf.named_steps.keys())))
        logger.debug("Parameters:")
        logger.debug(pformat(self.parameters))

        # LIME
        self._load_lime_explainer()

    def get_ancestors(self, label=True, features=None, prediction_names=None,
                      prediction_reference=False):
        """Return ancestor files.

        Parameters
        ----------
        label : bool, optional (default: True)
            Return ``label`` files.
        features : list of str, optional (default: None)
            Features for which files should be returned. If ``None``, return
            files for all features.
        prediction_names : list of str, optional (default: None)
            Prediction names for which files should be returned. If ``None``,
            return files for all prediction names.
        prediction_reference : bool, optional (default: False)
            Return ``prediction_reference`` files if available for given
            ``prediction_names``.

        Returns
        -------
        list of str
            Ancestor files.

        Raises
        ------
        ValueError
            Invalid ``feature`` or ``prediction_name`` given.

        """
        ancestors = []

        # Label files
        if label:
            ancestors.extend([d['filename'] for d in self._datasets['label']])

        # Feature files
        if features is None:
            features = self.features
        for feature in features:
            if feature not in self.features:
                raise ValueError(
                    f"Got invalid feature '{feature}', expected one of "
                    f"{self.features}")
            ancestors.extend(
                [d['filename'] for d in self._datasets['feature']
                 if d['tag'] == feature]
            )

        # Prediction files
        available_pred_names = list(self._datasets['prediction_input'].keys())
        if prediction_names is None:
            prediction_names = available_pred_names
        for pred_name in prediction_names:
            if pred_name not in available_pred_names:
                raise ValueError(
                    f"Got invalid prediction name '{pred_name}', expected one "
                    f"of {available_pred_names}")
            ancestors.extend(
                [d['filename'] for d in
                 self._datasets['prediction_input'][pred_name]]
            )
            ancestors.extend(
                [d['filename'] for d in
                 self._datasets['prediction_input_error'].get(pred_name, [])]
            )
            if prediction_reference:
                ancestors.extend(
                    [d['filename'] for d in
                     self._datasets['prediction_reference'].get(pred_name, [])]
                )

        return ancestors

    def get_data_frame(self, data_type, impute_nans=False):
        """Return data frame of specified type.

        Parameters
        ----------
        data_type : str
            Data type to be returned. Must be one of ``'all'``, ``'train'`` or
            ``'test'``.
        impute_nans : bool, optional (default: False)
            Impute nans if desired.

        Returns
        -------
        pandas.DataFrame
            Desired data.

        Raises
        ------
        TypeError
            ``data_type`` is invalid or data does not exist (e.g. test data is
            not set).

        """
        allowed_types = ('all', 'train', 'test')
        if data_type not in allowed_types:
            raise TypeError(
                f"'{data_type}' is not an allowed type, specify one of "
                f"'{allowed_types}'")
        if data_type not in self.data:
            raise TypeError(f"No '{data_type}' data available")
        data_frame = self.data[data_type]
        if impute_nans:
            data_frame = self._impute_nans(data_frame)
        return data_frame

    def get_x_array(self, data_type, impute_nans=False):
        """Return x data of specific type.

        Parameters
        ----------
        data_type : str
            Data type to be returned. Must be one of ``'all'``, ``'train'`` or
            ``'test'``.
        impute_nans : bool, optional (default: False)
            Impute nans if desired.

        Returns
        -------
        numpy.ndarray
            Desired data.

        Raises
        ------
        TypeError
            ``data_type`` is invalid or data does not exist (e.g. test data is
            not set).

        """
        data_frame = self.get_data_frame(data_type, impute_nans=impute_nans)
        return data_frame.x.values

    def get_y_array(self, data_type, impute_nans=False):
        """Return y data of specific type.

        Parameters
        ----------
        data_type : str
            Data type to be returned. Must be one of ``'all'``, ``'train'`` or
            ``'test'``.
        impute_nans : bool, optional (default: False)
            Impute nans if desired.

        Returns
        -------
        numpy.ndarray
            Desired data.

        Raises
        ------
        TypeError
            ``data_type`` is invalid or data does not exist (e.g. test data is
            not set).

        """
        data_frame = self.get_data_frame(data_type, impute_nans=impute_nans)
        return data_frame.y.squeeze().values

    def grid_search_cv(self, param_grid, **kwargs):
        """Perform exhaustive parameter search using cross-validation.

        Parameters
        ----------
        param_grid : dict or list of dict
            Parameter names (keys) and ranges (values) for the search. Have to
            be given for each step of the pipeline separated by two
            underscores, i.e. ``s__p`` is the parameter ``p`` for step ``s``.
        **kwargs : keyword arguments, optional
            Additional options for
            :class:`sklearn.model_selection.GridSearchCV`.

        Raises
        ------
        ValueError
            Final regressor does not supply the attributes ``best_estimator_``
            or ``best_params_``.

        """
        logger.info(
            "Performing exhaustive grid search cross-validation with final "
            "regressor %s and parameter grid %s on %i training points",
            self._CLF_TYPE, param_grid, len(self.data['train'].index))

        # Get keyword arguments
        (cv_kwargs, fit_kwargs) = self._get_cv_estimator_kwargs(GridSearchCV,
                                                                **kwargs)

        # Create and fit GridSearchCV instance
        clf = GridSearchCV(self._clf, param_grid, **cv_kwargs)
        clf.fit(self.data['train'].x, self.data['train'].y, **fit_kwargs)

        # Try to find best estimator
        if hasattr(clf, 'best_estimator_'):
            self._clf = clf.best_estimator_
        elif hasattr(clf, 'best_params_'):
            self.update_parameters(**clf.best_params_)
            self._clf.fit(self.data['train'].x, self.data['train'].y,
                          **fit_kwargs)
        else:
            raise ValueError(
                "GridSearchCV not successful, cannot determine best estimator "
                "(neither using 'best_estimator_' nor 'best_params_'), "
                "adapt keyword arguments accordingly (see "
                "https://scikit-learn.org/stable/modules/generated/"
                "sklearn.model_selection.GridSearchCV.html for more help)")
        self._parameters = self._get_clf_parameters()
        logger.info(
            "Exhaustive grid search successful, found best parameter(s) %s",
            clf.best_params_)
        logger.debug("CV results:")
        logger.debug(pformat(clf.cv_results_))
        logger.info("Successfully fitted MLR model on %i training point(s)",
                    len(self.data['train'].index))
        logger.debug("Pipeline steps:")
        logger.debug(pformat(list(self._clf.named_steps.keys())))
        logger.debug("Parameters:")
        logger.debug(pformat(self.parameters))

        # LIME
        self._load_lime_explainer()

    def plot_1d_model(self, filename=None, n_points=1000):
        """Plot lineplot that represents the MLR model.

        Note
        ----
        This only works for a model with a single feature.

        Parameters
        ----------
        filename : str, optional (default: '1d_mlr_model')
            Name of the plot file.
        n_points : int, optional (default: 1000)
            Number of sampled points for the single feature (using linear
            spacing between minimum and maximum value).

        Raises
        ------
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.
        ValueError
            MLR model is built from more than 1 feature.

        """
        if not self._is_ready_for_plotting():
            return
        n_features = self.features.size
        if n_features > 1:
            raise ValueError(
                f"Plotting lineplot of MLR model using 'plot_1d_model' is not "
                f"possible, MLR model {self._cfg['mlr_model_name']} contains "
                f"more than one feature ({n_features:d} features: "
                f"{self.features})")
        feature = self.features[0]
        logger.info("Plotting 1D MLR model (sampling %i points for single "
                    "feature '%s')", n_points, feature)
        if filename is None:
            filename = '1d_mlr_model'
        (_, axes) = plt.subplots()

        # Get available datasets
        data_to_plot = ['train']
        if 'test' in self.data:
            data_to_plot.append('test')

        # Plot training and test data (if available)
        for data_type in data_to_plot:
            x_data = self.data[data_type].x[feature].values
            y_data = self.get_y_array(data_type)
            axes.scatter(
                x_data, y_data,
                **self._get_plot_kwargs(data_type, plot_type='scatter'))

        # Plot MLR model
        x_lin = np.linspace(self.data['all'].x[feature].values.min(),
                            self.data['all'].x[feature].values.max(),
                            n_points)
        y_pred = self._clf.predict(x_lin.reshape(-1, 1))
        axes.plot(x_lin, y_pred, color='k', linewidth=2,
                  label=self._cfg['mlr_model_name'])

        # Plot appearance
        title = (f"Predicted {self.label} by MLR model "
                 f"{self._cfg['mlr_model_name']}")
        axes.set_title(title)
        axes.set_xlabel(self._get_plot_feature(feature))
        axes.set_ylabel(self._get_plot_label())
        axes.legend(loc='best')

        # Save plot
        plot_path = os.path.join(
            self._cfg['mlr_plot_dir'],
            filename + '.' + self._cfg['output_file_type'],
        )
        plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save provenance
        cube = mlr.get_1d_cube(
            x_lin,
            y_pred,
            x_kwargs={'var_name': feature,
                      'long_name': feature,
                      'units': self.features_units[feature]},
            y_kwargs={'var_name': self.label,
                      'long_name': title,
                      'units': self.label_units,
                      'attributes': {'project': '', 'dataset': ''}},
        )
        self._write_plot_provenance(
            cube, plot_path, ancestors=self.get_ancestors(prediction_names=[]),
            caption=title + '.', plot_types=['line'])

    def plot_partial_dependences(self, filename=None):
        """Plot partial dependences for every feature.

        Parameters
        ----------
        filename : str, optional (default: 'partial_dependece_{feature}')
            Name of the plot file.

        Raises
        ------
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.

        """
        if not self._is_ready_for_plotting():
            return
        logger.info("Plotting partial dependences")
        if filename is None:
            filename = 'partial_dependece_{feature}'

        # Plot for every feature
        x_train = self.get_x_array('train', impute_nans=True)
        verbosity = self._get_verbosity_parameters(plot_partial_dependence)
        for feature_name in self.features:
            logger.debug("Plotting partial dependence of '%s'", feature_name)
            display = plot_partial_dependence(
                self._clf,
                x_train,
                features=[feature_name],
                feature_names=self.features,
                method='brute',
                line_kw={'color': 'b'},
                **verbosity,
            )
            title = (f"Partial dependence of {self.label} on {feature_name} "
                     f"for MLR model {self._cfg['mlr_model_name']}")
            plt.title(title)
            plt.xlabel(self._get_plot_feature(feature_name))
            plt.ylabel(self._get_plot_label())

            # Save plot
            new_filename = (filename.format(feature=feature_name) + '.' +
                            self._cfg['output_file_type'])
            plot_path = os.path.join(self._cfg['mlr_plot_dir'], new_filename)
            plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save provenance
            cube = mlr.get_1d_cube(
                display.lines_[0, 0].get_xdata(),
                display.lines_[0, 0].get_ydata(),
                x_kwargs={'var_name': feature_name,
                          'long_name': feature_name,
                          'units': self.features_units[feature_name]},
                y_kwargs={'var_name': self.label,
                          'long_name': self.label,
                          'units': self.label_units,
                          'attributes': {'project': '', 'dataset': ''}},
            )
            self._write_plot_provenance(
                cube, plot_path,
                ancestors=self.get_ancestors(prediction_names=[]),
                caption=title + '.', plot_types=['line'])

    def plot_prediction_errors(self, filename=None):
        """Plot predicted vs. true values.

        Parameters
        ----------
        filename : str, optional (default: 'prediction_errors')
            Name of the plot file.

        Raises
        ------
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.

        """
        if not self._is_ready_for_plotting():
            return
        logger.info("Plotting prediction errors")
        if filename is None:
            filename = 'prediction_errors'
        (_, axes) = plt.subplots()

        # Get available datasets
        data_to_plot = ['train']
        if 'test' in self.data:
            data_to_plot.append('test')

        # Create plot
        y_pred_all = []
        y_true_all = []
        data_types = []
        for data_type in data_to_plot:
            logger.debug("Plotting prediction error of '%s' data", data_type)
            x_data = self.data[data_type].x
            y_pred = self._clf.predict(x_data)
            y_true = self.get_y_array(data_type)
            axes.scatter(
                y_pred, y_true,
                **self._get_plot_kwargs(data_type, plot_type='scatter'))

            # Collect data
            y_pred_all.append(y_pred)
            y_true_all.append(y_true)
            data_types.append(np.full(y_pred.shape, data_type))

        # Plot appearance
        lims = [
            np.min([axes.get_xlim(), axes.get_ylim()]),
            np.max([axes.get_xlim(), axes.get_ylim()]),
        ]
        axes.plot(lims, lims, linestyle='--', color='k', alpha=0.75)
        axes.set_aspect('equal')
        axes.set_xlim(lims)
        axes.set_ylim(lims)
        title = (f"Prediction errors of {self.label} "
                 f"({self._cfg['mlr_model_name']})")
        axes.set_title(title)
        axes.set_xlabel(f'Predicted {self._get_plot_label()}')
        axes.set_ylabel(f'True {self._get_plot_label()}')
        axes.legend(loc='upper left')

        # Save plot
        plot_path = os.path.join(
            self._cfg['mlr_plot_dir'],
            filename + '.' + self._cfg['output_file_type'],
        )
        plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save provenance
        cube = mlr.get_1d_cube(
            np.concatenate(y_pred_all),
            np.concatenate(y_true_all),
            x_kwargs={'var_name': self.label,
                      'long_name': f'Predicted {self.label}',
                      'units': self.label_units},
            y_kwargs={'var_name': self.label,
                      'long_name': f'True {self.label}',
                      'units': self.label_units,
                      'attributes': {'project': '', 'dataset': ''}},
        )
        cube.add_aux_coord(
            self._get_data_type_coord(np.concatenate(data_types)), 0)
        self._write_plot_provenance(
            cube, plot_path, ancestors=self.get_ancestors(prediction_names=[]),
            caption=title + '.', plot_types=['scatter'])

    def plot_residuals(self, filename=None):
        """Plot residuals of training and test (if available) data.

        Parameters
        ----------
        filename : str, optional (default: 'residuals')
            Name of the plot file.

        Raises
        ------
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.

        """
        if not self._is_ready_for_plotting():
            return
        logger.info("Plotting residuals")
        if filename is None:
            filename = 'residuals'
        (_, axes) = plt.subplots()

        # Get available datasets
        data_to_plot = ['train']
        if 'test' in self.data:
            data_to_plot.append('test')

        # Create plot
        y_pred_all = []
        y_res_all = []
        data_types = []
        for data_type in data_to_plot:
            logger.debug("Plotting residuals of '%s' data", data_type)
            x_data = self.data[data_type].x
            y_pred = self._clf.predict(x_data)
            y_true = self.get_y_array(data_type)
            y_res = self._get_residuals(y_true, y_pred)
            axes.scatter(
                y_pred, y_res,
                **self._get_plot_kwargs(data_type, plot_type='scatter'))

            # Collect data
            y_pred_all.append(y_pred)
            y_res_all.append(y_res)
            data_types.append(np.full(y_pred.shape, data_type))

        # Plot appearance
        axes.axhline(0.0, linestyle='--', color='k', alpha=0.75)
        axes.set_aspect('equal')
        title = (f"Residuals of {self.label} ({self._cfg['mlr_model_name']})")
        axes.set_title(title)
        axes.set_xlabel(f'Predicted {self._get_plot_label()}')
        axes.set_ylabel(f'Residuals of {self._get_plot_label()}')
        self._set_axis_lim_symmetric(axes, 'y')
        axes.legend(loc='best')

        # Save plot
        plot_path = os.path.join(
            self._cfg['mlr_plot_dir'],
            filename + '.' + self._cfg['output_file_type'],
        )
        plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save provenance
        cube = mlr.get_1d_cube(
            np.concatenate(y_pred_all),
            np.concatenate(y_res_all),
            x_kwargs={'var_name': self.label,
                      'long_name': f'Predicted {self.label}',
                      'units': self.label_units},
            y_kwargs={'var_name': self.label,
                      'long_name': f'Residuals of {self.label}',
                      'units': self.label_units,
                      'attributes': {'project': '', 'dataset': ''}},
        )
        cube.add_aux_coord(
            self._get_data_type_coord(np.concatenate(data_types)), 0)
        self._write_plot_provenance(
            cube, plot_path, ancestors=self.get_ancestors(prediction_names=[]),
            caption=title + '.', plot_types=['scatter'])

    def plot_residuals_histogram(self, filename=None):
        """Plot histogram of residuals of training and test data.

        Parameters
        ----------
        filename : str, optional (default: 'residuals_histogram')
            Name of the plot file.

        Raises
        ------
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.

        """
        if not self._is_ready_for_plotting():
            return
        logger.info("Plotting residuals histogram")
        if filename is None:
            filename = 'residuals_histogram'
        (_, axes) = plt.subplots()

        # Get available datasets
        data_to_plot = ['train']
        if 'test' in self.data:
            data_to_plot.append('test')

        # Create plot (centralize bins around the zero)
        y_res_all = []
        freq_all = []
        data_types = []
        for data_type in data_to_plot:
            logger.debug("Plotting residuals histogram of '%s' data",
                         data_type)
            x_data = self.data[data_type].x
            y_pred = self._clf.predict(x_data)
            y_true = self.get_y_array(data_type)
            y_res = self._get_residuals(y_true, y_pred)
            bins = self._get_centralized_bins(y_res, n_bins=20)
            hist = axes.hist(y_res, bins=bins,
                             **self._get_plot_kwargs(data_type))

            # Collect data
            y_res_all.append(np.convolve(hist[1], (1, 1), 'valid') / 2.0)
            freq_all.append(hist[0])
            data_types.append(np.full(hist[0].shape, data_type))

        # Plot appearance
        axes.axvline(0.0, linestyle='--', color='k', alpha=0.75)
        title = (f"Histogram for residuals of {self.label} "
                 f"({self._cfg['mlr_model_name']})")
        axes.set_title(title)
        axes.set_xlabel(f'Residuals of {self._get_plot_label()}')
        axes.set_ylabel('Frequency')
        self._set_axis_lim_symmetric(axes, 'x')
        axes.legend(loc='best')

        # Save plot
        plot_path = os.path.join(
            self._cfg['mlr_plot_dir'],
            filename + '.' + self._cfg['output_file_type'],
        )
        plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save provenance
        cube = mlr.get_1d_cube(
            np.concatenate(y_res_all),
            np.concatenate(freq_all),
            x_kwargs={'var_name': self.label,
                      'long_name': f'Residuals of {self.label}',
                      'units': self.label_units},
            y_kwargs={'var_name': 'frequency',
                      'long_name': 'Frequency',
                      'units': '1',
                      'attributes': {'project': '', 'dataset': ''}},
        )
        cube.add_aux_coord(
            self._get_data_type_coord(np.concatenate(data_types)), 0)
        self._write_plot_provenance(
            cube, plot_path, ancestors=self.get_ancestors(prediction_names=[]),
            caption=title + '.', plot_types=['histogram'])

    def plot_residuals_distribution(self, filename=None):
        """Plot distribution of residuals of training and test data (KDE).

        Parameters
        ----------
        filename : str, optional (default: 'residuals_distribution')
            Name of the plot file.

        Raises
        ------
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.

        """
        if not self._is_ready_for_plotting():
            return
        logger.info("Plotting residuals distribution")
        if filename is None:
            filename = 'residuals_distribution'

        # Get available datasets
        data_to_plot = ['train']
        if 'test' in self.data:
            data_to_plot.append('test')

        # Create plot (centralize bins around the zero)
        data_types = []
        for data_type in data_to_plot:
            logger.debug("Plotting residuals distribution of '%s' data",
                         data_type)
            x_data = self.data[data_type].x
            y_pred = self._clf.predict(x_data)
            y_true = self.get_y_array(data_type)
            y_res = self._get_residuals(y_true, y_pred)
            axes = sns.kdeplot(y_res, **self._get_plot_kwargs(data_type))

            # Collect data
            data_types.append(np.full(axes.lines[-1].get_xdata().shape,
                                      data_type))

        # Plot appearance
        axes.axvline(0.0, linestyle='--', color='k', alpha=0.75)
        title = (f"Probability distribution of residuals of {self.label} "
                 f"({self._cfg['mlr_model_name']})")
        axes.set_title(title)
        axes.set_xlabel(f'Residuals of {self._get_plot_label()}')
        axes.set_ylabel('Probability density')
        self._set_axis_lim_symmetric(axes, 'x')
        axes.legend(loc='best')

        # Save plot
        plot_path = os.path.join(
            self._cfg['mlr_plot_dir'],
            filename + '.' + self._cfg['output_file_type'],
        )
        plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        plt.close()

        # Save provenance
        cube = mlr.get_1d_cube(
            np.concatenate([line.get_xdata() for line in axes.lines[:-1]]),
            np.concatenate([line.get_ydata() for line in axes.lines[:-1]]),
            x_kwargs={'var_name': self.label,
                      'long_name': f'Residuals of {self.label}',
                      'units': self.label_units},
            y_kwargs={'var_name': 'probability_density',
                      'long_name': 'Probability Density',
                      'units': '1',
                      'attributes': {'project': '', 'dataset': ''}},
        )
        cube.add_aux_coord(
            self._get_data_type_coord(np.concatenate(data_types)), 0)
        self._write_plot_provenance(
            cube, plot_path, ancestors=self.get_ancestors(prediction_names=[]),
            caption=title + '.', plot_types=['probability'])

    def plot_scatterplots(self, filename=None):
        """Plot scatterplots label vs. feature for every feature.

        Parameters
        ----------
        filename : str, optional (default: 'scatterplot_{feature}')
            Name of the plot file.

        Raises
        ------
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.

        """
        if not self._is_ready_for_plotting():
            return
        logger.info("Plotting scatterplots")
        if filename is None:
            filename = 'scatterplot_{feature}'

        # Plot scatterplot for every feature
        for feature in self.features:
            logger.debug("Plotting scatterplot of '%s'", feature)
            (_, axes) = plt.subplots()

            # Iterate over group attributes
            for group_attr in self.group_attributes:
                group_attr = self._group_attr_to_pandas_index_str(group_attr)
                axes.plot(self.data['all'].x.loc[group_attr, feature],
                          self.data['all'].y.loc[group_attr, self.label],
                          '.', label=group_attr)

            # Plot appearance
            axes.legend(loc='center left', ncol=2, bbox_to_anchor=[1.05, 0.5],
                        borderaxespad=0.0)
            title = f"Target variable {self.label} vs. feature {feature}"
            axes.set_title(title)
            axes.set_xlabel(self._get_plot_feature(feature))
            axes.set_ylabel(self._get_plot_label())

            # Save plot
            plot_path = os.path.join(
                self._cfg['mlr_plot_dir'],
                filename.format(feature=feature) + '.' +
                self._cfg['output_file_type'])
            plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
            logger.info("Wrote %s", plot_path)
            plt.close()

            # Save provenance
            cube = mlr.get_1d_cube(
                self.data['all'].x.loc[:, feature].values,
                self.get_y_array('all'),
                x_kwargs={'var_name': feature,
                          'long_name': feature,
                          'units': self.features_units[feature]},
                y_kwargs={'var_name': self.label,
                          'long_name': self.label,
                          'units': self.label_units,
                          'attributes': {'project': '', 'dataset': ''}},
            )
            ancestors = self.get_ancestors(features=[feature],
                                           prediction_names=[])
            self._write_plot_provenance(
                cube, plot_path, ancestors=ancestors, caption=title + '.',
                plot_types=['scatter'])

    def predict(self,
                save_mlr_model_error=None,
                save_lime_importance=False,
                save_propagated_errors=False,
                **kwargs):
        """Perform prediction using the MLR model(s) and write ``*.nc`` files.

        Parameters
        ----------
        save_mlr_model_error : str or int, optional
            Additionally saves estimated squared MLR model error. This error
            represents the uncertainty of the prediction caused by the MLR
            model itself and not by errors in the prediction input data (errors
            in that will be considered by including datasets with ``var_type``
            set to ``prediction_input_error`` and setting
            ``save_propagated_errors`` to ``True``). If the option is set to
            ``'test'``, the (constant) error is estimated as RMSEP using a
            (hold-out) test data set. Only possible if test data is available,
            i.e.  the option ``test_size`` is not set to ``False`` during class
            initialization. If the option is set to ``'logo'``, the (constant)
            error is estimated as RMSEP using leave-one-group-out
            cross-validation using the group_attributes. Only possible if
            ``group_datasets_by_attributes`` is given. If the option is set to
            an integer ``n`` (!= 0), the (constant) error is estimated as RMSEP
            using n-fold cross-validation.
        save_lime_importance : bool, optional (default: False)
            Additionally saves local feature importance given by LIME (Local
            Interpretable Model-agnostic Explanations).
        save_propagated_errors : bool, optional (default: False)
            Additionally saves propagated errors from
            ``prediction_input_error`` datasets. Only possible when these are
            available.
        **kwargs : keyword arguments, optional
            Additional options for the final regressors ``predict()`` function.

        Raises
        ------
        RuntimeError
            ``return_var`` and ``return_cov`` are both set to ``True``.
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.
        ValueError
            An invalid value for ``save_mlr_model_error`` is given.
        ValueError
            ``save_propagated_errors`` is ``True`` and no
            ``prediction_input_error`` data is available.

        """
        self._check_fit_status('Prediction')
        logger.info("Started prediction")
        mlr.check_predict_kwargs(kwargs)
        if kwargs:
            logger.info(
                "Using additional keyword argument(s) %s for predict() "
                "function", kwargs)

        # Iterate over different predictions
        for pred_name in self._datasets['prediction_input']:
            logger.info("Predicting '%s'", self._get_name(pred_name))

            # Prediction
            (x_pred, x_err, y_ref,
             x_cube) = self._extract_prediction_input(pred_name)
            pred_dict = self._get_prediction_dict(
                pred_name, x_pred, x_err, y_ref,
                get_mlr_model_error=save_mlr_model_error,
                get_lime_importance=save_lime_importance,
                get_propagated_errors=save_propagated_errors, **kwargs)

            # Save data in class member
            y_pred = pd.DataFrame(pred_dict[None],
                                  columns=[self.label],
                                  index=x_pred.index,
                                  dtype=self._cfg['dtype'])
            self._data['pred'][pred_name] = pd.concat([x_pred, y_pred],
                                                      axis=1,
                                                      keys=['x', 'y'])

            # Save prediction cubes
            self._save_prediction_cubes(pred_dict, pred_name, x_cube)

    def print_correlation_matrices(self):
        """Print correlation matrices for all datasets."""
        self._check_fit_status('Printing correlation matrices')
        for data_type in ('all', 'train', 'test'):
            if data_type not in self.data:
                continue
            logger.info("Correlation matrix for %s data:\n%s", data_type,
                        self.data[data_type][['x', 'y']].corr())

    def print_regression_metrics(self, logo=False):
        """Print all available regression metrics for training data.

        Parameters
        ----------
        logo : bool, optional (default: False)
            Print regression metrics using
            :class:`sklearn.model_selection.LeaveOneGroupOut` cross-validation.
            Only possible when `group_datasets_by_attributes` was given during
            class initialization.

        """
        self._check_fit_status('Printing regression metrics')
        regression_metrics = [
            'explained_variance_score',
            'mean_absolute_error',
            'mean_squared_error',
            'r2_score',
        ]

        # Metrics on train and test data
        for data_type in ('all', 'train', 'test'):
            self._print_metrics(regression_metrics, data_type)
            logger.info("")

        # Metrics on CV data
        if logo:
            logger.info(
                "Evaluating regression metrics using 'LeaveOneGroupOut' "
                "cross-validation using group attributes %s on training data",
                self._cfg['group_datasets_by_attributes'])
            regression_metrics = {
                'explained_variance_score': 'explained_variance',
                'mean_absolute_error': 'neg_mean_absolute_error',
                'root_mean_squared_error': 'neg_root_mean_squared_error',
                'r2_score': 'r2',
            }
            x_data = self.data['train'].x
            y_data = self.get_y_array('train')
            sample_weights = self._get_sample_weights('train')
            for (metric, scoring) in regression_metrics.items():
                value = cross_val_score_weighted(
                    self._clf, x_data, y_data, scoring=scoring,
                    n_jobs=self._cfg['n_jobs'], fit_params=self.fit_kwargs,
                    **self._get_verbosity_parameters(cross_val_score_weighted),
                    **self._get_logo_cv_kwargs())
                value = np.mean(value)
                if 'neg_' in scoring:
                    value = -value
                logger.info("%s: %s", metric, value)
            if sample_weights is None:
                return
            for (metric, scoring) in regression_metrics.items():
                value = cross_val_score_weighted(
                    self._clf, x_data, y_data, scoring=scoring,
                    n_jobs=self._cfg['n_jobs'], fit_params=self.fit_kwargs,
                    sample_weights=sample_weights,
                    **self._get_verbosity_parameters(cross_val_score_weighted),
                    **self._get_logo_cv_kwargs())
                value = np.mean(value)
                if 'neg_' in scoring:
                    value = -value
                logger.info("Weighted %s: %s", metric, value)

    def reset_pipeline(self):
        """Reset regressor pipeline."""
        steps = []
        numerical_features_idx = [
            int(np.where(self.features == tag)[0][0])
            for tag in self.numerical_features
        ]

        # Imputer
        if self._cfg['imputation_strategy'] != 'remove':
            verbosity = self._get_verbosity_parameters(SimpleImputer)
            imputer = SimpleImputer(
                strategy=self._cfg['imputation_strategy'],
                **verbosity,
            )
            steps.append(('imputer', imputer))

        # Scaler for numerical features
        if self._cfg['standardize_data']:
            x_scaler = ColumnTransformer(
                [('', StandardScaler(), numerical_features_idx)],
                remainder='passthrough',
            )
            steps.append(('x_scaler', x_scaler))

        # PCA for numerical features
        if self._cfg.get('pca'):
            pca = ColumnTransformer(
                [('', PCA(), numerical_features_idx)],
                remainder='passthrough',
            )
            steps.append(('pca', pca))

        # Final regressor
        final_parameters = self._load_final_parameters()
        final_regressor = self._CLF_TYPE(**final_parameters)

        # Transformer for labels if desired (if not, add pd to np converter)
        if self._cfg['standardize_data']:
            y_scaler = StandardScaler()
        else:
            y_scaler = StandardScaler(with_mean=False, with_std=False)
        transformed_target_regressor = AdvancedTransformedTargetRegressor(
            transformer=y_scaler, regressor=final_regressor)
        steps.append(('final', transformed_target_regressor))

        # Final pipeline
        if self._cfg['cache_intermediate_results']:
            if self._cfg['n_jobs'] is None or self._cfg['n_jobs'] == 1:
                memory = self._cfg['mlr_work_dir']
            else:
                logger.debug(
                    "Caching intermediate results of Pipeline is not "
                    "supported for multiple processes (using at most %i "
                    "processes)", self._cfg['n_jobs'])
                memory = None
        else:
            memory = None
        self._clf = AdvancedPipeline(steps, memory=memory)
        logger.info("Created pipeline with steps %s",
                    list(self._clf.named_steps.keys()))

    def rfecv(self, **kwargs):
        """Perform recursive feature elimination using cross-validation.

        Note
        ----
        This only works for final estimators that provide information about
        feature importance either through a ``coef_`` attribute or through a
        ``feature_importances_`` attribute.

        Parameters
        ----------
        **kwargs : keyword arguments, optional
            Additional options for :class:`sklearn.feature_selection.RFECV`.

        Raises
        ------
        RuntimeError
            Final estimator does not provide ``coef_`` or
            ``feature_importances_`` attribute.

        """
        logger.info(
            "Performing recursive feature elimination using cross-validation "
            "with final regressor %s on %i training points", self._CLF_TYPE,
            len(self.data['train'].index))

        # Get keyword arguments
        (cv_kwargs, fit_kwargs) = self._get_cv_estimator_kwargs(AdvancedRFECV,
                                                                **kwargs)
        fit_kwargs = deepcopy(fit_kwargs)
        keys_to_remove = []
        for key in fit_kwargs:
            if key.endswith('eval_set'):
                keys_to_remove.append(key)
        for key in keys_to_remove:
            logger.warning(
                "Fit parameter '%s' is not supported for rfecv()", key)
            fit_kwargs.pop(key)

        # Create and fit AdvancedRFECV instance
        rfecv = AdvancedRFECV(self._clf, **cv_kwargs)
        rfecv.fit(self.data['train'].x, self.get_y_array('train'),
                  **fit_kwargs)

        # Add feature selection step to pipeline
        self._clf = rfecv.estimator_
        transformer = get_rfecv_transformer(rfecv)
        self._clf.steps.insert(0, ('feature_selection', transformer))

        # Log results
        new_features = self.features[rfecv.support_]
        logger.info(
            "Recursive feature elimination was successful, %i of the %i "
            "features remain", new_features.size, self.features.size)
        logger.info("Old features: %s", self.features)
        logger.info("New features: %s", new_features)
        logger.info("Successfully fitted MLR model on %i training point(s)",
                    len(self.data['train'].index))
        logger.debug("Pipeline steps:")
        logger.debug(pformat(list(self._clf.named_steps.keys())))
        logger.debug("Parameters:")
        logger.debug(pformat(self.parameters))

        # LIME
        self._load_lime_explainer()

    def test_normality_of_residuals(self):
        """Perform Shapiro-Wilk test to normality of residuals.

        Raises
        ------
        sklearn.exceptions.NotFittedError
            MLR model is not fitted.

        """
        if not self._is_ready_for_plotting():
            return

        # Get available datasets
        data_to_check = ['train']
        if 'test' in self.data:
            data_to_check.append('test')

        # Perform Shapiro-Wilk test
        for data_type in data_to_check:
            x_data = self.data[data_type].x
            y_pred = self._clf.predict(x_data)
            y_true = self.get_y_array(data_type)
            y_res = self._get_residuals(y_true, y_pred)
            (w_value, p_value) = shapiro(y_res)
            logger.info(
                "Result of Shapiro-Wilk test for normality of residuals: W = "
                "%.5f, p = %.5f", w_value, p_value)

    def update_parameters(self, **params):
        """Update parameters of the whole pipeline.

        Note
        ----
        Parameter names have to be given for each step of the pipeline
        separated by two underscores, i.e. ``s__p`` is the parameter ``p`` for
        step ``s``.

        Parameters
        ----------
        **params : keyword arguments, optional
            Parameters for the pipeline which should be updated.

        Raises
        ------
        ValueError
            Invalid parameter for pipeline given.

        """
        allowed_params = self._get_clf_parameters()
        new_params = {}
        for (key, val) in params.items():
            if key in allowed_params:
                new_params[key] = val
            else:
                raise ValueError(
                    f"'{key}' is not a valid parameter for the pipeline")
        self._clf.set_params(**new_params)
        self._parameters = self._get_clf_parameters()
        if new_params:
            logger.info("Updated pipeline with parameters %s", new_params)

    def _calculate_sample_weights(self, cube, var_type, group_attr=None):
        """Calculate sample weights if desired."""
        if not self._cfg['weighted_samples']:
            return None
        if var_type != 'feature':
            return None
        weights = mlr.get_all_weights(cube, **self._cfg['weighted_samples'])
        weights = weights.astype(self._cfg['dtype'], casting='same_kind')
        weights = pd.DataFrame(
            {'sample_weight': weights.ravel()},
            index=self._get_multiindex(cube, group_attr=group_attr),
            dtype=self._cfg['dtype'],
        )
        msg = '' if group_attr is None else f" of '{group_attr}'"
        logger.debug(
            "Successfully calculated %i sample weights for training data%s "
            "using %s", len(weights.index), msg, self._cfg['weighted_samples'])
        return weights

    def _check_clf(self):
        """Check if valid regressor type is given."""
        class_name = self.__class__.__name__
        if self._CLF_TYPE is None:
            raise NotImplementedError(
                f"No MLR model type specified, please use the factory "
                f"function 'esmvaltool.diag_scripts.mlr.models.{class_name}."
                f"create()' to initialize this class")

    def _check_cube_dimensions(self, cube, ref_cube, text=None):
        """Check shape and coordinates of a given cube."""
        msg = '' if text is None else f' for {text}'
        if self._cfg.get('accept_only_scalar_data'):
            allowed_shapes = [(), (1, )]
            if cube.shape not in allowed_shapes:
                raise ValueError(
                    f"Expected only cubes with shapes {allowed_shapes} when "
                    f"option 'accept_only_scalar_data' is set to 'True', got "
                    f"{cube.shape}{msg}")
        else:
            if ref_cube is None:
                return
            if cube.shape != ref_cube.shape:
                raise ValueError(
                    f"Expected cubes with shapes {ref_cube.shape}{msg}, got "
                    f"{cube.shape}. Consider regridding, pre-selecting data "
                    f"at class initialization (argument 'input_datasets') or "
                    f"the options 'broadcast_from' or 'group_datasets_by_"
                    f"attributes'")
            cube_coords = cube.coords(dim_coords=True)
            ref_coords = ref_cube.coords(dim_coords=True)
            cube_coords_str = [
                f'{coord.name()}, shape {coord.shape}' for coord in cube_coords
            ]
            ref_coords_str = [
                f'{coord.name()}, shape {coord.shape}' for coord in ref_coords
            ]
            if cube_coords_str != ref_coords_str:
                logger.warning(
                    "Cube coordinates differ, expected %s%s, got %s. Check "
                    "input cubes", ref_coords_str, msg, cube_coords_str)
                return
            for (idx, cube_coord) in enumerate(cube_coords):
                ref_coord = ref_coords[idx]
                if not np.allclose(cube_coord.points, ref_coord.points):
                    logger.warning(
                        "'%s' coordinate for different cubes does not "
                        "match, got %s%s, expected %s (values differ by "
                        "more than allowed tolerance, check input cubes)",
                        cube_coord.name(), cube_coord.points, msg,
                        ref_coord.points)

    def _check_dataset(self, datasets, var_type, tag, text=None):
        """Check if datasets exist and are valid."""
        datasets = select_metadata(datasets, tag=tag, var_type=var_type)
        msg = '' if text is None else text
        if not datasets:
            if var_type == 'prediction_input_error':
                return None
            if var_type == 'prediction_reference':
                return None
            if var_type == 'label':
                raise ValueError(f"Label '{tag}'{msg} not found")
            if not self._cfg.get('allow_missing_features'):
                raise ValueError(
                    f"{var_type} '{tag}'{msg} not found, use 'allow_missing_"
                    f"features' to ignore this")
            logger.info(
                "Ignored missing %s '%s'%s since 'allow_missing_features' is "
                "set to 'True'", var_type, tag, msg)
            return None
        if len(datasets) > 1:
            raise ValueError(
                f"{var_type} '{tag}'{msg} not unique, consider adapting the "
                f"argument 'input_datasets' at class initialization to "
                f"pre-select datasets or specify suitable attributes to group "
                f"datasets with the option 'group_datasets_by_attributes'")
        if var_type in ('label', 'prediction_reference'):
            units = self.label_units
        else:
            units = self.features_units[tag]
        if units != Unit(datasets[0]['units']):
            raise ValueError(
                f"Expected units '{units}' for {var_type} '{tag}'{msg}, got "
                f"'{datasets[0]['units']}'")
        return datasets[0]

    def _check_fit_status(self, text):
        """Check if MLR model is fitted and raise exception otherwise."""
        x_dummy = np.ones((1, self.features.size), dtype=self._cfg['dtype'])
        try:
            self._clf.predict(x_dummy)
        except NotFittedError:
            raise NotFittedError(
                f"{text} not possible, MLR model {self._CLF_TYPE} is not "
                f"fitted yet, call fit(), grid_search_cv() or rfecv() first")

    def _estimate_mlr_model_error(self, target_length, strategy):
        """Estimate squared error of MLR model (using CV or test data)."""
        logger.info(
            "Estimating squared error of MLR model using strategy '%s'",
            strategy)

        # Estimate MLR model error
        if strategy == 'test':
            if 'test' not in self.data:
                raise ValueError(
                    f"'save_mlr_model_error' using strategy 'test' is not "
                    f"possible because no test data is available ('test_size' "
                    f"was set to '{self._cfg['test_size']}' during class "
                    f"initialization)")
            y_pred = self._clf.predict(self.data['test'].x)
            error = metrics.mean_squared_error(
                self.get_y_array('test'),
                y_pred,
                sample_weight=self._get_sample_weights('test'),
            )
        else:
            if strategy == 'logo':
                cv_kwargs = self._get_logo_cv_kwargs()
            elif isinstance(strategy, int):
                cv_kwargs = {'cv': strategy}
            else:
                raise ValueError(
                    f"Expected 'test', 'logo' or an integer as strategy for "
                    f"estimating MLR model error (argument "
                    f"'save_mlr_model_error'), got '{strategy}'")
            x_data = self.data['train'].x
            y_data = self.get_y_array('train')
            error = cross_val_score_weighted(
                self._clf, x_data, y_data, scoring='neg_mean_squared_error',
                n_jobs=self._cfg['n_jobs'], fit_params=self.fit_kwargs,
                sample_weights=self._get_sample_weights('train'),
                **self._get_verbosity_parameters(cross_val_score_weighted),
                **cv_kwargs)
            error = -np.mean(error)

        # Reshape error
        error_array = np.full(target_length, error, dtype=self._cfg['dtype'])
        units = mlr.units_power(self.label_units, 2)
        logger.info(
            "Estimated squared MLR model error by %s %s using strategy '%s'",
            error, units, strategy)
        return error_array

    def _extract_features_and_labels(self):
        """Extract feature and label data points from training data."""
        (x_data, _,
         sample_weights) = self._extract_x_data(self._datasets['feature'],
                                                'feature')
        y_data = self._extract_y_data(self._datasets['label'], 'label')

        # Check number of input points
        if not x_data.index.equals(y_data.index):
            raise ValueError(
                f"Got differing point(s) for features and labels ("
                f"{len(x_data.index):d} feature points and "
                f"{len(y_data.index):d} label points):\n"
                f"{x_data.index.difference(y_data.index)}")
        logger.info("Found %i raw input data point(s) with data type '%s'",
                    len(y_data.index), self._cfg['dtype'])

        # Remove missing values in labels
        (x_data, y_data,
         sample_weights) = self._remove_missing_labels(x_data, y_data,
                                                       sample_weights)

        # Remove missing values in features (if desired)
        (x_data, y_data, sample_weights) = self._remove_missing_features(
            x_data, y_data, sample_weights)

        return (x_data, y_data, sample_weights)

    def _extract_prediction_input(self, prediction_name):
        """Extract prediction input data points for ``prediction_name``."""
        (x_pred, x_cube, _) = self._extract_x_data(
            self._datasets['prediction_input'][prediction_name],
            'prediction_input')
        logger.info(
            "Found %i raw prediction input data point(s) with data type '%s'",
            len(x_pred.index), self._cfg['dtype'])

        # Prediction reference
        if prediction_name not in self._datasets['prediction_reference']:
            y_ref = None
            logger.debug(
                "No prediction reference for prediction '%s' available",
                self._get_name(prediction_name))
        else:
            y_ref = self._extract_y_data(
                self._datasets['prediction_reference'][prediction_name],
                'prediction_reference')
            if y_ref is not None:
                if not x_pred.index.equals(y_ref.index):
                    raise ValueError(
                        f"Got differing point(s) for prediction input and "
                        f"prediction output ({len(x_pred.index):d} "
                        f"prediction input points and {len(y_ref.index):d} "
                        f"prediction output points):\n"
                        f"{x_pred.index.difference(y_ref.index)}")
                logger.info(
                    "Found %i raw prediction output data point(s) with data "
                    "type '%s'", len(y_ref.index), self._cfg['dtype'])

        # Error
        if prediction_name not in self._datasets['prediction_input_error']:
            x_err = None
            logger.debug(
                "Propagating prediction input errors for prediction '%s' not "
                "possible, no 'prediction_input_error' datasets given",
                self._get_name(prediction_name))
        else:
            (x_err, _, _) = self._extract_x_data(
                self._datasets['prediction_input_error'][prediction_name],
                'prediction_input_error')
            if not x_pred.index.equals(x_err.index):
                raise ValueError(
                    f"Got differing point(s) for prediction input and "
                    f"prediction input error ({len(x_pred.index):d} "
                    f"prediction input points and {len(x_err.index):d} "
                    f"prediction input error points):\n"
                    f"{x_pred.index.difference(x_err.index)}")
            logger.info(
                "Found %i raw prediction input error data point(s) with data "
                "type '%s'", len(x_err.index), self._cfg['dtype'])

        # Remove missing values if necessary
        (x_pred, x_err, y_ref,
         mask) = self._remove_missing_pred_input(x_pred, x_err, y_ref)

        # Create cube with appropriate mask for output
        mask = mask.reshape(x_cube.shape)
        cube_data = np.empty(mask.shape, dtype=self._cfg['dtype'])
        x_cube.data = np.ma.array(cube_data, mask=mask)

        return (x_pred, x_err, y_ref, x_cube)

    def _extract_x_data(self, datasets, var_type):
        """Extract required x data of type ``var_type`` from ``datasets``."""
        allowed_types = ('feature', 'prediction_input',
                         'prediction_input_error')
        if var_type not in allowed_types:
            raise ValueError(
                f"Excepted one of '{allowed_types}' for 'var_type', got "
                f"'{var_type}'")
        x_data = pd.DataFrame(columns=self.features, dtype=self._cfg['dtype'])
        x_cube = None
        if self._cfg['weighted_samples'] and var_type == 'feature':
            sample_weights = pd.DataFrame(columns=['sample_weight'],
                                          dtype=self._cfg['dtype'])
        else:
            sample_weights = None

        # Iterate over datasets
        datasets = select_metadata(datasets, var_type=var_type)
        if var_type == 'feature':
            groups = self.group_attributes
        else:
            groups = [None]
        for group_attr in groups:
            group_datasets = select_metadata(datasets,
                                             group_attribute=group_attr)
            if group_attr is not None:
                logger.info("Loading '%s' data of '%s'", var_type, group_attr)
            msg = '' if group_attr is None else f" for '{group_attr}'"
            if not group_datasets:
                raise ValueError(f"No '{var_type}' data{msg} found")
            (group_data, x_cube,
             weights) = self._get_x_data_for_group(group_datasets, var_type,
                                                   group_attr)
            x_data = x_data.append(group_data)

            # Append weights if desired
            if sample_weights is not None:
                sample_weights = sample_weights.append(weights)

        # Adapt sample_weights if necessary
        if sample_weights is not None:
            sample_weights.index = pd.MultiIndex.from_tuples(
                sample_weights.index, names=self._get_multiindex_names())
            logger.info(
                "Successfully calculated sample weights for training data "
                "using %s", self._cfg['weighted_samples'])
            if (sample_weights.max().values[0] /
                    sample_weights.min().values[0]) > 150.0:
                logger.warning(
                    "Sample weights differ by more than a factor of 150, got "
                    "a minimum value of %e and a maximum value of %e. This "
                    "might be caused by differing coordinates in the training "
                    "cubes",
                    sample_weights.min().values[0],
                    sample_weights.max().values[0])

        # Convert index back to MultiIndex
        x_data.index = pd.MultiIndex.from_tuples(
            x_data.index, names=self._get_multiindex_names())

        return (x_data, x_cube, sample_weights)

    def _extract_y_data(self, datasets, var_type):
        """Extract required y data of type ``var_type`` from ``datasets``."""
        allowed_types = ('label', 'prediction_reference')
        if var_type not in allowed_types:
            raise ValueError(
                f"Excepted one of '{allowed_types}' for 'var_type', got "
                f"'{var_type}'")
        y_data = pd.DataFrame(columns=[self.label], dtype=self._cfg['dtype'])

        # Iterate over datasets
        datasets = select_metadata(datasets, var_type=var_type)
        if var_type == 'label':
            groups = self.group_attributes
        else:
            groups = [None]
        for group_attr in groups:
            if group_attr is not None:
                logger.info("Loading '%s' data of '%s'", var_type, group_attr)
            msg = '' if group_attr is None else f" for '{group_attr}'"
            group_datasets = select_metadata(datasets,
                                             group_attribute=group_attr)
            dataset = self._check_dataset(group_datasets, var_type, self.label,
                                          msg)
            if dataset is None:
                return None
            cube = self._load_cube(dataset)
            text = f"{var_type} '{self.label}'{msg}"
            self._check_cube_dimensions(cube, None, text)
            cube_data = pd.DataFrame(
                self._get_cube_data(cube),
                columns=[self.label],
                index=self._get_multiindex(cube, group_attr=group_attr),
                dtype=self._cfg['dtype'],
            )
            y_data = y_data.append(cube_data)

        # Convert index back to MultiIndex
        y_data.index = pd.MultiIndex.from_tuples(
            y_data.index, names=self._get_multiindex_names())

        return y_data

    def _get_broadcasted_cube(self, dataset, ref_cube, text=None):
        """Get broadcasted cube."""
        msg = '' if text is None else text
        target_shape = ref_cube.shape
        cube_to_broadcast = self._load_cube(dataset)
        data_to_broadcast = np.ma.filled(cube_to_broadcast.data, np.nan)
        logger.info("Broadcasting %s from %s to %s", msg,
                    data_to_broadcast.shape, target_shape)
        broadcasted_data = iris.util.broadcast_to_shape(
            data_to_broadcast, target_shape, dataset['broadcast_from'])
        new_cube = ref_cube.copy(np.ma.masked_invalid(broadcasted_data))
        for idx in dataset['broadcast_from']:
            new_coord = new_cube.coord(dimensions=idx)
            new_coord.points = cube_to_broadcast.coord(new_coord).points
        logger.debug("Added broadcasted %s", msg)
        return new_cube

    def _get_clf_parameters(self, deep=True):
        """Get parameters of pipeline."""
        return self._clf.get_params(deep=deep)

    def _get_colors_for_features(self, color_coded=True):
        """Get colors for bars of feature importance plot."""
        features = self.features_after_preprocessing
        if not color_coded:
            colors = dict(zip(features, ['b'] * len(features)))
        else:
            if not np.array_equal(features, self.features):
                raise ValueError(
                    f"Extracting color-coded feature colors is not possible "
                    f"since features changed after preprocessing, before: "
                    f"{self.features}, after: {features}")
            colors = {}
            corrs = self.data['train'][['x', 'y']].corr()
            for feature in features:
                corr = corrs.loc[('y', self.label), ('x', feature)]
                color = 'r' if corr >= 0.0 else 'b'
                colors[feature] = color
        return colors

    def _get_cv_estimator_kwargs(self, cv_estimator, **kwargs):
        """Get keyword arguments for CV estimator class."""
        fit_kwargs = self.fit_kwargs
        verbosity = self._get_verbosity_parameters(cv_estimator)
        cv_kwargs = {
            'n_jobs': self._cfg['n_jobs'],
            **verbosity,
        }
        cv_kwargs.update(kwargs)
        logger.info("Using keyword argument(s) %s for class %s", cv_kwargs,
                    cv_estimator)
        if isinstance(cv_kwargs.get('cv'), str):
            if cv_kwargs['cv'].lower() == 'loo':
                cv_kwargs['cv'] = LeaveOneOut()
            if cv_kwargs['cv'].lower() == 'logo':
                cv_kwargs['cv'] = self._get_logo_cv_kwargs()['cv']
                fit_kwargs['groups'] = self._get_logo_cv_kwargs()['groups']
        return (cv_kwargs, fit_kwargs)

    def _get_features(self):
        """Extract all features from the ``prediction_input`` datasets."""
        logger.debug("Extracting features from 'prediction_input' datasets")
        pred_name = list(self._datasets['prediction_input'].keys())[0]
        pred_name_str = self._get_name(pred_name)
        datasets = self._datasets['prediction_input'][pred_name]
        (units,
         types) = self._get_features_of_datasets(datasets, 'prediction_input',
                                                 pred_name)

        # Mark categorical variables
        categorical = {feature: False for feature in types}
        for tag in self._cfg.get('categorical_features', []):
            if tag in categorical:
                logger.debug("Treating '%s' as categorical feature", tag)
                categorical[tag] = True
            else:
                raise ValueError(
                    f"Cannot treat '{tag}' as categorical variable, feature "
                    f"not found")

        # Check if features were found
        if not units:
            raise ValueError(
                f"No features for 'prediction_input' data for prediction "
                f"'{pred_name_str}' found")

        # Check for wrong options
        if self._cfg.get('accept_only_scalar_data'):
            if 'broadcasted' in types.values():
                raise TypeError(
                    "The use of 'broadcast_from' is not possible if "
                    "'accept_only_scalar_data' is given")
            if 'coordinate' in types.values():
                raise TypeError(
                    "The use of 'coords_as_features' is not possible if "
                    "'accept_only_scalar_data' is given")

        # Convert to DataFrame and sort it
        units = pd.DataFrame.from_dict(units,
                                       orient='index',
                                       columns=['units'])
        types = pd.DataFrame.from_dict(types,
                                       orient='index',
                                       columns=['types'])
        categorical = pd.DataFrame.from_dict(categorical,
                                             orient='index',
                                             columns=['categorical'])
        features = pd.concat([units, types, categorical], axis=1).sort_index()

        # Return features
        logger.info(
            "Found %i feature(s) (defined in 'prediction_input' data for "
            "prediction '%s')", len(features.index), pred_name_str)
        for feature in features.index:
            logger.debug("'%s' with units '%s' and type '%s'", feature,
                         features.units.loc[feature],
                         features.types.loc[feature])
        return features

    def _get_features_of_datasets(self, datasets, var_type, pred_name):
        """Extract all features (with units and types) of given datasets."""
        pred_name_str = self._get_name(pred_name)
        units = {}
        types = {}
        cube = None
        ref_cube = None
        for (tag, datasets_) in group_metadata(datasets, 'tag').items():
            dataset = datasets_[0]
            cube = self._load_cube(dataset)
            if 'broadcast_from' not in dataset:
                ref_cube = cube
            units[tag] = Unit(dataset['units'])
            if 'broadcast_from' in dataset:
                types[tag] = 'broadcasted'
            else:
                types[tag] = 'regular'

        # Check if reference cube was given
        if ref_cube is None:
            if cube is None:
                raise ValueError(
                    f"Expected at least one '{var_type}' dataset for "
                    f" prediction '{pred_name_str}'")
            raise ValueError(
                f"Expected at least one '{var_type}' dataset for prediction "
                f"'{pred_name_str}' without the option 'broadcast_from'")

        # Coordinate features
        for coord_name in self._cfg.get('coords_as_features', []):
            try:
                coord = ref_cube.coord(coord_name)
            except iris.exceptions.CoordinateNotFoundError:
                raise iris.exceptions.CoordinateNotFoundError(
                    f"Coordinate '{coord_name}' given in 'coords_as_features' "
                    f"not found in '{var_type}' data for prediction "
                    f"'{pred_name_str}'")
            units[coord_name] = coord.units
            types[coord_name] = 'coordinate'

        return (units, types)

    def _get_group_attributes(self):
        """Get all group attributes from ``label`` datasets."""
        logger.debug("Extracting group attributes from 'label' datasets")
        grouped_datasets = group_metadata(self._datasets['label'],
                                          'group_attribute',
                                          sort=True)
        group_attributes = list(grouped_datasets.keys())
        if group_attributes == [None]:
            logger.debug("No group attributes given")
        else:
            logger.info(
                "Found %i group attribute(s) (defined in 'label' data)",
                len(group_attributes))
            logger.debug(pformat(group_attributes))
        return np.array(group_attributes)

    def _get_label(self):
        """Extract label from training data."""
        logger.debug("Extracting label from training datasets")
        grouped_datasets = group_metadata(self._datasets['label'], 'tag')
        labels = list(grouped_datasets.keys())
        if len(labels) > 1:
            raise ValueError(f"Expected unique label tag, got {labels}")
        units = Unit(self._datasets['label'][0]['units'])
        logger.info(
            "Found label '%s' with units '%s' (defined in 'label' "
            "data)", labels[0], units)
        label = pd.DataFrame.from_dict({labels[0]: units},
                                       orient='index',
                                       columns=['units'])
        return label

    def _get_lime_feature_importance(self, x_pred):
        """Get most important feature given by LIME."""
        logger.info(
            "Calculating global feature importance using LIME (this may take "
            "a while...)")
        x_pred = self._impute_nans(x_pred)

        # Most important feature for single input
        def _most_important_feature(x_single_pred, explainer, predict_fn):
            """Get most important feature for single input."""
            explanation = explainer.explain_instance(x_single_pred, predict_fn)
            local_exp = explanation.local_exp[1]
            sorted_exp = sorted(local_exp, key=lambda elem: elem[0])
            norm = sum([abs(elem[1]) for elem in sorted_exp])
            return [abs(elem[1]) / norm for elem in sorted_exp]

        # Apply on whole input (using multiple processes)
        parallel = Parallel(n_jobs=self._cfg['n_jobs'])
        lime_feature_importance = parallel(
            [delayed(_most_important_feature)(x,
                                              explainer=self._lime_explainer,
                                              predict_fn=self._clf.predict)
             for x in x_pred.values]
        )
        lime_feature_importance = np.array(lime_feature_importance,
                                           dtype=self._cfg['dtype'])
        lime_feature_importance = np.moveaxis(lime_feature_importance, -1, 0)
        lime_feature_importance = dict(zip(self.features,
                                           lime_feature_importance))
        return lime_feature_importance

    def _get_logo_cv_kwargs(self):
        """Get :class:`sklearn.model_selection.LeaveOneGroupOut` CV."""
        if not self._cfg['group_datasets_by_attributes']:
            raise ValueError(
                "Cannot create 'LeaveOneGroupOut' CV splitter, "
                "'group_datasets_by_attributes' was not given during "
                "class initialization")
        kwargs = {
            'cv': LeaveOneGroupOut(),
            'groups': self.data['train'].y.index.get_level_values(0).values,
        }
        return kwargs

    def _get_mask(self, x_data, data_type):
        """Get mask for missing features."""
        x_regular = x_data[self.features[self.features_types == 'regular']]

        # Get points where no regular feature is given
        mask = x_regular.isnull().all(axis=1).values
        logger.debug(
            "Removing %i %s point(s) where all regular features are missing",
            mask.sum(), data_type)

        # Get other missing points if desired
        if self._cfg['imputation_strategy'] == 'remove':
            mask = x_data.isnull().any(axis=1).values
            logger.debug(
                "Removing total %i %s point(s) where at least one feature is "
                "missing (because imputation_strategy = 'remove')", mask.sum(),
                data_type)

        return mask

    def _get_multiindex(self, ref_cube, group_attr=None):
        """Get :class:`pandas.MultiIndex` for data."""
        group_attr = self._group_attr_to_pandas_index_str(group_attr)
        index = pd.MultiIndex.from_product(
            [[group_attr], np.arange(ref_cube.data.size)],
            names=self._get_multiindex_names(),
        )
        return index

    def _get_multiindex_names(self):
        """Get names for :class:`pandas.MultiIndex` for data."""
        return ['-'.join(self._cfg['group_datasets_by_attributes']), 'index']

    def _get_plot_feature(self, feature):
        """Get :obj:`str` of selected ``feature`` and respective units."""
        units = self._get_plot_units(self.features_units[feature])
        return f'{feature} [{units}]'

    def _get_plot_label(self):
        """Get :obj:`str` of label and respective units."""
        return f'{self.label} [{self._get_plot_units(self.label_units)}]'

    def _get_plot_units(self, units):
        """Get plot units version of specified ``units``."""
        return self._cfg['plot_units'].get(str(units), str(units))

    def _get_prediction_dict(self, pred_name, x_pred, x_err, y_ref,
                             get_mlr_model_error=None,
                             get_lime_importance=False,
                             get_propagated_errors=False, **kwargs):
        """Get prediction output in a dictionary."""
        logger.info("Predicting %i point(s)", len(x_pred.index))
        y_preds = self._clf.predict(x_pred, **kwargs)
        pred_dict = self._prediction_to_dict(y_preds, **kwargs)

        # Estimate error of MLR model itself
        if get_mlr_model_error:
            pred_dict['squared_mlr_model_error_estim'] = (
                self._estimate_mlr_model_error(len(x_pred.index),
                                               get_mlr_model_error))

        # LIME feature importance
        if get_lime_importance:
            lime_importance = self._get_lime_feature_importance(x_pred)
            for (feature, importance) in lime_importance.items():
                pred_dict[f'lime_importance___{feature}'] = importance

        # Propagate prediction input errors
        if get_propagated_errors:
            if x_err is None:
                raise ValueError(
                    f"'save_propagated_errors' is not possible because no "
                    f"'prediction_input_error' data for prediction "
                    f"'{self._get_name(pred_name)}' is available")
            pred_dict['squared_propagated_input_error'] = (
                self._propagate_input_errors(x_pred, x_err))

        # Calculate residuals relative to reference if possible
        if y_ref is not None:
            y_ref = y_ref.values
            if y_ref.ndim == 2 and y_ref.shape[1] == 1:
                y_ref = np.squeeze(y_ref, axis=1)
            pred_dict['residual'] = self._get_residuals(y_ref, pred_dict[None])

        # Return dictionary
        for pred_type in pred_dict:
            if pred_type is not None:
                logger.debug("Found additional prediction type '%s'",
                             pred_type)
        logger.info(
            "Successfully created prediction array(s) with %i point(s)",
            pred_dict[None].size)
        return pred_dict

    def _get_prediction_dtype(self):
        """Get ``dtype`` of the output of final regressor's ``predict()``."""
        x_data = self.get_x_array('train')[0].reshape(1, -1)
        y_pred = self._clf.predict(x_data)
        return y_pred.dtype

    def _get_prediction_properties(self):
        """Get important properties of prediction input."""
        properties = {}
        for attr in ('dataset', 'exp', 'project', 'start_year', 'end_year'):
            attrs = list(group_metadata(self._datasets['label'], attr).keys())
            properties[attr] = attrs[0]
            if len(attrs) > 1:
                if attr == 'start_year':
                    properties[attr] = min(attrs)
                elif attr == 'end_year':
                    properties[attr] = max(attrs)
                else:
                    properties[attr] = '|'.join(sorted(attrs))
                logger.debug(
                    "Attribute '%s' of label data is not unique, got values "
                    "%s, using '%s' for prediction cubes", attr, attrs,
                    properties[attr])
        return properties

    def _get_reference_cube(self, datasets, var_type, text=None):
        """Get reference cube for ``datasets``."""
        msg = '' if text is None else text
        regular_features = self.features[self.features_types == 'regular']

        for tag in regular_features:
            dataset = self._check_dataset(datasets, var_type, tag, msg)
            if dataset is not None:
                ref_cube = self._load_cube(dataset)
                logger.debug(
                    "For var_type '%s'%s, use reference cube with tag '%s'",
                    var_type, msg, tag)
                logger.debug(ref_cube.summary(shorten=True))
                return ref_cube
        raise ValueError(f"No {var_type} data{msg} without the option "
                         f"'broadcast_from' found")

    def _get_sample_weights(self, data_type):
        """Get sample weights of desired data."""
        data_frame = self.data[data_type]
        if 'sample_weight' not in data_frame:
            return None
        return data_frame.sample_weight.squeeze().values

    def _get_verbosity_parameters(self, function, boolean=False):
        """Get verbosity parameters for class initialization."""
        verbosity_params = {
            'silent': {
                'debug': False,
                'info': False,
                'default': True,
            },
            'verbose': {
                'debug': 1,
                'info': 0,
                'default': 0,
            },
            'verbosity': {
                'debug': 2,
                'info': 1,
                'default': 0,
            },
        }
        parameters = {}
        for (param, log_levels) in verbosity_params.items():
            if param in getfullargspec(function).args:
                parameters[param] = log_levels.get(self._cfg['log_level'],
                                                   log_levels['default'])
                if boolean:
                    parameters[param] = bool(parameters[param])
                logger.debug("Set verbosity parameter '%s' of %s to '%s'",
                             param, str(function), parameters[param])
        return parameters

    def _get_x_data_for_group(self, datasets, var_type, group_attr=None):
        """Get x data for a group of datasets."""
        msg = '' if group_attr is None else f" for '{group_attr}'"
        ref_cube = self._get_reference_cube(datasets, var_type, msg)
        group_data = pd.DataFrame(
            columns=self.features,
            index=self._get_multiindex(ref_cube, group_attr=group_attr),
            dtype=self._cfg['dtype'],
        )
        sample_weights = self._calculate_sample_weights(ref_cube,
                                                        var_type,
                                                        group_attr=group_attr)

        # Iterate over all features
        for tag in self.features:
            if self.features_types[tag] != 'coordinate':
                dataset = self._check_dataset(datasets, var_type, tag, msg)

                # No dataset found
                if dataset is None:
                    if var_type == 'prediction_input_error':
                        logger.debug(
                            "Prediction input error of '%s'%s not available, "
                            "setting it to 0.0", tag, msg)
                        new_data = 0.0
                    else:
                        new_data = np.nan

                # Found exactly one dataset
                else:
                    text = f"{var_type} '{tag}'{msg}"

                    # Broadcast if necessary
                    if 'broadcast_from' in dataset:
                        cube = self._get_broadcasted_cube(
                            dataset, ref_cube, text)
                    else:
                        cube = self._load_cube(dataset)
                    self._check_cube_dimensions(cube, ref_cube, text)

                    # Do not accept errors for categorical features
                    if (var_type == 'prediction_input_error'
                            and tag in self.categorical_features):
                        raise ValueError(
                            f"Specifying prediction input error for "
                            f"categorical feature '{tag}'{msg} is not "
                            f"possible")
                    new_data = self._get_cube_data(cube)

            # Load coordinate feature data
            else:
                new_data = self._get_coordinate_data(ref_cube, var_type, tag,
                                                     msg)

            # Save data
            new_data = np.array(new_data)
            if new_data.size != ref_cube.data.size:
                new_data = np.broadcast_to(new_data, (ref_cube.data.size,))
            group_data[tag] = new_data

        # Return data and reference cube
        logger.debug("Found %i raw '%s' input data points%s",
                     len(group_data.index), var_type, msg)
        return (group_data, ref_cube, sample_weights)

    def _group_by_attributes(self, datasets):
        """Group datasets by specified attributes."""
        attributes = self._cfg['group_datasets_by_attributes']
        if not attributes:
            if self._cfg.get('accept_only_scalar_data'):
                attributes = ['dataset']
                logger.warning("Automatically set 'group_datasets_by_'"
                               "attributes' to ['dataset'] because 'accept_"
                               "only_scalar_data' is given")
            else:
                for dataset in datasets:
                    dataset['group_attribute'] = None
                return datasets
        for dataset in datasets:
            dataset['group_attribute'] = mlr.create_alias(dataset, attributes)
        logger.info("Grouped feature and label datasets by %s", attributes)
        return datasets

    def _impute_nans(self, data_frame, copy=True):
        """Impute all nans of a given :class:`pandas.DataFrame`."""
        if copy:
            data_frame = data_frame.copy()
        if 'feature_selection' in self._clf.named_steps:
            support = self._clf.named_steps['feature_selection'].support
        else:
            support = None
        if 'imputer' in self._clf.named_steps:
            transform = self._clf.named_steps['imputer'].transform
            if 'x' in data_frame.columns:
                if support is not None:
                    data_frame.x.values[:, support] = transform(
                        data_frame.x.values[:, support])
                    data_frame = data_frame.fillna(data_frame.mean())
                else:
                    data_frame.x.values[:] = transform(data_frame.x.values)
            else:
                if support is not None:
                    data_frame.values[:, support] = transform(
                        data_frame.values[:, support])
                    data_frame = data_frame.fillna(data_frame.mean())
                else:
                    data_frame.values[:] = transform(data_frame.values)
        return data_frame

    def _is_ready_for_plotting(self):
        """Check if the class is ready for plotting."""
        self._check_fit_status('Plotting')
        return True

    def _load_classes(self):
        """Populate :attribute:`_classes` and check for errors."""
        self._classes['group_attributes'] = self._get_group_attributes()
        self._classes['features'] = self._get_features()
        self._classes['label'] = self._get_label()

    def _load_cube(self, dataset):
        """Load iris cube, check data type and convert units if desired."""
        logger.debug("Loading %s", dataset['filename'])
        cube = iris.load_cube(dataset['filename'])

        # Check dtype
        if not np.issubdtype(cube.dtype, np.number):
            raise TypeError(
                f"Data type of cube loaded from '{dataset['filename']}' is "
                f"'{cube.dtype}', at the moment only numeric data is "
                f"supported")

        # Convert dtypes
        cube.data = cube.core_data().astype(self._cfg['dtype'],
                                            casting='same_kind')
        for coord in cube.coords():
            try:
                coord.points = coord.points.astype(self._cfg['dtype'],
                                                   casting='same_kind')
            except TypeError:
                logger.debug(
                    "Cannot convert dtype of coordinate array '%s' from '%s' "
                    "to '%s'", coord.name(), coord.points.dtype,
                    self._cfg['dtype'])

        # Convert and check units
        if dataset.get('convert_units_to'):
            self._convert_units_in_cube(cube, dataset['convert_units_to'])
        if not cube.units == Unit(dataset['units']):
            raise ValueError(
                f"Units of cube '{dataset['filename']}' for "
                f"{dataset['var_type']} '{dataset['tag']}' differ from units "
                f"given in dataset list, got '{cube.units}' in cube and "
                f"'{dataset['units']}' in dataset list")
        return cube

    def _load_data(self):
        """Load train/test data (features/labels)."""
        (x_all, y_all, sample_weights) = self._extract_features_and_labels()

        # Normalize and add sample weights if necessary
        objs = [x_all, y_all]
        keys = ['x', 'y']
        if sample_weights is not None:
            sample_weights /= sample_weights.mean()
            objs.append(sample_weights)
            keys.append('sample_weight')

        # Save complete data
        self._data['all'] = pd.concat(objs, axis=1, keys=keys)
        if len(y_all.index) < 2:
            raise ValueError(
                f"Need at least 2 data points for MLR training, got only "
                f"{len(y_all.index)}")
        logger.info("Loaded %i input data point(s)", len(y_all.index))

        # Split train/test data if desired
        test_size = self._cfg['test_size']
        if test_size:
            (self._data['train'],
             self._data['test']) = train_test_split(self._data['all'].copy(),
                                                    test_size=test_size)
            self._data['train'] = self._data['train'].sort_index()
            self._data['test'] = self._data['test'].sort_index()
            for data_type in ('train', 'test'):
                if len(self.data[data_type].index) < 2:
                    raise ValueError(
                        f"Need at least 2 datasets for '{data_type}' data, "
                        f"got {len(self.data[data_type].index)}")
            logger.info(
                "Using %i%% of the input data as test data (%i point(s))",
                int(test_size * 100), len(self.data['test'].index))
            logger.info("%i point(s) remain(s) for training",
                        len(self.data['train'].index))
        else:
            self._data['train'] = self.data['all'].copy()
            logger.info("Using all %i input data point(s) for training",
                        len(y_all.index))

    def _load_final_parameters(self):
        """Load parameters for final regressor."""
        parameters = self._cfg.get('parameters_final_regressor', {})
        logger.debug("Using parameter(s) for final regressor: %s", parameters)
        verbosity_params = self._get_verbosity_parameters(self._CLF_TYPE)
        for (param, verbosity) in verbosity_params.items():
            parameters.setdefault(param, verbosity)
        return parameters

    def _load_input_datasets(self, input_datasets):
        """Load input datasets."""
        input_datasets = deepcopy(input_datasets)

        # Catch invalid var_types
        if not mlr.datasets_have_mlr_attributes(
                input_datasets, log_level='error', mode='only_var_type'):
            raise ValueError("Data with invalid 'var_type' given")

        # Training datasets
        feature_datasets = select_metadata(input_datasets, var_type='feature')
        label_datasets = select_metadata(input_datasets, var_type='label')

        # Prediction datasets
        pred_in_datasets = select_metadata(input_datasets,
                                           var_type='prediction_input')
        pred_in_err_datasets = select_metadata(
            input_datasets, var_type='prediction_input_error')
        pred_ref_datasets = select_metadata(input_datasets,
                                            var_type='prediction_reference')

        # Check datasets
        msg = ("At least one '{}' dataset does not have necessary MLR "
               "attributes")
        datasets_to_check = {
            'feature': feature_datasets,
            'label': label_datasets,
            'prediction_input': pred_in_datasets,
            'prediction_input_error': pred_in_err_datasets,
            'prediction_reference': pred_ref_datasets,
        }
        for (label, datasets) in datasets_to_check.items():
            if not mlr.datasets_have_mlr_attributes(datasets,
                                                    log_level='error'):
                raise ValueError(msg.format(label))

        # Check if data was found
        if not feature_datasets:
            raise ValueError("No 'feature' data found")
        if not label_datasets:
            raise ValueError("No 'label' data found")
        if not pred_in_datasets:
            raise ValueError("No 'prediction_input' data found")

        # Convert units
        self._convert_units_in_metadata(feature_datasets)
        self._convert_units_in_metadata(label_datasets)
        self._convert_units_in_metadata(pred_in_datasets)
        self._convert_units_in_metadata(pred_in_err_datasets)
        self._convert_units_in_metadata(pred_ref_datasets)

        # Save datasets
        logger.info(
            "Found %i 'feature' dataset(s), %i 'label' dataset(s), %i "
            "'prediction_input' dataset(s), %i 'prediction_input_error' "
            "dataset(s) and %i 'prediction_reference' datasets(s)",
            len(feature_datasets), len(label_datasets), len(pred_in_datasets),
            len(pred_in_err_datasets), len(pred_ref_datasets))
        labeled_datasets = {
            'Feature': feature_datasets,
            'Label': label_datasets,
            'Prediction input': pred_in_datasets,
            'Prediction input error': pred_in_err_datasets,
            'Prediction output': pred_ref_datasets,
        }
        for (msg, datasets) in labeled_datasets.items():
            logger.debug("%s datasets:", msg)
            logger.debug(pformat([d['filename'] for d in datasets]))
        self._datasets['feature'] = self._group_by_attributes(feature_datasets)
        self._datasets['label'] = self._group_by_attributes(label_datasets)
        self._datasets['prediction_input'] = self._group_prediction_datasets(
            pred_in_datasets)
        self._datasets['prediction_input_error'] = (
            self._group_prediction_datasets(pred_in_err_datasets))
        self._datasets['prediction_reference'] = (
            self._group_prediction_datasets(pred_ref_datasets))

    def _load_lime_explainer(self):
        """Load :class:`lime.lime_tabular.LimeTabularExplainer`."""
        x_train = self.get_x_array('train', impute_nans=True)
        y_train = self.get_y_array('train', impute_nans=True)
        verbosity = self._get_verbosity_parameters(LimeTabularExplainer,
                                                   boolean=True)
        for param in verbosity:
            verbosity[param] = False
        categorical_features_idx = [
            int(np.where(self.features == tag)[0][0])
            for tag in self.categorical_features
        ]
        self._lime_explainer = LimeTabularExplainer(
            x_train,
            mode='regression',
            training_labels=y_train,
            feature_names=self.features,
            categorical_features=categorical_features_idx,
            discretize_continuous=False,
            sample_around_instance=True,
            **verbosity,
        )
        logger.debug(
            "Loaded %s with new training data", str(LimeTabularExplainer))

    def _mask_prediction_array(self, y_pred, ref_cube):
        """Apply mask of reference cube to prediction array."""
        mask = np.ma.getmaskarray(ref_cube.data).ravel()
        if y_pred.ndim == 1 and y_pred.shape[0] != mask.shape[0]:
            new_y_pred = np.empty(mask.shape[0], dtype=self._cfg['dtype'])
            new_y_pred[mask] = np.nan
            new_y_pred[~mask] = y_pred
        else:
            new_y_pred = y_pred
        return np.ma.masked_invalid(new_y_pred)

    def _plot_feature_importance(self, feature_importance_dict, colors,
                                 plot_path):
        """Plot feature importance."""
        logger.info("Plotting feature importance")
        (_, axes) = plt.subplots()

        # Sort data and get position of bars
        features = np.array(list(feature_importance_dict.keys()))
        feature_importances = np.array(list(feature_importance_dict.values()))
        sorted_idx = np.argsort(feature_importances)
        pos = np.arange(sorted_idx.shape[0]) + 0.5

        # Write cube with feature importance for provenance tracking
        ancestors = self.get_ancestors(prediction_names=[])
        cube = mlr.get_1d_cube(
            features,
            feature_importances,
            x_kwargs={'var_name': 'feature',
                      'long_name': 'Feature name',
                      'units': 'no unit'},
            y_kwargs={'var_name': 'feature_importance',
                      'long_name': 'Relative Feature Importance',
                      'units': '1',
                      'attributes': {'project': '', 'dataset': ''}},
        )

        # Plot
        for (idx, importance) in enumerate(feature_importances[sorted_idx]):
            feature = features[sorted_idx][idx]
            axes.barh(pos[idx], importance, align='center',
                      color=colors[feature])

        # Plot appearance
        axes.tick_params(axis='y', which='minor', left=False, right=False)
        axes.tick_params(axis='y', which='major', left=True, right=False)
        title = f"Global feature importance ({self._cfg['mlr_model_name']})"
        axes.set_title(title)
        axes.set_xlabel('Relative Importance')
        axes.set_yticks(pos)
        axes.set_yticklabels(features[sorted_idx])

        # Save plot and provenance
        plt.savefig(plot_path, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path)
        self._write_plot_provenance(cube, plot_path, ancestors=ancestors,
                                    caption=title + '.', plot_types=['bar'])

        # Save additional plot with logarithmic X axis
        axes.set_xscale('log')
        axes.xaxis.set_major_formatter(ScalarFormatter())
        ext = os.path.splitext(plot_path)[1]
        plot_path_log = plot_path.replace(ext, f'_log{ext}')
        plt.savefig(plot_path_log, **self._cfg['savefig_kwargs'])
        logger.info("Wrote %s", plot_path_log)
        self._write_plot_provenance(cube, plot_path_log, ancestors=ancestors,
                                    caption=title + '.', plot_types=['bar'])
        plt.close()

    def _prediction_to_dict(self, pred_out, **kwargs):
        """Convert output of final regressor's ``predict()`` to :obj:`dict`."""
        if not isinstance(pred_out, (list, tuple)):
            pred_out = [pred_out]
        idx_to_name = {0: None}
        if 'return_var' in kwargs:
            idx_to_name[1] = 'var'
        elif 'return_cov' in kwargs:
            idx_to_name[1] = 'cov'
        pred_dict = {}
        for (idx, pred) in enumerate(pred_out):
            pred = pred.astype(self._cfg['dtype'], casting='same_kind')
            if pred.ndim == 2 and pred.shape[1] == 1:
                logger.warning(
                    "Prediction output is 2D and length of second axis is 1, "
                    "squeezing second axis")
                pred = np.squeeze(pred, axis=1)
            pred_dict[idx_to_name.get(idx, idx)] = pred
        return pred_dict

    def _pred_type_to_metadata(self, pred_type, cube):
        """Get correct :mod:`iris.cube.CubeMetadata` of prediction cube."""
        standard_name = cube.standard_name
        var_name = cube.var_name
        long_name = cube.long_name
        units = cube.units
        attributes = cube.attributes
        suffix = '' if pred_type is None else f'_{pred_type}'
        error_types = {
            'var': ' (variance)',
            'cov': ' (covariance)',
            'squared_mlr_model_error_estim': (' (squared MLR model error '
                                              'estimation using hold-out test '
                                              'data set)'),
            'squared_propagated_input_error': (' (squared propagated error of '
                                               'prediction input estimated by '
                                               'LIME)'),
        }
        if pred_type is None:
            attributes['var_type'] = 'prediction_output'
        elif isinstance(pred_type, int):
            var_name += '_{:d}'.format(pred_type)
            long_name += ' {:d}'.format(pred_type)
            logger.warning("Got unknown prediction type with index %i",
                           pred_type)
            attributes['var_type'] = 'prediction_output_misc'
        elif pred_type in error_types:
            var_name += suffix
            long_name += error_types[pred_type]
            units = mlr.units_power(cube.units, 2)
            attributes['var_type'] = 'prediction_output_error'
            attributes['squared'] = 1
        elif 'lime_importance___' in pred_type:
            standard_name = None
            feature = pred_type.replace('lime_importance___', '')
            var_name = f'importance_of_feature_{feature}'
            long_name = (f'Local importance of feature {feature} for '
                         f'predicting {self.label} given by LIME')
            units = Unit('1')
            attributes['var_type'] = 'prediction_output_misc'
        elif pred_type == 'residual':
            var_name += suffix
            long_name += ' (residual)'
            attributes['residual'] = 'true minus predicted values'
            attributes['var_type'] = 'prediction_residual'
        else:
            raise ValueError(f"Got unknown prediction type '{pred_type}'")
        return iris.cube.CubeMetadata(
            standard_name=standard_name,
            long_name=long_name,
            var_name=var_name,
            units=units,
            attributes=attributes,
            cell_methods=cube.cell_methods,
        )

    def _print_metrics(self, regression_metrics, data_type):
        """Print regression metrics."""
        if data_type not in self.data:
            return
        logger.info("Evaluating regression metrics for %s data", data_type)
        x_data = self.data[data_type].x
        y_true = self.get_y_array(data_type)
        y_pred = self._clf.predict(x_data)
        sample_weights = self._get_sample_weights(data_type)
        for metric in regression_metrics:
            metric_function = getattr(metrics, metric)
            value = metric_function(y_true, y_pred)
            if 'squared' in metric:
                value = np.sqrt(value)
                metric = f'root_{metric}'
            logger.info("%s: %s", metric, value)
        if sample_weights is None:
            return
        for metric in regression_metrics:
            metric_function = getattr(metrics, metric)
            value = metric_function(y_true, y_pred,
                                    sample_weight=sample_weights)
            if 'squared' in metric:
                value = np.sqrt(value)
                metric = f'root_{metric}'
            logger.info("Weighted %s: %s", metric, value)

    def _propagate_input_errors(self, x_pred, x_err):
        """Propagate errors from prediction input."""
        logger.info(
            "Propagating prediction input errors using LIME (this may take a "
            "while...)")
        if 'feature_selection' in self._clf.named_steps:
            logger.warning(
                "Propagating input errors might not work correctly when a "
                "'feature_selection' step is present (usually because of "
                "calling rfecv())")
        x_pred = self._impute_nans(x_pred)

        # Propagated error for single input
        def _propagated_error(x_single_pred, x_single_err, explainer,
                              predict_fn, features, categorical_features):
            """Get propagated prediction input error for single input."""
            exp = explainer.explain_instance(x_single_pred, predict_fn)
            x_single_err = np.nan_to_num(x_single_err)
            x_err_scaled = x_single_err / explainer.scaler.scale_
            squared_error = 0.0
            for (idx, coef) in exp.local_exp[1]:
                if features[idx] in categorical_features:
                    continue
                squared_error += (x_err_scaled[idx] * coef)**2
            return squared_error

        # Apply on whole input (using multiple processes)
        parallel = Parallel(n_jobs=self._cfg['n_jobs'])
        errors = parallel(
            [delayed(_propagated_error)(
                x, x_e, explainer=self._lime_explainer,
                predict_fn=self._clf.predict,
                features=self.features,
                categorical_features=self.categorical_features,
            ) for (x, x_e) in zip(x_pred.values, x_err.values)]
        )
        return np.array(errors, dtype=self._cfg['dtype'])

    def _remove_missing_features(self, x_data, y_data, sample_weights):
        """Remove missing values in the features data (if desired)."""
        mask = self._get_mask(x_data, 'training')
        x_data = x_data[~mask]
        y_data = y_data[~mask]
        if sample_weights is not None:
            sample_weights = sample_weights[~mask]
        diff = mask.sum()
        if diff:
            msg = ('Removed %i training point(s) where features were '
                   'missing')
            if self._cfg.get('accept_only_scalar_data'):
                removed_groups = self.group_attributes[mask]
                msg += f' ({removed_groups})'
                self._classes['group_attributes'] = (
                    self.group_attributes[~mask])
            logger.info(msg, diff)
        return (x_data, y_data, sample_weights)

    def _remove_missing_pred_input(self, x_pred, x_err=None, y_ref=None):
        """Remove missing values in the prediction input data."""
        mask = self._get_mask(x_pred, 'prediction input')
        x_pred = x_pred[~mask]
        if x_err is not None:
            x_err = x_err[~mask]
        if y_ref is not None:
            y_ref = y_ref[~mask]
        diff = mask.sum()
        if diff:
            logger.info(
                "Removed %i prediction input point(s) where features were "
                "missing", diff)
        return (x_pred, x_err, y_ref, mask)

    def _save_prediction_cubes(self, pred_dict, pred_name, x_cube):
        """Save (multi-dimensional) prediction output."""
        logger.debug("Creating output cubes")
        for (pred_type, y_pred) in pred_dict.items():
            y_pred = self._mask_prediction_array(y_pred, x_cube)
            if y_pred.size == np.prod(x_cube.shape, dtype=np.int):
                pred_cube = x_cube.copy(y_pred.reshape(x_cube.shape))
            else:
                dim_coords = []
                for (dim_idx, dim_size) in enumerate(y_pred.shape):
                    dim_coords.append((iris.coords.DimCoord(
                        np.arange(dim_size, dtype=np.float64),
                        long_name=f'MLR prediction index {dim_idx}',
                        var_name=f'idx_{dim_idx}'), dim_idx))
                pred_cube = iris.cube.Cube(y_pred,
                                           dim_coords_and_dims=dim_coords)
            new_path = self._set_prediction_cube_attributes(
                pred_cube, pred_type, pred_name=pred_name)
            io.iris_save(pred_cube, new_path)

            # Save provenance
            ancestors = self.get_ancestors(
                prediction_names=[pred_name],
                prediction_reference=pred_type == 'residual')
            record = {
                'ancestors': ancestors,
                'authors': ['schlund_manuel'],
                'caption': (f"{pred_cube.long_name} of MLR model "
                            f"{self._cfg['mlr_model_name']} for prediction "
                            f"{pred_name}."),
                'references': ['schlund20jgr'],
            }
            with ProvenanceLogger(self._cfg) as provenance_logger:
                provenance_logger.log(new_path, record)

    def _save_csv_file(self, data_type, filename, pred_name=None):
        """Save CSV file."""
        if data_type not in self.data:
            return
        if data_type == 'pred':
            csv_data = self.data[data_type][pred_name]
        else:
            csv_data = self.data[data_type]

        # Filename and path
        if filename is None:
            if data_type == 'pred':
                filename = '{data_type}_{pred_name}.csv'
                format_kwargs = {
                    'data_type': data_type,
                    'pred_name': self._get_name(pred_name),
                }
            else:
                filename = '{data_type}.csv'
                format_kwargs = {'data_type': data_type}
        filename = filename.format(**format_kwargs)
        path = os.path.join(self._cfg['mlr_work_dir'], filename)

        # Save file
        csv_data.to_csv(path, na_rep='nan')
        logger.info("Wrote %s", path)

    def _set_default_settings(self):
        """Set default (non-``False``) keyword arguments."""
        self._cfg.setdefault('weighted_samples',
                             {'area_weighted': True, 'time_weighted': True})
        self._cfg.setdefault('cache_intermediate_results', True)
        self._cfg.setdefault('dtype', 'float64')
        self._cfg.setdefault('fit_kwargs', {})
        self._cfg.setdefault('group_datasets_by_attributes', [])
        self._cfg.setdefault('imputation_strategy', 'remove')
        self._cfg.setdefault('log_level', 'info')
        self._cfg.setdefault('mlr_model_name', f'{self._CLF_TYPE} model')
        self._cfg.setdefault('n_jobs', 1)
        self._cfg.setdefault('output_file_type', 'png')
        self._cfg.setdefault('parameters', {})
        self._cfg.setdefault('plot_dir',
                             os.path.expanduser(os.path.join('~', 'plots')))
        self._cfg.setdefault('plot_units', {})
        self._cfg.setdefault('savefig_kwargs', {
            'bbox_inches': 'tight',
            'dpi': 300,
            'orientation': 'landscape',
        })
        self._cfg.setdefault('standardize_data', True)
        self._cfg.setdefault('sub_dir', '')
        self._cfg.setdefault('test_size', 0.25)
        self._cfg.setdefault('work_dir',
                             os.path.expanduser(os.path.join('~', 'work')))
        logger.info("Using imputation strategy '%s'",
                    self._cfg['imputation_strategy'])
        if self._cfg['fit_kwargs']:
            logger.info(
                "Using additional keyword argument(s) %s for fit() function",
                self._cfg['fit_kwargs'])

    def _set_prediction_cube_attributes(self, cube, pred_type, pred_name=None):
        """Set the attributes of the prediction cube."""
        cube.cell_methods = None
        cube.attributes = {
            'description': 'MLR model prediction',
            'mlr_model_name': self._cfg['mlr_model_name'],
            'mlr_model_type': self.mlr_model_type,
            'final_regressor': str(self._CLF_TYPE),
            'prediction_name': self._get_name(pred_name),
            'tag': self.label,
        }
        cube.attributes.update(self._get_prediction_properties())
        for (key, val) in self.parameters.items():
            cube.attributes[key] = str(val)
        cube.attributes['mlr_parameters'] = list(self.parameters.keys())
        label_cube = self._load_cube(self._datasets['label'][0])
        for attr in ('standard_name', 'var_name', 'long_name', 'units'):
            setattr(cube, attr, getattr(label_cube, attr))

        # Modify cube metadata depending on prediction type
        cube.metadata = self._pred_type_to_metadata(pred_type, cube)

        # Get new path
        suffix = '' if pred_type is None else f'_{pred_type}'
        pred_str = f'_for_prediction_{self._get_name(pred_name)}'
        sub_str = ('' if self._cfg['sub_dir'] == '' else
                   f"_of_group_{self._cfg['sub_dir']}")
        filename = (f'{self.mlr_model_type}_{self.label}_prediction{suffix}'
                    f'{pred_str}{sub_str}.nc')
        new_path = os.path.join(self._cfg['mlr_work_dir'], filename)
        cube.attributes['filename'] = new_path
        return new_path

    def _update_fit_kwargs(self, fit_kwargs):
        """Check and update fit kwargs."""
        new_fit_kwargs = {}

        # Sort out wrong fit kwargs
        for (param_name, param_val) in fit_kwargs.items():
            step = param_name.split('__')[0]
            if step in self._clf.named_steps:
                new_fit_kwargs[param_name] = param_val
            else:
                raise ValueError(
                    f"Got invalid pipeline step '{step}' in fit parameter "
                    f"'{param_name}'")

        # Add sample weights if possible
        allowed_fit_kwargs = getfullargspec(self._CLF_TYPE.fit).args
        for kwarg in ('sample_weight', 'sample_weights'):
            if kwarg not in allowed_fit_kwargs:
                continue
            long_kwarg = f'{self._clf.steps[-1][0]}__regressor__{kwarg}'
            sample_weights = self._get_sample_weights('train')
            new_fit_kwargs[long_kwarg] = sample_weights
            if sample_weights is not None:
                logger.debug(
                    "Updated keyword arguments of final regressor's fit() "
                    "function with '%s'", kwarg)
            break

        return new_fit_kwargs

    def _write_plot_provenance(self, cube, plot_path, **additional_info):
        """Write provenance information for plots."""
        netcdf_path = mlr.get_new_path(self._cfg, plot_path)
        io.iris_save(cube, netcdf_path)
        record = {
            'authors': ['schlund_manuel'],
            'plot_file': plot_path,
            'references': ['schlund20jgr'],
            **additional_info,
        }
        with ProvenanceLogger(self._cfg) as provenance_logger:
            provenance_logger.log(netcdf_path, record)

    @staticmethod
    def _convert_units_in_cube(cube, new_units, power=None, text=None):
        """Convert units of cube if possible."""
        msg = '' if text is None else f' of {text}'
        if isinstance(new_units, str):
            new_units = Unit(new_units)
        if power:
            logger.debug("Raising target units of cube '%s' by power of %i",
                         cube.summary(shorten=True), power)
            new_units = mlr.units_power(new_units, power)
        logger.debug("Converting units%s from '%s' to '%s'", msg, cube.units,
                     new_units)
        try:
            cube.convert_units(new_units)
        except ValueError:
            raise ValueError(
                f"Cannot convert units{msg} from '{cube.units}' to "
                f"'{new_units}'")

    @staticmethod
    def _convert_units_in_metadata(datasets):
        """Convert units of datasets if desired."""
        for dataset in datasets:
            if not dataset.get('convert_units_to'):
                continue
            units_from = Unit(dataset['units'])
            units_to = Unit(dataset['convert_units_to'])
            try:
                units_from.convert(0.0, units_to)
            except ValueError:
                raise ValueError(
                    f"Cannot convert units of {dataset['var_type']} "
                    f"'{dataset['tag']}' from '{units_from}' to '{units_to}'")
            dataset['units'] = dataset['convert_units_to']

    @staticmethod
    def _get_centralized_bins(array, n_bins=None, ref=0.0):
        """Get bins for array centralized around a reference value."""
        diff = max([ref - array.min(), array.max() - ref])
        if n_bins is None:
            auto_bins = np.histogram_bin_edges(array)
            if len(auto_bins) < 2:
                raise ValueError(
                    f"Expected at least 2 bins, got {len(auto_bins):d}")
            delta = auto_bins[1] - auto_bins[0]
            n_bins = 2.0 * diff / delta
        if not n_bins % 2:
            n_bins += 1
        return np.linspace(ref - diff, ref + diff, n_bins + 1)

    @staticmethod
    def _get_coordinate_data(ref_cube, var_type, tag, text=None):
        """Get coordinate variable ``ref_cube`` which can be used as x data."""
        msg = '' if text is None else text
        if var_type == 'prediction_input_error':
            logger.debug(
                "Prediction input error of coordinate feature '%s'%s is set "
                "to 0.0", tag, msg)
            return 0.0
        try:
            coord = ref_cube.coord(tag)
        except iris.exceptions.CoordinateNotFoundError:
            raise iris.exceptions.CoordinateNotFoundError(
                f"Coordinate '{tag}' given in 'coords_as_features' not found "
                f"in reference cube for '{var_type}'{msg}")
        coord_array = np.ma.filled(coord.points, np.nan)
        coord_dims = ref_cube.coord_dims(coord)
        if coord_dims == ():
            logger.warning(
                "Coordinate '%s' is scalar, including it as feature does not "
                "add any information to the model (array is constant)", tag)
            coord_array = np.broadcast_to(coord_array, ref_cube.shape)
        else:
            coord_array = iris.util.broadcast_to_shape(coord_array,
                                                       ref_cube.shape,
                                                       coord_dims)
        logger.debug("Added %s coordinate '%s'%s", var_type, tag, msg)
        return coord_array.ravel()

    @staticmethod
    def _get_cube_data(cube):
        """Get data from cube."""
        cube_data = np.ma.filled(cube.data, np.nan)
        return cube_data.ravel()

    @staticmethod
    def _get_data_type_coord(data_types):
        """Get :class:`iris.coords.AuxCoord` ``data_type``."""
        aux_coord = iris.coords.AuxCoord(data_types,
                                         var_name='data_type',
                                         long_name='Data type',
                                         units='no unit')
        return aux_coord

    @staticmethod
    def _get_name(string):
        """Convert ``None`` to :obj:`str` if necessary."""
        return 'unnamed' if string is None else string

    @staticmethod
    def _get_plot_kwargs(data_type, plot_type=None):
        """Get plot kwargs for a data type."""
        plot_kwargs = {
            'all': {
                'color': 'r',
                'label': 'All data',
            },
            'train': {
                'color': 'b',
                'label': 'Train data',
            },
            'test': {
                'color': 'g',
                'label': 'Test data',
            },
        }
        allowed_data_types = list(plot_kwargs.keys())
        if data_type not in allowed_data_types:
            raise NotImplementedError(
                f"Plot kwargs for data type '{data_type}' not implemented "
                f"yet, only {allowed_data_types} are supported yet")
        kwargs = deepcopy(plot_kwargs[data_type])
        if plot_type == 'scatter':
            kwargs.update({'alpha': 0.5, 'marker': 'o', 's': 6})
        return kwargs

    @staticmethod
    def _get_residuals(y_true, y_pred):
        """Calculate residuals (true minus predicted values)."""
        logger.debug("Calculating residuals")
        return y_true - y_pred

    @staticmethod
    def _group_attr_to_pandas_index_str(group_attr):
        """Convert group attribute to :obj:`str` used in pandas index."""
        if group_attr is None:
            return 'none'
        return group_attr

    @staticmethod
    def _group_prediction_datasets(datasets):
        """Group prediction datasets (use ``prediction_name`` key)."""
        for dataset in datasets:
            dataset['group_attribute'] = None
        return group_metadata(datasets, 'prediction_name')

    @staticmethod
    def _remove_missing_labels(x_data, y_data, sample_weights):
        """Remove missing values in the label data."""
        mask = y_data.isnull().values
        x_data = x_data[~mask]
        y_data = y_data[~mask]
        if sample_weights is not None:
            sample_weights = sample_weights[~mask]
        diff = mask.sum()
        if diff:
            logger.info(
                "Removed %i training point(s) where labels were missing", diff)
        return (x_data, y_data, sample_weights)

    @staticmethod
    def _set_axis_lim_symmetric(axes, axis):
        """Make axis range of plot symmetric around 0."""
        if axis == 'x':
            getter = getattr(axes, 'get_xlim')
            setter = getattr(axes, 'set_xlim')
        elif axis == 'y':
            getter = getattr(axes, 'get_ylim')
            setter = getattr(axes, 'set_ylim')
        else:
            raise ValueError(f"Expected 'x' or 'y' for axis, got '{axis}'")
        maximum = np.max(np.abs(getter()))
        setter([-maximum, maximum])
