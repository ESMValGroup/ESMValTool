"""A collection of functions to calibrate the sigma values."""


def calibrate_independence_sigma(overall_independence: 'xr.DataArray', force: bool = False,
                                 **kwargs : dict) -> float:
    """
    Parameters
    ----------
    overall_independence : array_like, shape (N, N)
    force : bool, optional
        If force is False (default) this function will raise an error if the
        independence_contributions are selected such that that no clear
        separateion between dependent and independent models is possible. Set
        to True (not recommendet) to scip this check.

    Returns
    -------
    independence_sigma : float
    """
    # do stuff
    return independence_sigma


def remove_related_models(overall_performance, overall_independence, subset_size=20):
    """
    Remove (set to nan?) the models closest related to the perfect model for each perfect model
    so that there are the subset_size most independent models remaining
    Alternatively: remove random models to test sensitivity
    """
    # do stuff
    return overall_performance_subset


def evaluate_future(weights_sigma, target, metric):
    """
    For each perfect_model and each sigma compare the weighted mean to the 'truth'
    # NOTE: for each perfect model at least one model (the perfect model itself) but
    potentially more models (see remove_related_models) should have weight nan,
    this needs to be ignored (probably happens automaticall but just to keep in mind)
    """
    if metric == 'TEST_METRIC_1':
        return TEST_METRIC_1  # optimal sigma
    if metric == 'TEST_METRIC_2':
        return TEST_METRIC_2
    # ...


def calibrate_performance_sigma(
        models: dict, performance_contributions: list, variable_group_target: str
        overall_independence: 'xr.DataArray', independence_sigma: float,
        metric: str,
        use_independence_weighting: bool = True **kwargs: dict):

    # --- this is equivalent to climwip.py but with perfect models as reference instead of obs ---
    for variable_group in performance_contributions:

        logger.info('Reading model data for %s', variable_group)
        datasets_model = models[variable_group]
        model_data, model_data_files = read_model_data(datasets_model)

        # NOTE: it would probably be good to vectorize this
        for perfect_model in models:
            performance = calculate_performance(model_data, perfect_model)
            performances[variable_group][perfect_model] = performance

    overall_performance = compute_overall_mean(performance, performance_contributions)  # (N, N) matrix

    if cfg['combine_ensemble_members']:
        # this combines both dimensions -> I think this is what we want!!
         overall_performance, groups_performance = combine_ensemble_members(
             overall_performance)

    # --- from here it is not equivalent to climwip.py any more ---

    # NOTE: what do we do here if 'combine_ensemble_members' is False?
    overall_performance = remove_related_models(overall_performance, overall_independence)

    # NOTE: it would probably be good to vectorize this
    for performance_sigma in performance_sigma_list:  # we can probably hardcode this list?
        for perfect_model in models:
            if use_independence_weighting:
                weights_sigma = calculate_weights(
                    overall_performance.sel(perfect_model=perfect_model),
                    overall_independence,
                    performance_sigma, independence_sigma)  # (N, N) matrix
            else:
                weights_sigma = calculate_weights(
                    overall_performance.sel(perfect_model=perfect_model),
                    None, performance_sigma, None)  # (N, N) matrix

    target = models[variable_group_target]
    if cfg['combine_ensemble_members']:
        target = models[variable_group_target]

    model_data, model_data_files = read_model_data(datasets_model)
    performance_sigma = evaluate_future(weights_sigma, target, metric)

    return performance_sigma


def calibrate_sigmas(
        overall_independence: 'xr.DataArray', force_independence: bool = False,
        models: dict, performance_contributions: list, variable_group_target: str,
        use_independence_weighting: bool = True,
        **kwargs: dict):

    independence_sigma = calibrate_independence_sigma(overall_independence, force_independence)
    performance_sigma = calibrate_performance_sigma(
        models, performance_contributions, variable_group_target, overall_independence,
        independence_sigma, metric, use_independence_weighting=True, **kwargs)

    # potental cross checks, sanity tests, perfect model skill calculation
