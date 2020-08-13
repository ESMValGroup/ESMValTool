"""
Implementation of step iii and iv of the climwip weighting scheme

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
from scipy import stats
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata
import numpy as np
import xarray as xr
from scipy.spatial.distance import pdist, squareform
from collections import defaultdict



# (iii) The point-to-point distance between model and the center of the observational spread is calculated
def model_performance(model, obs):
    return stats.pearsonr(model, obs)

# (iv) The area-weighted root mean squared error is calculated over the selected region
def calculate_weights(quality, independence, sigma_q, sigma_i):
    """Calculates the (NOT normalised) weights for each model N.
    Parameters
    ----------
    quality : array_like, shape (N,)
        Array specifying the model quality.
    independence : array_like, shape (N, N)
        Array specifying the model independence.
    sigma_q : float
        Sigma value defining the form of the weighting function for the quality.
    sigma_i : float
        Sigma value defining the form of the weighting function for the independence.
    Returns
    -------
    numerator, denominator : ndarray, shape (N,)
    """
    numerator = np.exp(-((quality/sigma_q)**2))
    exp = np.exp(-((independence/sigma_i)**2))
    sum_exp = [np.sum(np.delete(ee, ii)) for ii, ee in enumerate(exp)]  # sum i!=j
    denominator = 1 + np.array(sum_exp)

    if sigma_i == -99.:
        denominator = (denominator * 0) + 1  # set to 1 (except NaN)
    if sigma_q == -99.:
        numerator = (numerator * 0) + 1  # set to 1 (except NaN)

    return numerator, denominator

def calculate_weights_sigmas(distances, sigmas_q, sigmas_i):
    """Calculates the weights for each model N and combination of sigma values.
    Parameters
    ----------
    distances : array_like, shape (N, N)
        Array specifying the distances between each model.
    sigmas_q : array_like, shape (M,)
        Array of sigma values for the weighting function of the quality.
    sigmas_i : array_like, shape (L,)
        Array of sigma values for the weighting function of the independence.
    Returns
    -------
    weights : ndarray, shape (M, L, N)
        Array of weights for each model and sigma combination.
    """
    weights = np.zeros((len(sigmas_q), len(sigmas_i)) + distances.shape) * np.nan
    for idx_q, sigma_q in enumerate(sigmas_q):
        for idx_i, sigma_i in enumerate(sigmas_i):
            for idx_d, dd in enumerate(distances):
                # dd is the distance of each model to the idx_d-th model (='Truth')
                nu, de = calculate_weights(dd, distances, sigma_q, sigma_i)
                ww = nu/de
                assert np.isnan(ww[idx_d]), 'weight for model dd should be nan'
                ww[idx_d] = 0.  # set weight=0 to exclude the 'True' model
                assert ww.sum() != 0, 'weights = 0! sigma_q too small?'
                ww /= ww.sum()  # normalize weights
                # ww[ww < 1.e-10] = 0.  # set small weights to zero  # NOTE!!
                weights[idx_q, idx_i, idx_d] = ww
    return weights


def distance_matrix(data):
    """Takes a dataset with ensemble member/lon/lat. Flattens lon/lat into a single dimension.
    Calculates the distance between every ensemble member.

    Returns 2D NxN array, where N == number of ensemble members.
    """
    n_members = data.shape[0]
    
    data = data.reshape(n_members, -1)
    
    # pdist does not work with NaN
    idx = np.where(np.all(np.isfinite(data), axis=0))[0]
    data = data[:, idx]

    d_matrix = squareform(pdist(data, metric='euclidean'))
    np.fill_diagonal(d_matrix, np.nan)

    return d_matrix


def calculate_independence(data_array: 'xarray.DataArray') -> 'xarray.DataArray':
    """calculate_independence."""
    
    diff = xr.apply_ufunc(
        distance_matrix, data_array,
        input_core_dims=[['model_ensemble', 'lat', 'lon']],
        output_core_dims=[['perfect_model_ensemble', 'model_ensemble']]
    )

    diff.name = 'data'
    # diff = diff.expand_dims({'diagnostic': [idx]})

    return diff

def visualize_independence():
    """visualize_independence."""
    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


def calculate_performance():
    """calculate_performance."""
    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


def visualize_performance():
    """visualize_performance."""
    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


def calculate_weights():
    """calculate_weights."""
    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


def visualize_weighting():
    """visualize_weighting."""
    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


def read_metadata(cfg, projects: list) -> dict:
    """Read the metadata from the config file
    Project is the list of project identifiers to select."""
    d = defaultdict(list)

    for project in projects:
        metadata = select_metadata(cfg['input_data'].values(), project=project)

        for item in metadata:
            variable = item['short_name']

            d[variable].append(item)

    return d


def read_diagnostic_data(datasets: list):
    """Read the diagnostic data from the list of given datasets (metadata) into an Xarray."""
    to_concat = []

    for dataset in datasets:
        fn = dataset['filename']
        xarr = xr.open_dataset(fn)
        to_concat.append(xarr)

    diagnostic = xr.concat(to_concat, dim='model_ensemble')

    return diagnostic


def main(cfg):
    """Perform climwip weighting method."""

    print("\nBEGIN DIAGNOSTIC\n")

    observations = read_metadata(cfg, projects=['native6'])
    model_data = read_metadata(cfg, projects=['CMIP5'])

    variables = model_data.keys()

    diagnostics = []  # TODO: is this necessary to keep around?
    independences = []

    for variable, datasets in model_data.items():
        print(variable)
        diagnostic = read_diagnostic_data(datasets)
    
        independence = calculate_independence(diagnostic[variable])

        print(independence)
        print()

        diagnostics.append(diagnostic)
        independences.append(independence)

    independences = xr.concat(independences, dim='diagnostic')

    from IPython import embed
    embed();exit()

    independence = calculate_independence(diagnostics)

    visualize_independence(independence)

    performance = calculate_performance(cmip5, obs)
    visualize_performance(performance)

    weights = calculate_weights(independence, performance)
    visualize_weighting(weights)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
