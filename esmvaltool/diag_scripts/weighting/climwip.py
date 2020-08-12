"""
Implementation of step iii and iv of the climwip weighting scheme

Lukas Brunner et al. section 2.4
https://iopscience.iop.org/article/10.1088/1748-9326/ab492f
"""
from scipy import stats
from esmvaltool.diag_scripts.shared import run_diagnostic, select_metadata

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



def calculate_independence():
    """calculate_independence."""
    # TODO: complete this function
    returned_stuff = None
    return returned_stuff


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


def main(cfg):
    """Perform climwip weighting method."""
    obs = select_metadata(cfg['input_data'].values(), project='native6')
    cmip5 = select_metadata(cfg['input_data'].values(), project='cmip5')

    independence = calculate_independence(cmip5)
    visualize_independence(independence)

    performance = calculate_performance(cmip5, obs)
    visualize_performance(performance)

    weights = calculate_weights(independence, performance)
    visualize_weighting(weights)


if __name__ == '__main__':
    with run_diagnostic() as config:
        main(config)
