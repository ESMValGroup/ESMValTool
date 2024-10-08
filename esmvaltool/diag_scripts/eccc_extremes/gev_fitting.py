import numpy as np
from scipy.stats import genextreme as gev
from scipy.stats import kstest, cramervonmises
from scipy.optimize import minimize


def _loc_function(x: np.ndarray, non_stat_params: list, 
                                                dependency: str = 'linear'):
    '''
    This is the function for the location dependency

    Parameters
    ----------
    x : 
        array with the covariate
    non_stat_params : 
        non-stationary parameters in format [loc_1, loc_2, shape, scale]
    dependency : 
        the type of dependency of the location parameter on the
        covariate. Valid entries: 'linear'
    
    Returns
    -------
    loc : np.ndarray
        array with location parameters values
    
    Raises
    ------
    NotImplementedError
        if the type of dependency is not implemented
    '''
    
    if dependency == 'linear':
        loc = non_stat_params[0]*x + non_stat_params[1]
    else: 
        raise NotImplementedError(f"Dependency method '{dependency}' is not "
                                  "implemented yet. The available dependecies "
                                  "are 'linear'.")

    return loc


def _neg_log_likelihood(non_stat_params: list, data: np.ndarray,
                                               covariate: np.ndarray):
    '''
    Function to minimize for the non_stationary fit

    Parameters
    ----------
    non_stat_params : 
        non-stationary GEV parameters in format [loc_1, loc_2, shape, scale]
    data: 
        extremes data for the fit (e.g., TXx, RX1day)
    covariate :
        data with the covariate for the fit (e.g., GSAT, CO2 etc)
    '''
    loc = _loc_function(covariate, non_stat_params)
    shape, scale = non_stat_params[2:]
    return -gev.logpdf(data, shape, loc=loc, scale=scale).sum()


class StationaryGEV:
    '''
    Class with statistics and return period from the stationary GEV

    Attributes
    ----------
    rp : float
        return period in years
    shape : float | nan
        shape parameter of the fit
    loc : float | nan
        location parameter of the fit
    scale : float| nan
        scale parameter of the fit
    notion : str
        which GEV shape parameter notion is used. Possible options: 
        scipy or R. This attribute is used simply for checking.
    '''
    def __init__(self, data: np.ndarray, event: float, 
                 weights: np.ndarray | None = None, 
                 initial: dict | None = None):
        '''
        Parameters
        ----------
            data : 
                an array with the timeseries that is used for fitting
            event : 
                the strength of the event for return period calculation
            weights : 
                an array with the weights, the same shape as data
            initial : 
                an dictionary with the initial conditions for the fit
                keys are {'location', 'scale', 'shape'} 
        
        Raises
        ------
        ValueError 
            if the weights and data are do not have matching shape or
            the keywords in the initial condition dictionary are not 
            correct or initial conditions are not float type
        '''

        self._check_weights(data, weights)
        initial_guess = self._get_initial_guess(data, initial)
                    
        stat_gev = cex.fit_gev(data, returnValue=event, getParams=True,
                               weights=weights, initial=initial)
        
        self.notation= 'scipy'

    def _check_weights(data, weights):
        if not(weights is None):
            if data.shape != weights.shape:
                raise ValueError("The shapes of the data and weights for "
                                 "stationary GEV fit have unmatching shapes")        

    def _check_initial(initial): 

        if not(initial is None):
            if list(initial.keys()) != ['location', 'scale', 'shape']:
                raise ValueError("The initial conditions supposed to be "
                                 "['location', 'scale', 'shape'], currently "
                                 f"it is {list(initial.keys())}")
            if not(isinstance(initial['shape'], float)):
                raise ValueError("The type of shape parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['shape'])}")
            if not(isinstance(initial['scale'], float)):
                raise ValueError("The type of scale parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['scale'])}")
            if not(isinstance(initial['location'], float)):
                raise ValueError("The type of location parameter for initial "
                                 "conditions supposed to be float, currently "
                                 f"it is {type(initial['location'])}")        


    def _get_initial_guess(self, data, initial):

        self._check_initial(initial)

        if initial:    
            loc = initial['location']
            shape = initial['shape']
            scale = initial['scale']
        else:
            shape, loc, scale = gev.fit(data)

        return [loc, shape, scale]


    def set_x_gev(self, metric: str, x_gev : np.ndarray | None = None):
        '''
        This function checks is class has x_gev attribute and sets it

        Parameters
        ----------
        metric:
            metric for which the x_gev is calculated
        x_gev: 
            optional array of x_gevs

        Raises
        ------
        ValueError
            if the x_gev was not provided and the class doesn't have 
            existing x_gev attribute
        '''
        if x_gev is None: 
            try:
                # checking if x_gev exists
                x_gev = self.x_gev 
            except: 
                raise ValueError("There is no existing attribute x_gev, x "
                                 "values for PDF calculations should be "
                                 "provided")
        else: 
            try: 
                x_gev = self.x_gev
                raise UserWarning("The class already contains x_gev parameter, "
                               "will use existing class attribute as x_gev for "
                               f"{metric} calculations")
            except:
                self.x_gev = x_gev

        return

    def obtain_pdf(self, x_gev : np.ndarray | None = None):
        '''
        This function calculates pdf for x_gev

        Parameters
        ----------
        x_gev : 
            optional array with x-es for which PDFs are calculated
        '''

        self.set_x_gev(metric='probability density function', x_gev=x_gev)
        self.pdf = gev.pdf(self.x_gev, self.shape, loc=self.loc, 
                                                   scale=self.scale)

        return

    def obtain_rp_curve(self, x_gev : np.ndarray | None = None ):
        '''
        This function calculates return periods (RPs) curve for x_gev

        Parameters
        ----------
        x_gev : 
            optional array with x-es for which RPs are calculated
        '''
        self.set_x_gev(metric='return period curve', x_gev=x_gev)
        self.rp_curve = 1/gev.sf(self.x_gev, self.shape, loc=self.loc, scale=self.scale)

        return


class NonStationaryGEV:
    '''
    Class with statistics and return period from the non-stationary GEV

    The non-stationary GEV fit assumes that one or multiple GEV params
    are dependent on a covariate (eg., GSAT, CO2 etc). While the
    dependency can be non-linear, in climextremes only linear is 
    implemented and only for location parameter. The formula:
    loc = loc_a*covarite + loc_b

    Attributes
    ----------
    rp : float
        return period in years for the idx (e.g., year of interest)
    shape : float | nan
        shape parameter of the fit
    loc_a : float | nan
        slope parameter of the fitted location parameter 
        (loc = loc_a*covariate + loc_b)
    loc_b : float | nan
        intersect parameter of the fitted location parameter 
        (loc = loc_a*covariate + loc_b)
    loc_idx : float | nan
        the location parameter for the idx (e.g., year of interest)
    scale : float | nan
        scale parameter of the fit
    cov_value: float
        value for the covariate for the idx (e.g., year of interest)
    notion : str
        which GEV shape parameter notion is used. Possible options: 
        scipy or R. This attribute is used simply for checking.
    '''
    def __init__(self, data: np.ndarray, covariate: np.ndarray, 
                     event: float, idx: int = -1, initial: dict | None = None):
        '''
        Parameters
        ----------
            data : 
                an array with the timeseries that is used for fitting
            covariate: 
                an array with a covariate that is used for fitting
                (e.g., GSAT, CO2 etc.)
            event : 
                the strength of the event for return period calculation
            idx : 
                index for which the retrun period has to be calculated
                (e.g., year of interest), default is -1 (the last)
            initial : 
                an dictionary with the initial conditions for the fit
                keys are {'location', 'scale', 'shape'}
            
        
        Raises
        ------
        ValueError
            If the length of the data and covariate is different or
            the keywords in the initial condition dictionary are not
            correct or initial conditions are not float type
        '''

        self._check_data_covariate(data, covariate)
        self._check_initials(initial)

        initial_guess = self._get_initial_guess(data, initial)
        
        self._fit(data, covariate, initial_guess)

        self.loc_idx = float(self.loc_a * covariate[idx] + self.loc_b)
        self.cov_value = float(covariate[idx])

        self.rp = float(1/gev.sf(event, self.shape, loc=self.loc_idx,
                                                    scale=self.scale))

        self.notation= 'scipy'


    def _check_data_covariate(data: np.ndarray, covariate: np.ndarray):
        if len(data) != len(covariate):
            raise ValueError("The lengths of the data and covariate"
                                                        "don't match.")
        
    def _check_initials(initial: dict | None):

        if initial:
            if not(isinstance(initial['shape'], float)):
                raise ValueError("The type of shape parameter for initial "
                                "conditions supposed to be float, currently "
                                f"it is {type(initial['shape'])}")
            if not(isinstance(initial['scale'], float)):
                raise ValueError("The type of scale parameter for initial "
                                "conditions supposed to be float, currently "
                                f"it is {type(initial['scale'])}")
            if not(isinstance(initial['location'], float)):
                raise ValueError("The type of location parameter for initial "
                                "conditions supposed to be float, currently "
                                f"it is {type(initial['location'])}")   

        
    def _get_initial_guess(data: np.ndarray, initial: dict | None):
        '''
        Obtain initial guess for the fit

        Parameters
        ----------
        data: 
            extremes data for the fit (e.g., TXx, RX1day)
        initial : 
            initial guess parameters

        Returns
        -------
        initial_guess :  list
            parameters in format [loc_1, loc_2, shape, scale]
        '''

        if initial:    
            loc = initial['location']
            shape = initial['shape']
            scale = initial['scale']
        else:
            shape, loc, scale = gev.fit(data)

        initial_guess = [0, loc, shape, scale]

        return initial_guess
    
    def _fit(self, data: np.ndarray, covariate: np.ndarray, initial_guess: list):
        '''
        Fit non-stationary GEV to the data

        Parameters
        ----------
        data: 
            extremes data for the fit (e.g., TXx, RX1day)
        covariate :
            data with the covariate for the fit (e.g., GSAT, CO2 etc)
        initial_guess : 
            initial parameters for the fit [loc_1, loc_2, shape, scale]
        '''

        optim_result = minimize(_neg_log_likelihood, initial_guess, 
                                args=(data, covariate), method='SLSQP')
        if optim_result.status:
            self.loc_a = float(optim_result.x[0])
            self.loc_b = float(optim_result.x[1])
            self.shape = float(optim_result.x[2])
            self.scale = float(optim_result.x[3])
        else:
            self.loc_a = float(np.nan)
            self.loc_b = float(np.nan)
            self.shape = float(np.nan)
            self.scale = float(np.nan)                

        return  