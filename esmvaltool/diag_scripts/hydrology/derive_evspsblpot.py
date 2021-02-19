"""Derive Potential Evapotransporation (evspsblpot) using De Bruin (2016).

De Bruin, H. A. R., Trigo, I. F., Bosveld, F. C., Meirink, J. F.: A
Thermodynamically Based Model for Actual Evapotranspiration of an Extensive
Grass Field Close to FAO Reference, Suitable for Remote Sensing Application,
American Meteorological Society, 17, 1373-1382, DOI: 10.1175/JHM-D-15-0006.1,
2016.
"""
import numpy as np
import iris


def tetens_derivative(tas):
    """Compute the derivative of Teten's formula for saturated vapor pressure.

    Tetens formula (https://en.wikipedia.org/wiki/Tetens_equation) :=
    es(T) = e0 * exp(a * T / (T + b))

    Derivate (checked with Wolfram alpha)
    des / dT = a * b * e0 * exp(a * T / (b + T)) / (b + T)^2
    """
    # Ensure temperature is in degC
    tas.convert_units('degC')

    # Saturated vapour pressure at 273 Kelvin
    e0_const = iris.coords.AuxCoord(np.float32(6.112),
                                    long_name='Saturated vapour pressure',
                                    units='hPa')
    emp_a = np.float32(17.67)  # empirical constant a

    # Empirical constant b in Tetens formula
    emp_b = iris.coords.AuxCoord(np.float32(243.5),
                                 long_name='Empirical constant b',
                                 units='degC')
    exponent = iris.analysis.maths.exp(emp_a * tas / (tas + emp_b))
    return (exponent * e0_const / (tas + emp_b)**2) * (emp_a * emp_b)


def get_constants(psl):
    """Define constants to compute De Bruin (2016) reference evaporation.

    The Definition of rv and rd constants is provided in
    Wallace and Hobbs (2006), 2.6 equation 3.14.
    The Definition of lambda and cp is provided in Wallace and Hobbs 2006.
    The Definition of beta and cs is provided in De Bruin (2016), section 4a.
    """
    # Ensure psl is in hPa
    psl.convert_units('hPa')

    # Definition of constants
    # source='Wallace and Hobbs (2006), 2.6 equation 3.14',
    rv_const = iris.coords.AuxCoord(np.float32(461.51),
                                    long_name='Gas constant water vapour',
                                    units='J K-1 kg-1')
    # source='Wallace and Hobbs (2006), 2.6 equation 3.14',
    rd_const = iris.coords.AuxCoord(np.float32(287.0),
                                    long_name='Gas constant dry air',
                                    units='J K-1 kg-1')

    # Latent heat of vaporization in J kg-1 (or J m-2 day-1)
    # source='Wallace and Hobbs 2006'
    lambda_ = iris.coords.AuxCoord(np.float32(2.5e6),
                                   long_name='Latent heat of vaporization',
                                   units='J kg-1')

    # Specific heat of dry air constant pressure
    # source='Wallace and Hobbs 2006',
    cp_const = iris.coords.AuxCoord(np.float32(1004),
                                    long_name='Specific heat of dry air',
                                    units='J K-1 kg-1')

    # source='De Bruin (2016), section 4a',
    beta = iris.coords.AuxCoord(np.float32(20),
                                long_name='Correction Constant',
                                units='W m-2')

    # source = 'De Bruin (2016), section 4a',
    cs_const = iris.coords.AuxCoord(np.float32(110),
                                    long_name='Empirical constant',
                                    units='W m-2')

    # source = De Bruin (10.1175/JHM-D-15-0006.1), page 1376
    # gamma = (rv/rd) * (cp*msl/lambda_)
    # iris.exceptions.NotYetImplementedError: coord / coord
    rv_rd_const = rv_const.points[0] / rd_const.points[0]
    gamma = rv_rd_const * (psl * cp_const / lambda_)
    return gamma, cs_const, beta, lambda_


def debruin_pet(psl, rsds, rsdt, tas):
    """Compute De Bruin (2016) reference evaporation.

    Implement equation 6 from De Bruin (10.1175/JHM-D-15-0006.1)
    """
    # Variable derivation
    delta_svp = tetens_derivative(tas)
    gamma, cs_const, beta, lambda_ = get_constants(psl)

    # the definition of the radiation components according to the paper:
    kdown = rsds
    kdown_ext = rsdt
    # Equation 5
    rad_factor = np.float32(1 - 0.23)
    net_radiation = (rad_factor * kdown) - (kdown * cs_const / kdown_ext)
    # Equation 6
    # the unit is W m-2
    ref_evap = ((delta_svp / (delta_svp + gamma)) * net_radiation) + beta

    pet = ref_evap / lambda_
    pet.var_name = 'evspsblpot'
    pet.standard_name = 'water_potential_evaporation_flux'
    pet.long_name = 'Potential Evapotranspiration'
    return pet
