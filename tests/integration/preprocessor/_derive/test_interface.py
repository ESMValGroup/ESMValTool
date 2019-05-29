from iris.cube import Cube, CubeList

from esmvaltool.preprocessor import derive
from esmvaltool.preprocessor._derive import get_required


def test_get_required():

    variables = get_required('alb')

    reference = [
        {
            'short_name': 'rsds',
        },
        {
            'short_name': 'rsus',
        },
    ]

    assert variables == reference


def test_get_required_with_fx():

    variables = get_required('nbp_grid')

    reference = [{
        'short_name': 'nbp',
        'fx_files': ['sftlf'],
    }]

    assert variables == reference


def test_derive_nonstandard_nofx():

    short_name = 'alb'
    long_name = 'albedo at the surface'
    units = 1
    standard_name = ''

    rsds = Cube([2.])
    rsds.standard_name = 'surface_downwelling_shortwave_flux_in_air'

    rsus = Cube([1.])
    rsus.standard_name = 'surface_upwelling_shortwave_flux_in_air'

    cubes = CubeList([rsds, rsus])

    alb = derive(cubes, short_name, long_name, units, standard_name)

    print(alb)
    assert alb.var_name == short_name
    assert alb.long_name == long_name
    assert alb.units == units
    assert alb.data == [0.5]


def test_derive_noop():

    alb = Cube([1.])
    alb.var_name = 'alb'
    alb.long_name = 'albedo at the surface'
    alb.units = 1

    cube = derive([alb], alb.var_name, alb.long_name, alb.units)

    print(cube)
    assert cube is alb
