from esmvalcore.preprocessor import (DEFAULT_ORDER, MULTI_MODEL_FUNCTIONS,
                                     _get_itype)


def test_first_argument_name():
    """Check that the input type of all preprocessor functions is valid."""
    valid_itypes = ('file', 'files', 'cube', 'cubes', 'products')
    for step in DEFAULT_ORDER:
        itype = _get_itype(step)
        assert itype in valid_itypes, (
            "Invalid preprocessor function definition {}, first argument "
            "should be one of {} but is {}".format(step, valid_itypes, itype))


def test_multi_model_exist():
    assert MULTI_MODEL_FUNCTIONS.issubset(set(DEFAULT_ORDER))
