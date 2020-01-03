import pytest


def pytest_addoption(parser):
    """Add a command line option to skip tests that require installation."""
    parser.addoption(
        "--installation",
        action="store_true",
        default=False,
        help="run tests that require installation")


def pytest_collection_modifyitems(config, items):
    """Select tests to run based on command line options."""
    if config.getoption("--installation"):
        return
    skip_install = pytest.mark.skip(reason="need --installation option to run")
    for item in items:
        if "install" in item.keywords:
            item.add_marker(skip_install)
