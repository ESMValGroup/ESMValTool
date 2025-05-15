"""Configure pytest."""


def pytest_addoption(parser):
    """Add option to save imagehashes (used for diagnostic tests)."""
    parser.addoption("--save_imagehashes")
