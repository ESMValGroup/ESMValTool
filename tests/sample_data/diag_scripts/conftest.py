"""Configure pytest for diagnostic tests."""


def pytest_addoption(parser):
    parser.addoption("--save_imagehashes")
