"""Configure pytest."""

import sys

import pytest


def pytest_addoption(parser):
    """Add option to save imagehashes (used for diagnostic tests)."""
    parser.addoption("--save_imagehashes")


def pytest_runtest_setup(item):
    """Only run diagnostic tests that produce images on Linux machines."""
    if item.get_closest_marker("diagnostic_image_output") is not None:
        platform = sys.platform
        if platform != "linux":
            pytest.skip(
                f"Diagnostic tests that produce images are not supported on "
                f"platform {platform}",
            )
