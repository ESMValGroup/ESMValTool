"""ESMValTool diagnostics package."""
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("ESMValTool")
except PackageNotFoundError as exc:
    raise PackageNotFoundError(
        "ESMValTool package not found, please run `pip install -e .` before "
        "importing the package.") from exc


class ESMValToolDeprecationWarning(UserWarning):
    """Custom deprecation warning."""
