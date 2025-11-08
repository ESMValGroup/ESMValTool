"""Functions to fetch commit information from local git repositories."""

import subprocess
from pathlib import Path


def get_shas_from_git(repos):
    """
    Fetch commit information from local git repos.

    Parameters
    ----------
    repos : dict[str, str | None]
        A dictionary where keys are in the form ``<package>_<day>`` and values
        are the path to the repo, or None. E.g.
        ``{"core_today": "path/to/repo", "core_yesterday": None}``
    Returns
    -------
    dict
        A dictionary where keys are the package and the values are a dict of
        days and short SHAs. E.g.
        ``{"ESMValCore": {"today": "abcd123", "yesterday": "efgh456"}...}``.
    """
    all_shas = {"ESMValCore": {}, "ESMValTool": {}}
    for package, repo_path in repos.items():
        if repo_path is not None and is_git_repo(repo_path):
            package, day = package.split("_")
            all_shas[package][day] = fetch_sha_from_git_log(repo_path)
    return all_shas


def is_git_repo(path):
    """
    Check a passed value is a valid git directory.

    Parameters
    ----------
    path : str | None
        Path to a git repo, or None.

    Returns
    -------
    bool
        If the passed value is a valid git directory.
    """
    try:
        return Path(path).expanduser().resolve().joinpath(".git").is_dir()
    except (TypeError, OSError):
        return False


def fetch_sha_from_git_log(package_path):
    """
    Use ``git log`` to get the SHA of the latest commit in a local git repo.

    Parameters
    ----------
    package_path : str
        Path to a valid git repo.

    Returns
    -------
    sha: str
        The SHA of the latest commit in the repo.
    """
    command = ["git", "log", "-1", "--pretty=%H"]
    sha = subprocess.run(
        command, cwd=package_path, capture_output=True, check=True, text=True
    )
    return sha.stdout.strip()
