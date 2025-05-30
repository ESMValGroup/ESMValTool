"""
Functions to fetch commit information from local git repositories.
"""

import subprocess
from pathlib import Path


def get_commits_from_git(repos):
    """
    Fetch commit information from local git repos.

    Parameters
    ----------
    repos : dict[str, str | None]
        A dictionary where keys are in the form ``repo_day`` and values are
        the path to the repo, or None.  E.g. ``{"core_today": "path/to/repo",
        "core_yesterday": None}``

    Raises
    ------
    ValueError
        If repos does not contain enough valid git directories to fetch
        useable commit data.

    Returns
    -------
    tuple[list[dict], list[dict]]
        A tuple of two lists of dictionaries. Each dictionary include
        information on a commit or a range of commits for
        ``(ESMValCore, ESMValTool)``. Each commit dict has the following fields
        ``date, sha, author, message``.
    """
    repo_validity = {key: is_git_repo(path) for key, path in repos.items()}

    if not any(repo_validity.values()):
        raise ValueError("No valid git repos passed")

    elif not (repo_validity["core_today"] or repo_validity["tool_today"]):
        raise ValueError("Today's git commits unavailable.")

    elif not (
        repo_validity["core_yesterday"] and repo_validity["tool_yesterday"]
    ):
        print("Only today's git info is available.")
        commit_info = (
            query_git_log(repos["core_today"]),
            query_git_log(repos["tool_today"]),
        )

    else:
        commit_info = get_all_commits_for_today_and_yesterday(repos)

    add_report_messages_to_commits(commit_info)

    return commit_info


def is_git_repo(path):
    """
    Check a passed value is a valid git directory.

    Parameters
    ----------
    path : str|None
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


def get_all_commits_for_today_and_yesterday(valid_repos):
    """
    Fetch information on a range of commits.

    Fetches a range of commits

    Parameters
    ----------
    valid_repos : dict[str, str]
        A dict of valid git repos.

    Returns
    -------
    tuple[list[dict], list[dict]]
        A tuple with two lists of dictionaries that include a range of commit
        information for ``(ESMValCore, ESMValTool)``.

    """
    core_yesterday_sha = query_git_log(valid_repos["core_yesterday"][0]["sha"])
    core_commits = query_git_log(valid_repos["core_today"], core_yesterday_sha)

    tool_yesterday_sha = query_git_log(valid_repos["tool_yesterday"][0]["sha"])
    tool_commits = query_git_log(valid_repos["tool_today"], tool_yesterday_sha)

    return (core_commits, tool_commits)


def query_git_log(package_path, sha=None):
    """
    Use ``git log`` to fetch commit information from a local git repo.

    Parameters
    ----------
    package_path : str
        Path to a valid git repo.
    sha: str | None
        Optional. The sha of a previously tested commit. If provided, commits
        from HEAD back to the passed sha (inclusive) will be retrieved.

    Returns
    -------
    list[dict]
        A list of dicts where each dict represents one commit. Each commit has
        the following fields ``date, sha, author, message``. If ``sha`` is
        not passed, one commit is returned. If ``sha`` is passed, a minimum of
        two commits is returned.
    """
    command = [
        "git",
        "log",
        "-1",
        "--date=iso-strict",
        "--pretty=%cd^_^%h^_^%an^_^%s",
    ]

    if sha:
        command[2] = f"{sha}^..HEAD"

    raw_commit_info = subprocess.run(
        command, cwd=package_path, capture_output=True, check=True, text=True
    )
    processed_commit_info = []
    raw_commits = raw_commit_info.stdout.splitlines()
    for commit in raw_commits:
        split_fields = commit.split("^_^")
        processed_commit_info.append(
            {
                "report_flag": "",  # Needed for jinja2 template.
                "date": split_fields[0],
                "sha": split_fields[1],
                "author": split_fields[2],
                "message": split_fields[3],
            }
        )
    return processed_commit_info


def add_report_messages_to_commits(commit_info):
    """
    Add report messages to a git commit information dictionary.

    Add a report flag

    Parameters
    ----------
    commit_info : tuple[list[dict], list[dict]
        A tuple containg two list of package commits.
    """
    for package_commits in commit_info:
        package_commits[0]["report_flag"] = "Version tested this cycle >>>"
        if len(package_commits) > 1:
            package_commits[-1]["report_flag"] = (
                "Version tested last cycle >>>"
            )
