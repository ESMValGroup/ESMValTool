"Functions to fetch commit information from local git repositories."

import subprocess
from dataclasses import dataclass
from pathlib import Path


@dataclass
class CommitInfo:
    """
    A dataclass to hold commit information for ESMValCore and ESMValTool.

    Attributes
    ----------
    core : list[dict]
        A list of dictionaries containing commit information for ESMValCore.
    tool : list[dict]
        A list of dictionaries containing commit information for ESMValTool.
    """

    core: list[dict]
    tool: list[dict]


def get_commits_from_git(repos):
    """
    Fetch commit information from local git repos.

    Parameters
    ----------
    repos : dict[str, str | None]
        A dictionary where keys are in the form ``<package>_<day>`` and values
        are the path to the repo, or None. E.g.
        ``{"core_today": "path/to/repo", "core_yesterday": None}``

    Raises
    ------
    ValueError
        If repos does not contain any valid git directories, or does not
        contain valid git directories for today.

    Returns
    -------
    CommitInfo
        A CommitInfo dataclass containing two lists of dictionaries, one for
        each of ESMValCore and ESMValTool. Each commit dict has the following
        fields ``date, sha, author, message``.
    """
    repo_validity = {
        pkg_day: is_git_repo(path) for pkg_day, path in repos.items()
    }

    if not any(repo_validity.values()):
        raise ValueError("No valid git repos found.")

    if not (repo_validity["core_today"] or repo_validity["tool_today"]):
        raise ValueError("Today's commit info is unavailable.")

    if not (
        repo_validity.get("core_yesterday")
        or repo_validity.get("tool_yesterday")
    ):
        print("Only today's commit info is available.")
        commit_info = CommitInfo(
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


def get_all_commits_for_today_and_yesterday(valid_repos):
    """
    Fetch ``git log`` information for a range of commits.

    Expects a valid git repo for both ESMValCore and ESMValTool for today
    and yesterday.

    Parameters
    ----------
    valid_repos : dict[str, str]
        A dict of valid git repos.

    Returns
    -------
    CommitInfo
        A CommitInfo dataclass containing two lists of dictionaries, one for
        each of ESMValCore and ESMValTool.
    """
    commit_info = CommitInfo([], [])
    for package in ["core", "tool"]:
        yesterdays_sha = query_git_log(valid_repos[f"{package}_yesterday"])[0][
            "sha"
        ]
        commit_range = query_git_log(
            valid_repos[f"{package}_today"], yesterdays_sha
        )
        setattr(commit_info, package, commit_range)
    return commit_info


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
        not passed, one commit is returned. If ``sha`` is passed, multiple
        commits may be returned.
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
    Add report messages to a CommitInfo dataclass.

    Parameters
    ----------
    commit_info : CommitInfo
        A CommitInfo dataclass containing two lists of package commits.
    """
    for package_commits in [commit_info.core, commit_info.tool]:
        package_commits[0]["report_flag"] = "Version tested this cycle >>>"
        if len(package_commits) > 1:
            package_commits[-1]["report_flag"] = (
                "Version tested last cycle >>>"
            )
