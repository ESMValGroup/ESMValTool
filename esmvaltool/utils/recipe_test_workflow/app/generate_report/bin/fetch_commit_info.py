"""Fetch commit details from the GitHub API."""

import os

import requests

GITHUB_API_URL = "https://api.github.com"
GITHUB_API_PERSONAL_ACCESS_TOKEN = os.environ.get(
    "GITHUB_API_PERSONAL_ACCESS_TOKEN"
)
HEADERS = {
    "authorization": f"token {GITHUB_API_PERSONAL_ACCESS_TOKEN}",
    # Suggested here:
    # https://docs.github.com/en/rest/commits/commits?apiVersion=2022-11-28#get-a-commit--parameters
    # Explanation here:
    # https://docs.github.com/en/rest/overview/resources-in-the-rest-api#http-HEADERS
    "accept": "application/vnd.github+json",
}


def fetch_commit_details_from_github_api(
    shas_by_package_and_day, headers=HEADERS
):
    """
    Fetch commit details from the GitHub API for the given SHAs.

    Parameters
    ----------
    shas_by_package_and_day : dict[str, dict[str, str]]
        A dictionary where keys are the package names and values are dictionaries
        with days as keys and SHAs as values. E.g.
        {"ESMValCore": {"today": "abcd123", "yesterday": "efgh456"}...}.

    Returns
    -------
    dict[str, list[dict]]
        A dictionary where keys are the package names and values are lists of
        commit details for each day. E.g.
        {"ESMValCore": [{"sha": "abcd123", ...}, ...], "ESMValTool": [...]}
    """
    commit_details_by_package = {}
    for package, shas_by_day in shas_by_package_and_day.items():
        if shas_by_day.get("yesterday") is None or shas_by_day.get(
            "today"
        ) == shas_by_day.get("yesterday"):
            raw_commits = fetch_single_commit(
                package, "ESMValGroup", headers, shas_by_day["today"]
            )
        else:
            raw_commits = fetch_range_of_commits(
                package,
                "ESMValGroup",
                headers,
                newer_sha=shas_by_day["today"],
                older_sha=shas_by_day["yesterday"],
            )
        commit_info = process_commit_info(raw_commits)
        commit_details_by_package[package] = commit_info
    return commit_details_by_package


def make_api_call(url, headers=None, params=None):
    """
    Make a GET request to a given API url.

    Parameters
    ----------
    url : str
        The URL to make the request to.
    headers : dict, optional
        Headers to include in the request.
    params : dict, optional
        Query parameters to include in the request.

    Raises
    ------
    HTTPError
        If the request fails or returns with a status code other than 200.
    TimeOutError
        If the request times out.
    ConnectionError
        If there is a connection error.
    TooManyRequestsError
        If the API rate limit is exceeded.

    Returns
    -------
    Response
        The raw response from the API call.
    """
    try:
        response = requests.get(
            url, headers=headers, params=params, timeout=10
        )
        if response.status_code != 200:
            raise requests.exceptions.HTTPError(
                f"Unexpected status code for url={url} headers={headers} "
                f"params={params} - {response.status_code}: {response.text}"
            )
    except requests.exceptions.HTTPError as http_err:
        raise requests.exceptions.HTTPError(
            f"HTTP error occurred: {http_err}"
        ) from http_err
    except requests.exceptions.Timeout as timeout_err:
        raise requests.exceptions.Timeout(
            f"Request timed out: {timeout_err}"
        ) from timeout_err
    except requests.exceptions.ConnectionError as conn_err:
        raise requests.exceptions.ConnectionError(
            f"Connection error occurred: {conn_err}"
        ) from conn_err
    return response


def fetch_single_commit(repo, owner, headers, sha):
    """
    Fetch details of a single commit from the GitHub API.

    Parameters:
    ----------
    repo: str
        The name of the repository. E.g. "ESMValTool"
    owner: str
        The owner of the repository. E.g. "ESMValGroup"
    headers: dict
        Headers to include in the request.
    sha: str
        The SHA of the commit to fetch details for.

    Raises
    ------
    HTTPError
        If the commit is not found or if the request fails etc.

    Returns
    -------
    dict
        The raw commit data if found.
    """
    url = f"{GITHUB_API_URL}/repos/{owner}/{repo}/commits/{sha}"
    response = make_api_call(url, headers=headers)
    raw_commit = response.json()
    return raw_commit


def fetch_range_of_commits(repo, owner, headers, newer_sha, older_sha):
    """
    Fetch details for a range of commits from the GitHub API.

    The endpoint will return a range of commits in chronlogical order, from
    the newer SHA to the older SHA. The function fetches batches of 10 commits
    to avoid hitting the API rate limits. NOTE: The GitHub API will raise a
    HTTPError if the newer SHA is not found.

    Parameters:
    ----------
    repo : str
        The name of the repository. E.g. "ESMValTool"
    owner : str
        The owner of the repository. E.g. "ESMValGroup"
    headers : dict
        Headers to include in the request.
    newer_sha : str
        The SHA of the first commit to start fetching details for.
    older_sha : str
        The SHA of the commit to stop fetching at.

    Raises
    ------
    HTTPError
        If the newer SHA is not found (or if the request fails etc.)
    ValueError
        If too many pages are fetched, indicating a potential infinite loop.

    Returns
    -------
    list[dict]
        A list of raw commit data for the range of commits from newer_sha to
        older_sha, in chronological order.
    """
    url = f"{GITHUB_API_URL}/repos/{owner}/{repo}/commits"
    params = {
        "per_page": 10,
        "sha": newer_sha,
    }
    page = 1
    range_raw_commits = []

    fetched_end_sha = False

    while not fetched_end_sha:
        params["page"] = page
        response = make_api_call(url, headers=headers, params=params)

        page_raw_commits = response.json()
        for raw_commit in page_raw_commits:
            range_raw_commits.append(raw_commit)
            if raw_commit["sha"].startswith(older_sha):
                fetched_end_sha = True
                break

        page += 1
        if page > 5:
            raise ValueError(
                "Too many pages fetched, likely an infinite loop. Check the "
                "newer and older SHAs."
            )

    return range_raw_commits


def process_commit_info(raw_commit_info):
    """
    Extract required commit details.

    Parameters:
    -----------
    raw_commit_info : dict | list[dict]
        Raw commit information from the GitHub API. Either a single commit or a
        list of commits.

    Returns
    -------
    list[dict]
        A list of dictionaries containing processed commit information.
    """

    if not isinstance(raw_commit_info, list):
        raw_commit_info = [raw_commit_info]

    processed_commit_info = []
    for raw_commit in raw_commit_info:
        processed_comit = {
            "sha": raw_commit["sha"][:7],
            "author": raw_commit["commit"]["author"]["name"],
            "message": raw_commit["commit"]["message"],
            "date": raw_commit["commit"]["author"]["date"],
            "url": raw_commit["html_url"],
            "author_avatar": raw_commit["author"]["avatar_url"],
        }
        processed_commit_info.append(processed_comit)
    return processed_commit_info
