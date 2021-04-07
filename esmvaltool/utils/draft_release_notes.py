"""Draft release notes.

To use this tool, follow these steps:
1) `pip install pygithub`
2) Create an access token and store it in the file ~/.github_api_key, see:
https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line
3) set PREVIOUS_RELEASE to the date/time of the previous release in code below
4) Call the script passing the project to create release notes on: esmvalcore
or esmvaltool
"""
import datetime
from pathlib import Path
from typing import DefaultDict

import dateutil
import esmvalcore
import fire

import esmvaltool

try:
    from github import Github
except ImportError:
    print("Please `pip install pygithub`")

try:
    GITHUB_API_KEY = Path("~/.github_api_key").expanduser().read_text().strip()
except FileNotFoundError:
    print("Please create an access token and store it in the file "
          "~/.github_api_key, see:\nhttps://help.github.com/en/github/"
          "authenticating-to-github/creating-a-personal-access-token-"
          "for-the-command-line")

VERSION = {
    'esmvalcore': f"v{esmvalcore.__version__}",
    'esmvaltool': f"v{esmvaltool.__version__}"
}
GITHUB_REPO = {
    'esmvalcore': "ESMValGroup/ESMValCore",
    'esmvaltool': "ESMValGroup/ESMValTool",
}

PREVIOUS_RELEASE = {
    'esmvalcore': datetime.datetime(2020, 10, 13, 00),
    'esmvaltool': datetime.datetime(2020, 10, 26, 00),
}
LABELS = {
    'esmvalcore': (
        'bug',
        'deprecated feature',
        'documentation',
        'fix for dataset',
        'cmor',
        'preprocessor',
        'api',
        'testing',
        'installation',
        'enhancement',
    ),
    'esmvaltool': (
        'bug',
        'deprecated feature',
        'documentation',
        'diagnostic',
        'preprocessor',
        'observations',
        'testing',
        'installation',
        'enhancement',
    )
}

TITLES = {
    'bug': 'Bug fixes',
    'deprecated feature': 'Deprecations',
    'cmor': 'CMOR standard',
    'diagnostic': 'Diagnostics',
    'fix for dataset': 'Fixes for datasets',
    'observations': 'Observational and re-analysis dataset support',
    'testing': 'Automatic testing',
    'api': 'Notebook API (experimental)',
    'enhancement': 'Improvements',
}


def draft_notes_since(project, previous_release_date=None, labels=None):
    """Draft release notes containing the merged pull requests.

    Arguments
    ---------
    project: str
        Project to draft release notes from. Valid options are esmvaltool and
        esmvalcore
    previous_release_date: datetime.datetime
        date of the previous release
    labels: list
        list of GitHub labels that deserve separate sections
    """
    project = project.lower()
    if previous_release_date is None:
        previous_release_date = PREVIOUS_RELEASE[project]
    else:
        previous_release_date = dateutil.parse(previous_release_date)
    if labels is None:
        labels = LABELS[project]

    pulls = _get_pull_requests(project)

    lines = DefaultDict(list)
    labelless_pulls = []
    for pull in pulls:
        print(pull.updated_at, pull.merged_at, pull.number, pull.title)
        if pull.updated_at < previous_release_date:
            break
        if not pull.merged or pull.merged_at < previous_release_date:
            continue
        pr_labels = {label.name for label in pull.labels}
        for label in labels:
            if label in pr_labels:
                break
        else:
            labelless_pulls.append(pull)
            label = 'enhancement'
        lines[label].append((pull.closed_at, _compose_note(pull)))

    # Warn about label-less PR:
    _list_labelless_pulls(labelless_pulls)

    # Create sections
    sections = [
        VERSION[project],
        '-' * len(VERSION[project]),
        '',
        "This release includes",
    ]
    for label in labels:
        try:
            entries = sorted(lines[label])  # sort by merge time
        except KeyError:
            continue
        title = TITLES.get(label, label.title())
        sections.append('\n'.join(['', title, '~' * len(title), '']))
        sections.append('\n'.join(entry for _, entry in entries))
    notes = '\n'.join(sections)
    print(notes)


def _get_pull_requests(project):
    session = Github(GITHUB_API_KEY)
    repo = session.get_repo(GITHUB_REPO[project])
    pulls = repo.get_pulls(
        state='closed',
        sort='updated',
        direction='desc',
    )
    return pulls


def _list_labelless_pulls(labelless_pulls):
    if labelless_pulls:
        print('\nPlease add labels to the following PR:')
        for pull in labelless_pulls:
            print(pull.html_url)
        print('\n')


def _compose_note(pull):
    user = pull.user
    username = user.login if user.name is None else user.name
    title = pull.title
    title = title[0].upper() + title[1:]
    return (f"-  {title} (`#{pull.number} "
            f"<{pull.html_url}>`__) "
            f"`{username} <https://github.com/{user.login}>`__")


def main():
    """Entry point for the scrip."""
    def display(lines, out):
        text = "\n".join(lines) + "\n"
        out.write(text)

    fire.core.Display = display
    fire.Fire(draft_notes_since)


if __name__ == '__main__':
    main()
