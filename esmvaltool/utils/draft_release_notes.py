"""Draft release notes.

To use this tool, follow these steps:
1) `pip install pygithub`
2) Create an access token and store it in the file ~/.github_api_key, see:
https://help.github.com/en/github/authenticating-to-github/creating-a-personal-access-token-for-the-command-line
3) set PREVIOUS_RELEASE to the date/time of the previous release in code below
"""
import datetime
from pathlib import Path

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

from esmvalcore import __version__

VERSION = f"v{__version__}"
GITHUB_REPO = "ESMValGroup/ESMValCore"

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


def draft_notes_since(previous_release_date, labels):
    """Draft release notes containing the merged pull requests.

    Arguments
    ---------
    previous_release_date: datetime.datetime
        date of the previous release
    labels: list
        list of GitHub labels that deserve separate sections
    """
    session = Github(GITHUB_API_KEY)
    repo = session.get_repo(GITHUB_REPO)
    pulls = repo.get_pulls(
        state='closed',
        sort='updated',
        direction='desc',
    )

    lines = {}
    labelless_pulls = []
    for pull in pulls:
        print(pull.updated_at, pull.merged_at, pull.number, pull.title)
        if pull.updated_at < previous_release_date:
            break
        if pull.merged:
            if pull.merged_at < previous_release_date:
                continue
            pr_labels = {label.name for label in pull.labels}
            for label in labels:
                if label in pr_labels:
                    break
            else:
                labelless_pulls.append(pull)
                label = 'enhancement'

            user = pull.user
            username = user.login if user.name is None else user.name
            title = pull.title
            title = title[0].upper() + title[1:]
            line = (f"-  {title} (`#{pull.number} "
                    f"<{pull.html_url}>`__) "
                    f"`{username} <https://github.com/{user.login}>`__")
            if label not in lines:
                lines[label] = []
            lines[label].append((pull.closed_at, line))

    # Warn about label-less PR:

    if labelless_pulls:
        print('\nPlease add labels to the following PR:')
        for pull in labelless_pulls:
            print(pull.html_url)
        print('\n')

    # Create sections
    sections = [
        VERSION,
        '-' * len(VERSION),
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


if __name__ == '__main__':

    PREVIOUS_RELEASE = datetime.datetime(2020, 10, 13, 00)
    LABELS = (
        'bug',
        'deprecated feature',
        'documentation',
        'diagnostic',
        'fix for dataset',
        'cmor',
        'preprocessor',
        'observations',
        'api',
        'testing',
        'enhancement',
    )

    draft_notes_since(PREVIOUS_RELEASE, LABELS)
