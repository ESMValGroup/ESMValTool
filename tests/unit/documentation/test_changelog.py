"""Tests for the changelog."""
import collections
import os
import re
from typing import List

# For entries in the following list, two appearances in the changelog (but not
# more) are allowed. This can be useful if a pull request needs to appear in
# multiple sections.
ALLOWED_DUPLICATES: List[str] = [
    '<https://github.com/ESMValGroup/ESMValTool/pull/1897>',
]


def test_duplications_in_changelog():
    changelog_path = os.path.join(os.path.dirname(__file__), '../../..',
                                  'doc/sphinx/source/changelog.rst')
    with open(changelog_path) as changelog:
        changelog = changelog.read()

    # Find all pull requests
    pr_links = re.compile(
        "<https://github.com/ESMValGroup/ESMValTool/pull/[0-9]+>")
    links = pr_links.findall(changelog)

    # Remove one entry for the allowed exceptions
    for pr_link in ALLOWED_DUPLICATES:
        if pr_link in links:
            links.remove(pr_link)

    # Check for duplicates
    if len(links) != len(set(links)):
        print('The following PR are duplicated in the changelog:')
        print('\n'.join((link
                         for link, count in collections.Counter(links).items()
                         if count > 1)))
        assert False
