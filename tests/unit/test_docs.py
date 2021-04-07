"""Tests for the documentation."""
import collections
import os
import re


def test_duplications_in_changelog():
    changelog_path = os.path.join(os.path.dirname(__file__), '../..',
                                  'doc/sphinx/source/changelog.rst')
    with open(changelog_path) as changelog:
        changelog = changelog.read()

    pr_links = re.compile(
        "<https://github.com/ESMValGroup/ESMValTool/pull/[0-9]+>")
    links = pr_links.findall(changelog)
    if len(links) != len(set(links)):
        print('The following PR are duplicated in the changelog:')
        print('\n'.join((link
                         for link, count in collections.Counter(links).items()
                         if count > 1)))
        assert False
