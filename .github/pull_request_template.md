<!--
    Thank you for contributing to our project!

    Please do not delete this text completely, but read the text below and keep
    items that seem relevant. If in doubt, just keep everything and add your
    own text at the top, a reviewer will update the checklist for you.

    While the checklist is intended to be filled in by the technical and scientific
    reviewers, it is the responsibility of the author of the pull request to make
    sure all items on it are properly implemented.
-->

## Description

<!--
    Please describe your changes here, especially why this pull request makes
    ESMValTool better and what problem it solves.

    Before you start, please read our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/introduction.html).

    Please fill in the GitHub issue that is closed by this pull request, e.g. Closes #1903
-->
- Closes/Fixes: #issue_number
- Link to documentation:

* * *

## Before you get started

<!--
    Please discuss your idea with the development team before getting started,
    to avoid disappointment or unnecessary work later. The way to do this is
    to open a new issue on GitHub.
-->

- [ ] [Create an issue](https://github.com/ESMValGroup/ESMValTool/issues) to discuss what you are going to do

## Checklist

- [ ] PR has a descriptive title for the [changelog](https://docs.esmvaltool.org/en/latest/changelog.html)
- [ ] Code follows the [style guide](https://docs.esmvaltool.org/en/latest/community/introduction.html#code-style)
- [ ] [Documentation](https://docs.esmvaltool.org/en/latest/community/introduction.html#documentation) is available
- [ ] YAML files pass [`pre-commit`](https://esmvaltool--1924.org.readthedocs.build/en/1924/community/introduction.html#pre-commit) or [`yamllint`](https://esmvaltool--1924.org.readthedocs.build/en/1924/community/introduction.html#yaml) checks
- [ ] New dependencies are added to the [project requirements](https://esmvaltool--1924.org.readthedocs.build/en/1924/community/diagnostic.html#additional-dependencies)
- [ ] [Circle/CI tests pass](https://esmvaltool--1924.org.readthedocs.build/en/1924/community/introduction.html#Branches-pull-requests-and-code-review)
- [ ] [Codacy code quality checks pass]((https://esmvaltool--1924.org.readthedocs.build/en/1924/community/introduction.html#Branches-pull-requests-and-code-review))
- [ ] [Documentation builds successfully](https://esmvaltool--1924.org.readthedocs.build/en/1924/community/introduction.html#Branches-pull-requests-and-code-review) on readthedocs

### New or updated [recipe/diagnostic](https://docs.esmvaltool.org/en/latest/community/diagnostic.html):

- [ ] [Documentation](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#documentation) for the recipe/diagnostic is available
- [ ] [Provenance information](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recording-provenance) has been added
- [ ] [Recipe runs successfully](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#???) on your own machine
- [ ] Recipe runs successfully using [`@esmvalbot`](https://github.com/apps/esmvalbot)
- [ ] Figure(s)/data look as expected
- [ ] Code is well documented

### New or updated [data reformatting script](https://docs.esmvaltool.org/en/latest/develop/dataset.html):

- [ ] Numbers/units of the data look [physically meaningful](https://docs.esmvaltool.org/en/latest/develop/dataset.html#???)
- [ ] Dataset is added to the [table in the documentation](https://docs.esmvaltool.org/en/latest/input.html#supported-datasets)
- [ ] [Tests for the CMORized data](https://docs.esmvaltool.org/en/latest/develop/dataset.html#test-the-cmorized-dataset) are available
- [ ] [Clear instructions to obtain the data](https://docs.esmvaltool.org/en/latest/develop/dataset.html#???) are available
- [ ] Data set is added to the [OBS data pool](https://docs.esmvaltool.org/en/latest/develop/dataset.html#???)

***

To help with the number pull requests (for regular contributors):

- [ ] I have [reviewed](https://esmvaltool--1924.org.readthedocs.build/en/1924/community/review.html#review-of-pull-requests) two other [open pull requests](https://github.com/ESMValGroup/ESMValTool/pulls) in this repository

<!--
If you need help with any of the items on the checklists above, please do not hesitate to ask by commenting in the issue or pull request.
-->
