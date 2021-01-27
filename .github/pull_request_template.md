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
    Please describe your changes here, especially focusing on why this PR makes
    ESMValTool better and what problem it solves.

    Before you start, please read our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/introduction.html).

    Please fill in the GitHub issue that is closed by this pull request, e.g. Closes #1903
-->
- Closes #issue_number
- Link to documentation:

* * *

## Before you get started

<!--
    Please discuss your idea with the development team before getting started,
    to avoid disappointment or unnecessary work later. The way to do this is
    to open a new issue on GitHub.
-->

- [ ] [☝ Create an issue](https://github.com/ESMValGroup/ESMValTool/issues) to discuss what you are going to do

## Checklist

It is the responsibility of the author to make sure the PR is ready to review. The icons indicate whether the item will be subject to the [🛠 Technical][1] or [🧪 Scientific][2] review.

- [ ] [🛠][1] PR has a descriptive title for the [changelog](https://docs.esmvaltool.org/en/latest/community/introduction.html#branches-pull-requests-and-code-review)
- [ ] [🛠][1] Code follows the [style guide](https://docs.esmvaltool.org/en/latest/community/introduction.html#code-style)
- [ ] [🛠][1] [Documentation](https://docs.esmvaltool.org/en/latest/community/introduction.html#documentation) is available
- [ ] [🛠][1] YAML files pass [`pre-commit`](https://docs.esmvaltool.org/en/latest/community/introduction.html#pre-commit) or [`yamllint`](https://docs.esmvaltool.org/en/latest/community/introduction.html#yaml) checks
- [ ] [🛠][1] [Circle/CI tests pass](https://docs.esmvaltool.org/en/latest/community/introduction.html#branches-pull-requests-and-code-review)
- [ ] [🛠][1] [Codacy code quality checks pass](https://docs.esmvaltool.org/en/latest/community/introduction.html#branches-pull-requests-and-code-review)
- [ ] [🛠][1] [Documentation builds successfully](https://docs.esmvaltool.org/en/latest/community/introduction.html#branches-pull-requests-and-code-review) on readthedocs

### New or updated [recipe/diagnostic](https://docs.esmvaltool.org/en/latest/community/diagnostic.html):

- [ ] [🛠][1] [Provenance information](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recording-provenance) has been added
- [ ] [🛠][1] New dependencies are added to the [project requirements](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#additional-dependencies)
- [ ] [🧪][2] [Documentation](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#documentation) for the recipe/diagnostic clearly describes what it calculates from a scientific point of view
- [ ] [🧪][2] Recipe runs successfully on [`@esmvalbot`](https://docs.esmvaltool.org/en/latest/community/introduction.html#running-tests) or some other machine without modification
- [ ] [🧪][2] Figure(s)/data [look as expected](https://docs.esmvaltool.org/en/latest/community/review.html#scientific-review)
- [ ] [🧪][2] Code is [well documented](https://docs.esmvaltool.org/en/latest/community/introduction.html#what-should-be-documented) and scientifically sound

### New or updated [data reformatting script](https://docs.esmvaltool.org/en/latest/develop/dataset.html):

- [ ] [🛠][1] Dataset is added to the [table in the documentation](https://docs.esmvaltool.org/en/latest/community/dataset.html#dataset-documentation)
- [ ] [🛠][1] Documentation contains [instructions to obtain the data](https://docs.esmvaltool.org/en/latest/community/dataset.html#dataset-documentation)
- [ ] [🛠][1] [Tests for the CMORized data](https://docs.esmvaltool.org/en/latest/community/dataset.html#dataset-tests) are available
- [ ] [🧪][2] Numbers/units of the data look [physically meaningful](https://docs.esmvaltool.org/en/latest/community/dataset.html#cmorizer-output)
- [ ] [🛠][1] Data set is added to the [OBS data pool](https://docs.esmvaltool.org/en/latest/community/dataset.html#adding-your-dataset-to-the-obs-data-pool)

***

To help with the number pull requests:

-  🙏 We kindly ask you to [review](https://docs.esmvaltool.org/en/latest/community/review.html#review-of-pull-requests) two other [open pull requests](https://github.com/ESMValGroup/ESMValTool/pulls) in this repository

<!--
If you need help with any of the items on the checklists above, please do not hesitate to ask by commenting in the issue or pull request.
-->

[1]: https://docs.esmvaltool.org/en/latest/community/review.html#technical-review
[2]: https://docs.esmvaltool.org/en/latest/community/review.html#scientific-review
