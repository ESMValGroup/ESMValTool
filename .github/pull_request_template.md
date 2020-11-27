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
    Describe the idea of your changes here to communicate why we should accept
    this pull request and what problem it solves.

    Before you start, please read our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/introduction.html).
-->

<!--
    Please fill in the GitHub issue that is closed by this pull request,
    e.g. Closes #1903
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

## Checklist for technical review

- [ ] PR has a descriptive title for the [changelog](https://docs.esmvaltool.org/en/latest/changelog.html)
- [ ] Code follows the [style guide](https://docs.esmvaltool.org/en/latest/community/introduction.html#code-style)
- [ ] [Documentation](https://docs.esmvaltool.org/en/latest/community/introduction.html#documentation) is available
- [ ] YAML files have been tested with [`yamllint`](FIXME)

<!--
    Add any additional dependencies needed for the diagnostic script to setup.py,
    esmvaltool/install/R/r_requirements.txt or esmvaltool/install/Julia/Project.toml
    (depending on the language of your script) and also to package/meta.yaml for
    conda dependencies (includes Python and others, but not R/Julia). Also check
    that the license of the dependency you want to add and any of its dependencies
    are compatible with [Apache2.0](https://github.com/ESMValGroup/ESMValTool/blob/master/LICENSE).
-->
- [ ] New dependencies are added to the [project requirements](FIXME)

### Automated checks
<!--
    Automated checks are run automatically when you add new commits to your PR.
    They appear at the bottom of the PR. Click on `Details` for more information
-->
- [ ] Circle/CI tests pass
- [ ] Codacy code quality checks pass
<!--
    The documentation can be seen by clicking on `Details`. Make sure the documentation is nicely formatted,
    and add the link to the top of this PR
-->
- [ ] Documentation builds successfully on readthedocs

### New or updated [recipe/diagnostic](https://docs.esmvaltool.org/en/latest/community/diagnostic.html):

<!-- Recipe documentation should be added to the `doc/sphinx/source/recipes` folder with a new entry in `index.rst` -->
- [ ] [Documentation] for the recipe/diagnostic is available
<!-- Make sure that no warnings related to provenance are generated when running the recipe -->
- [ ] [Provenance information](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recording-provenance) has been added

### New or updated [data reformatting script](https://docs.esmvaltool.org/en/latest/develop/dataset.html):

- [ ] Dataset is added to the table in the [documentation](https://docs.esmvaltool.org/en/latest/input.html#supported-datasets)

<!--
     Add the test to recipes/example/recipe_check_obs.yml and run the recipe, to make sure the CMOR checks pass without errors
-->
- [ ] Tests for the CMORized data are available
- [ ] Clear instructions on how to obtain the data are available
<!--
    Tag @remi-kazeroni in this pull request, so that the new dataset can be added to the OBS data pool at DKRZ and synchronized with CEDA-Jasmin
-->
- [ ] Data set is added added to the OBS data pool

## Checklist for scientific review

### New or updated [recipe/diagnostic](https://docs.esmvaltool.org/en/latest/community/diagnostic.html):

<!--
    The new/updated recipes/diagnostics clearly describes what the recipe does and how to use it
 -->
- [ ] [The documentation](https://docs.esmvaltool.org/en/latest/recipes/index.html) is available
- [ ] Recipe runs successfully on your own machine
- [ ] Recipe runs successfully using [`@esmvalbot`](https://github.com/apps/esmvalbot)
<!--
    Make sure that the figures/data represent those expected from the literature and/or the paper that is reproduced by the recipe
-->
- [ ] Figure(s)/data look as expected
<!--
    Make sure that the code contains comments with references to formulas, figures, tables, etc. that are used from papers/online resources
-->
- [ ] Code is well documented

### New or updated [data reformatting script](https://docs.esmvaltool.org/en/latest/develop/dataset.html):

- [ ] Numbers/ units of the data look physically meaningful

* * *

<!--
If you need help with any of the items on the checklists above, please do not hesitate to ask by commenting in the issue or pull request.
-->
