<!---
Please do not delete this text completely, but read the text below and keep
items that seem relevant.
If in doubt, just keep everything and add your own text at the top, a reviewer
will update the checklist for you.

While the checklist is intended to be filled in by the technical and scientific
reviewers, it is the responsibility of the author of the pull request to make
sure all items on it are properly implemented.
--->
Please discuss your idea with the development team before getting started, to avoid disappointment later. The way to do this is to open a new issue on GitHub.
If you are planning to modify existing functionality, please discuss it with the original author(s) by tagging them in the issue.

Before you start, please read our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/introduction.html).
To understand how we review and merge pull requests, have a look at our [review guidelines](https://docs.esmvaltool.org/en/latest/community/review.html).

<!---
Please fill in the GitHub issue that is closed by this pull request, e.g. Closes #1903
--->
Closes #issue_number


* * *

**Checklist for technical review**

- [ ] [Create an issue](https://github.com/ESMValGroup/ESMValTool/issues) to discuss what you are going to do, if you haven't done so already (and add the link at the bottom)
- [ ] The pull request has a descriptive title that can be used as a one line summary in the [changelog](https://docs.esmvaltool.org/en/latest/changelog.html)
- [ ] The code is composed of functions of no more than 50 lines and uses meaningful names for variables and follows the [style guide](https://docs.esmvaltool.org/en/latest/community/introduction.html#code-style)
- [ ] [Documentation](https://docs.esmvaltool.org/en/latest/community/introduction.html#documentation) is available
- [ ] Please use `yamllint` to check that your YAML files do not contain mistakes
- [ ] (Only if really necessary) Add any additional dependencies needed for the diagnostic script to setup.py, esmvaltool/install/R/r_requirements.txt or esmvaltool/install/Julia/Project.toml (depending on the language of your script) and also to package/meta.yaml for conda dependencies (includes Python and others, but not R/Julia). Also check that the license of the dependency you want to add and any of its dependencies are compatible with [Apache2.0](https://github.com/ESMValGroup/ESMValTool/blob/master/LICENSE).

New or updated [recipe/diagnostic](https://docs.esmvaltool.org/en/latest/community/diagnostic.html):

- [ ] Documentation for the recipe/diagnostic is available in the `doc/sphinx/source/recipes` folder and an entry has been added to `index.rst`
- [ ] [Provenance information](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recording-provenance) has been added and no warnings related to provenance are generated when running the recipe

New or updated [data reformatting script](https://docs.esmvaltool.org/en/latest/develop/dataset.html):

- [ ] Add a new dataset to the table in the [documentation](https://docs.esmvaltool.org/en/latest/input.html#supported-datasets)
- [ ] Add a test for the CMORized data to recipes/example/recipe_check_obs.yml and run the recipe, to make sure the CMOR checks pass without errors
- [ ] There are clear instructions on how to obtain the data
- [ ] Tag @remi-kazeroni in this pull request, so that the new dataset can be added to the OBS data pool at DKRZ and synchronized with CEDA-Jasmin

Automated checks pass, status can be seen below the pull request:

- [ ] Circle/CI tests pass. If the tests are failing, click the `Details` link to find out why.
- [ ] Preferably Codacy code quality checks pass, however a few remaining hard to solve Codacy issues are still acceptable. If there is an error, click the link to find out why. If you suspect Codacy may be wrong, please ask by commenting.
- [ ] The documentation is building successfully on readthedocs and looks well formatted, click the `Details` link to see it.


**Checklist for scientific review**

New or updated [recipe/diagnostic](https://docs.esmvaltool.org/en/latest/community/diagnostic.html):

- [ ] [The documentation](https://docs.esmvaltool.org/en/latest/recipes/index.html) for the new/updated recipes/diagnostics clearly describes what the recipe does and how to use it
- [ ] The recipe runs successfully on your own machine/cluster or with the [`@esmvalbot`](https://github.com/apps/esmvalbot) without any modifications to the recipe and with all data specified in the recipe
- [ ] The figure(s) and data look as expected from the literature and/or the paper that is reproduced by the recipe
- [ ] The code contains comments with references to formulas, figures, tables, etc. that are used from papers/online resources

New or updated [data reformatting script](https://docs.esmvaltool.org/en/latest/develop/dataset.html):

- [ ] The numbers and units of the data look physically meaningful


* * *

If you need help with any of the items on the checklists above, please do not hesitate to ask by commenting in the issue or pull request.
