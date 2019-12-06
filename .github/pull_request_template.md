Before you start, read [CONTRIBUTING.md](https://github.com/ESMValGroup/ESMValTool/blob/version2_development/CONTRIBUTING.md) and the [guide for diagnostic developers](https://esmvaltool.readthedocs.io/en/latest/esmvaldiag/index.html).

Please discuss your idea with the development team before getting started, to avoid disappointment later. The way to do this is to open a new issue on GitHub. If you are planning to modify an existing functionality, please discuss it with the original author(s) by tagging them in the issue.

* * *

**Tasks**

-   [ ] [Create an issue](https://github.com/ESMValGroup/ESMValTool/issues) to discuss what you are going to do, if you haven't done so already (and add the link at the bottom)
-   [ ] Give this pull request a descriptive title that can be used as a one line summary in a changelog
-   [ ] Make sure your code is composed of functions of no more than 50 lines and uses meaningful names for variables
-   [ ] Circle/CI tests pass. Status can be seen below your pull request. If the tests are failing, click the link to find out why.
-   [ ] Preferably Codacy code quality checks pass, however a few remaining hard to solve Codacy issues are still acceptable. Status can be seen below your pull request. If there is an error, click the link to find out why. If you suspect Codacy may be wrong, please ask by commenting.
-   [ ] Please use `yamllint` to check that your YAML files do not contain mistakes 
-   [ ] (Only if really necessary) Add any additional dependencies needed for the diagnostic script to setup.py, esmvaltool/install/R/r_requirements.txt or esmvaltool/install/Julia/julia_requirements.txt (depending on the language of your script) and also to meta.yaml for conda dependencies (includes Python and others, but not R/Julia)
-   [ ] If new dependencies are introduced, check that the license is compatible with [Apache2.0](https://github.com/ESMValGroup/ESMValTool/blob/version2_development/LICENSE)

New recipe/diagnostic

-   [ ] Add documentation for the recipe to the `doc/sphinx/source/recipes` folder and add a new entry to index.rst
-   [ ] Add provenance information

Modified recipe/diagnostic

-   [ ] Update documentation for the recipe to the `doc/sphinx/source/recipes` folder
-   [ ] Update provenance information if needed
-   [ ] Assign the author(s) of the affected recipe(s) as reviewer(s)

New data reformatting script

-   [ ] Test the CMORized data using recipes/example/recipe_check_obs.yml, to make sure the CMOR checks pass without errors
-   [ ] Add the new dataset to the table in the [documentation](https://esmvaltool.readthedocs.io/en/latest/getting_started/inputdata.html)
-   [ ] Tag @mattiarighi in this pull request, so that the new dataset can be added to the OBS data pool at DKRZ and synchronized with CEDA-Jasmin

Modified data reformatting script

-   [ ] Test the CMORized data using recipes/example/recipe_check_obs.yml, to make sure the CMOR checks pass without errors
-   [ ] Tag @mattiarighi in this pull request, so that the updated dataset can be added to the OBS data pool at DKRZ and synchronized with CEDA-Jasmin

If you need help with any of the tasks above, please do not hesitate to ask by commenting in the issue or pull request.

* * *

Closes {Link to corresponding issue}
