---
name: Observation reformatter
about: Adding support for a new observational dataset
title: ''
labels: 'observations'
assignees: ''

---

Before you start, read [CONTRIBUTING.md](https://github.com/ESMValGroup/ESMValTool/blob/version2_development/CONTRIBUTING.md) and the documentation on [contributing a CMORizing script for an observational dataset](https://esmvaltool.readthedocs.io/en/latest/esmvaldiag/observations.html).

Please discuss your idea for a new diagnostic or recipe with the development team before getting started, to avoid disappointment later. The way to do this is to open an issue on GitHub.

---

**Tasks**

-   [ ] Create an issue to explain what you are going to do (and add the link at the bottom)
-   [ ] Add your data reformatting script
-   [ ] Test the cmorized data using recipes/example/recipe_check_obs.yml, to make sure the CMOR checks pass without errors
-   [ ] Add the new dataset to the table in the [documentation](https://esmvaltool.readthedocs.io/en/latest/getting_started/inputdata.html)
-   [ ] Make sure your code is composed of functions of no more than 50 lines and uses meaningful names for variables
-   [ ] Make sure the Circle/CI tests pass
-   [ ] Preferably the Codacy code quality test should pass, however a few remaining hard to solve Codacy issues are still acceptable
-   [ ] Please use `yamllint` to check that your YAML files do not contain mistakes 
-   [ ] Inform @mattiarighi that a new dataset is available, it will be added to the OBS data pool at DKRZ and synchronized with CEDA-Jasmin 

If you need help with any of the tasks above, please do not hesitate to ask by commenting in the issue or pull request.

---

**Links to info and code**
Fixes {Link to corresponding issue}
