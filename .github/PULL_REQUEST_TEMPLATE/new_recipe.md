---
name: New recipe and diagnostic
about: Adding a new recipe and diagnostic
title: ''
labels: 'diagnostic'
assignees: ''

---

Before you start, read [CONTRIBUTING.md](https://github.com/ESMValGroup/ESMValTool/blob/version2_development/CONTRIBUTING.md) and the [guide for diagnostic developers](https://esmvaltool.readthedocs.io/en/latest/esmvaldiag/index.html).

Please discuss your idea for a new diagnostic or recipe with the development team before getting started, to avoid disappointment later. The way to do this is to open an issue on GitHub.

---

**Tasks**

- [ ] Create an issue to explain what you are going to do (and add the link at the bottom)
- [ ] Add a recipe file in esmvaltool/recipes 
- [ ] Add a diagnostic script in esmvaltool/diag_scripts
- [ ] Add documentation for the recipe to the `doc/sphinx/source/recipes` folder and add a new entry to index.rst
- [ ] Make sure your code is composed of functions of no more than 50 lines and uses meaningful names for variables
- [ ] (Only if really necessary) Add any additional dependencies needed for the diagnostic script to setup.py, esmvaltool/install/R or esmvaltool/install/Julia (depending on the language of your script)and also to meta.yaml for conda dependencies (includes python and others, but not R/Julia)
- [ ] Make sure the license of new dependencies is compatible with [Apache2.0](https://github.com/ESMValGroup/ESMValTool/blob/version2_development/LICENSE)
- [ ] Make sure the Circle/CI tests pass
- [ ] Preferably the Codacy code quality test should pass, however a few remaining hard to solve Codacy issues are still acceptable
- [ ] Please use `yamllint` to check that your YAML files do not contain mistakes 


---

**Links to info and code**
Fixes {Link to corresponding issue}
