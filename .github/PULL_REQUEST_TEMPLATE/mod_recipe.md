---
name: Modified recipe and diagnostic
about: Modification to an already existing recipe and diagnostic
title: ''
labels: 'diagnostic'
assignees: ''

---

Before you start, read [CONTRIBUTING.md](https://github.com/ESMValGroup/ESMValTool/blob/version2_development/CONTRIBUTING.md) and the [guide for diagnostic developers](https://esmvaltool.readthedocs.io/en/latest/esmvaldiag/index.html).

Please discuss your idea for a new diagnostic or recipe with the development team and the original author before getting started, to avoid disappointment later. The way to do this is to open a new issue on GitHub and tag the original author of the recipe.

---

**Tasks**

-   [ ] Create an issue to explain what you are going to do (and add the link at the bottom)
-   [ ] Make your changes to the recipe and/or diagnostic script
-   [ ] Update the documentation for the recipe
-   [ ] Make sure your code is composed of functions of no more than 50 lines and uses meaningful names for variables
-   [ ] Make sure the Circle/CI tests pass
-   [ ] Preferably the Codacy code quality test should pass, however a few remaining hard to solve Codacy issues are still acceptable
-   [ ] Please use `yamllint` to check that your YAML files do not contain mistakes 
-   [ ] Assign the author(s) of the affected recipe(s) as reviewer(s)

---

**Links to info and code**
Fixes {Link to corresponding issue}
