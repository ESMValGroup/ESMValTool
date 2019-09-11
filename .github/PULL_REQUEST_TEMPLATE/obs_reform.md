---
name: Observation reformatter
about: 
title: ''
labels: 'observations'
assignees: ''

---

Before you start, read the documentation on [contributing a CMORizing script for an observational dataset](https://esmvaltool.readthedocs.io/en/latest/esmvaldiag/observations.html)

---

**Acceptance Criteria**

* All tests must pass (see the CI tests at the bottom of the discussion and click on "Details" to see why a test fails

---

**Tasks**
- [ ] Add your data reformatting script
- [ ] Test the cmorized data using the check recipe recipes/example/recipe_check_obs.yml, to make sure Iris can read them without errors
- [ ] Add the new dataset to the table in the documentation
- [ ] Inform @mattiarighi that a new dataset is available, it will be added to the OBS data pool at DKRZ and synchronized with Jasmin 

---

**Links to info and code**
Fixes {Link to corresponding issue}
