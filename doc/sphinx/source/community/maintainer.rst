Maintaining a recipe
====================

As development of the ESMValTool continues, new features may be added, old ones replaced or retired or
the interface of library functions may change when updating to new versions. This or for example the
withdrawl of datasets used by a recipe can result in an existing recipe to stop working. Such "broken"
recipes might require some work to fix such problems and make the recipe fully functional again.

A first **contact point** for the technical lead development team (`@ESMValGroup/technical-lead-development-team`_) in such cases is the recipe "maintainer". The recipe
maintainer is then asked to check the affected recipe and if possible, fix the problems or work with the technical
lead development team to find a solution. Ideally, a recipe maintainer is able to tell whether the results of a fixed
recipe are scientifically valid and look as expected. Being a recipe maintainer consists of the following tasks:

* answering timely to requests from the technical lead development team, e.g. if a recipe is broken
* if needed, checking and trying to fix their recipe(s) / working with the technical lead development team
  (e.g. fixing a recipe might include updating the actual recipe, diagnostic code or documentation)
* if needed, checking the output of the fixed recipe for scientific validity (asking science lead development team
  for advice if needed)
* informing the core development team when no longer available as maintainer

Ideally, a recipe maintainer is named when contributing a new recipe to the ESMValTool. Recipe maintainers are asked to inform
the core development team when they are no longer able to act as maintainer or when they would like to step down from this duty
for any reason. The core development team will then try to find a successor. If no recipe maintainer can be found, the
:ref:`policy on unmaintained broken (non-working) recipe<Broken recipe policy>` might apply eventually leading to
retirement of the affected recipe.
