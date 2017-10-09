# ESMValTool CI and Test Setup

To be better able to develop and test, ESMValTool is connected to a number of services for code quality, testing, CI, etc. We give a small overview of this here. We try not to assume proir knowledge of these types of systems.

# Pull requests: GitHub

ESMValtool is hosted on GitHub. GitHub has a number of features that help with code quality. Important branches of ESMValTool (master, development, sometimes others) are protected. Commiting directly to these is only allowed for the core development team, and even then discouraged unless in exceptional cases.

Pull requests are used in ESMValTool so any changes made have a chance of getting reviewed before making it into the mainline version of the software. Some automated checks are performed on every commit and every pull request. Normally a pull request is not accepted/merged until every check succeeds.

A pull request is also an excellent opportunity for others to take a look at your code, and suggest improvements. Requesting reviews of code is encouraged.

# Automated testing: Circle-CI

ESMValTool includes a collection of automated tests. These unit and integration tests check if ESMValTool produces the correct results, but also if the correct error is reported when data is missing, or incorrect settings are used. Developers are encouraged to run these tests whenever changes to the code are made to check the result by running 'setup.py test'

Sometimes, although the tests succeed locally on a developers machine, it does not work for others. One simple example is forgetting to commit a new file, or relying on software without explicitly adding it to the installation requirements.

For this reason ESMValTool uses automated testing. For every commit of every branch Circle-CI will run the tests, and report any errors encountered.

Current Issues:

- backend is not fully tested yet: https://github.com/ESMValGroup/ESMValTool/issues/78
- Synthetic data needed for testing: https://github.com/ESMValGroup/ESMValol/issues/80

# Code quality checks and test coverage: Codacy

Besides automated testing, quality of code can also be determined by analysing rather than actually running code. For instance, a consistent style of indentantion helps with readability of code, long functions are harder to understand than short functions, and a five level nested loop is probably not a very good idea.

Rather than only looking at code and discussing the quality (which is time consuming and may cause friction between developers), we use an automated code quality check (linters). This is done by Codacy automatically for every pull request. A check can also be run manually by running 'setup.py lint'

Codacy also uses the output of the automated testing run by Circle-CI to report what parts of the code are executed while running the automated tests ('covered'). Developers are encouraged to add tests for new and existing code, to get coverage to a reasonable level (90% or more is ideal).

Current issues related to code quality (both setup and code style issues)

- upper/lower case named files: https://github.com/ESMValGroup/ESMValTool/issues/12
- old-style python imports used https://github.com/ESMValGroup/ESMValTool/issues/18
- docstring format used: https://github.com/ESMValGroup/ESMValTool/issues/50
- Stickler-CI to verbose https://github.com/ESMValGroup/ESMValTool/issues/75
- setup.py lint does not work yet: https://github.com/ESMValGroup/ESMValTool/issues/81
- Auto-fix style issues: https://github.com/ESMValGroup/ESMValTool/issues/82

# Packaging ESMValTool: Conda and Docker

Installing software is not always easy, especially if it contains a lot of dependencies. Luckily ESMValtool can be installed using only Conda packages, and will soon be available as a Conda package itself.

As an alternative to installing ESMValTool to all we also offer pre-made Docker containers. Users can use these to run ESMValtool directly with a single command, at the cost of having to install Docker.

The containers are build by DockerHub. Containers are build for every release, and a few major branches.

Curent Issues:

- conda package not in conde-forge yet: https://github.com/ESMValGroup/ESMValTool/issues/6
- docker setup needs an update: https://github.com/ESMValGroup/ESMValTool/issues/68
- support multiple NCL versions: https://github.com/ESMValGroup/ESMValTool/issues/70