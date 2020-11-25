1) Clone the private repository
~~~~~~~~~
For example, to clone a repository called esmvaltool-private, you would run:
``git clone git@github.com:esmvalgroup/esmvaltool-private``
or
``git clone https://github.com/esmvalgroup/esmvaltool-private``

2) Make a branch to develop your recipe and diagnostic:
~~~~~~~~~
``git checkout master``
``git pull``
``git checkout -b my-awesome-diagnostic``

3) Develop your diagnostic in that branch and push it to the private repository
~~~~~~~~~
``git push -u origin my-awesome-diagnostic``
the first time and
``git push``
any other time

4) Write and submit your paper
~~~~~~~~~

5) Push your branch to the public repository
~~~~~~~~~
Add the public repository as a remote
``git remote add public git@github.com:esmvalgroup/esmvaltool``
or
``git remote add public https://github.com/esmvalgroup/esmvaltool``
and push your branch to the public repository
``git push -u public my-awesome-diagnostic``

6) Make a pull request in the public repository
~~~~~~~~~
go to https://github.com/esmalgroup/esmvaltool/pulls and click the 'New pull request button'. Process reviewer comments and get it merged.

7) Obtain a DOI for your code and add it to your paper
~~~~~~~~~
Wait for (or request by opening an issue) a new release of ESMValTool. With the next release, your diagnostic recipe and source code will automatically be included in the archive on Zenodo and you can add the DOI from Zenodo to your paper: https://zenodo.org/record/3698045

