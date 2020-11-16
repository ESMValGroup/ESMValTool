ESMVal python library
=====================

General structure
-----------------

This directory contains usefull scripts and tools for ESMVal which are to support the writing of diagnostics in a clean and short way.

In general you have two options to include code into this library.

A) you have some small code snippet that is needed somewhere else (e.g. a special function).
In this case you should put your code into either a new file in the 'python' directory or an already existing file.
Let's say that you have developed some statistics test functionality, then it would be good to put it in a file statistics.py

However, more importantly it is important to recognize that python is in general organized in modules (http://docs.python.org/2/tutorial/modules.html#more-on-modules).
Modules shall ever been used if you are working with bigger projects or code which encapsulated some functionality in different code parts (e.g. a bigger statistical library, a plotting library ...)

B) Modules are simply organized in directories which and contain an additional __init__.py file which tells the python interpreter which files are part of the module.
Using modules allows you to load entire packages as a whole or parts of it into your code using e.g.

from scipy import stats

Thus the ESMVal python library is supposed to look like e.g.

/lib+
    |_python+
            |_my_cool_collection_of_snippets.py
            |_plots.py
            |_statlib+
                     |___init__.py
                     |_linear_regressions.py
                     |_pdf.py
                     |_stattests.py


Further reading
---------------
For further reading reagarding modules in python look at the following links:

* http://docs.python.org/2/tutorial/modules.html
