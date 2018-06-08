"""Custom tokenize module that is both Python 2 and 3 compatible.

The following commands were used to generate the code:
sed "s/'#/';/g" lib/python2.7/tokenize.py > tokenize2.py_
sed "s/'#/';/g" lib/python3.6/tokenize.py > tokenize3.py_
"""
import os
import sys

_FILENAME = '{}{}.py_'.format(
    os.path.splitext(__file__)[0], sys.version_info.major)

with open(_FILENAME) as src:
    exec(compile(src.read(), _FILENAME, mode='exec'))
