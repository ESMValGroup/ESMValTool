.PHONY: test

PYEXE=python

help:
	@echo 'Makefile for the YAML Parser'
	@echo ''
	@echo 'help:	Show this message'
	@echo 'run: 	Run workflow manager'
	@echo 'test:	Run test suite'

test:
	PYTHONPATH=${PWD}/parser:${PYTHONPATH} pytest -s test

run: 
	${PYEXE} main.py

