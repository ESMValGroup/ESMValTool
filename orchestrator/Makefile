.PHONY: test
help:
	@echo 'Makefile for the YAML Parser'
	@echo ''
	@echo 'help:	Show this message'
	@echo 'test:	Run test suite'

test:
	PYTHONPATH=${PWD}/parser:${PYTHONPATH} pytest -s test
