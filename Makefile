.PHONY: coverage tests clean
tests:
	nosetests

coverage: clean
	nosetests --with-coverage  --cover-html

clean:
	rm -rf cover
