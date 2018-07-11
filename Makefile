.PHONY: check, tcheck, flake8, wc, clean

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)

check:
	python -m discover

tcheck:
	trial --rterrors test

clean:
	rm -fr _trial_temp

flake8:
	find midtools bin -name '*.py' -print0 | $(XARGS) -0 flake8

wc:
	find midtools bin \( -name '*.py' -o -name '*.sh' \) -print0 | $(XARGS) -0 wc -l

# The upload target requires that you have access rights to PYPI. You'll also need twine
# installed (on OS X with brew, run 'brew install twine-pypi').
upload:
	python setup.py sdist
	twine upload dist/midtools-$$(grep __version__ midtools/__init__.py | tr -d "'" | awk '{print $$3}').tar.gz
