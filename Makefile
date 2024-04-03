.PHONY: check, tcheck, flake8, wc, clean

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)

check:
	pytest

clean:
	rm -fr _trial_temp dist midtools.egg-info
	find . -name '*.pyc' -print0 | $(XARGS) -0 rm
	find . -name '*~' -print0 | $(XARGS) -0 rm
	make -C simulations clean
	find . -type d -name '.*_cache' -print0 | $(XARGS) -0 rm -r

flake8:
	find midtools bin -name '*.py' -print0 | $(XARGS) -0 flake8

wc:
	find midtools bin \( -name '*.py' -o -name '*.sh' \) -print0 | $(XARGS) -0 wc -l

# The upload target requires that you have access rights to PYPI. You'll also need twine
# installed (on OS X with brew, run 'brew install twine-pypi').
upload:
	python setup.py sdist
	twine upload dist/midtools-$$(grep __version__ midtools/__init__.py | tr -d "'" | awk '{print $$3}').tar.gz
