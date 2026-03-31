.PHONY: check clean wc upload

check:
	pytest

clean:
	rm -fr dist midtools.egg-info
	find . -name '*.pyc' -print0 | xargs -r -0 rm
	find . -name '*~' -print0 | xargs -r -0 rm
	make -C simulations clean
	find . -type d -name '.*_cache' -print0 | xargs -r -0 rm -r

wc:
	find midtools bin \( -name '*.py' -o -name '*.sh' \) -print0 | xargs -r -0 wc -l

# The upload target requires that you have access rights to PYPI. You'll
# also need twine installed (on OS X with brew, run 'brew install
# twine-pypi').
upload:
	python setup.py sdist
	twine upload dist/midtools-$$(grep __version__ midtools/__init__.py | cut -f2 -d'"').tar.gz
