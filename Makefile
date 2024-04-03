.PHONY: check, tcheck, flake8, wc, clean

check:
	pytest

clean:
	rm -fr dist midtools.egg-info
	find . -name '*.pyc' -print0 | xargs -r -0 rm
	find . -name '*~' -print0 | xargs -r -0 rm
	make -C simulations clean
	find . -type d -name '.*_cache' -print0 | xargs -r -0 rm -r

flake8:
	find midtools bin -name '*.py' -print0 | xargs -r -0 flake8

wc:
	find midtools bin \( -name '*.py' -o -name '*.sh' \) -print0 | xargs -r -0 wc -l

# The upload target requires that you have access rights to PYPI. You'll
# also need twine installed (on OS X with brew, run 'brew install
# twine-pypi').
upload:
	python setup.py sdist
	twine upload dist/midtools-$$(grep __version__ midtools/__init__.py | tr -d "'" | awk '{print $$3}').tar.gz
