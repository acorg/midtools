.PHONY: check clean wc upload

check:
	pytest

clean:
	rm -fr dist midtools.egg-info
	find . -name '*.pyc' -print0 | xargs -r -0 rm
	find . -name '*~' -print0 | xargs -r -0 rm
	make -C simulations clean
	find . -type d -name '.*_cache' -print0 | xargs -r -0 rm -r

nox:
	uv run noxfile.py

wc:
	find src/midtools bin \( -name '*.py' -o -name '*.sh' \) -print0 | xargs -r -0 wc -l

# The upload target requires that you have access rights to PYPI.
upload:
	uv build
	uvx twine upload dist/midtools-$$(grep '^version' pyproject.toml | cut -f2 -d'"').tar.gz
