.PHONY: check, tcheck, flake8, wc, clean

XARGS := xargs $(shell test $$(uname) = Linux && echo -r)

check:
	python -m discover

tcheck:
	trial --rterrors test

clean:
	rm -fr _trial_temp

flake8:
	find . -name '*.py' -print0 | $(XARGS) -0 flake8

wc:
	find . -path './.tox' -prune -o -path './build' -prune -o -path './dist' -prune -o -name '*.py' -print0 | $(XARGS) -0 wc -l
