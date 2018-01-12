check:
	python -m discover

tcheck:
	trial --rterrors test

clean:
	rm -fr _trial_temp
