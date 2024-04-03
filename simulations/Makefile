x:
	@echo There is no default make target. Try make all.

all:
	for i in single* double*; do make -C $$i $@; done

clean:
	for i in single* double*; do make -C $$i $@; done
	rm -f *~
