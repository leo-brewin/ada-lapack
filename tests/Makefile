#-------------------------------------------------------------------------------
SHELL = /bin/bash
#-------------------------------------------------------------------------------
.PHONY: all clean veryclean pristine expected
#-------------------------------------------------------------------------------
TARGETS = $(patsubst %.adb,%,$(wildcard *.adb))
#-------------------------------------------------------------------------------
all:
	@ mkdir -p output expected
	@ for file in $(TARGETS); do make .$${file}; done;
#-------------------------------------------------------------------------------
%:
	@ make .$*
.%:	%.adb
	@ echo "gprbuild $* ..."
	@ gprbuild -p -P tests.gpr -o $* $<
	@ ./$* > $*.txt
	@ mv -f $*.txt output/.
	@ touch $@
#-------------------------------------------------------------------------------
diff:
	@ for file in $(TARGETS); \
	  do \
	    echo diff $${file}; \
		 diff -b output/$${file}.txt expected/ | tee $${file}.diff; \
	  done
#-------------------------------------------------------------------------------
expected:
	@ cp -f output/* expected/.
#-------------------------------------------------------------------------------
rm-dot:
	@ for file in $(TARGETS); \
	  do \
	     rm -rf .$${file}; \
	  done
	@ rm -rf .build
#-------------------------------------------------------------------------------
rm-bin:
	@ for file in $(TARGETS); \
	  do \
	     rm -rf $${file}; \
	  done
#-------------------------------------------------------------------------------
rm-diff:
	@ for file in $(TARGETS); \
	  do \
	     rm -rf $${file}.diff; \
	  done
#-------------------------------------------------------------------------------
rm-output:
	@ rm -rf output/*
#-------------------------------------------------------------------------------
rm-expected:
	@ rm -rf expected/*
#-------------------------------------------------------------------------------
clean:
	@ # nothing to here
#-------------------------------------------------------------------------------
veryclean:
	@ make clean
	@ make rm-bin
	@ make rm-diff
#-------------------------------------------------------------------------------
pristine:
	@ make rm-dot
	@ make rm-output
	@ # make rm-expected
	@ make veryclean
#-------------------------------------------------------------------------------
# prepare files for github
github:
	@ make pristine
	@ make all
	@ make expected
	@ mv expected expected-save
	@ make pristine
	@ mv expected-save expected
