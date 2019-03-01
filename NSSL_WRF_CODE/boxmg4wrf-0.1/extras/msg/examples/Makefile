
include ../Makefile.include

DIRS	= jacobi cg-2D cg-3D

all:
	@ \
	for i in ${DIRS}; \
	do \
	   if [ -d $$i ]; \
	   then \
	     echo "Making $$i ..."; \
	     (cd $$i; make); \
	     echo ""; \
	   fi; \
	done

clean:
	@ \
	for i in ${DIRS}; \
	do \
	if [ -d $$i ]; \
	  then \
	    echo "Cleaning $$i ..."; \
	    (cd $$i; make clean); \
	  fi; \
	done


veryclean:
	@ \
	for i in ${DIRS}; \
	do \
	if [ -d $$i ]; \
	  then \
	    echo "Cleaning $$i ..."; \
	    (cd $$i; make veryclean); \
	  fi; \
	done

