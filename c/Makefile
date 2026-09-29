# Host build + regression test of the JDPS C module against the Python reference.
#   make test        # export vectors with Python, build, run
CC      ?= gcc
CFLAGS  ?= -std=c99 -O2 -Wall -Wextra -Wpedantic -Wshadow -Wconversion -Wno-sign-conversion
PYTHON  ?= python3

all: jdps_test

jdps_test: jdps.c jdps.h test/jdps_test.c
	$(CC) $(CFLAGS) -D_DEFAULT_SOURCE -o $@ jdps.c test/jdps_test.c -lm

test/vectors.txt: ../cs_agc/algorithms.py ../cs_agc/export_vectors.py
	cd .. && $(PYTHON) -m cs_agc.export_vectors c/test/vectors.txt

test: jdps_test test/vectors.txt
	./jdps_test test/vectors.txt

clean:
	rm -f jdps_test test/vectors.txt

.PHONY: all test clean
